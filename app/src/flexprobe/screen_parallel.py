# flexiprobes — probe design for 10x Genomics Flex v1 / Flex v2 / Visium.
# Copyright (C) 2026 Adrian Ilich and contributors
#
# This file is part of flexiprobes.
#
# flexiprobes is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# flexiprobes is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Process-parallel off-target screen: search AND score a chunk in slices.

WHY THIS EXISTS. On a whole-transcriptome run the off-target screen dominates
the wall clock, far ahead of tiling and filtering. ``RunConfig.workers`` does
not touch it: it feeds the tiling ProcessPoolExecutor and nothing else, and
``KmerScreen`` has no thread of its own. Without this module the screen runs
single-threaded however many cores the node has.

WHAT CROSSES THE PROCESS BOUNDARY, AND WHY IT IS THE VERDICTS. The obvious
split — workers search, parent scores — sends raw hits back, and on a human
transcriptome one chunk's hits can run to gigabytes (see
``screen_cache.save_verdicts``). Pickling that per chunk would cost more than
the parallelism returns, and peak memory would be worker-slice + pickle buffer
+ parent merge, i.e. roughly tripled. Scoring INSIDE the worker instead sends
back only the ~11 ``VERDICT_FIELDS`` per probe, some three orders of magnitude
less, so the pipe carries about a megabyte per chunk.

That is also what makes this memory-NEUTRAL rather than memory-hungry. Each of
W workers searches 1/W of the chunk and therefore holds 1/W of its hits; the
sum across workers is what one serial process would hold. The marginal cost of
a worker is the ``kmer_index.verify`` transient — an ``(n_cand, 25)`` int64
index array at 200 bytes per proposed position, ~6 MB for a typical query half
and a few hundred MB for one sitting in a repeat.

WHY FORK IS REQUIRED, NOT MERELY PREFERRED. The reference (~450 MB uint8) and
the index offsets (134 MB) are inherited copy-on-write, and ``index.pos`` (1.42
GB) is mmap'd and therefore shared through the page cache. Their data buffers
are never written during a search, and CPython's refcounting dirties only the
small array header, not the payload — so W workers cost ~6 MB each for the
search, not ~2 GB each. Under ``spawn`` every worker would rebuild all of it
and the memory model above would be false, so ``spawn`` is refused rather than
silently paid for: :func:`available` returns False and the caller stays serial.

SCOPE. k-mer engine only. BLAST-engine scoring can consult ``SubjectWindows``,
whose ``prefetch`` warms a whole chunk's subjects in ONE ``blastdbcmd`` call;
splitting a chunk into slices would fragment that into one spawn per slice,
which is precisely the cost that prefetch exists to avoid.

FORKING A PROCESS THAT MAY HAVE THREADS. CPython warns that ``fork()`` in a
multi-threaded process can deadlock, and numpy's backing libraries do start
thread pools. It is survivable here because a worker's whole life is numpy
array arithmetic over inherited read-only buffers plus a pure-Python scoring
loop — it acquires no lock the parent could have been holding, and it never
returns to the parent's own machinery. What it does NOT survive is
oversubscription: W workers each running a multi-threaded BLAS would put W*N
threads on N cores. Pin the library thread counts in the job script
(``OMP_NUM_THREADS=1``, ``OPENBLAS_NUM_THREADS=1``) so ``screen_workers`` is
the only source of parallelism and the memory model above stays true.

The index build is serial and paid once, before the pool starts, so it caps
the speed-up on small references. ``results.csv`` is expected to be
byte-identical to a serial run (see CORRECTNESS).

CORRECTNESS. ``specificity.batch_evaluate`` is a pure per-candidate loop with
no cross-candidate state, so scoring a contiguous slice returns exactly what
scoring the whole chunk returns for those candidates. The parent then applies
the merged verdicts through ``screen_cache.apply_verdicts`` — the SAME call the
resume path already makes for every chunk reloaded from cache — so the parallel
path rejoins the serial one at a junction that is exercised on every resumed
run rather than at a new one.
"""
from __future__ import annotations

import logging
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from typing import Any, Dict, List, Optional, Sequence

from .screen_cache import VERDICT_FIELDS

logger = logging.getLogger(__name__)

# Set in the PARENT before the pool is created; workers read it as inherited
# state. Deliberately not passed as ``initargs``: that would make the transfer
# look like it works under spawn, where the screen's mmap'd index cannot be
# pickled cheaply and the whole memory argument in the module docstring
# collapses. Requiring fork keeps the sharing explicit.
_CTX: Dict[str, Any] = {}


def available() -> bool:
    """True when the fork start method this module depends on exists."""
    return "fork" in multiprocessing.get_all_start_methods()


def _slice_bounds(n: int, w: int) -> List[range]:
    """Contiguous, near-equal slices covering ``range(n)``.

    Contiguous rather than strided so each worker's candidates come from
    neighbouring transcripts, which keeps its index lookups clustered in the
    same region of the mmap'd position array.
    """
    if n <= 0 or w <= 1:
        return [range(n)] if n > 0 else []
    w = min(w, n)
    step, extra = divmod(n, w)
    out: List[range] = []
    lo = 0
    for i in range(w):
        hi = lo + step + (1 if i < extra else 0)
        out.append(range(lo, hi))
        lo = hi
    return out


def _search_and_score(args: Any) -> Dict[str, Dict[str, Any]]:
    """Worker body: search a candidate slice, score it, return its verdicts."""
    from .specificity import batch_evaluate

    ci, cands = args
    screen = _CTX["screen"]
    hits = screen.search_chunk(cands, ci, _CTX["output_dir"])
    batch_evaluate(
        cands, hits, _CTX["target_ids"], subjects=None, **_CTX["spec_kwargs"]
    )
    # Only the verdicts travel home. `hits` dies with this frame.
    return {c.probe_id: {f: getattr(c, f) for f in VERDICT_FIELDS}
            for c in cands}


class ParallelScreen:
    """Pool that searches+scores one chunk at a time across ``workers`` slices.

    The pool is created once per run, not once per chunk: forking a process
    that holds a 450 MB reference is cheap but not free, and a whole-transcriptome
    run with hundreds of chunks would pay it hundreds of times.
    """

    def __init__(
        self, screen: Any, *, workers: int, output_dir: str,
        target_ids: Any, spec_kwargs: Dict[str, Any],
    ) -> None:
        self.workers = workers
        _CTX.update(
            screen=screen, output_dir=output_dir,
            target_ids=target_ids, spec_kwargs=spec_kwargs,
        )
        # Forked EAGERLY, here, while _CTX is fresh and before the first chunk
        # mutates any candidate: a worker forked later would inherit verdicts
        # written by apply_verdicts and carry state it has no business seeing.
        self._pool: Optional[ProcessPoolExecutor] = ProcessPoolExecutor(
            max_workers=workers,
            mp_context=multiprocessing.get_context("fork"),
        )

    def close(self) -> None:
        """Shut the pool down. Called beside ``screen.close()``."""
        if self._pool is not None:
            self._pool.shutdown(wait=True)
            self._pool = None
        _CTX.clear()

    def __enter__(self) -> "ParallelScreen":
        return self

    def __exit__(self, *exc: Any) -> None:
        self.close()

    def verdicts_for(
        self, ci: int, chunk: Sequence[Any],
    ) -> Dict[str, Dict[str, Any]]:
        """Search + score ``chunk`` in parallel; return every probe's verdict.

        Raises whatever a worker raised. That is deliberate: a screen failure is
        deterministic by construction, so retrying it would only reproduce it.
        A batch scheduler should treat the non-zero exit as final.
        """
        assert self._pool is not None, "ParallelScreen is closed"
        slices = _slice_bounds(len(chunk), self.workers)
        futures = [
            self._pool.submit(_search_and_score, (ci, list(chunk[r.start:r.stop])))
            for r in slices
        ]
        merged: Dict[str, Dict[str, Any]] = {}
        for fut in futures:
            merged.update(fut.result())
        return merged
