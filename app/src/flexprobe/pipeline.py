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

"""Pipeline orchestrator: ingest -> tile -> filter -> BLAST -> specificity -> rank -> chemistry -> export."""
from __future__ import annotations
import copy
import dataclasses
import hashlib
import functools
import json
import logging
import os
import platform
import re
import shutil
import sys
import threading
import time
from collections import Counter
from datetime import datetime, timezone
try:                                  # POSIX only; None on Windows
    import resource
except ImportError:                   # pragma: no cover
    resource = None                   # type: ignore[assignment]
from typing import (
    Any, Callable, Dict, FrozenSet, Iterator, List, NamedTuple, NoReturn,
    Optional, Sequence, Set, Tuple,
)

from . import orientation
from . import reference_qc
from . import engine_identity
from . import screen_cache
from . import screen_parallel
from .reference_qc import EXCLUDE_ALT_SCAFFOLDS_DEFAULT
from .models import (
    RunConfig, ProbeCandidate, ProbePair, RunManifest,
)
from .errors import FlexprobeError, PipelineError, CancellationError
from .io_fasta import parse_fasta
from .tiling import tile_sequence
from . import blast
from .blast import (
    SubjectWindows, ensure_db_for_path, is_db_prefix, resolve_db_prefix,
)
# check_version / discover_blast / run_blastn / parse_hits_tsv are deliberately
# NOT imported by name. They are called by screen.BlastScreen through the
# ``blast`` module, and a name bound here as well would be a trap: a test could
# patch ``pipeline.run_blastn``, see nothing take effect, and conclude the fake
# was wired up. Patch ``flexprobe.blast.<name>`` instead — where it is defined.
from .screen import (ENGINE_BLAST, ENGINE_KMER, Screen, make_screen,
                     screen_labels)
from .chemistry import (
    FLEX_V1_PROBE_BARCODES,
    assemble_pair, get_chemistry, resolve_tn_mode,
)
from .constants import (
    BLASTN_SHORT_MISMATCH_COST, OFF_TARGET_PUSHDOWN_GAP_MARGIN,
    PREPARED_MARKER,
    RECOUNT_EXACT_WARN_HALVES,
    REFERENCE_STRAND_DEFAULT,
    TN_MODE_STRICT_A,
)
from .consensus import (
    StrainAlignment, build_consensus,
    constitutive_start_positions, degenerate_conserved_windows,
    resolve_to_reference, resolve_window_per_record,
    score_probe_against_all_strains, strain_unique_start_positions,
)
from .constants import (
    PROBE_ROLE_CONSENSUS, PROBE_ROLE_VARIANT,
)
from .errors import BlastError
from .ranking import RankingWeights, rank_candidates
from .specificity import batch_evaluate, StrainRule


from ._utils import safe_tm as _safe_tm  # noqa: E402  (shared NaN-honest Tm)
from ._utils import gc_fraction as _gc_fraction, is_nan as _is_nan  # noqa: E402


# Ensembl-style cDNA headers carry ``gene:ENSG00000223997.1`` — grab the token
# after ``gene:`` up to the next whitespace. Used to group isoforms by gene.
_GENE_FIELD_RE = re.compile(r"(?:^|\s)gene:(\S+)")


def _build_gene_groups(
    records: List[Tuple[str, str, str]],
) -> Tuple[Dict[str, str], Dict[str, Set[str]]]:
    """Group input transcripts by gene for isoform-aware specificity.

    Parses each record's FASTA description for an Ensembl-style ``gene:`` field.
    Records without one fall back to their own id (each becomes its own group,
    so coverage is trivially 1.0 and grouping has no effect). Returns
    ``(gene_of, group_members)`` where ``gene_of`` maps transcript id -> gene
    key and ``group_members`` maps gene key -> set of its transcript ids.
    """
    gene_of: Dict[str, str] = {}
    group_members: Dict[str, Set[str]] = {}
    for gid, desc, _seq in records:
        m = _GENE_FIELD_RE.search(desc or "")
        gkey = m.group(1) if m else gid
        gene_of[gid] = gkey
        group_members.setdefault(gkey, set()).add(gid)
    return gene_of, group_members


def _usable_cores() -> int:
    """CPUs this PROCESS may actually use — not the machine's core count.

    `os.cpu_count()` reports the physical node and ignores cgroup/cpuset limits,
    so under a scheduler it over-reports badly: a PBS job with ncpus=16 on a
    64-core node would give blastn 63 threads on 16 allocated cores. That is
    heavy oversubscription (slower, not faster) and can get the job killed for
    exceeding its allocation. A desktop has no cpuset, so the problem shows
    only under a scheduler.

    Order: process_cpu_count (3.13+, honours affinity) -> sched_getaffinity
    (Linux, honours cpuset/taskset) -> cpu_count (macOS/Windows fallback).
    """
    for probe in ("process_cpu_count",):          # Python 3.13+
        fn = getattr(os, probe, None)
        if fn is not None:
            n = fn()
            if n:
                return n
    getaffinity = getattr(os, "sched_getaffinity", None)   # Linux only
    if getaffinity is not None:
        try:
            n = len(getaffinity(0))
            if n:
                return n
        except OSError:
            pass
    return os.cpu_count() or 2


def _resolve_db_prefixes(
    db_paths: List[str], makeblastdb_path: Optional[str] = None,
) -> List[str]:
    """The BLAST DB prefixes ``db_paths`` refers to, building any FASTAs.

    If the user gave a raw FASTA, a DB is built next to it (cached by SHA-256).
    If they gave a DB prefix whose sidecars are missing, raise with a clear
    message.

    Factored out of ``_run_screen`` so the subject-window retrieval used by
    the full-half recount resolves EXACTLY the same prefixes the search ran
    against. Resolving them twice by different routes can leave the manifest
    describing a different (even 0-byte) database from the one searched (see
    blast.resolve_db_prefix).
    """
    resolved: List[str] = []
    missing: List[str] = []
    for db in db_paths:
        if is_db_prefix(db):
            resolved.append(db)
            continue
        try:
            resolved.append(ensure_db_for_path(
                db, makeblastdb_path=makeblastdb_path,
            ))
        except BlastError:
            missing.append(db)
    if missing:
        raise BlastError(
            "BLAST DB sidecars missing for: " + ", ".join(missing)
            + ". Pass a FASTA (will be auto-indexed) or a real DB prefix, or "
            + "build the DB explicitly with `flexprobe build-db`."
        )
    return resolved


def _make_parallel_factory(
    config: Any, target_ids: Any, spec_kwargs: Dict[str, Any], *,
    needs_subjects: bool, output_dir: str,
    emit: Optional[Callable[[str, float], None]] = None,
) -> Optional[Callable[[Any], Any]]:
    """A factory that wraps a prepared screen in a :class:`ParallelScreen`.

    Returns None — meaning "run the serial path" — whenever the parallel screen
    would be wrong or unavailable. Every rejection SAYS SO rather than silently
    running serially: the whole reason to set ``screen_workers`` is a multi-day
    run, and discovering at the end that it was ignored is the expensive
    failure.
    """
    def _say(msg: str) -> None:
        logger.warning("screen_workers=%d ignored: %s", config.screen_workers,
                       msg)
        if emit:
            emit(f"Parallel screen unavailable ({msg}); running serially", 45)

    if config.screen_workers <= 1:
        return None
    if config.screen.engine != ENGINE_KMER:
        _say(f"only the {ENGINE_KMER!r} engine supports it, this run uses "
             f"{config.screen.engine!r}")
        return None
    if needs_subjects:
        # Defensive: implied by the engine check today, stated because the two
        # are independent facts and the engine check is the one likely to be
        # relaxed first.
        _say("this run retrieves subject windows, whose per-chunk blastdbcmd "
             "prefetch a per-slice split would fragment")
        return None
    if not screen_parallel.available():
        _say("the 'fork' start method is unavailable, and the shared reference "
             "and mmap'd index this depends on are inherited, not pickled")
        return None

    # Never more workers than slices are worth, nor than the node can run. The
    # screen is a random gather over a large reference and so is bound by
    # memory bandwidth well before core count — oversubscribing costs cache,
    # not just scheduling.
    workers = min(config.screen_workers, max(1, _usable_cores()))
    if workers < config.screen_workers:
        logger.info("screen_workers capped %d -> %d by usable cores",
                    config.screen_workers, workers)
    if workers <= 1:
        return None

    def _factory(screen: Any) -> Any:
        return screen_parallel.ParallelScreen(
            screen, workers=workers, output_dir=output_dir,
            target_ids=target_ids, spec_kwargs=spec_kwargs,
        )
    return _factory


def _run_screen(
    passing: List[ProbeCandidate],
    db_paths: List[str],
    output_dir: str,
    *,
    min_pident: float = 0.0,
    min_aln: int = 0,
    query_half_len: int = 25,
    resume: bool = True,
    chunk_label: str = "BLAST chunk",
    emit: Optional[Callable[[str, float], None]] = None,
    cancel: Optional[Callable[[], None]] = None,
    per_chunk: Optional[Callable[
        [int, List[ProbeCandidate], Optional[Dict[str, List[Dict[str, Any]]]]],
        None]] = None,
    spec_params: Optional[Dict[str, Any]] = None,
    parallel_factory: Optional[Callable[[Any], Any]] = None,
    cache_mode: str = screen_cache.MODE_VERDICTS,
    chunk_size: Optional[int] = None,
    max_target_seqs: Optional[int] = None,
    engine: str = ENGINE_BLAST,
    exclude_non_primary: bool = EXCLUDE_ALT_SCAFFOLDS_DEFAULT,
    reference_strand: str = REFERENCE_STRAND_DEFAULT,
    rebuild_kmer_index: bool = False,
) -> Optional[Dict[str, List[Dict[str, Any]]]]:
    """Run the off-target search across one or more DBs; hits keyed by probe_id.

    If ``per_chunk`` is given it is called as ``per_chunk(ci, chunk_candidates,
    chunk_hits)`` for every chunk — whether freshly BLASTed or reloaded from the
    resume cache — and the function returns ``None`` instead of accumulating all
    hits. This lets the caller score specificity chunk-by-chunk so a whole-
    transcriptome screen never holds every probe's hits in RAM at once (which
    can exhaust memory at the end of the screen). Without ``per_chunk``, every
    chunk's hits are merged and the full ``hits_by_probe`` dict is returned.

    ``per_chunk`` is called with ``chunk_hits=None`` for a chunk RESUMED from the
    verdict cache: the candidates already carry their scored fields, so the
    callback must only tally them (re-scoring is impossible — the hits are gone
    by design). ``spec_params`` describes the scoring config and is folded into
    the cache fingerprint whenever ``per_chunk`` is given, so verdicts can never
    be reused across differing thresholds.

    Implementation notes:
      • blastn accepts multiple DBs via a single space-separated ``-db``
        argument. We hand it the user's list directly; missing DBs raise
        BlastError with a helpful path-by-path check.
      • Probe queries are split into chunks of ``BLAST_QUERY_CHUNK`` so a
        whole-transcriptome run reports incremental progress (each chunk
        bumps the progress bar) and peak memory stays bounded.
      • Hits are grouped by qseqid (== probe_id) so the specificity
        engine can score them per probe.
    """
    # Resolved HERE rather than as parameter defaults, which bind at def time:
    # a default of BLAST_QUERY_CHUNK would freeze the value at import and make
    # the module constant unpatchable, silently breaking any caller (or test)
    # that expects to change it by setting the constant.
    # Re-imported from the module (not read from a name bound at import time)
    # so the constant stays patchable and an explicit argument still wins.
    from . import constants as _k
    if chunk_size is None:
        chunk_size = _k.BLAST_QUERY_CHUNK
    if max_target_seqs is None:
        max_target_seqs = _k.BLAST_MAX_TARGET_SEQS


    screen: Screen = make_screen(
        engine,
        min_pident=min_pident, min_aln=min_aln,
        query_half_len=query_half_len,
        num_threads=max(1, _usable_cores() - 1),
        max_target_seqs=max_target_seqs,
        # k-mer engine only, and BlastScreen ignores it: BLAST searches a
        # database built by makeblastdb, so there is no record set for the
        # pipeline to filter at search time. A BLAST run's manifest still
        # RECORDS the reference's composition (see _reference_compositions), it
        # simply cannot act on it.
        exclude_non_primary=exclude_non_primary,
        # Both engines honour it (BlastScreen.__init__ and KmerScreen.__init__
        # both take it); "both" (the default) searches both strands.
        reference_strand=reference_strand,
        # k-mer engine only; BlastScreen accepts and ignores it.
        rebuild_kmer_index=rebuild_kmer_index,
        derive_qcov=_derive_qcov_hsp_perc,
        # DB resolution is INJECTED, not performed here: it builds databases
        # with makeblastdb and its result is what the manifest describes (see
        # _blast_db_sizes / _blast_db_sha256s), so it stays owned by the
        # pipeline — but an engine that needs no BLAST database never calls it,
        # and must not pay for makeblastdb just to start.
        resolve_db=lambda paths: _resolve_db_prefixes(
            paths, blast.discover_blast()["makeblastdb"]),
        # NOT cache_dir: an index is a property of the REFERENCE and belongs
        # beside it, so every run against the same transcriptome shares one
        # (KmerIndex.cache_prefix). This is only where it lands when the
        # reference directory is read-only.
        fallback_cache_dir=os.path.join(output_dir, "_kmer_index"),
    )
    screen.prepare(db_paths)

    # When streaming (per_chunk given) we never accumulate the full hit set:
    # on a whole transcriptome that dict alone can exceed the memory cap.
    streaming = per_chunk is not None
    # Verdict caching only applies to the streaming path: the non-streaming path
    # merges raw hits into one dict and genuinely needs them, so it always
    # caches hits regardless of the knob.
    if cache_mode not in screen_cache.CACHE_MODES:
        raise ValueError(
            f"screen_cache_mode must be one of {screen_cache.CACHE_MODES!r}, "
            f"got {cache_mode!r}"
        )
    cache_verdicts = streaming and cache_mode == screen_cache.MODE_VERDICTS
    if cache_verdicts and resume and spec_params is None:
        # Enforced, not just documented. A verdict cache keyed WITHOUT the
        # scoring config would produce a digest identical to a non-streaming
        # (hits-on-disk) run at the same defaults — two different on-disk formats
        # under one fingerprint. Unreachable from run_pipeline (which always
        # passes spec_params), but a ValueError is cheaper than the debugging.
        raise ValueError(
            "spec_params is required when caching verdicts: they are only valid "
            "for the scoring config that produced them, so it must key the cache"
        )

    # The parallel screen scores in its workers and hands back verdicts, so it
    # is only meaningful where verdicts are what this run deals in. Under
    # `hits` mode the caller wants the raw hits kept (to re-score without
    # re-searching), and those are exactly what the workers throw away — so
    # rather than quietly returning a different thing, the parallel path is
    # simply not taken and the serial one runs.
    parallel = None
    if parallel_factory is not None and cache_verdicts:
        parallel = parallel_factory(screen)
    hits_by_probe: Dict[str, List[Dict[str, Any]]] = (
        {} if streaming else {c.probe_id: [] for c in passing}
    )
    probe_ids = None if streaming else set(hits_by_probe)
    n = len(passing)
    n_chunks = (n + chunk_size - 1) // chunk_size

    # Resumable BLAST: checkpoint each chunk's hits to disk so a job killed by an
    # HPC walltime can be resubmitted and pick up where it left off. The
    # fingerprint guards against reusing a cache built from a different candidate
    # set (e.g. a changed tile_step).
    cdir: Optional[str] = None
    if resume and n_chunks > 0:
        # spec_params is load-bearing for verdicts: they are only valid for the
        # exact scoring config that produced them, so omitting it would let a run
        # with different thresholds silently reuse another run's verdicts. For
        # `hits` we deliberately leave it OUT — raw hits are scoring-agnostic, and
        # that omission is what preserves re-score-without-re-BLAST.
        # Hash the EFFECTIVE mode, not the requested one: the non-streaming path
        # always caches hits even when `verdicts` was asked for, and a
        # fingerprint that misdescribes its own on-disk format is exactly the
        # hazard CACHE_SCHEMA exists to prevent.
        # The engine's own parameters key the cache: an engine that changes
        # which hits come back must invalidate, never silently reuse another
        # engine's answers. `engine=` is in the digest for the same reason.
        fp = screen_cache.fingerprint(
            (c.probe_id for c in passing),
            chunk_size=chunk_size,
            spec_params=(spec_params if cache_verdicts else None),
            cache_mode=(screen_cache.MODE_VERDICTS if cache_verdicts
                        else screen_cache.MODE_HITS),
            engine=screen.name,
            **screen.fingerprint_parts(),
        )
        cdir = screen_cache.open_cache(
            output_dir, fp,
            stamp=engine_identity.current_stamp(
                blastn=(screen.name == ENGINE_BLAST)),
        )

    def _merge(chunk_hits: Dict[str, List[Dict[str, Any]]]) -> None:
        for base, hlist in chunk_hits.items():
            if base in hits_by_probe:
                hits_by_probe[base].extend(hlist)

    for ci in range(n_chunks):
        if cancel:
            cancel()

        chunk = passing[ci * chunk_size:(ci + 1) * chunk_size]

        if cdir is not None:
            resumed = False
            if cache_verdicts:
                # The cache holds this chunk's VERDICTS, not its hits. Restore
                # them onto the candidates and let per_chunk tally.
                # apply_verdicts returns None on a partial/corrupt chunk, which
                # falls through to a fresh BLAST rather than half-applying.
                cached_v = screen_cache.load_verdicts(cdir, ci)
                if cached_v is not None and \
                        screen_cache.apply_verdicts(chunk, cached_v) is not None:
                    # cache_verdicts is defined as `streaming and ...`, and
                    # streaming is `per_chunk is not None` — so this holds by
                    # construction. Stated once per chunk, not per candidate.
                    assert per_chunk is not None
                    per_chunk(ci, chunk, None)   # hits=None => tally only
                    resumed = True
            else:
                cached = screen_cache.load_chunk(cdir, ci)
                if cached is not None:
                    # `hits` mode re-scores from the cached hits — which is the
                    # point of the mode: a threshold change is free here.
                    if per_chunk is not None:    # i.e. streaming
                        per_chunk(ci, chunk, cached)
                    else:
                        _merge(cached)
                    resumed = True
            if resumed:
                if emit and n_chunks > 1:
                    pct = 45 + 5 * (ci + 1) / n_chunks
                    emit(f"{chunk_label} {ci + 1}/{n_chunks} (resumed)", pct)
                continue

        if parallel is not None:
            # Parallel path: worker slices search AND score, so only the
            # ~11 VERDICT_FIELDS per probe come back, never the raw hits
            # (gigabytes per chunk on a whole transcriptome; see screen_parallel for why that choice is what keeps this
            # memory-neutral). apply_verdicts is the SAME call the resume
            # branch above makes, so the parallel and serial paths rejoin at a
            # junction every resumed run already exercises.
            verdicts = parallel.verdicts_for(ci, chunk)
            if screen_cache.apply_verdicts(chunk, verdicts) is None:
                # All-or-nothing by construction: apply_verdicts refuses a
                # partial chunk rather than leaving some probes at their
                # model defaults (specificity_pass=False), which would drop
                # real probes silently. A worker that returned less than its
                # slice is a bug, not a retryable fault.
                raise RuntimeError(
                    f"parallel screen returned verdicts for {len(verdicts)} of "
                    f"{len(chunk)} candidates in chunk {ci}; refusing to apply "
                    f"a partial chunk"
                )
            assert per_chunk is not None   # implied by cache_verdicts
            per_chunk(ci, chunk, None)     # already scored => tally only
            if cdir is not None:
                screen_cache.save_verdicts(cdir, ci, chunk)
            if emit and n_chunks > 1:
                pct = 45 + 5 * (ci + 1) / n_chunks
                emit(f"{chunk_label} {ci + 1}/{n_chunks} done", pct)
            continue

        # Query each 25-bp HALF separately (not the 50-bp tile): off-/on-target
        # hits are judged per 25-mer (e.g.
        # off-target = >=20/25 at >=80%). Query ids are suffixed ``__L`` /
        # ``__R`` and mapped back — see screen.half_queries.
        chunk_hits: Dict[str, List[Dict[str, Any]]] = screen.search_chunk(
            chunk, ci, output_dir,
        )
        # When accumulating (not streaming) the caller owns a global probe-id
        # set; drop anything outside it. The screen already restricts to this
        # chunk's ids, so this only matters for the accumulating path.
        if probe_ids is not None:
            chunk_hits = {k: v for k, v in chunk_hits.items() if k in probe_ids}
        if cache_verdicts:
            # Score FIRST, then checkpoint the verdicts. The `hits` order below
            # (persist, then score) survives a crash mid-scoring, but writes the
            # entire hit dict to disk: gigabytes per chunk, which over a
            # whole-transcriptome run can exceed the disk it runs on. Verdicts
            # are orders of magnitude smaller and are all resume needs.
            # Trade-off: a crash between BLAST and the save loses THIS chunk's
            # BLAST, but scoring takes seconds, so the window is tiny and
            # bounded to one chunk.
            assert per_chunk is not None   # implied by cache_verdicts
            per_chunk(ci, chunk, chunk_hits)   # score now, hold nothing
            if cdir is not None:
                screen_cache.save_verdicts(cdir, ci, chunk)
        else:
            # `hits`: persist BEFORE consuming, so a crash mid-scoring still
            # leaves a reusable checkpoint. Also the non-streaming path, which
            # merges raw hits and genuinely needs them.
            if cdir is not None:
                screen_cache.save_chunk(cdir, ci, chunk_hits)
            if per_chunk is not None:    # i.e. streaming
                per_chunk(ci, chunk, chunk_hits)
            else:
                _merge(chunk_hits)
        # Per-chunk temp files are the engine's own business; search_chunk
        # cleans up whatever it wrote.
        if emit and n_chunks > 1:
            pct = 45 + 5 * (ci + 1) / n_chunks
            emit(f"{chunk_label} {ci + 1}/{n_chunks} done", pct)
    if parallel is not None:
        parallel.close()
    screen.close()
    return None if streaming else hits_by_probe


# ---------------------------------------------------------------------------
# Degenerate consensus probes: the parent / variant family
# ---------------------------------------------------------------------------
#
# A consensus PARENT is a mixed-base string that is never synthesised and never
# orderable. Its VARIANTS are the concrete ACGT oligos — one per distinct
# per-record resolution — and they are what actually goes into the tube. Every
# variant of a site is synthesised into the SAME pool, so the pool behaves like
# its worst member: the parent is judged on its variants' worst case, and a
# family that is not intact (a parent whose variants were not all carried
# forward together) is dropped whole rather than scored on a subset.


def _degenerate_stability(
    parent: ProbeCandidate, variant_seqs: List[str],
) -> None:
    """Populate a parent's Tm/GC ranges from its variants' 50-mer halves.

    Computed at TILE time, before the filter funnel, because the per-half Tm
    gate is the pipeline's duplex-stability check and a parent has no Tm of its
    own — primer3 rejects mixed-base input. Without these numbers the filter has
    nothing to gate a parent on (see ``filters._filter_body``).

    ``*_lhs_*`` / ``*_rhs_*`` are the mRNA window's left and right halves, the
    same convention as ``gc_left`` / ``gc_right`` and ``results.csv``'s
    ``gc_lhs_half``. Note this is NOT the convention of ``tm_target_lhs``, which
    is the synthesised LHS oligo = RC(mRNA RIGHT half); Tm is invariant under
    reverse complement in the nearest-neighbour model, so the two conventions
    give the same pair of numbers with the labels swapped.
    """
    if not variant_seqs:
        return
    half = len(parent.sequence) // 2
    tms_l, tms_r, gcs_l, gcs_r = [], [], [], []
    for seq in variant_seqs:
        left, right = seq[:half], seq[half:]
        tm_l, tm_r = _safe_tm(left), _safe_tm(right)
        if not _is_nan(tm_l):
            tms_l.append(tm_l)
        if not _is_nan(tm_r):
            tms_r.append(tm_r)
        gcs_l.append(_gc_fraction(left))
        gcs_r.append(_gc_fraction(right))
    if tms_l:
        parent.tm_lhs_min, parent.tm_lhs_max = min(tms_l), max(tms_l)
    if tms_r:
        parent.tm_rhs_min, parent.tm_rhs_max = min(tms_r), max(tms_r)
    if gcs_l:
        parent.gc_lhs_min, parent.gc_lhs_max = min(gcs_l), max(gcs_l)
    if gcs_r:
        parent.gc_rhs_min, parent.gc_rhs_max = min(gcs_r), max(gcs_r)


def _reconcile_degenerate_families(
    candidates: List[ProbeCandidate],
) -> List[ProbeCandidate]:
    """Drop any parent/variant family that is no longer intact.

    A family survives only if the parent AND every one of its ``n_variants``
    children are still present. Anything else is dropped whole: a parent scored
    on a subset of its pool would understate the pool's worst case (which is the
    only thing the aggregate means), and an orphaned variant has lost the parent
    that records what it is a resolution OF.

    Called after the filter funnel and again after the pre-screen cap — the two
    places that can remove one family member without the other.

    This is also why ``filter_report.json``'s ``counts.passed`` and
    ``design_summary.json``'s ``stats.filter_passed`` legitimately disagree: the
    former is tallied before either call, the latter after both. See the
    comments at the call sites in ``run_pipeline`` for the exact accounting.
    """
    parents = {
        c.probe_id: c for c in candidates
        if c.probe_role == PROBE_ROLE_CONSENSUS
    }
    if not parents:
        return candidates
    kids: Dict[str, List[ProbeCandidate]] = {}
    for c in candidates:
        if c.probe_role == PROBE_ROLE_VARIANT and c.parent_probe_id:
            kids.setdefault(c.parent_probe_id, []).append(c)
    intact = {
        pid for pid, p in parents.items()
        if len(kids.get(pid, [])) == p.n_variants and p.n_variants > 0
    }
    out = []
    for c in candidates:
        if c.probe_role == PROBE_ROLE_CONSENSUS:
            if c.probe_id in intact:
                out.append(c)
        elif c.probe_role == PROBE_ROLE_VARIANT:
            if c.parent_probe_id in intact:
                out.append(c)
        else:
            out.append(c)
    return out


def _rollup_degenerate_specificity(candidates: List[ProbeCandidate]) -> int:
    """Roll each variant's specificity verdict up onto its consensus parent.

    Everything a child's verdict is made of travels up, not just the boolean:
    the status string (so ``failed_junction_duplex`` stays distinguishable from
    ``failed_off_target`` on the parent row), the margin, and the junction-duplex
    evidence. A parent that reported -1 in the clean-bp columns while its
    variants reported 9+9 would describe an unmeasured screen for a probe that
    was measured, and would hide the very rejection it inherited.

    **All children must pass.** The parent's string is never BLASTed — it is not
    a real sequence and blastn would not take it — so the screen runs on the
    concrete variants, separately. That separation is the point: a variant
    specific to strain B can hit a human transcript where strain A's variant
    does not, and screening a single "representative" would miss it. Since every
    variant is ordered into the same pool, one off-target variant contaminates
    the pool, so the parent inherits the WORST verdict, and its status names
    which variant failed and why rather than just saying "off_target".

    Returns the number of parents that passed.
    """
    kids: Dict[str, List[ProbeCandidate]] = {}
    for c in candidates:
        if c.probe_role == PROBE_ROLE_VARIANT and c.parent_probe_id:
            kids.setdefault(c.parent_probe_id, []).append(c)
    n_passed = 0
    for parent in candidates:
        if parent.probe_role != PROBE_ROLE_CONSENSUS:
            continue
        children = kids.get(parent.probe_id, [])
        if not children:
            parent.specificity_pass = False
            parent.specificity_status = "variant_missing"
            parent.specificity_margin = 0.0
            continue
        failed = [c for c in children if not c.specificity_pass]
        if failed:
            worst = failed[0]
            parent.specificity_pass = False
            parent.specificity_status = (
                f"variant_failed:{worst.probe_id}:{worst.specificity_status}"
            )
        else:
            parent.specificity_pass = True
            # Every child agrees on a status when they all pass; take the first
            # so an unscreened run reports "not_screened", not a bare "pass".
            parent.specificity_status = children[0].specificity_status
            n_passed += 1
        # The pool is only as specific as its least specific member.
        parent.specificity_margin = min(c.specificity_margin for c in children)
        # Same rule for the junction-duplex evidence: the parent reports the
        # WORST child's site, i.e. the largest min(left, right), because that is
        # the one that would ligate on the wrong gene. -1 children (a variant
        # the scan could not measure) are skipped rather than treated as 0; if
        # no child was measured the parent stays -1, which is the truth.
        measured = [
            (c.offtarget_clean_left_bp, c.offtarget_clean_right_bp)
            for c in children
            if c.offtarget_clean_left_bp >= 0 and c.offtarget_clean_right_bp >= 0
        ]
        if measured:
            left, right = max(measured, key=lambda lr: (min(lr), sum(lr)))
            parent.offtarget_clean_left_bp = left
            parent.offtarget_clean_right_bp = right
    return n_passed


class _IngestResult(NamedTuple):
    """The Stage 1 locals that the rest of :func:`run_pipeline` reads.

    ``records`` is deliberately absent: the parsed record list is consumed
    entirely inside the stage (duplicate-id check, gene grouping) and
    nothing downstream reads it.
    """

    gene_seqs: Dict[str, str]
    target_ids: Set[str]
    gene_aware: bool
    gene_of: Optional[Dict[str, str]]
    group_members: Optional[Dict[str, Set[str]]]


def _stage_ingest(
    config: RunConfig,
    target_ids: Optional[Set[str]],
    _emit: Callable[..., None],
    _check_cancel: Callable[[], None],
) -> _IngestResult:
    """Stage 1: parse the input FASTA and set up gene-aware (isoform) groups.

    Called from :func:`run_pipeline`. ``_emit`` and ``_check_cancel`` are that
    function's own closures, passed in so the progress events, and the stage
    timings they drive (:class:`_StageTimer`), fire from one place.
    """
    _emit("ingest", "Parsing FASTA", 5)
    _check_cancel()
    records = parse_fasta(
        config.fasta_path, split_header_on_pipe=config.split_header_on_pipe,
    )
    gene_seqs = {gid: seq for gid, _desc, seq in records}
    # Belt-and-braces: parse_fasta already rejects duplicate ids, so this dict
    # cannot lose a record. Asserting it here keeps the "N records in => N genes
    # designed against" invariant explicit: a parser change that merges ids
    # fails loudly instead of silently designing against fewer targets than
    # the user supplied.
    if len(gene_seqs) != len(records):
        raise PipelineError(
            f"{len(records)} FASTA records collapsed to {len(gene_seqs)} unique "
            f"ids — record ids must be unique."
        )
    _emit("ingest", f"Loaded {len(records)} genes", 10, {"gene_count": len(records)})

    # Is this the sense strand? Nothing downstream can tell: every filter, the
    # chemistry and the off-target screen are all computed against whichever
    # strand was supplied, so a genome/vRNA-orientation input yields a complete,
    # well-scored, orderable panel that detects nothing. Warn only — the check
    # is a heuristic and some legitimate targets are non-coding. See
    # flexprobe.orientation.
    try:
        orientation.check_orientation(
            [(gid, seq) for gid, _desc, seq in records]
        )
    except Exception:                      # pragma: no cover - advisory only
        # A sanity check must never be the thing that kills a run.
        logger.debug("Input-orientation check failed", exc_info=True)

    if target_ids is None:
        target_ids = set(gene_seqs.keys())

    # Gene-aware (isoform) mode: group transcripts by gene so sibling isoforms
    # aren't off-targets and each probe is scored on isoform coverage. Ignored
    # in multi-strain mode, which has its own consensus logic.
    gene_aware = config.group_by_gene and not config.multi_strain
    if gene_aware and not any(_GENE_FIELD_RE.search(desc or "")
                              for _gid, desc, _seq in records):
        # Grouping is on by default. Without a single gene: label every
        # record would be its own gene, and a probe that also hits ANOTHER
        # input record would turn from on-target into off-target: a silent
        # change for viral segments or a custom multi-gene panel. So it is
        # switched off here, visibly, and the run is exactly an ungrouped one.
        gene_aware = False
        _msg = ("Gene grouping is on, but no input record carries an "
                "Ensembl-style gene: header field, so there is nothing to "
                "group; running without gene grouping (unchanged behaviour).")
        logger.info(_msg)
        _emit("ingest", _msg, 11)
    gene_of: Optional[Dict[str, str]] = None
    group_members: Optional[Dict[str, Set[str]]] = None
    if gene_aware:
        gene_of, group_members = _build_gene_groups(records)
        _emit("ingest",
              f"Gene-aware: {len(group_members)} gene(s) across "
              f"{len(records)} transcripts", 11,
              {"gene_count": len(group_members)})
    return _IngestResult(
        gene_seqs, target_ids, gene_aware, gene_of, group_members,
    )


class _MultiStrainTiling(NamedTuple):
    """What the multi-strain arm of Stage 2 hands back to :func:`run_pipeline`.

    ``candidate_iter`` is a LAZY generator that is not consumed until the
    filter stage, hundreds of lines later. It closes over the window lists,
    the alignment and the per-window sequence overrides built here, so those
    stay alive inside the generator's closure. Nothing here may force it.
    """

    alignment: StrainAlignment
    gene_seqs: Dict[str, str]
    target_ids: Set[str]
    on_target_records: Optional[FrozenSet[str]]
    consensus_strain_rule: StrainRule
    candidate_iter: Iterator[ProbeCandidate]


def _consensus_provenance(
    alignment: Optional["StrainAlignment"],
) -> Dict[str, Any]:
    """Durable record of WHICH ALIGNER built the consensus, for the manifest.

    ``_emit`` forwards only ``(stage, msg, pct)`` to the timeline writer, whose
    record is a fixed ten keys, so ``ConsensusResult.aligner`` passed as
    ``_emit``'s ``data`` reaches only a log line and the GUI callback. This
    function puts it in the manifest instead. The bundled MAFFT and the
    Biopython pairwise fallback do not produce the same alignment, so a
    finished run has to be able to say which one ran.

    Returns an EMPTY dict for a single-species run, where no consensus is built
    and there is therefore no aligner to name. ``params['multi_strain']`` is
    False in that case and distinguishes "nothing to record" from "failed to
    record", so an absent key here must not be read as the latter.

    Lives at module level, called from both manifest sites, so the full-run and
    prepare-only halves of a split run cannot record different things.
    """
    if alignment is None:
        return {}
    rec: Dict[str, Any] = {
        "aligner": alignment.aligner,
        "used_pairwise": alignment.used_pairwise,
        "n_strains": len(alignment.strain_ids),
        "strain_ids": list(alignment.strain_ids),
        "variant_positions": len(alignment.variant_positions),
        "fragment_name": alignment.fragment_name,
    }
    # Present ONLY when an aligner that was found failed and the built-in
    # fallback ran instead. A logger.warning alone is easy to miss (a bundled
    # MAFFT that fails on Windows would fall back unnoticed), so the manifest
    # records it. Absent rather than empty otherwise, so a run where nothing
    # failed carries no extra key.
    if alignment.aligner_fallback:
        rec["aligner_fallback"] = alignment.aligner_fallback
    return rec


def _stage_multi_strain_tiling(
    config: RunConfig,
    gene_seqs: Dict[str, str],
    output_dir: str,
    _emit: Callable[..., None],
) -> _MultiStrainTiling:
    """Stage 2, multi-strain arm: build the consensus and the candidate stream.

    Called from :func:`run_pipeline`'s ``if config.multi_strain:`` branch.
    Its two nested helpers (``_annotate_strain_binding`` and
    the ``_multi_strain_iter`` generator) stay nested here because they
    close over this function's locals and are meaningless without them.

    ``_emit`` is run_pipeline's own closure, passed in rather than
    reimplemented: it drives :class:`_StageTimer`, so the stage timings stay
    consistent with the rest of the run.

    ``gene_seqs`` goes in as the parsed input records and comes back REPLACED
    by the single virtual consensus gene -- that substitution is what makes
    the downstream screen / specificity / chemistry stages see one gene, and
    it is why the value is returned rather than mutated in place.
    """
    # Treat all FASTA records as strains of one fragment. Build the
    # consensus once, then either restrict tiling to conserved windows
    # (probes work for every strain) or to strain-unique windows (probes
    # distinguish one strain from the rest).
    # Fragment name = the input FASTA's basename (which in this
    # workflow is the fragment label, e.g. "OROV_L").
    fragment_name = (
        os.path.splitext(os.path.basename(config.fasta_path))[0]
        or "consensus"
    )
    _records = [(gid, seq) for gid, seq in gene_seqs.items()]
    # Rotate so config.consensus_reference leads. build_consensus treats the
    # first record as both the alignment coordinate system and the strain
    # that ambiguous consensus bases resolve to, so rotating is exactly
    # equivalent to reordering the input FASTA — and it is what makes the
    # OTHER records' flavour of each pan-strain window reachable at all.
    # Skipped for strain_unique, where strain_unique_target indexes these
    # same records and rotating underneath it would silently retarget the run.
    if (config.consensus_strategy != "strain_unique"
            and 0 < config.consensus_reference < len(_records)):
        _r = config.consensus_reference
        _records = _records[_r:] + _records[:_r]
        _emit("tile",
              f"Consensus template: {_records[0][0]} "
              f"(record {_r} of {len(_records)})", 15)
    alignment = build_consensus(
        _records,
        fragment_name=fragment_name,
    )
    _emit("tile",
          f"Built consensus from {alignment.n_strains} strain(s); "
          f"{len(alignment.variant_positions)} variant positions",
          18,
          {"variant_positions": len(alignment.variant_positions),
           "used_pairwise": alignment.used_pairwise,
           "aligner": alignment.aligner})
    # Replace gene_seqs with the consensus so the downstream BLAST /
    # specificity / chemistry stages all reference one virtual gene.
    consensus_id = alignment.fragment_name or "consensus"
    gene_seqs = {consensus_id: alignment.consensus}
    # BLAST on-target classification matches the DB record names, not
    # the synthetic consensus id. The strains are what the user put in
    # the DB, so target on-target hits should be judged against strain ids.
    target_ids = set(alignment.strain_ids)
    # Same ids, a different question: which of them must a probe BIND?
    # Every one of them under "conserved", only the target under
    # "strain_unique" (where mismatching the others is the design goal),
    # and see filters.on_target_record_ids for "constitutive". This is what
    # scopes the junction-mismatch rule; getting it wrong in the
    # strain-unique direction would reject the run's whole point, so it is
    # derived in exactly one place rather than re-decided at the filter.
    from .filters import on_target_record_ids as _on_target_record_ids
    on_target_records = _on_target_record_ids(config, alignment.strain_ids)
    # Persist the consensus so the user can inspect / re-align it.
    try:
        cons_path = os.path.join(output_dir, "consensus.fa")
        with open(cons_path, "w") as _cf:
            _cf.write(f">{consensus_id}_consensus\n{alignment.consensus}\n")
            for sid, aseq in zip(alignment.strain_ids, alignment.aligned_seqs):
                _cf.write(f">{sid}_aligned\n{aseq}\n")
    except OSError:
        logger.warning("Could not write consensus.fa", exc_info=True)

    # Per-window probe sequence override (used by constitutive/isoform mode
    # where each window's probe is the shared sequence of its isoform group,
    # not a slice of a single reference sequence). Empty otherwise.
    window_seq_by_start: Dict[int, str] = {}
    # Per-window covered RECORD IDS, the companion to the sequence override:
    # which isoforms the window's shared sequence actually belongs to. Only
    # constitutive mode populates it, and it becomes
    # ``ProbeCandidate.covered_record_ids`` — the per-candidate scope for the
    # junction-mismatch rule. Without it a partial-quorum run has no scope at
    # all (filters.on_target_record_ids returns None for that case) and the
    # rule is skipped on the one mode where window selection does NOT
    # already keep mismatches out of the band.
    covered_by_start: Dict[int, Tuple[str, ...]] = {}
    # Windows that carry IUPAC ambiguity, keyed by start. Only the
    # "conserved" strategy at a non-zero tolerance ever populates this;
    # empty everywhere else, which keeps the other two strategies (and every
    # tolerance-0 run) on the plain one-candidate-per-window path.
    degenerate_by_start: Dict[int, Any] = {}
    # Specificity rule for this consensus run. Constitutive (isoform) probes
    # need only hit at least one of their isoforms (pan), not every record.
    consensus_strain_rule = StrainRule.EXACT_SET
    if config.consensus_strategy == "strain_unique":
        tgt = config.strain_unique_target
        if not (0 <= tgt < alignment.n_strains):
            raise PipelineError(
                f"strain_unique_target index {tgt} out of range "
                f"(found {alignment.n_strains} strains)"
            )
        starts = strain_unique_start_positions(
            alignment,
            target_strain_idx=tgt,
            tile_length=config.tile_length,
            tile_step=config.tile_step,
            min_mismatches_vs_others=config.strain_unique_min_mismatches,
            per_half=config.strain_unique_per_half,
            junction_duplex_bp=config.strain_unique_junction_duplex_bp,
        )
        start_list = [s for s, _mm in starts]
        _unit = (
            f"a single {config.tile_length // 2} bp half"
            if config.strain_unique_per_half
            else f"the whole {config.tile_length} bp window"
        )
        _duplex_note = (
            f", none sealing a >= {config.strain_unique_junction_duplex_bp} bp "
            f"junction duplex against a non-target strain"
            if (config.strain_unique_per_half
                and config.strain_unique_junction_duplex_bp > 0)
            else ""
        )
        _emit("tile",
              f"Strain-unique mode: {len(start_list)} windows differ "
              f"from other strains by >= {config.strain_unique_min_mismatches} bp "
              f"within {_unit}{_duplex_note}",
              20)
        # Strain-unique probes use the target strain's actual sequence
        # (not the IUPAC consensus) so the probe matches that strain
        # perfectly. The consensus's variant positions would be
        # ambiguity codes — orderable only as a mixed-base oligo and
        # not what the user wants for diagnostic specificity.
        seq = alignment.aligned_seqs[tgt]
    elif config.consensus_strategy == "constitutive":
        # Isoform-collapse: windows shared by a quorum of isoforms. Each
        # probe is designed from the shared sequence of the largest isoform
        # group at that window.
        min_iso = config.isoform_min_count or None
        triples = constitutive_start_positions(
            alignment,
            tile_length=config.tile_length,
            tile_step=config.tile_step,
            min_isoforms=min_iso,
        )
        start_list = [s for s, _w, _idx in triples]
        window_seq_by_start = {s: w for s, w, _idx in triples}
        # Map the covered isoform INDICES onto record ids while the
        # alignment is in scope. At the default quorum every window's group
        # is every record, so this is the same set on_target_record_ids
        # already returns and the gate is unchanged; below the quorum it is
        # the only place the group is known.
        _sids = alignment.strain_ids
        covered_by_start = {
            s: tuple(_sids[i] for i in idx if 0 <= i < len(_sids))
            for s, _w, idx in triples
        }
        consensus_strain_rule = StrainRule.PAN_STRAIN
        _need = config.isoform_min_count or alignment.n_strains
        _emit("tile",
              f"Constitutive (isoform) mode: {len(start_list)} windows "
              f"shared by >= {_need} of {alignment.n_strains} isoforms",
              20)
        seq = alignment.consensus
    else:
        # "conserved": windows every strain can be probed from. At the
        # default tolerance of 0 that means every strain agrees exactly;
        # above 0 it admits a bounded number of variant positions per half
        # with the ligation junction kept clear, and those windows fan out
        # into a degenerate parent plus one concrete variant per record.
        _degen_windows = degenerate_conserved_windows(
            alignment,
            tile_length=config.tile_length,
            tile_step=config.tile_step,
            max_degenerate_per_half=config.conserved_max_degenerate_per_half,
            junction_clearance=config.filters.degenerate_junction_clearance,
        )
        start_list = [w.start for w in _degen_windows]
        degenerate_by_start = {
            w.start: w for w in _degen_windows if w.degenerate_positions
        }
        if config.conserved_max_degenerate_per_half > 0:
            _emit("tile",
                  f"Conserved mode (tolerance "
                  f"{config.conserved_max_degenerate_per_half}/half): "
                  f"{len(start_list)} windows across all "
                  f"{alignment.n_strains} strains, {len(degenerate_by_start)} "
                  f"of them degenerate", 20,
                  {"degenerate_windows": len(degenerate_by_start)})
        else:
            _emit("tile",
                  f"Conserved mode: {len(start_list)} windows identical "
                  f"across all {alignment.n_strains} strains",
                  20)
        # Exact-identity windows are the same in every strain by
        # construction, so consensus == every strain's window; degenerate
        # windows carry the IUPAC consensus, which is the parent string.
        seq = alignment.consensus
    # First input sequence = the default/reference strain. Any ambiguous
    # (N / IUPAC) base in a designed probe window defaults to this strain
    # so the exported oligos are concrete, orderable ACGT with no N's.
    reference_seq = alignment.aligned_seqs[0] if alignment.aligned_seqs else ""
    assert alignment is not None      # assigned at the top of this branch
    _aln = alignment

    def _annotate_strain_binding(pc: ProbeCandidate, start: int) -> None:
        """Attach per-strain mismatch/ΔTm columns for ONE candidate.

        Scored against the candidate's own sequence, so a variant reports
        how IT binds each strain (perfect on the strain it resolves from,
        mismatched on the others) while the parent — whose IUPAC codes read
        as "any of these bases" — reports the pan-strain zero that is the
        whole claim being made for it.
        """
        scores = score_probe_against_all_strains(
            _aln, start, tile_length=config.tile_length,
            probe_seq=pc.sequence,
            junction_window=config.filters.junction_window,
        )
        pc.strain_binding = [
            {
                "strain_id": s.strain_id,
                "total_mismatches": s.total_mismatches,
                "lhs_mismatches": s.lhs_mismatches,
                "rhs_mismatches": s.rhs_mismatches,
                "junction_mismatches": s.junction_mismatches,
                "mismatch_positions": s.mismatch_positions,
                "estimated_dTm": s.estimated_dTm,
                "likely_binds": s.likely_binds,
            }
            for s in scores
        ]
        if scores:   # guard: empty reduction would raise ValueError
            pc.min_strain_mismatches = min(s.total_mismatches for s in scores)
            pc.max_strain_mismatches = max(s.total_mismatches for s in scores)
            pc.worst_strain_dTm = min(s.estimated_dTm for s in scores)

    def _multi_strain_iter():
        for start in start_list:
            end = start + config.tile_length
            dw = degenerate_by_start.get(start)
            if dw is None:
                # Plain path: one candidate, no parent
                # and no children. Probe sequence is the per-window shared
                # sequence in constitutive mode, else a slice of the
                # reference/consensus; resolve_to_reference is a no-op on a
                # window that is already pure ACGT.
                probe_seq = window_seq_by_start.get(start) or seq[start:end]
                probe_seq = resolve_to_reference(
                    probe_seq, reference_seq[start:end])
                pc = ProbeCandidate(
                    gene_id=consensus_id,
                    probe_id=f"{consensus_id}_t{start:06d}",
                    sequence=probe_seq,
                    start=start,
                    end=end,
                    # Empty outside constitutive mode, where it leaves the
                    # junction-mismatch scope exactly as the run-level one.
                    covered_record_ids=list(
                        covered_by_start.get(start, ())),
                )
                _annotate_strain_binding(pc, start)
                yield pc
                continue

            # Degenerate window: one never-orderable consensus parent plus
            # one concrete variant per distinct per-record resolution.
            parent_id = f"{consensus_id}_t{start:06d}"
            resolutions = resolve_window_per_record(
                dw.sequence, _aln.aligned_seqs, _aln.strain_ids, start,
            )
            parent = ProbeCandidate(
                gene_id=consensus_id,
                probe_id=parent_id,
                sequence=dw.sequence,
                start=start,
                end=end,
                probe_role=PROBE_ROLE_CONSENSUS,
                degenerate_positions=list(dw.degenerate_positions),
                iupac_codes=list(dw.iupac_codes),
                n_variants=len(resolutions),
            )
            # Populate the parent's Tm/GC ranges BEFORE it reaches the
            # filter funnel — they are the only thing the funnel can gate a
            # mixed-base string on (filters._filter_body).
            _degenerate_stability(parent, [rs for rs, _ids in resolutions])
            _annotate_strain_binding(parent, start)
            yield parent
            for vi, (resolved, sids) in enumerate(resolutions):
                # Ids are derived from the parent's and from the resolution
                # ORDER, which follows the input record order — so they are
                # stable across runs. The pipeline is resumable and cached;
                # an id that moved between runs would silently invalidate
                # every checkpoint keyed on it.
                child = ProbeCandidate(
                    gene_id=consensus_id,
                    probe_id=f"{parent_id}_v{vi:02d}",
                    sequence=resolved,
                    start=start,
                    end=end,
                    probe_role=PROBE_ROLE_VARIANT,
                    parent_probe_id=parent_id,
                    resolves_to=list(sids),
                    degenerate_positions=list(dw.degenerate_positions),
                    iupac_codes=list(dw.iupac_codes),
                    n_variants=len(resolutions),
                )
                _annotate_strain_binding(child, start)
                yield child
    candidate_iter = _multi_strain_iter()
    return _MultiStrainTiling(
        alignment=alignment,
        gene_seqs=gene_seqs,
        target_ids=target_ids,
        on_target_records=on_target_records,
        consensus_strain_rule=consensus_strain_rule,
        candidate_iter=candidate_iter,
    )


def _stream_filter(
    candidate_iter,
    thresholds,
    *,
    reject_cap: int = 25,
    cancel: Optional[Callable[[], None]] = None,
    sink: Optional[Callable[[ProbeCandidate], None]] = None,
    on_target_records: Optional[FrozenSet[str]] = None,
):
    """Stream candidates through the filter funnel without holding them all.

    Returns ``(passing, counters, rejects)`` where ``passing`` is the list of
    candidates that cleared every filter, ``counters`` is the per-reason
    tally (incl. 'total'/'passed'), and ``rejects`` is a *capped* sample of
    rejected candidates — up to ``reject_cap`` per failure reason — for the
    filter report. Peak memory is O(passing + reject_cap·#reasons), so a
    whole-transcriptome run never materializes every tile.

    If ``sink`` is given, each passing candidate is handed to it and NOT kept in
    ``passing`` (which is returned empty). The pre-screen cap uses this to feed a
    bounded per-gene accumulator, so even ``passing`` itself is never fully
    materialized on a whole-transcriptome run.

    ``on_target_records`` is the multi-strain junction-mismatch scope (the set
    of input records every candidate must bind), forwarded verbatim to the
    filter. ``None`` — the default, and what every single-sequence run passes —
    leaves that rule inert; see ``filters.on_target_record_ids``.
    """
    from .filters import evaluate_candidate, new_filter_counters

    counters = new_filter_counters()
    passing: List[ProbeCandidate] = []
    rejects: List[ProbeCandidate] = []
    per_reason_kept: Dict[str, int] = {}
    for i, pc in enumerate(candidate_iter):
        if cancel and (i & 0x3FFF) == 0:   # check cancel every ~16k tiles
            cancel()
        if evaluate_candidate(pc, thresholds, counters,
                              on_target_records=on_target_records):
            if sink is not None:
                sink(pc)
            else:
                passing.append(pc)
        else:
            keep = False
            for r in pc.filter_fail_reasons:
                if per_reason_kept.get(r, 0) < reject_cap:
                    per_reason_kept[r] = per_reason_kept.get(r, 0) + 1
                    keep = True
            if keep:
                rejects.append(pc)
    return passing, counters, rejects


def _tile_filter_one_gene(args):
    """Tile + filter a single gene. Top-level (picklable) so it can run in a
    worker process for whole-transcriptome parallelism. Returns
    ``(gene_id, passing, counters, rejects)`` where ``rejects`` is already
    capped to ``reject_cap`` per failure reason.

    No junction-mismatch scope is threaded here, and none is needed: the
    caller only takes this path when ``not config.multi_strain`` (see
    ``_use_parallel``), and a single-sequence candidate carries no
    ``strain_binding`` for that rule to read. If the parallel path is ever
    opened to multi-strain runs, the scope has to come with it.
    """
    gid, gseq, tile_length, tile_step, thresholds, reject_cap = args

    def _iter():
        for start, end, subseq in tile_sequence(gseq, length=tile_length, step=tile_step):
            yield ProbeCandidate(
                gene_id=gid, probe_id=f"{gid}_t{start:06d}",
                sequence=subseq, start=start, end=end,
            )

    passing, counters, rejects = _stream_filter(
        _iter(), thresholds, reject_cap=reject_cap,
    )
    return gid, passing, counters, rejects


def _merge_filter_results(results, reject_cap: int):
    """Merge per-gene ``_tile_filter_one_gene`` outputs in a FIXED order:
    concatenate passing (preserving the given gene order), sum counters, and
    re-cap the reject sample per reason so the merged report matches the serial
    path's bounds."""
    from .filters import new_filter_counters

    passing: List[ProbeCandidate] = []
    counters = new_filter_counters()
    rejects: List[ProbeCandidate] = []
    per_reason_kept: Dict[str, int] = {}
    for _gid, gene_passing, gene_counters, gene_rejects in results:
        passing.extend(gene_passing)
        for k, v in gene_counters.items():
            counters[k] = counters.get(k, 0) + v
        for pc in gene_rejects:
            keep = False
            for r in pc.filter_fail_reasons:
                if per_reason_kept.get(r, 0) < reject_cap:
                    per_reason_kept[r] = per_reason_kept.get(r, 0) + 1
                    keep = True
            if keep:
                rejects.append(pc)
    return passing, counters, rejects


def _candidates_per_gene(pairs: Sequence[Any],
                         gene_seqs: Dict[str, Any]) -> Dict[str, int]:
    """Pairs per gene — the count the review dialog labels its column with.

    ONE pass over ``pairs``. The obvious comprehension,

        {gid: sum(1 for p in pairs if p.lhs.gene_id == gid) for gid in gene_seqs}

    rescans every pair once per gene: O(genes x pairs). On a whole
    transcriptome (tens of thousands of genes, millions of pairs) that costs
    hours, spent building a diagnostic dict. It does not show in the stage
    timings, because ``stage_timer`` has already snapshotted into the manifest
    by the time this runs.

    Every gene in ``gene_seqs`` stays a key, including genes with no surviving
    pair: consumers index by gene and a missing key is a KeyError, not a zero.
    """
    counts = Counter(p.lhs.gene_id for p in pairs)
    return {gid: counts.get(gid, 0) for gid in gene_seqs}


def _write_filter_report(
    output_dir: str,
    rejected_candidates: List[ProbeCandidate],
    filter_stats: Dict[str, int],
    config: RunConfig,
    *,
    max_examples: int = 25,
) -> None:
    """Write ``filter_report.json`` summarizing why tiles were rejected.

    ``rejected_candidates`` is the capped reject sample from
    ``_stream_filter`` (already only failing candidates). Includes the
    headline counts plus up to ``max_examples`` examples per failure reason.

    ``filter_stats["passed"]`` here counts funnel survivors row-by-row
    (parent and variant independently) and is captured BEFORE
    ``_reconcile_degenerate_families`` runs. ``design_summary.json``'s
    ``stats.filter_passed`` is measured after it (twice, in fact — see the
    call sites in ``run_pipeline``), so the two numbers can legitimately
    differ by the count of families reconciliation dropped. See the comments
    at both call sites for the full explanation.
    """
    import json as _json

    examples: Dict[str, List[Dict[str, Any]]] = {}
    for c in rejected_candidates:
        if c.passed_filters:
            continue
        for reason in c.filter_fail_reasons:
            bucket = examples.setdefault(reason, [])
            if len(bucket) < max_examples:
                bucket.append({
                    "probe_id": c.probe_id,
                    "gene_id": c.gene_id,
                    "start": c.start,
                    "sequence": c.sequence,
                    "gc_content": round(c.gc_content, 4),
                    "gc_left": round(c.gc_left, 4),
                    "gc_right": round(c.gc_right, 4),
                    "homopolymer_max": c.homopolymer_max,
                    "reasons": list(c.filter_fail_reasons),
                })
    import dataclasses as _dc

    report = {
        # asdict, not a hand-listed subset: a hand-listed subset can omit a
        # threshold (such as reject_junction_base_n), and the report would then
        # understate the funnel that produced its own counts. A threshold added
        # to FilterThresholds later appears here automatically.
        "thresholds": _dc.asdict(config.filters),
        "counts": filter_stats,
        "examples": examples,
    }
    tmp = os.path.join(output_dir, "filter_report.json.tmp")
    final = os.path.join(output_dir, "filter_report.json")
    with open(tmp, "w") as fh:
        _json.dump(report, fh, indent=2)
    os.replace(tmp, final)


def _run_null_model(
    config: RunConfig,
    passing: List[ProbeCandidate],
    output_dir: str,
    *,
    spec_kwargs: Dict[str, Any],
    target_ids: Set[str],
    half_len: int,
    emit,
    cancel,
    subjects: Any = None,
) -> Optional[Dict[str, Any]]:
    """Screen shuffled copies of this run's windows and report the FP rate.

    Uses the run's own ``spec_kwargs`` / ``target_ids`` / databases, so the null
    population goes through the identical gate — a calibration against a
    slightly different screen would measure nothing. Writes
    ``null_model.json`` and returns its payload.

    The shuffled candidates are screened WITHOUT resume: they are synthetic and
    differ every seed, so a candidate-fingerprinted cache would either miss on
    every chunk or, worse, collide with the real run's cache.
    """
    import random as _random

    from . import null_model as _nm
    from .specificity import batch_evaluate as _batch_evaluate

    if not passing:
        return None
    rng = _random.Random(config.null_model_seed)
    sample = _nm.sample_windows(
        passing, max_windows=config.null_model_max_windows, rng=rng,
    )
    nulls, degraded = _nm.build_null_candidates(
        sample, shuffles=config.null_model_shuffles, rng=rng,
    )
    if not nulls:
        return None
    emit("specificity",
         f"Null model: screening {len(nulls)} shuffled windows "
         f"({len(sample)} sampled x {config.null_model_shuffles})", 67)

    null_hits: Dict[str, List[Dict[str, Any]]] = {}

    def _collect(_ci, _chunk, chunk_hits):
        if chunk_hits:
            null_hits.update(chunk_hits)

    off_pident_pd, off_aln_pd = _off_target_pushdown(config.specificity, half_len)
    # Written to a scratch subdirectory so the null run's temporary BLAST query
    # and TSV files never land beside the design outputs.
    null_dir = os.path.join(output_dir, "_null_model")
    os.makedirs(null_dir, exist_ok=True)
    _run_screen(
        nulls, config.screen.db_paths, null_dir,
        chunk_size=config.screen.query_chunk,
        max_target_seqs=config.screen.max_target_seqs,
        min_pident=min(off_pident_pd, config.specificity.on_target_pident_min),
        min_aln=min(off_aln_pd, config.specificity.on_target_aln_min),
        query_half_len=half_len,
        resume=False,
        emit=None,
        cancel=cancel,
        per_chunk=_collect,
        # The null calibration must go through the SAME engine and see the SAME
        # reference records the run did. A rate measured by a different engine,
        # or against a different record set, is not comparable to the run's --
        # which is the only thing it is for.
        #
        # Omitting `engine` here would leave this on the BLAST default while
        # the run used k-mer, and that is not a cosmetic discrepancy. The BLAST
        # path recounts off-target mismatches by reading subject bases back
        # through `blastdbcmd -entry`, so against a k-mer run's plain-FASTA
        # references (which carry no accession index) every chance core would
        # score as a 0-mismatch off-target and the null would report a
        # fabricated 100% false-positive rate. The CLI refuses to start a BLAST
        # run in that state; for this path the guard is the measurability
        # check after the evaluation below.
        engine=config.screen.engine,
        exclude_non_primary=config.screen.exclude_alt_scaffolds,
        # Same reasoning as `engine` just above: the null calibration must
        # count off-targets by the SAME rule the run did, and which strand(s)
        # count is part of that rule.
        reference_strand=config.screen.reference_strand,
    )
    # Same `subjects` handle as the real run: the calibration compares the null's
    # rejection rate against the run's, and counting the two by different rules
    # (exact vs conservative) would make the comparison meaningless.
    _batch_evaluate(nulls, null_hits, target_ids, subjects=subjects,
                    **spec_kwargs)

    payload = _nm.summarise(
        nulls,
        n_real_windows=len(passing),
        n_sampled=len(sample),
        shuffles=config.null_model_shuffles,
        degraded=degraded,
        seed=config.null_model_seed,
    )
    # The run's own rejection rate, alongside the null's, so the comparison the
    # calibration exists to support is in one place rather than assembled by the
    # reader from two files.
    real_rejected = sum(1 for c in passing if not c.specificity_pass)
    payload["real_windows_screened"] = len(passing)
    payload["real_rejected"] = real_rejected
    payload["real_rejection_rate"] = round(
        real_rejected / len(passing), 6) if passing else 0.0
    try:
        shutil.rmtree(null_dir, ignore_errors=True)
    except OSError:
        pass
    path = _nm.write_report(output_dir, payload)
    # An unmeasurable calibration carries no `false_positive_rate` key, so the
    # line that reports one cannot be formatted. Say what happened instead of
    # printing a rate the payload deliberately does not contain -- and WARN,
    # because a null model that ran but measured nothing is a silent no-op
    # otherwise, and the user asked for the number.
    if not payload.get("measurable", True):
        logger.warning(
            "Null model: NOT MEASURED. %s Report: %s",
            payload.get("unmeasurable_reason", ""), path,
        )
        emit("specificity",
             f"Null model: not measured ({payload['null_windows_screened']} "
             f"windows could not be recounted exactly) — {path}", 68,
             {"null_model": payload})
        return payload
    emit("specificity",
         f"Null model: gate rejects {payload['false_positive_rate']:.1%} of "
         f"composition-matched noise (vs {payload['real_rejection_rate']:.1%} "
         f"of real windows) — {path}", 68,
         {"null_model": payload})
    return payload


def _off_target_pushdown(spec, half_len: int) -> Tuple[float, int]:
    """Server-side (identity, length) floors for the off-target screen.

    BLAST is given the loosest thresholds that still cannot drop a hit the
    Python classifier would have counted — the pushdown exists to cut hit
    volume, never to make a decision.

    In full-half mode the Python rule is "at least ``off_target_min_mismatch``
    mismatches across the whole half"; a hit disqualifies at
    ``m <= min_mismatch - 1`` mismatches (4 at the defaults). The floors are
    therefore the WEAKEST alignment BLAST could report for such a hit.

    THE LENGTH FLOOR IS NOT ``keep``.
    --------------------------------
    ``keep = half_len - min_mismatch + 1`` (21/25) assumes that a disqualifying
    hit aligns over at least that many query positions. It does not: BLAST
    reports the MAXIMUM-SCORING local segment, not the full-length alignment.
    Under blastn-short's default reward +1 / penalty -3, a 25-mer sitting 4
    mismatches from a subject scores 21 - 12 = 9 over its full length, so any
    exact core of 9 nt or more outscores it and is reported instead. A 21 nt
    floor would throw those cores away server-side, before
    ``specificity.full_half_mismatches`` ever saw them, and make the recount
    inert. run_blastn passes -penalty -1, so the same hit scores 21 - 4 = 17
    and the full-length alignment wins. The bound below is written in terms of
    the constant so it stays correct if the scoring changes.

    The correct bound comes from that same scoring arithmetic. For a hit with
    ``m`` mismatches over the half, the maximum-scoring segment scores at least
    the full-length alignment, ``half_len - (reward + |penalty|) * m``, and its
    LENGTH is at least its score (every column contributes at most +1). So over
    all disqualifying hits (``m <= min_mismatch - 1``)::

        length >= half_len - COST * (min_mismatch - 1)

    where COST is ``constants.BLASTN_SHORT_MISMATCH_COST``, being reward plus
    the absolute penalty. At -penalty -1 that is 2, giving 25 - 2*4 = 17. At
    blastn-short's default -penalty -3 it would be 4, giving 25 - 4*4 = 9. A
    cost that is too small only lowers the floor and keeps hits that need not
    be kept, so the safe direction is down.

    This is the floor used here. Gaps cannot break it: a gap column adds length
    and subtracts score, so it only widens the margin. The search is
    -ungapped in any case.

    The IDENTITY floor is ``keep``-based and safe. The worst case is a hit
    reported at full length with ``min_mismatch - 1`` mismatches, i.e.
    ``keep / half_len`` = 84%; dividing by ``half_len + GAP_MARGIN`` instead
    keeps it below that with slack for gapped HSPs.

    NOTE: it is tempting to push ``-qcov_hsp_perc 100`` here, since the rule is
    "over the full 25 nt". That would be wrong for the same reason a 21 nt
    floor is: BLAST trims, and demanding full query coverage hides exactly the
    hits the screen exists to catch. The full-25-nt accounting belongs in Python
    (specificity.full_half_mismatches), which retrieves the subject bases
    opposite the unaligned positions and compares them.

    COST: this widens what BLAST emits, by more than an order of magnitude
    against a whole reference transcriptome compared with a 21 nt floor. That
    is intrinsic, not tunable: most of the extra hits are the short exact
    cores, and whether one of them is a 4-mismatch off-target or noise cannot
    be decided without the subject. :func:`_warn_recount_cost` warns before a
    large screen.
    """
    if not getattr(spec, "off_target_full_half", False) or half_len <= 0:
        return float(spec.off_target_pident_min), int(spec.off_target_aln_min)
    keep = max(1, half_len - max(0, spec.off_target_min_mismatch - 1))
    pident = 100.0 * keep / (half_len + OFF_TARGET_PUSHDOWN_GAP_MARGIN)
    # run_blastn OVERRIDES blastn-short's defaults with -reward 1 -penalty -1,
    # so a mismatch costs 2 relative to a match, not 4. The two numbers are one
    # fact in two files: see constants.BLASTN_SHORT_MISMATCH_COST.
    min_len = max(1, half_len - BLASTN_SHORT_MISMATCH_COST
                  * max(0, spec.off_target_min_mismatch - 1))
    return pident, min_len


def _warn_recount_cost(
    spec,
    n_halves: int,
    emit: Optional[Callable[[str, str, Optional[float]], None]] = None,
    stage: str = ENGINE_BLAST,
    noun: str = "BLAST",
) -> Optional[str]:
    """Warn once, before the search, when the exact screen will be expensive.

    The exact recount is the only screen, and its cost is simply large: it
    multiplies the hit volume many times over against a reference
    transcriptome, and a whole-transcriptome run screens a very large number
    of halves. Hit volume at that scale can exhaust disk or walltime. So the
    cost is surfaced BEFORE the search commits to it, rather than discovered at
    90% of walltime.

    Returns the message (for tests / the manifest) or None.
    """
    if n_halves < RECOUNT_EXACT_WARN_HALVES:
        return None
    msg = (
        f"Off-target recount mode is 'exact' and this run screens {n_halves:,} "
        f"probe halves. The exact screen drops the BLAST length floor so trimmed "
        f"HSPs are still reported, which raises the hit volume by a large "
        f"factor on a reference transcriptome. Expect a correspondingly larger "
        f"BLAST stage. There is no cheaper mode: a stricter length floor "
        f"discards the trimmed hits the screen exists to catch. Budget for the "
        f"{noun} stage accordingly."
    )
    logger.warning(msg)
    if emit is not None:
        emit(stage, msg, None)
    return msg


def _derive_qcov_hsp_perc(min_aln: int, query_half_len: int) -> Optional[float]:
    """Server-side query-coverage floor for the BLAST off-target pushdown.

    Derived from the alignment-length floor with a 3-bp safety margin so a
    *gapped* boundary HSP (which loses a few query positions to gaps) is never
    dropped by BLAST before the Python classifier sees it. Kept no stricter than
    the post-filter. Returns None when there's no length floor to push down."""
    if min_aln <= 0 or query_half_len <= 0:
        return None
    target_len = max(1, min_aln - 3)
    return max(0.0, min(100.0, target_len / query_half_len * 100.0))


def _offtarget_recount_summary(
    subjects: Any, config: RunConfig,
) -> Dict[str, Any]:
    """Run-level provenance for the full-half off-target mismatch counts.

    Mirrors ``ProbeCandidate.offtarget_recount_basis`` at run scope. The basis is
    the WEAKEST one any hit used: a single fallback makes the whole run's counts
    lower bounds, and a manifest that averaged instead of taking the worst case
    would let a mostly-exact run be read as an exact one.
    """
    if not config.specificity.off_target_full_half:
        return {"basis": "legacy_local",
                "note": "off_target_full_half=False; no full-half count made"}
    if config.screen.engine != ENGINE_BLAST:
        # An end-to-end engine compares all 25 positions of every half itself,
        # so there is no retrieval to succeed or fail and no SubjectWindows to
        # report on. Falling through to `not_screened` below would say the run
        # was never screened, flatly contradicting the per-probe
        # `offtarget_recount_basis` column (which reads `exact`) and inviting a
        # reader to discard a perfectly good run.
        return {"basis": "exact",
                "note": (f"{config.screen.engine} engine: every hit spans the "
                         "full probe half by construction, so no subject "
                         "retrieval was needed"),
                "blastdbcmd_available": False,
                "blastdbcmd_calls": 0}
    if subjects is None:
        return {"basis": "not_screened"}
    total = subjects.n_exact + subjects.n_fallback
    return {
        "basis": ("conservative" if subjects.n_fallback else "exact"),
        "blastdbcmd_available": bool(subjects.available),
        "hits_counted": total,
        "hits_exact": subjects.n_exact,
        "hits_conservative": subjects.n_fallback,
        "blastdbcmd_calls": subjects.n_calls,
    }


def _blast_db_sizes(db_paths: List[str]) -> Dict[str, int]:
    """Total on-disk bytes of each BLAST DB prefix (sum of its index files).

    A DB-size proxy recorded in the manifest so a specificity "pass" is
    interpretable against the database it was screened over."""
    import glob
    sizes: Dict[str, int] = {}
    for configured in db_paths or []:
        # Key by what the user configured, but measure the prefix blastn read.
        prefix = resolve_db_prefix(configured)
        if prefix is None:
            continue          # absent means "unknown", never a bogus 0
        total = 0
        for f in glob.glob(prefix + ".*"):
            try:
                total += os.path.getsize(f)
            except OSError:
                pass
        sizes[configured] = total
    return sizes


# Above this the DB is not hashed (an absent key means "unknown", the same
# honest convention tool_versions uses). Nothing in practice comes close — the
# human transcriptome DB is ~1 GB — but a run must never be delayed by hours
# because someone pointed it at a genome-scale database.
_DB_HASH_MAX_BYTES = 20 * 1024 ** 3   # 20 GiB


def _blast_db_sha256s(db_paths: List[str]) -> Dict[str, str]:
    """SHA-256 of each BLAST DB actually used, keyed by DB prefix.

    This fills ``RunManifest.db_sha256s``, the field that answers "which
    database was this screened against?". ``blast_db_sizes_bytes`` alone is
    not enough: any two different databases can share a size.

    The digest covers the DB's on-disk index files (``<prefix>.*``, sorted, with
    each file's name and size folded in so a rename or a split can't collide) —
    NOT the source FASTA. That is deliberate: the index files are what blastn
    actually read. Never raises; a DB that can't be hashed is simply absent from
    the result.
    """
    import glob

    out: Dict[str, str] = {}
    for configured in db_paths or []:
        # Same resolution the search used; a configured FASTA path has no
        # sidecars of its own, so globbing it directly would find nothing.
        prefix = resolve_db_prefix(configured)
        if prefix is None:
            logger.debug("No BLAST DB found for %s", configured)
            continue
        files = sorted(f for f in glob.glob(prefix + ".*") if os.path.isfile(f))
        if not files:
            logger.debug("No BLAST DB files found for prefix %s", prefix)
            continue
        try:
            total = sum(os.path.getsize(f) for f in files)
            if total > _DB_HASH_MAX_BYTES:
                logger.warning(
                    "BLAST DB %s is %.1f GB — skipping the manifest checksum "
                    "(over the %.0f GB cap).",
                    prefix, total / 1024 ** 3, _DB_HASH_MAX_BYTES / 1024 ** 3,
                )
                continue
            h = hashlib.sha256()
            for f in files:
                h.update(f"{os.path.basename(f)}:{os.path.getsize(f)}\0".encode())
                with open(f, "rb") as fh:
                    for chunk in iter(lambda: fh.read(1 << 20), b""):
                        h.update(chunk)
            out[configured] = h.hexdigest()
        except OSError:
            logger.warning("Could not checksum BLAST DB %s", prefix, exc_info=True)
    return out


def _reference_sha256s(screen: Any,
                       compositions: Dict[str, Dict[str, Any]]) -> Dict[str, str]:
    """``RunManifest.db_sha256s``: the identity of every reference screened.

    What is digested is what the engine READ. The BLAST engine read a
    database's index files, so those are hashed (``_blast_db_sha256s``). The
    k-mer engine reads the FASTA itself, so the FASTA's bytes are hashed; the
    composition pass already computed that digest while loading, so it is
    reused rather than read twice.

    A screened run that ends up missing a reference's digest says so in the
    log, naming the reference. A --skip-screen run screened nothing and
    records ``{}``.
    """
    if screen.skip_screen:
        return {}
    if screen.engine == ENGINE_BLAST:
        out = _blast_db_sha256s(screen.db_paths)
    else:
        out = {}
        for configured in screen.db_paths or []:
            digest = (compositions.get(configured) or {}).get("sha256")
            if not digest:
                fasta = configured
                if not (os.path.exists(fasta) and not os.path.isdir(fasta)):
                    fasta = next(
                        (configured + ext for ext in
                         (".fa", ".fasta", ".fna", ".ffn", ".frn")
                         if os.path.exists(configured + ext)), "")
                try:
                    if fasta and os.path.getsize(fasta) <= _DB_HASH_MAX_BYTES:
                        digest = blast._sha256_file(fasta)
                except OSError:
                    logger.warning("Could not hash reference %s", fasta,
                                   exc_info=True)
            if digest:
                out[configured] = digest
    missing = [p for p in screen.db_paths or [] if p not in out]
    if missing:
        logger.warning(
            "The manifest records no sha256 for %d screened reference(s): %s. "
            "This run cannot prove which reference it used.",
            len(missing), ", ".join(missing))
    return out


def _reference_compositions(
    db_paths: List[str], *, skip_screen: bool,
) -> Dict[str, Dict[str, Any]]:
    """What each off-target reference CONTAINED, keyed by configured path.

    ``db_sha256s`` identifies the reference; this characterises it. They answer
    different questions and a Methods section needs both: two runs can cite the
    same hash and still be screened against references whose composition makes
    their rejection rates incomparable. For example, GRCh38 cdna.all includes
    ALT haplotype and PATCH scaffold records, a small share of its bases that
    can carry near-copies of primary-assembly genes and so produce spurious
    perfect-50bp off-target calls.

    Read from what the loader already measured during the screen
    (``reference_qc.measured``). The header-only scan is the fallback for a run
    that never built a k-mer reference — a BLAST-engine run, where ``blastn``
    searched a database instead — and it costs one sequential read with no
    per-base allocation.

    Never raises: a provenance record must not be the thing that fails a
    finished run.
    """
    out: Dict[str, Dict[str, Any]] = {}
    # Nothing was screened, so describing a reference the run never opened would
    # be a fabricated claim — the same rule _blast_db_sha256s follows.
    if skip_screen:
        return out
    for configured in db_paths or []:
        try:
            comp = reference_qc.measured(configured)
            if comp is None:
                # A configured path can be a BLAST DB prefix rather than a
                # FASTA; resolve the same way the k-mer engine does and give up
                # quietly if there is no readable FASTA behind it.
                fasta = configured
                if not (os.path.exists(fasta) and not os.path.isdir(fasta)):
                    fasta = next(
                        (configured + ext for ext in
                         (".fa", ".fasta", ".fna", ".ffn", ".frn")
                         if os.path.exists(configured + ext)),
                        "",
                    )
                if not fasta:
                    continue
                comp = reference_qc.measured(fasta) or reference_qc.scan_fasta(
                    fasta)
            out[configured] = comp.to_manifest()
        except Exception:       # pragma: no cover - provenance is never fatal
            logger.warning("Could not summarise composition of reference %s",
                           configured, exc_info=True)
    return out


def _portable_inputs(inputs: Dict[str, str], output_dir: str) -> Dict[str, str]:
    """The manifest's ``inputs`` with a portable path and the file's name.

    ``fasta`` becomes relative to the output directory. ``fasta_name``
    keeps the basename beside ``fasta_sha256``, so the record still says which
    file it was when the run tree has been moved or the path means nothing.
    """
    from .config_io import portable_path
    out = dict(inputs)
    fasta = out.get("fasta")
    if isinstance(fasta, str) and fasta:
        out["fasta_name"] = os.path.basename(fasta)
        out["fasta"] = portable_path(fasta, output_dir) or fasta
    return out


def _portable_keys(by_path: Dict[str, Any], output_dir: str) -> Dict[str, Any]:
    """A path-keyed manifest block rekeyed by portable paths."""
    from .config_io import portable_path
    return {(portable_path(k, output_dir) or k): v for k, v in by_path.items()}


def _hash_input_fasta(path: str) -> Dict[str, str]:
    """``{path, sha256, bytes}`` for the input FASTA (never raises).

    A bare path identifies a location rather than a file. Months later that
    path may hold a different (or no) sequence, so the hash and size are what
    make the run reproducible.
    """
    info: Dict[str, str] = {"fasta": path}
    try:
        from .blast import _sha256_file
        info["fasta_sha256"] = _sha256_file(path)
        info["fasta_bytes"] = str(os.path.getsize(path))
    except OSError:
        logger.warning("Could not hash input FASTA %s", path, exc_info=True)
    return info


def _code_version() -> Dict[str, str]:
    """Which build of the engine ran. Never empty, never raises.

    ``source`` says which kind of answer this is:

    * ``"git"``: the package was imported from a checkout that TRACKS it. Then
      ``commit``, ``dirty`` and (when git can describe it) ``describe``. Many
      commits can ship under one version string, so the commit is the precise
      answer.
    * ``"installed"``: no such checkout (an installed wheel, the container, the
      packaged app). Then ``installed_version``, which is ``__version__``.

    Outside a checkout an empty record would leave a run untraceable to any
    build (the usual case for an installed package in a container). A version
    string is coarser than a commit, but it is a real claim, not a placeholder.

    The checkout must track this package's own ``__init__.py``. Otherwise a
    package installed into a venv that happens to sit inside some unrelated git
    repository would stamp that repository's commit.
    """
    import subprocess

    from . import __version__

    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    installed = {"source": "installed", "installed_version": __version__}

    def _git(*args: str) -> Optional[str]:
        try:
            r = subprocess.run(
                ["git", "-C", pkg_dir, *args],
                capture_output=True, text=True, timeout=10,
            )
        except Exception:
            return None
        return r.stdout.strip() if r.returncode == 0 else None

    if _git("ls-files", "--error-unmatch", "__init__.py") is None:
        return installed
    commit = _git("rev-parse", "HEAD")
    if not commit:
        return installed
    info = {"source": "git", "commit": commit}
    status = _git("status", "--porcelain")
    # A run made from a modified working tree is NOT reproducible from the
    # commit alone; say so rather than implying the commit is the whole story.
    info["dirty"] = "true" if status else "false"
    describe = _git("describe", "--tags", "--always", "--dirty")
    if describe:
        info["describe"] = describe
    return info


def _collect_tool_versions(screen_used: bool = True) -> Dict[str, str]:
    """External tool + library versions for the run manifest (never raises).

    Recorded so a run's provenance states what actually produced a result
    set. A tool that can't be found or probed is simply left out — an
    absent key means "unknown", which is honest, whereas a placeholder string
    would be a fabricated provenance claim.

    ``screen_used`` means "blastn ran": True only when the BLAST engine screened.
    It is False for a ``--skip-screen`` run and for a k-mer engine run, neither
    of which invokes blastn: recording the version of a binary that never ran
    would be the same kind of fabricated claim as a placeholder.
    """
    versions: Dict[str, str] = {}

    if screen_used:
        from .blast import blastn_version
        bv = blastn_version()
        if bv:
            versions["blastn"] = bv

    # primer3 is the primer3-py binding, so this is the wrapper version (what
    # pip pins), not the bundled libprimer3 it wraps.
    try:
        import primer3
        pv = getattr(primer3, "__version__", None)
        if not pv:
            from importlib.metadata import version as _pkg_version
            pv = _pkg_version("primer3-py")
        versions["primer3"] = str(pv)
    except Exception:
        logger.debug("Could not determine primer3 version", exc_info=True)

    versions["python"] = platform.python_version()
    return versions


def _is_windows() -> bool:
    """True on Windows. A function, not a constant, so tests can simulate it."""
    return os.name == "nt"


class _WinProcessProbe:
    """CPU time, peak memory and process-tree memory on Windows, via ctypes.

    WHY. ``resource`` does not exist on Windows and ``ps`` is not a Windows
    command, so without this a Windows run could not measure ``cpu_s``,
    ``peak_rss_mb`` or the tree peak. ``psutil`` would answer all three but is deliberately not a dependency, so
    the handful of kernel32 calls it would make are made here directly.

    WHAT EACH FIGURE IS, since none is the POSIX quantity of the same name.
    CPU is ``GetProcessTimes`` user plus kernel time for THIS process only:
    Windows keeps no total for exited children, so a worker pool's CPU is not
    included. Peak memory is ``PeakWorkingSetSize`` for this process only. The
    tree figure sums ``WorkingSetSize`` over this process and every descendant
    found in a ``CreateToolhelp32Snapshot`` walk. Working set is the Windows
    analogue of RSS.

    Constructed only through :func:`_windows_probe`, which returns None off
    Windows or when kernel32 cannot be loaded. Every method returns None rather
    than raising: this is instrumentation and must never fail a run.
    """

    _PROCESS_QUERY_LIMITED_INFORMATION = 0x1000
    # GetProcessMemoryInfo's documentation asks for PROCESS_VM_READ as well.
    _PROCESS_VM_READ = 0x0010
    _TH32CS_SNAPPROCESS = 0x00000002

    def __init__(self) -> None:
        import ctypes as c

        dword, long_, handle, bool_ = c.c_uint32, c.c_int32, c.c_void_p, c.c_int32

        class FILETIME(c.Structure):
            _fields_ = [("lo", dword), ("hi", dword)]

        class PROCESS_MEMORY_COUNTERS(c.Structure):
            _fields_ = [("cb", dword), ("PageFaultCount", dword),
                        ("PeakWorkingSetSize", c.c_size_t),
                        ("WorkingSetSize", c.c_size_t),
                        ("QuotaPeakPagedPoolUsage", c.c_size_t),
                        ("QuotaPagedPoolUsage", c.c_size_t),
                        ("QuotaPeakNonPagedPoolUsage", c.c_size_t),
                        ("QuotaNonPagedPoolUsage", c.c_size_t),
                        ("PagefileUsage", c.c_size_t),
                        ("PeakPagefileUsage", c.c_size_t)]

        class PROCESSENTRY32W(c.Structure):
            _fields_ = [("dwSize", dword), ("cntUsage", dword),
                        ("th32ProcessID", dword),
                        ("th32DefaultHeapID", c.c_size_t),     # ULONG_PTR
                        ("th32ModuleID", dword), ("cntThreads", dword),
                        ("th32ParentProcessID", dword),
                        ("pcPriClassBase", long_), ("dwFlags", dword),
                        ("szExeFile", c.c_wchar * 260)]

        # getattr: WinDLL exists only on Windows, and a plain attribute access
        # would fail type-checking everywhere else.
        k32 = getattr(c, "WinDLL")("kernel32", use_last_error=True)
        pft = c.POINTER(FILETIME)
        k32.GetCurrentProcess.argtypes = []
        k32.GetCurrentProcess.restype = handle
        k32.GetProcessTimes.argtypes = [handle, pft, pft, pft, pft]
        k32.GetProcessTimes.restype = bool_
        k32.K32GetProcessMemoryInfo.argtypes = [
            handle, c.POINTER(PROCESS_MEMORY_COUNTERS), dword]
        k32.K32GetProcessMemoryInfo.restype = bool_
        k32.CreateToolhelp32Snapshot.argtypes = [dword, dword]
        k32.CreateToolhelp32Snapshot.restype = handle
        k32.Process32FirstW.argtypes = [handle, c.POINTER(PROCESSENTRY32W)]
        k32.Process32FirstW.restype = bool_
        k32.Process32NextW.argtypes = [handle, c.POINTER(PROCESSENTRY32W)]
        k32.Process32NextW.restype = bool_
        k32.OpenProcess.argtypes = [dword, bool_, dword]
        k32.OpenProcess.restype = handle
        k32.CloseHandle.argtypes = [handle]
        k32.CloseHandle.restype = bool_

        self._c = c
        self._k32 = k32
        self._filetime = FILETIME
        self._pmc = PROCESS_MEMORY_COUNTERS
        self._entry = PROCESSENTRY32W
        self._invalid_handle = (1 << (8 * c.sizeof(c.c_void_p))) - 1

    # -- primitives ----------------------------------------------------------
    def _times(self, h: Any) -> Optional[Tuple[int, int]]:
        """``(creation time, user+kernel time)`` in 100 ns units, or None."""
        ft = [self._filetime() for _ in range(4)]
        if not self._k32.GetProcessTimes(h, *[self._c.byref(f) for f in ft]):
            return None
        creation, _exit, kernel, user = [(f.hi << 32) | f.lo for f in ft]
        return creation, kernel + user

    def _memory(self, h: Any) -> Optional[Any]:
        pmc = self._pmc()
        pmc.cb = self._c.sizeof(pmc)
        if not self._k32.K32GetProcessMemoryInfo(h, self._c.byref(pmc), pmc.cb):
            return None
        return pmc

    def _processes(self) -> Optional[List[Tuple[int, int]]]:
        """``(pid, parent pid)`` for every process on the machine, or None."""
        snap = self._k32.CreateToolhelp32Snapshot(self._TH32CS_SNAPPROCESS, 0)
        if not snap or snap == self._invalid_handle:
            return None
        try:
            entry = self._entry()
            entry.dwSize = self._c.sizeof(entry)
            out: List[Tuple[int, int]] = []
            ok = self._k32.Process32FirstW(snap, self._c.byref(entry))
            while ok:
                out.append((int(entry.th32ProcessID),
                            int(entry.th32ParentProcessID)))
                ok = self._k32.Process32NextW(snap, self._c.byref(entry))
            return out
        finally:
            self._k32.CloseHandle(snap)

    def _probe(self, pid: int) -> Optional[Tuple[int, int]]:
        """``(creation time, working set bytes)`` for ``pid``, or None."""
        h = self._k32.OpenProcess(
            self._PROCESS_QUERY_LIMITED_INFORMATION | self._PROCESS_VM_READ,
            False, pid)
        if not h:
            return None
        try:
            times = self._times(h)
            mem = self._memory(h)
            if times is None or mem is None:
                return None
            return times[0], int(mem.WorkingSetSize)
        finally:
            self._k32.CloseHandle(h)

    # -- the three figures ---------------------------------------------------
    def cpu_seconds(self) -> Optional[float]:
        try:
            times = self._times(self._k32.GetCurrentProcess())
        except Exception:                   # pragma: no cover - Windows only
            return None
        return None if times is None else times[1] / 1e7

    def peak_rss_mb(self) -> Optional[float]:
        try:
            mem = self._memory(self._k32.GetCurrentProcess())
        except Exception:                   # pragma: no cover - Windows only
            return None
        if mem is None:
            return None
        return round(mem.PeakWorkingSetSize / (1024.0 * 1024.0), 1)

    def tree_rss(self, root_pid: int) -> Optional[Tuple[float, int]]:
        """``(summed working set in MiB, process count)`` for the tree, or None.

        A Windows parent pid is never cleared when the parent exits, and pids
        are reused, so a process can name as its parent a pid that now belongs
        to a different, later process. A child created BEFORE its supposed
        parent is therefore skipped. A child that exits between the snapshot
        and being opened is skipped too, with its own descendants, which is
        the lower-bound direction the report already states.
        """
        try:
            procs = self._processes()
            if procs is None:
                return None
            root = self._probe(root_pid)
            if root is None:
                return None
            kids: Dict[int, List[int]] = {}
            for pid, ppid in procs:
                if pid != ppid:
                    kids.setdefault(ppid, []).append(pid)
            total, n = root[1], 1
            seen, stack = {root_pid}, [(root_pid, root[0])]
            while stack:
                pid, created = stack.pop()
                for child in kids.get(pid, ()):
                    if child in seen:
                        continue
                    seen.add(child)
                    info = self._probe(child)
                    if info is None or info[0] < created:
                        continue
                    total += info[1]
                    n += 1
                    stack.append((child, info[0]))
            return total / (1024.0 * 1024.0), n
        except Exception:                   # pragma: no cover - Windows only
            return None


@functools.lru_cache(maxsize=1)
def _load_windows_probe() -> Optional[_WinProcessProbe]:
    try:
        return _WinProcessProbe()
    except Exception:
        logger.debug("Windows resource back end unavailable", exc_info=True)
        return None


def _windows_probe() -> Optional[_WinProcessProbe]:
    """The Windows back end, or None off Windows or if it cannot load."""
    if not _is_windows():
        return None
    return _load_windows_probe()


class _TreeRssSampler:
    """Peak RSS of the WHOLE PROCESS TREE, sampled, in MiB.

    WHY THIS EXISTS, AND WHY ``_peak_rss_mb`` CANNOT ANSWER IT.
    ``_StageTimer._peak_rss_mb`` returns
    ``max(getrusage(RUSAGE_SELF).ru_maxrss, getrusage(RUSAGE_CHILDREN).ru_maxrss)``
    and ``RUSAGE_CHILDREN.ru_maxrss`` is the high-water mark of the **largest
    single child**, never the sum of children alive at once. A screen running 24
    workers holds the parent plus 24 worker resident sets simultaneously and that
    field cannot see it.

    Reading the per-process figure as the job's footprint is a real hazard: on
    a multi-worker whole-transcriptome screen the scheduler's own
    ``resources_used.mem`` can be several times the per-process peak, so a
    memory request sized from the per-process number would get the screen
    OOM-killed. This sampler measures the quantity a request has to cover.

    HOW. One ``ps`` call per sample, which works identically on macOS and Linux
    and needs no ``/proc`` parsing and no third-party dependency — ``psutil`` is
    deliberately not a dependency of this project. The parent and every
    descendant are summed, so the figure is what the cgroup would charge. A
    sample is cheap: at the default interval a 12 h run costs a few thousand
    subprocess calls against a screen that runs for hours.

    LIMITS, STATED. It is a SAMPLE, so a spike shorter than the interval is
    missed and the reported peak is a lower bound. Shared pages between parent
    and children are counted once per process, so the sum OVERSTATES unique
    memory — which is the safe direction for sizing a request, but it is not the
    same quantity as PSS. Compare against ``qstat -xf`` ``resources_used.mem``
    rather than treating this as ground truth.

    WINDOWS. There is no ``ps``. The same walk is made through
    :class:`_WinProcessProbe`, summing working set rather than RSS.

    NOTHING MEASURED IS NULL. When no sample succeeded the peak, its GiB form
    and the process count are ``None``, ``available`` is False, ``reason`` says
    why, and ``method`` is ``None``: a 0.0 reads as a measurement.
    """

    #: What each back end actually sums, for the report's ``method``.
    _METHODS = {
        "ps": "sum of ps RSS over the parent process and every "
              "descendant, sampled; a lower bound on the true peak and "
              "an overstatement of unique memory (shared pages counted "
              "per process). Compare with qstat resources_used.mem.",
        "windows": "sum of WorkingSetSize (K32GetProcessMemoryInfo) over the "
                   "parent process and every descendant found in a "
                   "CreateToolhelp32Snapshot walk, sampled; a lower bound on "
                   "the true peak and an overstatement of unique memory "
                   "(shared pages counted per process).",
    }

    #: Steady-state seconds between samples, reached by backing off from
    #: :data:`FIRST_INTERVAL_S`. A FIXED interval of this length would give a
    #: run lasting seconds a single sample at t=0, before memory has grown,
    #: and a "peak" far below the parent's own high-water mark.
    INTERVAL_S = 5.0
    #: First interval. Short, so a run measured in seconds still gets a usable
    #: series, then doubled until it reaches INTERVAL_S — so a multi-hour screen
    #: still pays a negligible cost.
    FIRST_INTERVAL_S = 0.2
    #: Below this many samples the peak is not reported as a peak. One sample is
    #: an observation, not a maximum, and labelling it a peak is how a number
    #: with no support reaches a table.
    MIN_SAMPLES_FOR_PEAK = 5

    def __init__(self, interval_s: Optional[float] = None) -> None:
        self.interval_s = self.INTERVAL_S if interval_s is None else interval_s
        self.peak_mib: Optional[float] = None
        self.peak_procs: Optional[int] = None
        self.samples = 0
        self.failures = 0
        self.backend = "windows" if _is_windows() else "ps"
        self.last_error: Optional[str] = None
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None

    # -- measurement -------------------------------------------------------
    @staticmethod
    def _sum_tree(root_pid: int,
                  rows: Sequence[Tuple[int, int, int]]) -> Tuple[float, int]:
        """``(summed MiB, process count)`` from ``(pid, ppid, rss KiB)`` rows."""
        kids: Dict[int, List[int]] = {}
        rss: Dict[int, int] = {}
        for pid, ppid, kb in rows:
            rss[pid] = kb
            kids.setdefault(ppid, []).append(pid)
        # Walk the tree from root_pid. Iterative, so a deep tree cannot recurse
        # past the limit, and `seen` guards against a cycle from a reused pid.
        total_kb, n, seen, stack = 0, 0, {root_pid}, [root_pid]
        while stack:
            pid = stack.pop()
            if pid in rss:
                total_kb += rss[pid]
                n += 1
            for child in kids.get(pid, ()):
                if child not in seen:
                    seen.add(child)
                    stack.append(child)
        return total_kb / 1024.0, n

    def _snapshot(self, root_pid: int) -> Optional[Tuple[float, int]]:
        """``(summed RSS in MiB, process count)`` for ``root_pid`` and its tree.

        Returns None on any failure, with the cause in ``last_error``, rather
        than raising: this is instrumentation and must never be able to fail a
        multi-day run. A failure must never appear as a zero in the report.
        """
        if self.backend == "windows":
            probe = _windows_probe()
            if probe is None:
                self.last_error = "the Windows process back end did not load"
                return None
            snap = probe.tree_rss(root_pid)
            if snap is None:
                self.last_error = "the Windows process snapshot failed"
            return snap
        # Imported locally, matching this module's existing convention (see
        # _code_version): the GUI pays the import cost of this module at
        # startup, and subprocess is only needed on the sampling path.
        import subprocess
        try:
            out = subprocess.run(
                ["ps", "-axo", "pid=,ppid=,rss="],
                capture_output=True, text=True, timeout=30,
            ).stdout
        except (OSError, subprocess.SubprocessError) as exc:
            self.last_error = (f"ps could not be run "
                               f"({type(exc).__name__}: {exc})")
            return None
        rows: List[Tuple[int, int, int]] = []
        for line in out.splitlines():
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                rows.append((int(parts[0]), int(parts[1]), int(parts[2])))
            except ValueError:
                continue
        mib, n = self._sum_tree(root_pid, rows)
        if not n:
            self.last_error = "ps output did not list this process"
            return None
        return mib, n

    # -- lifecycle ---------------------------------------------------------
    def _loop(self) -> None:
        me = os.getpid()
        wait = min(self.FIRST_INTERVAL_S, self.interval_s)
        while not self._stop.is_set():
            snap = self._snapshot(me)
            if snap is not None:
                mib, n = snap
                self.samples += 1
                if self.peak_mib is None or mib > self.peak_mib:
                    self.peak_mib, self.peak_procs = mib, n
            else:
                self.failures += 1
            self._stop.wait(wait)
            wait = min(wait * 2.0, self.interval_s)

    def start(self) -> None:
        if self._thread is not None:
            return
        self._thread = threading.Thread(
            target=self._loop, name="flexprobe-tree-rss", daemon=True)
        self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        t, self._thread = self._thread, None
        if t is not None:
            t.join(timeout=self.interval_s + 5.0)

    def report(self) -> Dict[str, Any]:
        enough = self.samples >= self.MIN_SAMPLES_FOR_PEAK
        peak = self.peak_mib
        measured = peak is not None and self.samples > 0
        return {
            # FALSE means nothing below was measured: the figures are null and
            # `reason` says why. Never read a null here as zero.
            "available": measured,
            "reason": None if measured else (
                self.last_error or "no sample was taken before the sampler "
                                   "stopped"),
            "backend": self.backend,
            "peak_tree_rss_mib": round(peak, 1) if peak is not None else None,
            "peak_tree_rss_gib": (round(peak / 1024.0, 2)
                                  if peak is not None else None),
            "processes_at_peak": self.peak_procs if measured else None,
            "samples": self.samples,
            # FALSE means the figure above is an observation, not a maximum.
            # Read it as such or not at all; do not put it in a table.
            "sufficient_samples": enough,
            "min_samples_for_peak": self.MIN_SAMPLES_FOR_PEAK,
            "steady_interval_s": self.interval_s,
            "first_interval_s": self.FIRST_INTERVAL_S,
            "failed_samples": self.failures,
            # Names the back end that produced the figure, and only when one did.
            "method": self._METHODS[self.backend] if measured else None,
        }


def _cores(cpu_s: Optional[float], wall_s: float) -> Optional[float]:
    """Mean cores over a span: None when CPU was not measured."""
    if cpu_s is None:
        return None
    return round(cpu_s / wall_s, 2) if wall_s > 0 else 0.0


class _StageTimer:
    """Per-stage wall clock, CPU and peak memory — recorded AS IT GOES.

    Driven off the existing ``_emit`` stage boundaries rather than a second set
    of hand-placed timing calls: whenever a progress event names a different
    stage than the one currently running, the elapsed time since the last
    boundary is banked against the stage that *was* running. So a segment of
    wall clock is always attributed to the stage the pipeline last announced,
    and adding or reordering a stage keeps the accounting correct for free.

    Stages accumulate (rather than being recorded once) because the streaming
    paths interleave: tile/filter alternate per gene in the streamed filter
    loop, and blast/specificity alternate per chunk when specificity is scored
    inside the BLAST loop. The per-stage totals are therefore the true share of
    wall clock spent in each, not a single contiguous span.

    WHY IT APPENDS TO DISK INSTEAD OF ONLY SUMMARISING AT THE END. A snapshot
    taken when the run finishes tells you nothing about a run that never
    finishes — and the expensive failures here are all kills. A slow step that
    runs after the last timing snapshot shows up only as an unexplained gap
    between the scheduler's walltime and the manifest's own total. Every
    boundary is therefore flushed to ``run_timeline.jsonl`` as it happens, so a
    job killed at its walltime still says exactly which step it died in and how
    long that step had been running.

    Each record also carries CPU seconds and peak RSS, which answers the other
    recurring question — how much memory and how many cores a step actually
    needs — without re-running anything. ``ru_maxrss`` is a high-water mark, so
    the per-stage value is "peak so far", not "peak within this stage".

    CPU covers this process plus REAPED children, so a worker pool's time lands
    in whichever stage was open when the pool shut down. That is right for
    tile/filter (the pool is scoped to it) and is why the screen's ~1.0 core
    reads honestly. On Windows CPU and peak memory cover this process only; see
    :class:`_WinProcessProbe`.

    A figure that could not be measured is ``None`` (null in JSON), never 0.0,
    and ``run_timings.json``'s ``measurement`` block says which were available,
    by what method, and why not.
    """

    #: ru_maxrss is kilobytes on Linux, bytes on macOS/BSD.
    _RSS_DIVISOR = 1024.0 * 1024.0 if sys.platform == "darwin" else 1024.0

    def __init__(self) -> None:
        self._start = time.monotonic()
        self._mark = self._start
        self._stage: Optional[str] = None
        self._totals: Dict[str, float] = {}
        self._cpu_totals: Dict[str, float] = {}
        self._peak_rss: Dict[str, float] = {}
        self._events: Dict[str, int] = {}
        self._cpu_mark: Optional[float] = self._cpu_now()
        self._sink: Optional[str] = None
        self._seq = 0
        # Created lazily in open_sink, so helper timers do not spawn a thread.
        self._tree_rss: Optional["_TreeRssSampler"] = None

    # -- resource sampling -------------------------------------------------
    @staticmethod
    def _cpu_now() -> Optional[float]:
        """User+sys seconds for this process and its reaped children.

        None when nothing can measure it, so an unmeasured stage never reads
        as having used no CPU at all.
        """
        if resource is None:
            probe = _windows_probe()
            return probe.cpu_seconds() if probe is not None else None
        try:
            me = resource.getrusage(resource.RUSAGE_SELF)
            kids = resource.getrusage(resource.RUSAGE_CHILDREN)
        except (OSError, ValueError):       # pragma: no cover - platform guard
            return None
        return (me.ru_utime + me.ru_stime + kids.ru_utime + kids.ru_stime)

    @classmethod
    def _peak_rss_mb(cls) -> Optional[float]:
        """Peak RSS in MB for this process, and the largest single child.

        None when nothing can measure it, for the same reason as _cpu_now.
        """
        if resource is None:
            probe = _windows_probe()
            return probe.peak_rss_mb() if probe is not None else None
        try:
            me = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            kids = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        except (OSError, ValueError):       # pragma: no cover - platform guard
            return None
        return round(max(me, kids) / cls._RSS_DIVISOR, 1)

    @classmethod
    def measurement(cls) -> Dict[str, Dict[str, Any]]:
        """Which resource figures this platform measured, how, and why not.

        Written into ``run_timings.json`` so a null in ``stages`` carries its
        explanation next to it rather than in a log line.
        """
        if resource is not None:
            cpu_method: Optional[str] = (
                "getrusage user+sys time, RUSAGE_SELF plus RUSAGE_CHILDREN "
                "(reaped children only)")
            rss_method: Optional[str] = (
                "max of getrusage ru_maxrss over RUSAGE_SELF and "
                "RUSAGE_CHILDREN (the largest single child, not a sum)")
        elif _windows_probe() is not None:
            cpu_method = ("GetProcessTimes user+kernel time for this process "
                          "only; Windows keeps no total for exited children")
            rss_method = ("K32GetProcessMemoryInfo PeakWorkingSetSize for this "
                          "process only")
        else:
            cpu_method = rss_method = None
        none_reason = ("no source on this platform: the resource module is "
                       "unavailable and no Windows back end loaded")

        def _block(value: Optional[float],
                   method: Optional[str]) -> Dict[str, Any]:
            ok = value is not None
            return {
                "available": ok,
                "method": method if ok else None,
                "reason": None if ok else (
                    none_reason if method is None
                    else f"{method} returned no value"),
            }

        return {"cpu": _block(cls._cpu_now(), cpu_method),
                "peak_rss": _block(cls._peak_rss_mb(), rss_method)}

    # -- sink --------------------------------------------------------------
    def open_sink(self, path: str) -> None:
        """Start appending one JSON record per boundary to ``path``.

        Called once ``output_dir`` exists. A fresh run truncates any timeline
        from a previous attempt in the same directory; a RESUMED run appends,
        so the chain's links read as one continuous timeline.
        """
        self._sink = path
        # Start the concurrent-tree sampler here, not in __init__: __init__ runs
        # for short-lived helper timers too, and a thread per helper would be
        # noise. open_sink is called once, by a real run, in its output dir.
        if self._tree_rss is None:
            self._tree_rss = _TreeRssSampler()
            self._tree_rss.start()

    def _write(self, msg: str, pct: Optional[float], wall: float,
               cpu: Optional[float]) -> None:
        if not self._sink:
            return
        self._seq += 1
        rec = {
            "seq": self._seq,
            "t": datetime.now(timezone.utc).isoformat(),
            "elapsed_s": round(time.monotonic() - self._start, 3),
            "stage": self._stage,
            "msg": msg,
            "pct": pct,
            "step_s": round(wall, 3),
            "step_cpu_s": None if cpu is None else round(cpu, 3),
            "stage_total_s": round(self._totals.get(self._stage or "", 0.0), 3),
            "peak_rss_mb": self._peak_rss_mb(),
        }
        try:
            with open(self._sink, "a", encoding="utf-8") as fh:
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                os.fsync(fh.fileno())
        except OSError:
            # Never let instrumentation fail a multi-day run.
            logger.debug("Could not append to %s", self._sink, exc_info=True)

    # -- boundaries --------------------------------------------------------
    def enter(self, stage: str, msg: str = "",
              pct: Optional[float] = None) -> None:
        now = time.monotonic()
        cpu = self._cpu_now()
        wall_delta = now - self._mark
        cpu_delta = (max(0.0, cpu - self._cpu_mark)
                     if cpu is not None and self._cpu_mark is not None
                     else None)
        if self._stage is not None:
            self._totals[self._stage] = (
                self._totals.get(self._stage, 0.0) + wall_delta
            )
            # Unmeasured CPU and memory are left unbanked, so the stage reports
            # null for them instead of an accumulated 0.0.
            if cpu_delta is not None:
                self._cpu_totals[self._stage] = (
                    self._cpu_totals.get(self._stage, 0.0) + cpu_delta
                )
            self._events[self._stage] = self._events.get(self._stage, 0) + 1
            rss = self._peak_rss_mb()
            if rss is not None:
                self._peak_rss[self._stage] = max(
                    self._peak_rss.get(self._stage, 0.0), rss
                )
        prev_stage = self._stage
        self._stage = stage
        # Attribute the record to the stage that just ENDED when there was one,
        # so a boundary line reads as "this is what that step cost".
        if prev_stage is not None:
            self._stage = prev_stage
            self._write(msg, pct, wall_delta, cpu_delta)
            self._stage = stage
        self._mark = now
        self._cpu_mark = cpu

    def snapshot(self) -> Dict[str, float]:
        """Per-stage wall totals plus ``total``, without ending the running stage.

        This is the manifest's ``elapsed_seconds``, and `record_stage_elapsed` folds
        later stages into the same dict. Richer per-stage detail lives in
        ``resources()`` and ``run_timeline.jsonl`` instead of changing this.
        """
        now = time.monotonic()
        totals = dict(self._totals)
        if self._stage is not None:
            totals[self._stage] = (
                totals.get(self._stage, 0.0) + (now - self._mark)
            )
        out = {k: round(v, 3) for k, v in totals.items()}
        out["total"] = round(now - self._start, 3)
        return out

    def resources(self) -> Dict[str, Dict[str, Optional[float]]]:
        """Per-stage wall, CPU, mean cores and peak RSS — the Methods table.

        ``cores_mean`` is cpu/wall for the stage: ~1.0 says the stage is
        serial, and is the number that decides whether parallelising it can
        help at all.

        Read it with the children caveat above: a worker pool's whole CPU is
        credited at the boundary where the pool is reaped, so a SHORT stage
        that happens to close a pool can report more cores than the node has.
        On the multi-hour stages this exists to measure, it does not arise.

        ``cpu_s``, ``cores_mean`` and ``peak_rss_mb`` are None when they could
        not be measured; ``wall_s`` always can be.
        """
        now = time.monotonic()
        cpu = self._cpu_now()
        walls = dict(self._totals)
        cpus = dict(self._cpu_totals)
        if self._stage is not None:
            walls[self._stage] = walls.get(self._stage, 0.0) + (now - self._mark)
            if cpu is not None and self._cpu_mark is not None:
                cpus[self._stage] = (
                    cpus.get(self._stage, 0.0) + max(0.0, cpu - self._cpu_mark)
                )
        live = self._peak_rss_mb()
        out: Dict[str, Dict[str, Optional[float]]] = {}
        for stage, wall in walls.items():
            cpu_s = cpus.get(stage)
            # The still-open stage has no banked high-water yet; sample now.
            peak = self._peak_rss.get(stage)
            if peak is None and stage == self._stage:
                peak = live
            out[stage] = {
                "wall_s": round(wall, 3),
                "cpu_s": None if cpu_s is None else round(cpu_s, 3),
                "cores_mean": _cores(cpu_s, wall),
                "peak_rss_mb": peak,
                "events": float(self._events.get(stage, 0)),
            }
        total_wall = now - self._start
        out["total"] = {
            "wall_s": round(total_wall, 3),
            "cpu_s": None if cpu is None else round(cpu, 3),
            "cores_mean": _cores(cpu, total_wall),
            "peak_rss_mb": self._peak_rss_mb(),
            "events": float(self._seq),
        }
        return out

    def _tree_rss_report(self) -> Optional[Dict[str, Any]]:
        """The tree sampler's report, plus the comparison that makes it usable.

        NEITHER FIGURE DOMINATES THE OTHER. ``ru_maxrss`` is an EXACT high-water mark, but
        only of the parent and the largest single child. The tree figure SUMS
        every concurrent process, but it is SAMPLED, so a spike between samples
        is missed and it is a lower bound.

        So on a long screen with many workers the tree figure wins, because
        summing 24 workers beats any one of them. On a short single-process run
        the exact figure wins, because there is nothing to sum and sampling can
        miss the peak.

        ``recommended_floor_mib`` is therefore the MAXIMUM of the two: the best
        available lower bound on what the job really needed. It is a floor, not
        an estimate. Size a scheduler request above it, and check it against
        ``qstat -xf`` ``resources_used.mem`` once the job has run.
        """
        if self._tree_rss is None:
            return None
        rep = self._tree_rss.report()
        exact = self._peak_rss_mb()          # parent and largest single child
        tree = rep["peak_tree_rss_mib"]
        rep["parent_exact_rss_mib"] = exact
        # Either figure can be null (not measured). The floor is the larger of
        # those that exist, and null only when neither does.
        known = [v for v in (tree, exact) if v is not None]
        floor = round(max(known), 1) if known else None
        rep["recommended_floor_mib"] = floor
        rep["recommended_floor_gib"] = (
            round(floor / 1024.0, 2) if floor is not None else None)
        if tree is None and exact is None:
            rep["which_is_higher"] = None
        elif exact is None or (tree is not None and tree > exact):
            rep["which_is_higher"] = "tree_sampled"
        else:
            rep["which_is_higher"] = "parent_exact"
        return rep

    def write_report(self, path: str) -> None:
        """Write ``run_timings.json`` — the whole run's cost, one file.

        Written LAST, after every stage including the summary build, so unlike
        the manifest's ``elapsed_seconds`` it has no blind spot at the tail.
        """
        try:
            with open(path, "w", encoding="utf-8") as fh:
                if self._tree_rss is not None:
                    self._tree_rss.stop()
                json.dump({
                    # 3: an unmeasured cpu_s / cores_mean / peak_rss_mb is null
                    # rather than 0.0, `measurement` says which were available,
                    # and peak_tree_rss carries `available` and `reason`.
                    "schema": 3,
                    "stages": self.resources(),
                    "measurement": self.measurement(),
                    # The quantity `stages[*].peak_rss_mb` CANNOT report: the
                    # summed resident set of the parent and every concurrent
                    # worker. See _TreeRssSampler for why the per-process
                    # figure must not be read as a job footprint.
                    "peak_tree_rss": self._tree_rss_report(),
                    "timeline": os.path.basename(self._sink) if self._sink
                    else None,
                }, fh, indent=2, sort_keys=True)
                fh.write("\n")
        except OSError:
            logger.debug("Could not write %s", path, exc_info=True)


def _resolve_probe_barcode_seq(chemistry, label: str) -> str:
    if not chemistry.has_probe_barcode:
        return ""
    if not label:
        label = chemistry.default_probe_barcode
    table = FLEX_V1_PROBE_BARCODES
    return table.get(label, label)


def _pair_half(source: ProbeCandidate, oligo: str) -> ProbeCandidate:
    """Build one half of a ProbePair from the tile candidate it came from.

    The half is a COPY of ``source`` with a single field overridden: ``sequence``
    becomes the assembled oligo instead of the 50-bp mRNA window. Everything
    else — the specificity verdict and the off-target evidence behind it, the
    recount basis, isoform coverage, the score components, the degenerate
    parent/variant role — is carried over verbatim, because the half IS that
    candidate; only the strand it reports differs.

    Copying rather than enumerating the fields is the point. Halves built by
    hand from a curated list silently default every field the list misses: a
    fully screened run would hand callers ``pair.lhs.offtarget_recount_basis ==
    "not_screened"`` while results.csv says ``exact``, and
    ``pair.rhs.probe_role`` would read "plain" for a degenerate variant,
    quietly disarming the order writer's last-line consensus check. A copy
    cannot drift that way: a new field is provenance on the halves the day it
    is added.

    ``dataclasses.replace`` is shallow, so the list/dict fields would otherwise
    be SHARED with the scored candidate that results.csv is written from; the
    containers are deep-copied so mutating a half can never rewrite the artefact
    (or the other half). Found by walking ``fields()`` rather than by naming
    them, for the same no-drift reason.
    """
    half = dataclasses.replace(source, sequence=oligo)
    for f in dataclasses.fields(half):
        value = getattr(half, f.name)
        if isinstance(value, (list, dict, set)):
            setattr(half, f.name, copy.deepcopy(value))
    return half


logger = logging.getLogger("flexprobe.pipeline")


def effective_ranking_weights(config: RunConfig, *,
                              gene_aware: bool) -> RankingWeights:
    """The ``RankingWeights`` a run ranks with (and records).

    ``config.ranking_weights`` holds only the user's overrides; the run then
    adjusts two weights by mode, below. One function, so the ranking stage and
    the run records (``manifest.json`` and ``run_config.json``, key
    ``effective_ranking_weights``) cannot disagree.
    """
    # Build the ranking weights once: start from the defaults and apply any
    # per-run overrides from config.ranking_weights (validated in
    # RunConfig.validate, so unknown keys never reach here).
    from .ranking import (
        DEFAULT_ISOFORM_COVERAGE_WEIGHT, DEFAULT_WEIGHTS,
    )
    if config.ranking_weights:
        from dataclasses import asdict as _asdict
        weights = RankingWeights(**{**_asdict(DEFAULT_WEIGHTS), **config.ranking_weights})
    else:
        weights = DEFAULT_WEIGHTS
    # In gene-aware mode, give isoform coverage a sensible default weight unless
    # the user explicitly set one — otherwise the feature wouldn't affect
    # ranking. ``replace`` makes a fresh instance, never mutating DEFAULT_WEIGHTS.
    if gene_aware and not (config.ranking_weights and
                           "isoform_coverage" in config.ranking_weights):
        from dataclasses import replace as _replace
        weights = _replace(weights, isoform_coverage=DEFAULT_ISOFORM_COVERAGE_WEIGHT)
    # On-target mismatch placement assumes every input record is a record the
    # probe is SUPPOSED to bind, which is only true of a pan-strain
    # ("conserved") design. In strain_unique mode the mismatches against the
    # other strains ARE the deliverable, and in constitutive (isoform) mode a
    # probe need only hit its own isoform group at all (StrainRule.PAN_STRAIN);
    # in both, penalising those mismatches would rank the design backwards —
    # rewarding the least discriminating probe. So the term is switched off
    # there unless the user asked for it by name. Single-sequence runs need no
    # such guard: they carry no strain_binding, so the term is inapplicable and
    # drops out on its own.
    if (config.multi_strain and config.consensus_strategy != "conserved"
            and not (config.ranking_weights
                     and "mismatch_placement" in config.ranking_weights)):
        from dataclasses import replace as _replace
        weights = _replace(weights, mismatch_placement=0.0)
    return weights


def _stage_rank(
    specific: List[ProbeCandidate],
    config: RunConfig,
    gene_seqs: Dict[str, str],
    *,
    gene_aware: bool,
    _emit,
    _check_cancel,
) -> List[ProbeCandidate]:
    """Stage 6: score and order the specificity survivors, gene by gene.

    Returns the flat ranked candidate list (``final_ranked``) — every gene's
    ranked block concatenated in gene-encounter order, which is what the
    results table, the chemistry stage and the run summary all consume.

    ``ranked_by_gene`` and the assembled ``RankingWeights`` stay local: nothing
    downstream reads them, and the weights in particular are per-run derived
    state that must not be mistaken for the module-level ``DEFAULT_WEIGHTS``.

    ``_emit`` / ``_check_cancel`` are ``run_pipeline``'s own closures, passed in
    rather than rebuilt: ``_emit`` drives the ``_StageTimer``, so both WHICH
    calls happen and WHERE they sit in this body are load-bearing.
    """
    _emit("rank", "Ranking candidates", 70)
    _check_cancel()
    ranked_by_gene: Dict[str, List[ProbeCandidate]] = {}
    for c in specific:
        ranked_by_gene.setdefault(c.gene_id, []).append(c)

    weights = effective_ranking_weights(config, gene_aware=gene_aware)

    final_ranked: List[ProbeCandidate] = []
    for gid, cands in ranked_by_gene.items():
        seq_len = len(gene_seqs.get(gid, ''))
        ranked = rank_candidates(
            cands, seq_length=seq_len,
            weights=weights,
            primer3_scoring=config.primer3_scoring,
        )
        final_ranked.extend(ranked)
    _emit("rank", f"Ranked {len(final_ranked)} candidates across {len(ranked_by_gene)} genes", 80)

    return final_ranked


@dataclasses.dataclass(frozen=True)
class _ChemistryResult:
    """What Stage 7 hands forward to export and to the manifest.

    Six values rather than one because they are consumed by three different
    readers: ``pairs`` by the order writers, ``multi_barcode`` /
    ``barcode_pairs`` by the IDT writer alone, and ``bc_labels`` /
    ``probe_bc_label`` / ``probe_bc_seq`` by the provenance block. Named
    fields, not a bare tuple — ``probe_bc_label`` and ``probe_bc_seq`` are both
    strings and transposing them would silently write the wrong manifest.
    """
    pairs: List[ProbePair]
    bc_labels: List[str]
    multi_barcode: bool
    probe_bc_label: str
    probe_bc_seq: str
    barcode_pairs: List[Tuple[str, str]]


def _stage_chemistry_and_pair(
    final_ranked: List[ProbeCandidate],
    config: RunConfig,
    *,
    _emit,
    _check_cancel,
) -> _ChemistryResult:
    """Stage 7: apply the chemistry handles and assemble the LHS/RHS pairs.

    MUTATES ``final_ranked`` in place: every paired candidate is annotated with
    its chemistry-derived fields (probe_rc, the target halves, the oligos, the
    Tms, tn_junction_ok) so results.csv can show them. That is deliberate and
    load-bearing — the caller keeps writing the SAME candidate objects — so the
    annotation is a side effect of this stage, not something the return value
    carries.

    ``chemistry`` itself and the ``_bc_table`` alias stay local: nothing
    downstream reads either.

    ``_emit`` / ``_check_cancel`` are ``run_pipeline``'s own closures, passed in
    rather than rebuilt: ``_emit`` drives the ``_StageTimer``, so both WHICH
    calls happen and WHERE they sit in this body are load-bearing.
    """
    _emit("chemistry", "Applying chemistry handles", 82)
    _check_cancel()
    chemistry = get_chemistry(config.modality)
    # ``config.probe_barcodes`` is a list. Single-element = either singleplex
    # or single-sample multiplex (one barcode embedded in the oligo). Multi-
    # element = Flex v1 multi-sample multiplex: we leave a {BC} placeholder
    # in the RHS string so each export expands one row per barcode.
    bc_labels = (
        list(config.probe_barcodes) if chemistry.has_probe_barcode else []
    )
    multi_barcode = len(bc_labels) > 1
    probe_bc_label = bc_labels[0] if bc_labels else ""
    probe_bc_seq = _resolve_probe_barcode_seq(chemistry, probe_bc_label)
    # (label, 8-mer) for every selected barcode — used to expand the one-shot
    # IDT export per barcode (the discovery path expands later in select_export).
    _bc_table = FLEX_V1_PROBE_BARCODES
    barcode_pairs = (
        [(lbl, _bc_table.get(lbl, lbl)) for lbl in bc_labels]
        if chemistry.has_probe_barcode else []
    )

    pairs: List[ProbePair] = []
    for c in final_ranked:
        # Skip tiles that aren't exactly the chemistry's expected length —
        # mostly relevant when a user picks a non-50 tile_size in the GUI.
        if len(c.sequence) != chemistry.tile_length:
            continue
        # A degenerate parent is never synthesised: it has no reverse
        # complement worth the name, no target halves and no oligos. Its
        # variants are the pair-forming, orderable probes and they are in this
        # same list. The parent stays in results.csv (it is what the reviewer
        # picks a SITE by) but it never becomes a ProbePair, so it can never
        # reach the IDT or Cell Ranger writers.
        if c.probe_role == PROBE_ROLE_CONSENSUS:
            continue
        ap = assemble_pair(
            c.sequence, chemistry,
            multiplex=config.multiplex,
            probe_barcode=probe_bc_label or None,
            leave_barcode_placeholder=multi_barcode,
        )
        # Annotate the source candidate with chemistry-derived fields so
        # results.csv can show them.
        c.probe_rc = ap.probe_rc
        c.target_lhs = ap.target_lhs
        c.target_rhs = ap.target_rhs
        c.lhs_oligo = ap.lhs_oligo
        c.rhs_oligo = ap.rhs_oligo
        c.tn_junction_ok = ap.tn_junction_ok
        c.tm_target_lhs = _safe_tm(ap.target_lhs)
        c.tm_target_rhs = _safe_tm(ap.target_rhs)
        c.tm_full_50 = _safe_tm(ap.probe_rc)

        # Both halves are copies of the scored candidate (see _pair_half), so
        # their provenance matches the candidate's own results.csv row field for
        # field. Only `sequence` differs, and only because a half is ordered as
        # its own oligo: LHS on one, RHS on the other. The two share a probe_id
        # — the _LHS/_RHS suffix is added exactly once, by the order writer.
        lhs_probe = _pair_half(c, ap.lhs_oligo)
        rhs_probe = _pair_half(c, ap.rhs_oligo)
        pairs.append(ProbePair(
            lhs=lhs_probe, rhs=rhs_probe,
            probe_barcode_label=probe_bc_label,
            probe_barcode_seq=probe_bc_seq,
            pair_score=c.composite_score,
        ))
    _emit("chemistry", f"Created {len(pairs)} probe pairs", 85)

    return _ChemistryResult(
        pairs=pairs,
        bc_labels=bc_labels,
        multi_barcode=multi_barcode,
        probe_bc_label=probe_bc_label,
        probe_bc_seq=probe_bc_seq,
        barcode_pairs=barcode_pairs,
    )


def _stage_export(
    final_ranked: List[ProbeCandidate],
    pairs: List[ProbePair],
    config: RunConfig,
    output_dir: str,
    *,
    barcode_pairs: List[Tuple[str, str]],
    multi_barcode: bool,
    _emit,
    _check_cancel,
) -> Tuple[Any, bool]:
    """Stage 8: write results.csv and, outside discovery mode, the order CSVs.

    Returns ``(export, overwrite)`` — the ``ExportSettings`` block and its
    overwrite flag. They are read straight off ``config`` here, but the
    provenance block that follows this stage reads them again for
    ``manifest.json``, and re-deriving them there would be a second source of
    truth for "is this run allowed to clobber existing files".

    ``_emit`` / ``_check_cancel`` are ``run_pipeline``'s own closures, passed in
    rather than rebuilt: ``_emit`` drives the ``_StageTimer``, so both WHICH
    calls happen and WHERE they sit in this body are load-bearing.
    """
    _emit("export", "Writing outputs", 88)
    _check_cancel()
    export = config.export
    overwrite = export.overwrite

    from .writers.results_table import write_results_table

    if export.write_results_table:
        rt_path = os.path.join(output_dir, 'results.csv')
        write_results_table(final_ranked, rt_path, overwrite=overwrite,
                            screen_engine=config.screen.engine)

    # Order-format CSVs are written immediately in one-shot mode. In
    # discovery mode they're deferred until the user reviews + selects
    # probes in the GUI (see flexprobe.select_export).
    if not config.discovery_only:
        if export.write_idt_opools:
            from .writers.idt_opools import write_idt_opools
            idt_path = os.path.join(output_dir, 'idt_opools.csv')
            write_idt_opools(
                pairs, idt_path, overwrite=overwrite,
                barcodes=barcode_pairs if multi_barcode else None,
            )

        if export.write_cellranger_probe_set:
            from .writers.cellranger_probe_set import write_cellranger_probe_set
            cr_path = os.path.join(output_dir, 'cellranger_probe_set.csv')
            # The FASTA basename is the closest flexiprobes concept to a 10x
            # panel identity / reference genome — this run designed directly
            # against that input, not against a curated whole-genome
            # reference (see cellranger_probe_set.py's header docstring).
            fasta_ident = (
                os.path.splitext(os.path.basename(config.fasta_path))[0]
                or "input"
            )
            write_cellranger_probe_set(
                pairs, cr_path,
                version=config.modality.version,
                spatial=config.modality.spatial,
                overwrite=overwrite,
                panel_name=f"flexiprobes custom probe set: {fasta_ident}",
                reference_genome=fasta_ident,
            )

    return export, overwrite


def _stage_provenance_and_summary(
    config: RunConfig,
    output_dir: str,
    *,
    pairs: List[ProbePair],
    final_ranked: List[ProbeCandidate],
    passing: List[ProbeCandidate],
    specific: List[ProbeCandidate],
    gene_seqs: Dict[str, str],
    gene_aware: bool,
    group_members: Optional[Dict[str, Set[str]]],
    alignment: Optional[StrainAlignment],
    bc_labels: List[str],
    probe_bc_label: str,
    probe_bc_seq: str,
    export: Any,
    overwrite: bool,
    stage_timer: "_StageTimer",
    tiles_total: int,
    filter_stats: Dict[str, int],
    gated_low_iso: int,
    screen_ran: bool,
    null_model_stats: Optional[Dict[str, Any]],
    _prescreen_capture: str,
    _subjects: Optional[SubjectWindows],
    _tally: Any,
    _emit,
) -> Dict[str, Any]:
    """Provenance + summary: manifest.json, run_config.json, cache cleanup,
    run_timings.json, and the dict ``run_pipeline`` returns.

    The parameter list is long because the block genuinely reads that much: the
    manifest is a deliberate DUMP of the whole run (see the comment on
    ``params`` below), and the stats dict reports every stage's counters
    side by side so a Methods table can be built from one file. Every argument
    is keyword-only past ``output_dir`` — a dozen of them are ``int``/``bool``
    and a positional call would transpose two of them without raising.

    Returns run_pipeline's own result dict verbatim (``pairs`` / ``manifest`` /
    ``stats``) rather than the parts, so the tail of run_pipeline is a single
    ``return``.

    ``_emit`` is ``run_pipeline``'s own closure, passed in rather than rebuilt:
    it drives the ``_StageTimer``, so both WHICH calls happen and WHERE they sit
    in this body are load-bearing — the ``done`` event in particular is emitted
    BEFORE the summary is built, and ``run_timings.json`` is written after it.
    """
    # params is the WHOLE RunConfig plus derived run context — deliberately a
    # dump, not a hand-curated allowlist. An allowlist silently omits whatever
    # it does not name (filter thresholds, specificity thresholds, ranking
    # weights, the screened databases), so "what GC range and off-target
    # identity produced these probes?" and "which database was this screened
    # against?" become unanswerable from the outputs. A dump also means a
    # config field added later is captured by DEFAULT instead of being dropped
    # until someone remembers to extend a list here.
    # Imported here, not with the export-stage writers above: the manifest is
    # built and written in THIS block, and nothing in the export stage calls it.
    from .writers.manifest import write_manifest as _write_manifest
    from .config_io import filled_from_defaults, portable_config_dict, save_config

    # Paths relative to the output directory the manifest lives in.
    params: Dict[str, Any] = portable_config_dict(config, output_dir)
    # What ranked the probes, all seven weights, beside what was asked
    # for (params['ranking_weights'], the overrides only).
    _effective_weights = dataclasses.asdict(
        effective_ranking_weights(config, gene_aware=gene_aware))
    params.update({
        'effective_ranking_weights': _effective_weights,
        # ---- derived run context (computed, not configured) ----
        # 'modality' is the structured ModalityDefinition (from the config);
        # 'modality_label' is the flat "flex_v1" string older readers expect.
        'modality_label': f"{config.modality.name}_{config.modality.version}",
        # The barcode actually applied (config.probe_barcodes is the request;
        # chemistry defaults can differ, e.g. Flex v2).
        'probe_barcode': probe_bc_label,
        'probe_barcode_seq': probe_bc_seq,
        # The barcodes ACTUALLY applied — chemistry can drop them entirely
        # (Flex v2 has no probe barcode), so this differs from the requested
        # params['probe_barcodes'] that came from the config. >1 means the RHS
        # oligo in results.csv carries the {BC} placeholder and the export step
        # expands one row per barcode at order time.
        'probe_barcodes_applied': bc_labels,
        # Recorded so the review dialog can draw to-scale visualizations.
        'gene_lengths': {gid: len(seq) for gid, seq in gene_seqs.items()},
        # Gene -> all its isoform (transcript) ids, for gene-aware runs. Lets
        # the per-gene coverage selection know each gene's complete isoform set
        # (so it can name uncovered isoforms) without re-reading the input
        # FASTA. Empty when not grouping by gene.
        # What ran, beside what was asked for (params['group_by_gene']):
        # grouping is switched off on input with no gene: field.
        'group_by_gene_applied': bool(gene_aware),
        'gene_isoforms': (
            {g: sorted(ts) for g, ts in group_members.items()}
            if gene_aware and group_members is not None else {}
        ),
        # Multi-strain mode metadata so the writers / GUI know to expose
        # per-strain columns and the review dialog can label things.
        'strain_ids': alignment.strain_ids if alignment is not None else [],
        'variant_positions': (
            len(alignment.variant_positions) if alignment is not None else 0
        ),
        # Actual variant positions for the review-dialog track ticks.
        # Capped so the manifest stays small on highly divergent inputs.
        'variant_position_list': (
            alignment.variant_positions[:10000] if alignment is not None else []
        ),
        # Legacy top-level alias; the authoritative value is params['filters'].
        'require_tn_junction': config.filters.require_tn_junction,
        # When screen_top_per_gene > 0, only the top-N candidates per gene (by
        # the BLAST-free prescreen_score: GC + homopolymer + placement) were
        # screened for off-targets — so specificity could not rescue a window
        # with slightly worse GC but far better specificity. Recorded so
        # downstream analysis knows BLAST saw a pre-filtered subset.
        'blast_saw_prefiltered_subset': bool(config.screen_top_per_gene),
        # Total on-disk size (bytes) of each off-target BLAST DB. A bigger DB
        # means more chances for a spurious near-match, so "specificity passed"
        # is only interpretable alongside the DB it was screened against
        # (the multiple-testing context of the screen).
        'blast_db_sizes_bytes': _blast_db_sizes(config.screen.db_paths),
        # How the off-target mismatch counts were made across the WHOLE run, so
        # the manifest answers "were these verdicts exact?" without the reader
        # having to aggregate a results.csv column themselves. `basis` is
        # 'exact' only when every hit was recounted against real subject bases;
        # one fallback anywhere makes it 'conservative', i.e. the counts are
        # lower bounds and the screen may have over-rejected. See
        # specificity.full_half_mismatches.
        'offtarget_recount': _offtarget_recount_summary(_subjects, config),
        # Which search engine produced the hits behind every verdict above.
        # Present in `params` too (it is a config field), but promoted here
        # because it changes how the NEIGHBOURING fields must be read:
        # `offtarget_recount` is 'exact' under the k-mer engine by construction
        # rather than by successful retrieval, `best_offtarget_bitscore` /
        # `_evalue` are undefined (NaN), and `offtarget_nearmiss_count` covers a
        # different mismatch band. A reader comparing two runs needs this before
        # comparing anything else.
        'screen_engine': config.screen.engine,
        # Whether the junction-duplex gate (`off_target_junction_duplex_bp`,
        # `params['specificity']['off_target_junction_duplex_bp']` above) was
        # actually able to run. Only `KmerScreen._junction_scan` produces a
        # `junction_duplex` record; `BlastScreen` has none, so under
        # `--screen-engine blast` (or `--skip-screen`) the gate never fires no
        # matter what threshold is configured, and `offtarget_clean_left_bp` /
        # `_right_bp` stay at their -1 "not measured" default (see
        # specificity._worst_junction_duplex / models.ProbeCandidate). Mirrors
        # `null_model.summarise`'s `measurable` field: a boolean the manifest
        # carries so a reader can tell "gate configured" from "gate ran"
        # without re-deriving it from `screen_engine` themselves.
        'junction_duplex_measured': (
            config.screen.engine == ENGINE_KMER
            and not config.screen.skip_screen
        ),
    })

    # Same scoping rule as db_sha256s, and the same reason: a --skip-screen
    # run screened nothing.
    compositions = _reference_compositions(
        config.screen.db_paths,
        skip_screen=config.screen.skip_screen,
    )
    manifest = RunManifest(
        inputs=_portable_inputs(_hash_input_fasta(config.fasta_path), output_dir),
        code_version=_code_version(),
        # blastn only when the BLAST engine actually screened. The k-mer engine
        # (the default) never invokes it, and a --skip-screen run invokes
        # nothing, so either would otherwise cite a binary that never ran.
        tool_versions=_collect_tool_versions(
            screen_used=(config.screen.engine == ENGINE_BLAST
                         and not config.screen.skip_screen)),
        # Which aligner built the consensus (see _consensus_provenance).
        # Empty on a single-species run.
        consensus=_consensus_provenance(alignment),
        # Snapshotted (not finalised) — the export stage is still running, since
        # writing this manifest is itself part of it.
        elapsed_seconds=stage_timer.snapshot(),
        params=params,
        # Only the DBs actually screened against: a --skip-screen run screened
        # nothing, and hashing a DB it never opened would be a fabricated claim.
        db_sha256s=_portable_keys(
            _reference_sha256s(config.screen, compositions), output_dir),
        reference_composition=_portable_keys(compositions, output_dir),
        seed=config.seed,
        config_filled_from_defaults=filled_from_defaults(config),
    )
    if export.write_manifest:
        man_path = os.path.join(output_dir, 'manifest.json')
        _write_manifest(manifest, man_path, overwrite=overwrite)

    # A machine-loadable copy of the config, next to the results it produced.
    # The CLI writes only design_summary.json and the HPC export tars only
    # output/, so without this an imported HPC archive would contain no
    # configuration (the bundle's config.json stays on the cluster). This one
    # rides along in the tarball and reloads verbatim via config_io.load_config, so a
    # run can be re-executed rather than merely read about.
    save_config(config, os.path.join(output_dir, 'run_config.json'),
                portable=True,
                extra={'effective_ranking_weights': _effective_weights})

    # The run finished and exported, so the resumable checkpoints are no longer
    # needed — drop both the per-chunk BLAST cache and the pre-screen candidate
    # checkpoint. (They survive any earlier crash/timeout so a resubmit can pick
    # up; a clean finish leaves nothing behind, so they never reach the tarball.)
    if config.resume_screen:
        from . import screen_cache, candidate_cache
        screen_cache.clear(output_dir)
        if config.keep_candidate_cache:
            # Belt-and-braces for the pre-screen artefacts: a second, independent
            # copy of the same candidate set in a different format, written by
            # different code. If the artefact writer has a bug, this is the
            # fallback — and a whole-transcriptome rerun to recover it costs
            # days. Still dot-prefixed, so `tar -C output .` picks it up but the
            # GUI's results view never shows it. The BLAST cache is NOT kept:
            # it is orders of magnitude larger and carries nothing the results
            # don't already have.
            logger.info("Keeping %s (keep_candidate_cache).",
                        candidate_cache.cache_dir(output_dir))
        else:
            candidate_cache.clear(output_dir)

    _emit("done", "Pipeline complete", 100, {
        "genes": len(gene_seqs),
        "tiles_total": tiles_total,
        "filter_passed": len(passing),
        "specificity_passed": len(specific),
        "pairs_output": len(pairs),
    })

    # Breakdown of the specificity verdict across all output candidates so
    # the GUI/home stats can warn "these probes were never off-target
    # screened" rather than implying a clean pass.
    _emit("summary", "Building run summary", 99)
    spec_status_counts: Dict[str, int] = {}
    for c in final_ranked:
        spec_status_counts[c.specificity_status] = (
            spec_status_counts.get(c.specificity_status, 0) + 1
        )

    stats: Dict[str, Any] = {
        "genes": len(gene_seqs),
        "tiles_total": tiles_total,
        # NOT the same quantity as filter_report.json's counts.passed (or
        # filter_stats["passed"] a few lines below, embedded verbatim in this
        # same dict for exactly this comparison). `passing` here has been
        # through _reconcile_degenerate_families TWICE since filter_stats was
        # tallied (once right after the filter funnel, once after the
        # Stage 3b pre-screen cap — see the comment at the first call), each
        # time dropping any consensus parent/variant family that lost a
        # member independently of its siblings. filter_stats["passed"] counts
        # funnel survivors row-by-row; this counts what is left of them once
        # every family is whole-or-absent. The gap between the two is real,
        # legitimate, and exactly the dropped-family count — not a bug, and
        # not double-counted (see filter_stats["passed"] embedded below for
        # the pre-reconciliation number this can be compared against).
        "filter_passed": len(passing),
        "specificity_passed": len(specific),
        "specificity_screened": screen_ran,
        "specificity_status_counts": spec_status_counts,
        "pairs_output": len(pairs),
        # The PRE-reconciliation tally (see "filter_passed" above) — kept
        # verbatim so filter_stats["passed"] here and "filter_passed" above
        # can be diffed from this one file without also opening
        # filter_report.json.
        "filter_stats": filter_stats,
        "discovery_only": config.discovery_only,
        "auto_select_per_gene": config.auto_select_per_gene,
        "multi_strain": config.multi_strain,
        # Gene-aware / isoform-coverage telemetry (None/0 when off).
        "group_by_gene": gene_aware,
        "min_isoform_coverage": config.min_isoform_coverage if gene_aware else 0.0,
        "gated_low_isoform_coverage": gated_low_iso,
        "mean_isoform_coverage": (
            round(sum(c.isoform_coverage for c in final_ranked)
                  / len(final_ranked), 4)
            if (gene_aware and final_ranked) else None
        ),
        # Which pre-screen candidate set this run preserved: "none",
        # "post_cap" (what BLAST screened) or "pre_cap" (every
        # filter-passing window). It determines what a null model can
        # legitimately be built over, so it belongs in the run's own summary
        # and not only in prescreen_index.json alongside the files.
        "prescreen_capture": _prescreen_capture,
        # Null-model calibration of the off-target gate: the fraction of
        # composition-matched shuffled windows the same gate rejects. None
        # when the calibration was not requested or could not run. Lives in
        # the run's own stats, not only in null_model.json, so a Methods
        # table can be built from the summary alone.
        "null_model": null_model_stats,
        "prescreen_precap_rows": (
            _tally.rows_written if (_tally is not None
                                    and config.emit_prescreen_precap_rows)
            else None
        ),
        # The review dialog uses this to label the candidate-count column.
        "candidates_per_gene": _candidates_per_gene(pairs, gene_seqs),
    }

    # Written LAST — after the summary build — so unlike the manifest's
    # elapsed_seconds it has no blind spot at the tail. run_timeline.jsonl next
    # to it carries the same run step by step, flushed as it went.
    stage_timer.write_report(os.path.join(output_dir, "run_timings.json"))

    return {
        "pairs": pairs,
        "manifest": manifest,
        "stats": stats,
    }


class ProgressEvent:
    """Structured progress event emitted by the pipeline.

    ``pct`` is None for messages that carry no progress information — a warning
    raised before a stage commits, say. That is a normal event, not a broken
    one: consumers must render it without a percentage rather than assume a
    float. Formatting None as a float raises and would kill the run.
    """
    def __init__(self, stage: str, message: str, pct: Optional[float] = 0.0,
                 data: Optional[Dict] = None):
        self.stage = stage
        self.message = message
        self.pct = pct
        self.data = data or {}
        self.timestamp = time.time()

    def __repr__(self) -> str:
        pct = "n/a" if self.pct is None else f"{self.pct:.0f}%"
        return f"ProgressEvent({self.stage!r}, {self.message!r}, pct={pct})"


@dataclasses.dataclass
class _ScreenResult:
    """What Stage 4 (the off-target screen) leaves for Stages 5, 5b and the
    manifest.

    Eight values because the screen decides three separate things that later
    stages each need a different part of: WHETHER a screen happened
    (``screen_ran`` / ``blast_skip_status``), whether it already scored
    specificity as it went (``streamed`` / ``spec_passed``), and the scoring
    context it built (``spec_kwargs`` / ``half_len`` / ``subjects``) — which
    Stage 5, the null model and the provenance block all re-use, and which
    must be the SAME objects, not rebuilt, or a cached verdict would be keyed
    to one rule and scored under another.

    ``screen_hits_by_probe`` comes back because the fallback arms rebind it to
    an empty-hit dict: Stage 5 distinguishes "screened, no hits" from "never
    screened" by ``screen_ran``, but still needs a dict to batch-score.
    """

    screen_ran: bool
    blast_skip_status: str
    streamed: bool
    spec_passed: int
    spec_kwargs: Dict[str, Any]
    half_len: int
    subjects: Optional[SubjectWindows]
    screen_hits_by_probe: Optional[Dict[str, List[Dict[str, Any]]]]


def _stage_screen(
    config: RunConfig,
    _screen: List[ProbeCandidate],
    output_dir: str,
    *,
    target_ids: Set[str],
    gene_of: Optional[Dict[str, str]],
    group_members: Optional[Dict[str, Set[str]]],
    consensus_strain_rule: StrainRule,
    screen_hits_by_probe: Optional[Dict[str, List[Dict[str, Any]]]],
    n_prescreen: int,
    _scr: Dict[str, str],
    _emit,
    _check_cancel,
) -> _ScreenResult:
    """Stage 4: run the off-target screen (or record why it did not run).

    Called from :func:`run_pipeline`. ``_screen`` is the concrete
    variants only — a degenerate consensus parent is never queried and inherits
    its variants' verdicts in Stage 5 — while ``n_prescreen`` is the size of the
    FULL post-cap set, which is what the recount cost warning is scaled from
    (two half-queries per candidate). Keeping those two counts distinct is
    load-bearing: swapping them would understate the cost warning on any run
    with degenerate parents.

    ``_emit`` / ``_check_cancel`` are ``run_pipeline``'s own closures, passed in
    rather than rebuilt: ``_emit`` drives the :class:`_StageTimer`, so both
    WHICH events fire and WHERE they sit in this body are load-bearing.

    The candidate objects in ``_screen`` are MUTATED in place by the streamed
    scoring sink (specificity verdicts, mismatch counts, recount basis). That is
    the point — the caller keeps holding the same objects — so the scoring is a
    side effect of this stage, not something the return value carries.
    """
    _scr_stage = _scr["stage"]
    _check_cancel()
    # ``screen_ran`` distinguishes "we actually screened these probes" from
    # "we couldn't, so they're unscreened". It gates the specificity fallback
    # so a genuine off-target rejection is never silently overwritten with a
    # pass.
    screen_ran = False
    blast_skip_status = "not_screened"   # status to stamp on unscreened probes
    # What actually gets BLASTed. A consensus parent is a mixed-base string:
    # blastn has no sensible interpretation of it, and screening one
    # "representative" resolution would be worse than not screening at all,
    # because a variant specific to strain B can hit a human transcript where
    # strain A's variant does not — which is exactly the risk the screen exists
    # to catch. So the concrete variants are screened SEPARATELY and their
    # verdicts roll up onto the parent afterwards (all children must pass).
    # Specificity parameters shared by the streamed per-chunk path (when we run
    # BLAST ourselves) and the batch path (caller-provided cached hits), so both
    # score identically.
    _half_len = max(1, config.tile_length // 2)
    _spec_kwargs: Dict[str, Any] = dict(
        strain_rule=consensus_strain_rule,
        off_target_pident_min=config.specificity.off_target_pident_min,
        off_target_aln_min=config.specificity.off_target_aln_min,
        on_target_pident_min=config.specificity.on_target_pident_min,
        on_target_aln_min=config.specificity.on_target_aln_min,
        # 10x off-target rule as written; see SpecificityThresholds. These are
        # part of _spec_kwargs (not just passed at the call site) so they also
        # key the verdict cache — a verdict is only valid for the rule that
        # produced it.
        off_target_full_half=config.specificity.off_target_full_half,
        off_target_min_mismatch=config.specificity.off_target_min_mismatch,
        off_target_half_len=_half_len,
        off_target_require_both_halves=(
            config.specificity.off_target_require_both_halves
        ),
        off_target_same_subject=config.specificity.off_target_same_subject,
        # OURS, not 10x's — and in _spec_kwargs for the same reason as the
        # switches above: it changes verdicts, so a cached verdict is only valid
        # for the threshold that produced it.
        off_target_junction_duplex_bp=(
            config.specificity.off_target_junction_duplex_bp
        ),
        # --junction-window. It shapes which off-target is
        # profiled as closest, so it keys the verdict cache like the rest.
        junction_window=config.filters.junction_window,
        extra_on_target_ids=set(config.extra_on_target_ids) or None,
        gene_of=gene_of,
        group_members=group_members,
    )
    # Subject-window retrieval for the full-half mismatch recount. Built once
    # per run and shared by every chunk, because its whole value is the cache:
    # a screen hits the same abundant subject thousands of times and one
    # blastdbcmd process per hit would cost more than the BLAST search.
    # Deliberately NOT part of _spec_kwargs — that dict keys the resume cache,
    # and an object's repr changes every process, which would make the
    # fingerprint miss on every restart. Its scoring-relevant *state* is folded
    # in as a stable string below instead.
    _subjects: Optional[SubjectWindows] = None
    streamed = False        # True once we've scored specificity inside the BLAST loop
    spec_passed = 0
    if screen_hits_by_probe is not None:
        # Pre-computed (cached) hits provided by the caller.
        _emit(_scr_stage, "Using pre-computed hits", 50)
        screen_ran = True
    elif not config.screen.skip_screen and config.screen.db_paths:
        n_dbs = len(config.screen.db_paths)
        _emit(_scr_stage, _scr["search"].format(n=n_dbs), 45)
        # Surface the exact screen's cost BEFORE committing to it, not at 90%
        # of walltime. Two halves are queried per candidate.
        _warn_recount_cost(config.specificity, 2 * n_prescreen, _emit,
                           stage=_scr_stage, noun=_scr["noun"])

        def _spec_sink(ci, chunk_cands, chunk_hits):
            # Score this chunk's probes the instant its hits are ready, then drop
            # the hits — so the whole-transcriptome screen never accumulates
            # every probe's hits in RAM.
            # chunk_hits is None for a chunk resumed from the verdict cache: the
            # candidates already carry their scored fields, so only tally.
            nonlocal spec_passed
            if chunk_hits is None:
                spec_passed += sum(1 for c in chunk_cands if c.specificity_pass)
                return
            if _subjects is not None and _subjects.available:
                # Warm the whole chunk's subjects in ONE blastdbcmd call before
                # scoring. Without this the recount pays a process spawn per
                # distinct subject, which on a human-vs-human screen dominates
                # the runtime (see blast.SubjectWindows.prefetch).
                _subjects.prefetch(
                    h.get("sseqid", "")
                    for hits in chunk_hits.values() for h in hits
                )
            spec_passed += batch_evaluate(
                chunk_cands, chunk_hits, target_ids,
                subjects=_subjects, **_spec_kwargs,
            )

        try:
            # Resolve (and, for a raw FASTA, build) the DBs BEFORE the search so
            # the recount reads exactly the prefixes blastn will search.
            #
            # Failures here are swallowed on purpose. Retrieval is an ACCURACY
            # aid, not a prerequisite: a missing DB, a missing blastn, anything
            # — the recount degrades to the conservative count and the search
            # itself reports the real error a few lines below, with its own
            # message and its own blast_failed handling. Letting this raise would
            # convert "counts are bounded" into "the run died", and would report
            # the failure from the wrong place.
            #
            # Skipped entirely unless the ENGINE needs it. Subject retrieval
            # exists to reconstruct the probe positions a LOCAL alignment never
            # covered; an end-to-end engine compares all 25 of them itself, so
            # building a SubjectWindows for it would resolve databases, shell
            # out to makeblastdb, and hard-fail on a DB lacking -parse_seqids —
            # reintroducing every dependency that engine was chosen to remove,
            # in service of a recount that is never consulted.
            _needs_subjects = config.screen.engine == ENGINE_BLAST
            _db_prefixes: List[str] = []
            if _needs_subjects:
                try:
                    _db_prefixes = _resolve_db_prefixes(
                        config.screen.db_paths,
                        blast.discover_blast().get("makeblastdb"),
                    )
                except (BlastError, OSError):
                    logger.debug("Could not resolve DBs for subject retrieval",
                                 exc_info=True)
                    _db_prefixes = []
            if _needs_subjects:
                _subjects = SubjectWindows(_db_prefixes)
                # Ask each DB, once, whether it can serve -entry at all — BEFORE
                # a single hit is scored. Exact mode is not merely degraded
                # without retrieval, it is guaranteed to return nothing: its
                # pushdown floor invites every chance short exact core, and the
                # conservative fallback then counts the unaligned positions as
                # matches, so noise reads as a 0-mismatch off-target and every
                # probe is rejected, while the same candidates pass under
                # bounded mode.
                #
                # So this is a hard error, not a warning. A run that completes
                # and writes an empty panel is worse than one that refuses: the
                # first looks like a biological result.
                # Gate on there BEING a resolvable database. An unresolvable
                # path is a different failure — the search itself reports it,
                # with its own message — and stubbed/cached-hit callers never
                # have one at all. What must not pass silently is a database
                # that exists and cannot answer -entry, because that is the one
                # that produces an empty panel and looks like biology.
                if (config.specificity.off_target_full_half
                        and _db_prefixes and not _subjects.verify()):
                    _how = (
                        "rebuild the database(s) with accessions:\n"
                        "      makeblastdb -in <fasta> -dbtype nucl "
                        "-parse_seqids -out <prefix>"
                    )
                    raise PipelineError(
                        "None of the configured BLAST databases can serve "
                        "sequence retrieval — they were built without "
                        "-parse_seqids, so they carry no accession index.\n\n"
                        "The off-target screen reads the subject bases opposite "
                        "the unaligned part of each probe half through "
                        "`blastdbcmd -entry`. Without it there is nothing to "
                        "count against: every chance short core would score as a "
                        "0-mismatch off-target and the panel would come out "
                        "empty.\n"
                        f"Fix: {_how}"
                    )
            if (_subjects is not None and not _subjects.available
                    and config.specificity.off_target_full_half):
                logger.warning(
                    "%s, so off-target mismatch counts cannot be recounted "
                    "against the subject sequence. Every count will be "
                    "CONSERVATIVE (unaligned probe positions counted as "
                    "matches), which over-rejects but never over-passes.",
                    ("blastdbcmd was not found (install ncbi-blast+; the "
                     "packaged desktop app and the container bundle it)"
                     if not _subjects.blastdbcmd_path
                     else "no BLAST DB could be resolved for retrieval"),
                )
            # Computed HERE, not before the DB capability check above: that
            # check can downgrade an impossible 'exact' to 'bounded', and the
            # server-side floor has to follow the mode that will actually run.
            # Computed first, the downgrade would be inert: the warning would
            # fire and the mode change, but BLAST would still be asked for every
            # chance short core, and the panel would still come out empty.
            _off_pident_pd, _off_aln_pd = _off_target_pushdown(
                config.specificity, _half_len,
            )
            _parallel_factory = _make_parallel_factory(
                config, target_ids, _spec_kwargs,
                needs_subjects=_needs_subjects, output_dir=output_dir,
                emit=lambda msg, pct: _emit(_scr_stage, msg, pct),
            )
            _run_screen(
                _screen, config.screen.db_paths, output_dir,
                # Lowest identity/length floor across on- and off-target
                # classification, so BLAST never drops a hit either path needs.
                min_pident=min(_off_pident_pd,
                               config.specificity.on_target_pident_min),
                min_aln=min(_off_aln_pd,
                            config.specificity.on_target_aln_min),
                query_half_len=_half_len,
                resume=config.resume_screen,
                chunk_size=config.screen.query_chunk,
                max_target_seqs=config.screen.max_target_seqs,
                chunk_label=_scr["chunk"],
                emit=lambda msg, pct: _emit(_scr_stage, msg, pct),
                cancel=_check_cancel,
                per_chunk=_spec_sink,
                # Keys the verdict cache to this exact scoring config — target
                # set included, since gene-aware on-target membership changes
                # every verdict.
                # ⚠️ If `require_on_target` or `required_strains` are ever
                # plumbed through from RunConfig, they MUST be added here too:
                # both change a verdict, and neither is currently in _spec_kwargs
                # (they're pinned at False/None, which is the only reason
                # their absence from the fingerprint is safe).
                # `offtarget_recount` is a STABLE description of whether the
                # exact recount was possible at all. It belongs in the digest
                # because it changes verdicts (exact vs conservative counts), but
                # the SubjectWindows object itself must not — its repr differs
                # every process, which would defeat resume entirely.
                spec_params={
                    **_spec_kwargs,
                    "target_ids": target_ids,
                    # `_subjects` is None in bounded mode, which retrieves
                    # nothing by design — the same cache-key meaning as a run
                    # whose DB cannot serve retrieval, since both produce
                    # conservative counts.
                    "offtarget_recount": (
                        "blastdbcmd"
                        if (_subjects is not None and _subjects.available)
                        else "conservative"
                    ),
                },
                cache_mode=config.screen_cache_mode,
                engine=config.screen.engine,
                exclude_non_primary=config.screen.exclude_alt_scaffolds,
                reference_strand=config.screen.reference_strand,
                # Main screen only. The null model screens after it and finds
                # the index this call rebuilt, so it never needs the flag.
                rebuild_kmer_index=config.screen.rebuild_kmer_index,
                parallel_factory=_parallel_factory,
            )
            screen_ran = True
            streamed = True
            _emit(_scr_stage,
                  f"Screened {len(_screen)} probes (streamed specificity)", 64)
        except BlastError as exc:
            logger.error("BLAST failed: %s", exc)
            _emit(_scr_stage, f"{_scr['noun']} failed ({exc}); probes left "
                          f"unscreened", 50)
            screen_hits_by_probe = {c.probe_id: [] for c in _screen}
            blast_skip_status = "blast_failed"
    elif not config.screen.skip_screen:
        # No DB configured and not explicitly skipped — warn loudly so users
        # don't silently accept un-screened probes.
        _emit(_scr_stage, f"No {_scr['target']} configured; skipping "
                          f"(specificity ungated)", 50)
        logger.warning(
            "No off-target reference configured (config.screen.db_paths is "
            "empty), so the screen is skipped and probes are not "
            "specificity-gated. Set ScreenResources.db_paths or pass "
            "--off-target-ref to enable real off-target screening."
        )
        screen_hits_by_probe = {c.probe_id: [] for c in _screen}
        blast_skip_status = "not_screened"
    else:
        _emit(_scr_stage, f"{_scr['noun']} skipped per config", 50)
        screen_hits_by_probe = {c.probe_id: [] for c in _screen}
        blast_skip_status = "not_screened"

    return _ScreenResult(
        screen_ran=screen_ran,
        blast_skip_status=blast_skip_status,
        streamed=streamed,
        spec_passed=spec_passed,
        spec_kwargs=_spec_kwargs,
        half_len=_half_len,
        subjects=_subjects,
        screen_hits_by_probe=screen_hits_by_probe,
    )


@dataclasses.dataclass
class _PrescreenPlan:
    """Everything the pre-screen cap and the pre-screen artefacts need, decided
    once, before a single tile is generated.

    These ten values are read by four later blocks (the streamed filter, the
    batch cap, the artefact writer and the manifest). They belong together because they answer ONE question — what, if anything,
    narrows and records the candidate set before the screen — and because
    several are only meaningful in combination: ``gene_key`` and ``weights``
    are None unless ``cap_active``, and ``tally`` exists only when the pre-cap
    counts are still observable (i.e. not on a resumed link).

    ``seq_len_of`` closes over ``gene_seqs``. On a multi-strain run that dict is
    REPLACED by the consensus in Stage 2, which is why the plan is built after
    tiling is set up and not before it.
    """

    cap_n: int
    cap_active: bool
    stream_cap: bool
    weights: Optional[Any]
    gene_key: Optional[Callable[[ProbeCandidate], str]]
    seq_len_of: Callable[[str], int]
    emit_prescreen: bool
    capture: str
    prescreen_gene_key: Callable[[ProbeCandidate], str]
    tally: Optional[Any]


def _plan_prescreen(
    config: RunConfig,
    gene_seqs: Dict[str, str],
    output_dir: str,
    *,
    gene_aware: bool,
    gene_of: Optional[Dict[str, str]],
    screen_hits_by_probe: Optional[Dict[str, List[Dict[str, Any]]]],
    cand_loaded: bool,
) -> _PrescreenPlan:
    """Decide the pre-screen cap and the pre-screen artefact capture.

    Called from :func:`run_pipeline`. Emits no progress event and
    touches no candidate: it only reads the config and the (possibly
    consensus-replaced) gene set. It MUST be called before the first tile is
    generated, because the :class:`PrecapTally` it may build has to see every
    candidate the filter passes, including the first.
    """
    # ---- Pre-screen cap setup (used to bound the filter-stage memory) ----
    # The per-gene cap keeps only the best N windows per gene before BLAST. It is
    # meaningful only when BLAST will actually run. On a normal (non-multi-strain)
    # run we apply it *while streaming the filter* — passing candidates flow into
    # a bounded per-gene accumulator — so the full pre-cap set is never held in
    # RAM. That is the filter-stage memory ceiling for a whole-transcriptome run.
    will_blast = (
        screen_hits_by_probe is None
        and not config.screen.skip_screen
        and bool(config.screen.db_paths)
    )
    _cap_n = config.screen_top_per_gene
    _cap_active = _cap_n > 0 and will_blast
    cap_weights = None
    _gene_key: Optional[Callable[[ProbeCandidate], str]] = None
    _seq_len_of = lambda gid: len(gene_seqs.get(gid, ''))      # noqa: E731
    if _cap_active:
        from .ranking import DEFAULT_WEIGHTS, RankingWeights
        if config.ranking_weights:
            from dataclasses import asdict as _asdict
            cap_weights = RankingWeights(
                **{**_asdict(DEFAULT_WEIGHTS), **config.ranking_weights})
        else:
            cap_weights = DEFAULT_WEIGHTS
        # Group by biological gene (Ensembl gene for isoform-aware runs, else the
        # record id). seq length for central placement comes from the candidate's
        # own transcript.
        if gene_aware and gene_of is not None:
            _gene_key = lambda c: gene_of.get(c.gene_id, c.gene_id)  # noqa: E731
        else:
            _gene_key = lambda c: c.gene_id                          # noqa: E731
    # Stream the cap on the normal path only; multi-strain carries heavy
    # per-strain candidate state and is not a whole-transcriptome memory case, so
    # it keeps the batch cap in Stage 3b below.
    _stream_cap = _cap_active and not config.multi_strain

    # ---- Pre-screen artefact setup ----
    # The pre-cap candidate count per gene can only be observed HERE. The cap
    # runs inside the filter stream via a bounded accumulator, so the full
    # pre-cap set never exists in one place; if we don't tally it as it goes by,
    # it is gone. (The post-cap set is materialised and can be written later.)
    _emit_prescreen = bool(config.emit_prescreen_candidates)
    # Named once, here, and carried through to prescreen_index.json and the run
    # stats. Whether a run captured the pre-cap or only the post-cap set changes
    # what can be concluded from it, so it is stated rather than inferred from
    # which files happen to be present.
    from .prescreen_export import (
        CAPTURE_NONE, CAPTURE_POSTCAP, CAPTURE_PRECAP,
    )
    _prescreen_capture = (
        CAPTURE_NONE if not _emit_prescreen
        else CAPTURE_PRECAP if config.emit_prescreen_precap_rows
        else CAPTURE_POSTCAP
    )
    _prescreen_gene_key: Callable[[ProbeCandidate], str]
    if gene_aware and gene_of is not None:
        _prescreen_gene_key = lambda c: gene_of.get(c.gene_id, c.gene_id)  # noqa: E731
    else:
        _prescreen_gene_key = lambda c: c.gene_id                          # noqa: E731
    _tally = None
    if _emit_prescreen and not cand_loaded:
        from .prescreen_export import PrecapTally
        _tally = PrecapTally(
            _prescreen_gene_key,
            output_dir=output_dir,
            seq_length_of=_seq_len_of,
            weights=cap_weights,
            write_rows=bool(config.emit_prescreen_precap_rows),
        )

    return _PrescreenPlan(
        cap_n=_cap_n,
        cap_active=_cap_active,
        stream_cap=_stream_cap,
        weights=cap_weights,
        gene_key=_gene_key,
        seq_len_of=_seq_len_of,
        emit_prescreen=_emit_prescreen,
        capture=_prescreen_capture,
        prescreen_gene_key=_prescreen_gene_key,
        tally=_tally,
    )


class _FilterResult(NamedTuple):
    """What Stage 3 hands back to :func:`run_pipeline`.

    ``reject_samples`` is deliberately absent: the capped per-reason examples
    exist only to be written into filter_report.json, which happens inside the
    stage, and nothing downstream reads them.

    ``stream_capped`` records whether the per-gene cap ALREADY ran inside the
    filter stream. Stage 3b's batch cap is the fallback for the paths that
    could not stream it (chiefly multi-strain), and running both would cap a
    capped set — so this flag is what keeps the cap applied exactly once.
    """

    passing: List[ProbeCandidate]
    filter_stats: Dict[str, int]
    tiles_total: int
    stream_capped: bool


def _stage_filter(
    config: RunConfig,
    gene_seqs: Dict[str, str],
    candidate_iter: Iterator[ProbeCandidate],
    output_dir: str,
    *,
    plan: _PrescreenPlan,
    on_target_records: Optional[FrozenSet[str]],
    _scr_stage: str,
    _emit,
    _check_cancel,
) -> _FilterResult:
    """Stage 3: consume the lazy tile stream and run the filter funnel.

    Called from :func:`run_pipeline`'s ``if not cand_loaded:``
    branch; the caller keeps that guard, because a resumed link reloads
    ``passing`` from the checkpoint instead and this stage must not run at all.

    Three paths converge on the same result — a process pool feeding a bounded
    per-gene accumulator, a process pool merged in gene order, and the serial
    stream — and all three are here rather than spread across the caller, which
    is the point: they must produce identical output, and that is only checkable
    if they are readable side by side.

    ``candidate_iter`` is consumed. It is the lazy generator built in Stage 2
    and is only forced here, so the whole-transcriptome memory ceiling is set by
    what this function chooses to retain, not by the tile count.

    ``_emit`` / ``_check_cancel`` are ``run_pipeline``'s own closures: ``_emit``
    drives the :class:`_StageTimer` and ``_check_cancel`` is what makes a
    multi-hour filter interruptible, so both the calls and their positions in
    this body are load-bearing.
    """
    # Bound as plain locals for readability; see _PrescreenPlan for what each
    # one means.
    _stream_cap = plan.stream_cap
    _cap_n = plan.cap_n
    _gene_key = plan.gene_key
    cap_weights = plan.weights
    _seq_len_of = plan.seq_len_of
    _tally = plan.tally
    _stream_capped = False
    _emit("filter", "Tiling + filtering", 30)
    _check_cancel()
    # Parallel path: multi-gene, non-multi-strain runs can fan genes out
    # across worker processes (the tile+filter work is independent per gene
    # and CPU-bound, so processes side-step the GIL). Results are merged in a
    # fixed gene order so the output is identical to the serial path. Default
    # workers=1 uses the serial streaming path.
    _use_parallel = (
        config.workers > 1 and not config.multi_strain and len(gene_seqs) > 1
    )
    reject_cap = 25
    _acc = None
    if _stream_cap:
        from .ranking import TopPerGeneAccumulator
        # _stream_cap implies _cap_active, which is the branch that assigns
        # both of these; they are only Optional because the no-cap path
        # never builds them.
        assert _gene_key is not None and cap_weights is not None
        _acc = TopPerGeneAccumulator(
            gene_key=_gene_key, seq_length_of=_seq_len_of,
            n=_cap_n, weights=cap_weights,
        )
    if _use_parallel:
        from concurrent.futures import ProcessPoolExecutor
        from .filters import new_filter_counters
        tasks = [
            (gid, gseq, config.tile_length, config.tile_step, config.filters, reject_cap)
            for gid, gseq in gene_seqs.items()
        ]
        if _acc is not None:
            # Feed each gene's passing into the bounded accumulator as futures
            # complete — never concatenate the full passing set.
            filter_stats = new_filter_counters()
            reject_samples = []
            _per_reason: Dict[str, int] = {}
            with ProcessPoolExecutor(max_workers=config.workers) as ex:
                for _gid, p, c, rj in ex.map(_tile_filter_one_gene, tasks):
                    _check_cancel()
                    if _tally is not None:
                        # Before the accumulator, which is where the cap
                        # happens — this is the only sight of the pre-cap set.
                        _tally.add_many(p)
                    _acc.add_many(p)
                    for k, v in c.items():
                        filter_stats[k] = filter_stats.get(k, 0) + v
                    for pc in rj:
                        keep = False
                        for r in pc.filter_fail_reasons:
                            if _per_reason.get(r, 0) < reject_cap:
                                _per_reason[r] = _per_reason.get(r, 0) + 1
                                keep = True
                        if keep:
                            reject_samples.append(pc)
        else:
            ordered: Dict[str, Any] = {}
            with ProcessPoolExecutor(max_workers=config.workers) as ex:
                for gid, p, c, rj in ex.map(_tile_filter_one_gene, tasks):
                    _check_cancel()
                    ordered[gid] = (gid, p, c, rj)
            # Re-assemble in the original gene order for deterministic output.
            results = [ordered[gid] for gid in gene_seqs if gid in ordered]
            passing, filter_stats, reject_samples = _merge_filter_results(
                results, reject_cap)
    else:
        # Declared up front: the three branches below produce a closure, a
        # bound method and None, and without the annotation the first one
        # to be seen fixes the variable's type for the other two.
        _sink: Optional[Callable[[ProbeCandidate], None]]
        if _acc is not None and _tally is not None:
            def _sink_tally_then_cap(c, _a=_acc, _t=_tally) -> None:
                _t.add(c)      # pre-cap: counted (and optionally written)
                _a.add(c)      # then the cap's bounded accumulator
            _sink = _sink_tally_then_cap
        elif _acc is not None:
            _sink = _acc.add
        else:
            # No cap: _stream_filter materialises `passing`, so the pre-cap
            # set is the post-cap set and is tallied from it below.
            _sink = None
        passing, filter_stats, reject_samples = _stream_filter(
            candidate_iter, config.filters, reject_cap=reject_cap,
            cancel=_check_cancel,
            sink=_sink,
            on_target_records=on_target_records,
        )
    tiles_total = filter_stats["total"]
    _emit("tile", f"Generated {tiles_total} tiles", 25,
          {"tile_count": tiles_total})
    if _acc is not None:
        n_passed = filter_stats.get("passed", 0)
        passing = _acc.result()
        _stream_capped = True
        _emit("filter", f"{n_passed} candidates passed filters", 40,
              {"filter_stats": filter_stats})
        _emit(_scr_stage,
              f"Pre-screen cap: top {_cap_n}/gene kept {len(passing)} of "
              f"{n_passed} candidates for screening", 44)
    else:
        _emit("filter", f"{len(passing)} candidates passed filters", 40,
              {"filter_stats": filter_stats})

    if _tally is not None:
        if _acc is None:
            # No streaming cap, so `passing` IS the pre-cap set (Stage 3b's
            # batch cap, if any, has not run yet). Tally it here.
            _tally.add_many(passing)
        _tally.close()

    # A parent and its variants are filtered independently, so the funnel
    # can keep one and drop another. Put the families back to intact-or-
    # absent before anything downstream scores a parent on a partial pool.
    #
    # THIS is where filter_report.json's counts.passed and
    # design_summary.json's stats.filter_passed part ways, deliberately.
    # filter_stats["passed"] below was tallied by _stream_filter /
    # _merge_filter_results BEFORE this line and is never revised
    # afterward, so it counts every row that individually cleared the
    # filter funnel — a consensus parent and each of its variants counted
    # independently. stats["filter_passed"] (written into design_summary
    # much later, from len(passing) at that point) is measured AFTER this
    # reconciliation — and after the Stage 3b pre-screen cap's own
    # reconciliation below — so it excludes any family this call drops
    # whole for having lost a sibling (parent kept, a variant filtered
    # out, or vice versa). filter_stats["passed"] - stats["filter_passed"]
    # is exactly the number of rows dropped by family reconciliation
    # (across both calls), not a bug and not double-counted — see
    # _reconcile_degenerate_families's docstring. Both numbers are real
    # and both are worth citing; a report quoting one should say
    # which, since they answer different questions ("how many candidates
    # cleared the funnel" vs "how many make it into a scoreable family").
    passing = _reconcile_degenerate_families(passing)

    # Filter explorer: persist the rejection counts + capped examples.
    # filter_stats here is intentionally the PRE-reconciliation tally (see
    # the comment above) — this file describes the filter funnel itself,
    # not what survived the later family bookkeeping.
    try:
        _write_filter_report(
            output_dir, reject_samples, filter_stats, config, max_examples=25,
        )
    except OSError:
        logger.warning("Could not write filter_report.json", exc_info=True)

    return _FilterResult(
        passing=passing,
        filter_stats=filter_stats,
        tiles_total=tiles_total,
        stream_capped=_stream_capped,
    )


class _SpecificityResult(NamedTuple):
    """What Stage 5 hands forward.

    ``spec_passed`` is deliberately absent: it is reported in this stage's own
    progress event and read nowhere afterwards — the manifest counts
    ``len(specific)``, which is a different number whenever the isoform gate
    fired.
    """

    specific: List[ProbeCandidate]
    gated_low_iso: int


def _stage_specificity(
    config: RunConfig,
    prescreen_set: List[ProbeCandidate],
    _screen: List[ProbeCandidate],
    *,
    target_ids: Set[str],
    sc: _ScreenResult,
    n_parents: int,
    gene_aware: bool,
    _emit,
    _check_cancel,
) -> _SpecificityResult:
    """Stage 5: turn screen hits into per-probe verdicts, or stamp "unscreened".

    Called from :func:`run_pipeline`. The whole stage forks on
    ``sc.screen_ran``, and that fork is the safety property: if
    a screen genuinely ran, a rejected probe STAYS rejected and the run honestly
    reports zero pairs; only when no screen ran at all are probes passed forward,
    and then they are stamped with a status saying so rather than quietly reading
    as clean.

    MUTATES the candidates in ``prescreen_set`` on the unscreened path (status,
    margin, recount basis). That is deliberate — those fields are what results.csv
    reads — and the recount basis is stamped explicitly rather than left at the
    model default, because a probe reloaded from a candidate checkpoint could
    otherwise carry an earlier link's basis and read as recounted when it was not.

    ``_emit`` / ``_check_cancel`` are ``run_pipeline``'s own closures; ``_emit``
    drives the :class:`_StageTimer`.
    """
    # Bound as plain locals for readability.
    screen_ran = sc.screen_ran
    streamed = sc.streamed
    spec_passed = sc.spec_passed
    blast_skip_status = sc.blast_skip_status
    screen_hits_by_probe = sc.screen_hits_by_probe
    _spec_kwargs = sc.spec_kwargs
    _subjects = sc.subjects
    _n_parents = n_parents
    _check_cancel()
    gated_low_iso = 0
    if screen_ran:
        if not streamed:
            # Caller-provided cached hits: score the whole set in one batch (the
            # streamed path already scored chunk-by-chunk inside the BLAST loop).
            _emit("specificity", "Evaluating specificity", 55)
            # Every path that reaches Stage 5 without having streamed either was
            # handed hits by the caller or assigned an empty-hit dict above.
            assert screen_hits_by_probe is not None
            spec_passed = batch_evaluate(
                # `_screen`, not `passing`: a degenerate consensus parent is
                # never queried, so it has no hits of its own and inherits its
                # variants' verdicts in the rollup below. `subjects` is what
                # makes the count a real recount rather than HSP geometry.
                _screen, screen_hits_by_probe, target_ids,
                subjects=_subjects, **_spec_kwargs,
            )
        # ONE aggregated line per run, not one per hit: a whole-transcriptome
        # screen has tens of millions of hits and per-hit logging would itself be
        # the outage. Reported whenever ANY hit fell back, because the run's
        # verdicts are then bounded rather than exact — the same distinction
        # results.csv carries per probe in `offtarget_recount_basis`.
        if _subjects is not None and _subjects.n_fallback:
            _total = _subjects.n_exact + _subjects.n_fallback
            logger.warning(
                "%d of %d off-target hits (%.1f%%) could not be recounted "
                "against the subject sequence; their unaligned probe positions "
                "were counted as MATCHES. Those probes' mismatch counts are "
                "LOWER BOUNDS, so the screen may have over-rejected — it cannot "
                "have over-passed. Common causes: a BLAST DB built without "
                "-parse_seqids (rebuild it), or subjects absent from the DB. "
                "Per-probe provenance is in the results.csv "
                "`offtarget_recount_basis` column.",
                _subjects.n_fallback, _total,
                100.0 * _subjects.n_fallback / max(1, _total),
            )
        # Parents were never screened themselves; they inherit their variants'
        # verdicts now, and only pass if every one of them did.
        if _n_parents:
            _parents_passed = _rollup_degenerate_specificity(prescreen_set)
            _emit("specificity",
                  f"{_parents_passed} of {_n_parents} degenerate parent(s) "
                  f"passed (every variant clean)", 65,
                  {"degenerate_parents_passed": _parents_passed})
        specific = [c for c in prescreen_set if c.specificity_pass]
        _emit("specificity", f"{spec_passed} probes passed specificity", 65,
              {"specificity_passed": spec_passed})
        # Optional isoform-coverage gate: drop probes that bind too few of their
        # gene's isoforms. Soft scoring already down-ranks them; this is the
        # hard cut for users who don't want low-coverage probes at all.
        if gene_aware and config.min_isoform_coverage > 0.0:
            before_iso_gate = len(specific)
            specific = [
                c for c in specific
                if c.isoform_coverage >= config.min_isoform_coverage
            ]
            gated_low_iso = before_iso_gate - len(specific)
            if gated_low_iso:
                _emit("specificity",
                      f"Gated {gated_low_iso} probe(s) below "
                      f"{config.min_isoform_coverage:.0%} isoform coverage", 66,
                      {"gated_low_isoform_coverage": gated_low_iso})
        # When BLAST genuinely ran, we do NOT resurrect rejected probes. An
        # empty ``specific`` here means "everything was off-target / had no
        # on-target hit" and the run honestly produces 0 pairs.
        if not specific:
            logger.warning(
                "BLAST ran but no probe passed specificity (all off-target or "
                "no on-target hit). Producing 0 pairs."
            )
    else:
        # BLAST did not run: probes are unscreened. Pass them all through so
        # the user still gets candidates, but stamp the honest status so the
        # CSV and stats make clear these were NOT off-target screened.
        for c in prescreen_set:
            c.specificity_pass = True
            c.specificity_status = blast_skip_status
            c.specificity_margin = 0.0
            # Stamped explicitly rather than left at the model default: a probe
            # reloaded from a candidate checkpoint could otherwise carry a
            # basis from an earlier link and read as recounted when it was not.
            c.offtarget_recount_basis = "not_screened"
        specific = prescreen_set
        spec_passed = 0
        _emit("specificity",
              f"BLAST not run — {len(prescreen_set)} probes passed forward UNSCREENED",
              65, {"specificity_passed": 0, "unscreened": len(prescreen_set)})

    return _SpecificityResult(specific=specific, gated_low_iso=gated_low_iso)


def _write_prescreen_artefacts(
    config: RunConfig,
    output_dir: str,
    _screen: List[ProbeCandidate],
    *,
    plan: _PrescreenPlan,
    _scr_stage: str,
    _emit,
) -> None:
    """Write the opt-in pre-screen candidate artefacts.

    Called from :func:`run_pipeline`. Must be called at exactly this
    point (after the cap, before the screen), because these files
    document themselves as the set the screen was given, which is the
    denominator any statement about the specificity screen needs. Called
    earlier they would describe a set that is about to be narrowed; called
    later, one that has already been scored.

    Best-effort but LOUD: a failure here must not lose a multi-day run, so it
    is caught, but it is logged at error level, because a silently missing
    artefact is indistinguishable from a run that was never asked for one, and
    the flag exists so that this set is kept rather than discarded.
    """
    _emit_prescreen = plan.emit_prescreen
    _prescreen_gene_key = plan.prescreen_gene_key
    _seq_len_of = plan.seq_len_of
    cap_weights = plan.weights
    _tally = plan.tally
    _prescreen_capture = plan.capture

    # Written HERE, after the cap and before BLAST, because `passing` at this
    # point is exactly the set BLAST will screen — the denominator any statement
    # about the specificity screen needs. Deliberately NOT derived from
    # .candidate_cache/: that is an internal, schema-versioned resume format
    # that is deleted on success and free to change shape.
    #
    # Best-effort like the checkpoint above: a full whole-transcriptome run must
    # not be lost to a failure writing a reporting artefact. But unlike the
    # checkpoint this one is LOUD, because a silently missing artefact is
    # indistinguishable from a run that was never asked for one, and the flag
    # exists so that this set is kept rather than discarded.
    if _emit_prescreen and _screen:
        from . import prescreen_export
        try:
            cpath = prescreen_export.write_candidates(
                output_dir, _screen,
                gene_key=_prescreen_gene_key, seq_length_of=_seq_len_of,
                weights=cap_weights,
            )
            _emit(_scr_stage,
                  f"Wrote pre-screen candidate set ({len(_screen)} rows) -> "
                  f"{os.path.basename(cpath)}", 44)
        except OSError:
            logger.error("Could not write the pre-screen candidate artefact "
                         "(%s) — the run continues, but the candidate set that "
                         "BLAST screened will NOT be recoverable from the "
                         "outputs.", prescreen_export.CANDIDATES_FILE,
                         exc_info=True)
        if _tally is not None:
            try:
                prescreen_export.write_cap_summary(
                    output_dir, _tally.counts, _screen,
                    gene_key=_prescreen_gene_key,
                )
                prescreen_export.write_index(
                    output_dir,
                    capture=_prescreen_capture,
                    tile_length=config.tile_length,
                    screen_top_per_gene=config.screen_top_per_gene,
                    screened_rows=len(_screen),
                    precap_rows=(_tally.rows_written
                                 if config.emit_prescreen_precap_rows else None),
                )
            except OSError:
                logger.error("Could not write the pre-screen cap summary (%s).",
                             prescreen_export.SUMMARY_FILE, exc_info=True)
        else:
            # Resumed from a checkpoint: the pre-cap counts were observed by the
            # link that did the filtering, and its summary file is already in
            # output/. Say so rather than leaving a gap someone has to guess at.
            logger.info(
                "Resumed from a candidate checkpoint, so the pre-cap per-gene "
                "counts were not re-observed; %s from the earlier link is the "
                "authoritative copy.", prescreen_export.SUMMARY_FILE)


def _prepare_config(config: RunConfig) -> None:
    """Validate the run configuration and resolve its chemistry-dependent gates.

    Called at the head of :func:`run_pipeline`. Runs before any
    output directory exists, so a bad configuration costs nothing.

    MUTATES ``config.filters.tn_junction_mode``: the "auto" TN-junction setting
    is resolved against the chemistry HERE, once, and the concrete mode is baked
    in. That is load-bearing rather than tidy — the filter funnel is deliberately
    chemistry-agnostic (it also runs inside worker processes and against the
    candidate cache, neither of which can see the modality), so "auto" must never
    reach it.

    Both warnings describe a gate that is configured but weaker than it looks,
    and both exist to make that LOUD rather than silent. Neither refuses the run.
    """
    # Fail fast on an invalid configuration with a single, clear message
    # rather than a confusing error deep in a stage (or silent nonsense).
    config.validate()

    # Resolve the chemistry-dependent TN-junction gate ONCE, here, and bake the
    # concrete mode into the thresholds the filter funnel sees. The funnel is
    # deliberately chemistry-agnostic (it also runs in worker processes and from
    # the candidate cache), so "auto" must never reach it.
    _chem = get_chemistry(config.modality)
    _tn_mode = resolve_tn_mode(config.filters, _chem)
    config.filters.tn_junction_mode = _tn_mode
    # Warn whenever the effective gate is weaker than what this chemistry's
    # tech note states as a requirement, so a relaxation is never silent.
    if _chem.requires_tn_junction and _tn_mode != TN_MODE_STRICT_A:
        logger.warning(
            "require_tn_junction is OFF (effective TN gate: %r) but %s %s "
            "expects the TN junction adenine (target position 25 == 'A'). "
            "Probes without it may ligate poorly. Pass --require-tn-junction "
            "to enforce it.",
            _tn_mode, config.modality.name, config.modality.version,
        )

    # The BLAST-under-junction-gate case is refused by RunConfig.validate,
    # called above, so no warning is needed here (it could never be reached).


def _raise_no_candidates(
    gene_seqs: Dict[str, str], tile_length: int,
) -> NoReturn:
    """Explain an empty candidate set and raise :class:`PipelineError`.

    Called from :func:`run_pipeline`. Three messages rather than one
    because "no candidates passed filters" is the right diagnosis only when
    there were candidates to filter: an input shorter than one tile produces no
    tiles at all, and telling that user to loosen the GC range sends them the
    wrong way. Which of the three fires is decided purely by how many input
    sequences are under ``tile_length``.
    """
    # Distinguish "filters rejected everything" from "nothing could even be
    # tiled because the sequences are shorter than one tile" \u2014 the latter
    # is a common input mistake that the generic message hides.
    tile_len = tile_length
    try:
        short = {
            gid: len(seq)
            for gid, seq in gene_seqs.items()
            if len(seq) < tile_len
        }
        n_genes = len(gene_seqs)
    except (AttributeError, TypeError):
        short, n_genes = {}, 0
    if short and len(short) == n_genes:
        detail = ", ".join(f"{gid} ({n} bp)" for gid, n in short.items())
        raise PipelineError(
            f"Every input sequence is shorter than the {tile_len}-bp tile "
            f"length, so no probes could be tiled: {detail}. Provide longer "
            f"target sequences or reduce tile_length."
        )
    if short:
        detail = ", ".join(f"{gid} ({n} bp)" for gid, n in short.items())
        raise PipelineError(
            f"No candidates passed filters. Note: these sequences are "
            f"shorter than the {tile_len}-bp tile length and produced no "
            f"tiles: {detail}. Check filter thresholds for the rest."
        )
    raise PipelineError(
        "No candidates passed filters \u2014 cannot continue. Loosen the "
        "filter thresholds (GC range, homopolymer, N-count) or verify the "
        "input sequence."
    )


class _CheckpointLoad(NamedTuple):
    """The candidate-cache lookup's answer, and the state it restores.

    ``use_cand_cache`` and ``fingerprint`` are returned even when nothing was
    loaded: they are what the SAVE side needs at the end of the stage, and
    recomputing the fingerprint there would risk the two halves disagreeing
    about which run this checkpoint belongs to.

    ``loaded`` False leaves ``passing`` empty and ``filter_stats`` at {} — the
    same starting state the filter stage assumes.
    """

    use_cand_cache: bool
    fingerprint: Optional[str]
    loaded: bool
    passing: List[ProbeCandidate]
    tiles_total: int
    filter_stats: Dict[str, int]


def _load_candidate_checkpoint(
    config: RunConfig,
    gene_seqs: Dict[str, str],
    output_dir: str,
    *,
    screen_hits_by_probe: Optional[Dict[str, List[Dict[str, Any]]]],
    _emit,
) -> _CheckpointLoad:
    """Reload the post-cap candidate set from an earlier link, if there is one.

    Called from :func:`run_pipeline`. A whole-transcriptome run
    generates hundreds of millions of tiles before the screen even starts; on an
    HPC resubmission with an unchanged input and config, this skips straight to
    the screen, which then resumes from its own chunk cache.

    Gated to the paths where that is both correct and worth it: not multi-strain
    (which rebuilds a consensus that must run), and only when a screen will
    actually happen (otherwise there is no long stage to resume into).
    """

    cand_loaded = False
    cand_fp: Optional[str] = None
    passing: List[ProbeCandidate] = []
    tiles_total = 0
    filter_stats: Dict[str, int] = {}
    _use_cand_cache = (
        config.resume_screen
        and not config.multi_strain
        and screen_hits_by_probe is None
        and not config.screen.skip_screen
        and bool(gene_seqs)
    )
    if _use_cand_cache:
        from . import candidate_cache
        cand_fp = candidate_cache.fingerprint(config)
        _hit = candidate_cache.load(output_dir, cand_fp)
        if _hit is not None:
            passing, _cmeta = _hit
            tiles_total = int(_cmeta.get("tiles_total", 0))
            filter_stats = {
                "total": tiles_total,
                "passed": int(_cmeta.get("filter_passed", len(passing))),
            }
            cand_loaded = True
            _emit("filter",
                  f"Resumed {len(passing)} candidates from checkpoint "
                  f"(skipped tiling + filtering)", 44,
                  {"resumed_candidates": len(passing)})

    return _CheckpointLoad(
        use_cand_cache=_use_cand_cache,
        fingerprint=cand_fp,
        loaded=cand_loaded,
        passing=passing,
        tiles_total=tiles_total,
        filter_stats=filter_stats,
    )


def _save_candidate_checkpoint(
    output_dir: str,
    cand_fp: Optional[str],
    prescreen_set: List[ProbeCandidate],
    *,
    use_cand_cache: bool,
    cand_loaded: bool,
    tiles_total: int,
    filter_stats: Dict[str, int],
    _scr_stage: str,
    _emit,
) -> None:
    """Persist the settled pre-screen candidate set for a resumed job.

    Called from :func:`run_pipeline`; the counterpart of
    :func:`_load_candidate_checkpoint`. Skipped when the set was itself
    reloaded from a checkpoint (``cand_loaded``), which would only rewrite
    what is already on disk.
    """

    _use_cand_cache = use_cand_cache

    # Checkpoint the final pre-screen candidate set so a resubmitted job skips the
    # tile+filter+cap work and goes straight to (resumable) BLAST. Best-effort:
    # a write failure must never abort an otherwise-good run.
    if _use_cand_cache and not cand_loaded and prescreen_set:
        from . import candidate_cache
        # _use_cand_cache is in the guard above, and it is the branch that
        # computes the fingerprint.
        assert cand_fp is not None
        try:
            candidate_cache.save(
                output_dir, cand_fp, prescreen_set,
                tiles_total=tiles_total,
                filter_passed=filter_stats.get("passed", len(prescreen_set)),
            )
            _emit(_scr_stage,
                  f"Checkpointed {len(prescreen_set)} candidates for resume", 44)
        except OSError:
            logger.warning("Could not write candidate checkpoint", exc_info=True)


def _plain_tile_iter(
    gene_seqs: Dict[str, str], *, length: int, step: int,
) -> Iterator[ProbeCandidate]:
    """Tile every input record, one window at a time — the single-sequence arm
    of Stage 2, and the counterpart of :func:`_stage_multi_strain_tiling`.

    Called from :func:`run_pipeline`. LAZY, and that is the whole
    point: a whole-transcriptome run yields hundreds of millions of windows and
    is only ever consumed inside the filter stage, which keeps a bounded subset.
    Forcing this generator into a list would reintroduce exactly the memory
    ceiling the streaming filter exists to avoid.

    ``probe_id`` is ``<record>_t<start padded to 6>``, which is what every
    downstream artefact — results.csv, the order CSVs, the resume caches — keys
    a probe by, so the format is part of the output contract, not a detail.
    """
    for gid, gseq in gene_seqs.items():
        for start, end, subseq in tile_sequence(gseq, length=length, step=step):
            yield ProbeCandidate(
                gene_id=gid,
                probe_id=f"{gid}_t{start:06d}",
                sequence=subseq,
                start=start,
                end=end,
            )


def _apply_batch_cap(
    passing: List[ProbeCandidate],
    *,
    plan: _PrescreenPlan,
    _scr_stage: str,
    _emit,
) -> List[ProbeCandidate]:
    """Stage 3b: the per-gene pre-screen cap, applied to a materialised set.

    Called from :func:`run_pipeline`. The normal path streams this
    cap inside the filter (see ``_PrescreenPlan.stream_cap``); this batch form
    covers the cases that path skips — chiefly multi-strain runs, which keep
    heavy per-strain state and are not a whole-transcriptome memory concern.
    Same result, just materialised first.

    Returns ``passing`` unchanged when no cap is configured, so the caller only
    has to decide whether the cap has ALREADY run, not whether it applies.
    """
    _cap_n = plan.cap_n
    if not (plan.cap_active and passing):
        return passing
    from .ranking import select_top_per_gene
    before_cap = len(passing)
    # cap_active is checked above, and it is the branch that assigns both of
    # these.
    assert plan.gene_key is not None and plan.weights is not None
    passing = select_top_per_gene(
        passing, gene_key=plan.gene_key, seq_length_of=plan.seq_len_of,
        n=_cap_n, weights=plan.weights,
    )
    # The cap ranks candidates individually, so it too can split a family.
    passing = _reconcile_degenerate_families(passing)
    _emit(_scr_stage,
          f"Pre-screen cap: top {_cap_n}/gene kept "
          f"{len(passing)} of {before_cap} candidates for screening", 44)
    return passing


PREPARE_MANIFEST_FILE = "prepare_manifest.json"

# Stage name a --prepare-only run files its cap / checkpoint / artefact events
# under, in place of the screen engine's own. See run_pipeline's `_prep_stage`.
STAGE_PREPARE = "prepare"


def _flexprobe_version() -> str:
    """This engine's version string, imported lazily.

    Local so the module keeps its existing import graph — `pipeline` is
    imported by the package's own `__init__`, and a top-level
    `from . import __version__` there would be a cycle.
    """
    from . import __version__
    return __version__


def _stage_prepare_only(
    config: RunConfig,
    output_dir: str,
    *,
    prescreen_set: List[ProbeCandidate],
    screen_set: List[ProbeCandidate],
    gene_seqs: Dict[str, str],
    gene_aware: bool,
    group_members: Optional[Dict[str, Set[str]]],
    tiles_total: int,
    filter_stats: Dict[str, int],
    stage_timer: "_StageTimer",
    cand_loaded: bool,
    _prescreen_capture: str,
    _tally: Any,
    alignment: Optional[StrainAlignment] = None,
    _emit,
) -> Dict[str, Any]:
    """Finish a ``--prepare-only`` run: provenance, marker, and an honest stop.

    Called from :func:`run_pipeline` at exactly one point — after the candidate
    checkpoint is saved and the pre-screen artefacts are written, and before
    :func:`_stage_screen`. That position is the whole contract: everything a
    later screen job needs is on disk, and nothing that needs the screen has
    been attempted.

    THREE THINGS IT MUST NOT DO, each of which would turn a working split into
    a silent data-loss bug:

    * It must not write ``manifest.json``. That file is the finished-run
      contract — a full funnel ran, results.csv sits beside it — and the GUI
      importer reads it that way. The prepare half writes ``prepare_manifest.json``
      instead, which also keeps it out of the way of the screen job's real
      manifest and its ``overwrite=False`` guard.
    * It must not clear the caches. ``_stage_provenance_and_summary`` deletes
      ``.candidate_cache/`` and ``.screen_cache/`` because a finished run has
      no use for them. Here the candidate checkpoint IS the deliverable.
    * It must not write ``.flexiprobes_DONE``. See ``constants.PREPARED_MARKER``.

    The manifest it does write is deliberately narrower than the full one. It
    omits ``probe_barcode``, ``offtarget_recount`` and
    ``junction_duplex_measured``, and its ``db_sha256s`` and
    ``reference_composition`` are empty, because chemistry, the recount and the
    reference were not touched by this process. Recording them from the config
    would be a claim about work that did not happen.
    """
    from .writers.manifest import write_manifest as _write_manifest
    from .config_io import filled_from_defaults, portable_config_dict, save_config

    # Paths relative to the output directory the manifest lives in.
    params: Dict[str, Any] = portable_config_dict(config, output_dir)
    params.update({
        'modality_label': f"{config.modality.name}_{config.modality.version}",
        'gene_lengths': {gid: len(seq) for gid, seq in gene_seqs.items()},
        # What ran, beside what was asked for (params['group_by_gene']):
        # grouping is switched off on input with no gene: field.
        'group_by_gene_applied': bool(gene_aware),
        'gene_isoforms': (
            {g: sorted(ts) for g, ts in group_members.items()}
            if gene_aware and group_members is not None else {}
        ),
        'require_tn_junction': config.filters.require_tn_junction,
        'blast_saw_prefiltered_subset': bool(config.screen_top_per_gene),
        'screen_engine': config.screen.engine,
        # THE two fields that stop this being mistaken for a finished run.
        # `prepare_only` says what was asked for; `run_complete` says what came
        # out. Both, rather than one, because a reader checking either alone
        # must get the right answer.
        'prepare_only': True,
        'run_complete': False,
        # What the screen job will pick up, so the handover is auditable from
        # the prepared directory alone rather than by re-deriving it.
        'prepared_candidates': len(prescreen_set),
        'prepared_screen_queue': len(screen_set),
        'prepared_marker': PREPARED_MARKER,
        # A prepare job that RESUMED from an existing checkpoint did no tiling
        # of its own. Recorded because it makes the prepare-stage timings in
        # this same file uninterpretable otherwise — they would read as a
        # transcriptome tiled in 60 seconds.
        'prepared_from_existing_checkpoint': bool(cand_loaded),
    })

    manifest = RunManifest(
        inputs=_portable_inputs(_hash_input_fasta(config.fasta_path),
                                output_dir),
        code_version=_code_version(),
        # screen_used=False: this process never invoked a search engine, so it
        # must not stamp a blastn version it did not run.
        tool_versions=_collect_tool_versions(screen_used=False),
        # A SPLIT multi-strain run builds its consensus in THIS half (tiling is
        # Stage 2, which prepare owns) and then exits, so without this the
        # screen half could never learn which aligner ran. Empty for a
        # single-species run.
        consensus=_consensus_provenance(alignment),
        elapsed_seconds=stage_timer.snapshot(),
        params=params,
        db_sha256s={},
        reference_composition={},
        seed=config.seed,
        config_filled_from_defaults=filled_from_defaults(config),
    )
    man_path = os.path.join(output_dir, PREPARE_MANIFEST_FILE)
    _write_manifest(manifest, man_path, overwrite=True)

    # Same file, same name, same loader as a full run writes — the screen job
    # reads config.json from the bundle, but a human comparing the two halves
    # needs the prepared directory to say what config produced its checkpoint.
    save_config(config, os.path.join(output_dir, 'run_config.json'),
                portable=True)

    stats: Dict[str, Any] = {
        "genes": len(gene_seqs),
        "tiles_total": tiles_total,
        "filter_passed": len(prescreen_set),
        "filter_stats": filter_stats,
        # Zero because nothing was screened, ranked, paired or exported — not
        # because those stages rejected everything. `prepare_only` below is
        # what distinguishes the two, and it is why cli.cmd_design must not
        # read a 0 here as the "designed nothing" failure it is for a full run.
        "specificity_passed": 0,
        "specificity_screened": False,
        "specificity_status_counts": {},
        "pairs_output": 0,
        "prepare_only": True,
        "run_complete": False,
        "prepared_candidates": len(prescreen_set),
        "prepared_screen_queue": len(screen_set),
        "prepared_from_existing_checkpoint": bool(cand_loaded),
        "discovery_only": config.discovery_only,
        "auto_select_per_gene": config.auto_select_per_gene,
        "multi_strain": config.multi_strain,
        "group_by_gene": gene_aware,
        "prescreen_capture": _prescreen_capture,
        "prescreen_precap_rows": (
            _tally.rows_written if (_tally is not None
                                    and config.emit_prescreen_precap_rows)
            else None
        ),
    }

    # LAST, after everything above is on disk. A reader (or the screen job's
    # shell guard) that finds this file may assume the checkpoint, the
    # artefacts, the prepare manifest and run_config.json are all complete, so
    # it cannot be written before they are.
    marker_path = os.path.join(output_dir, PREPARED_MARKER)
    with open(marker_path, "w") as fh:
        json.dump({
            "prepared_at": datetime.now(timezone.utc).isoformat(),
            "app_version": manifest.app_version or _flexprobe_version(),
            "candidates": len(prescreen_set),
            "screen_queue": len(screen_set),
            "tiles_total": tiles_total,
            "resumed_from_existing_checkpoint": bool(cand_loaded),
            "note": (
                "Prepared, NOT finished. Ingest, tiling and filtering are done "
                "and the candidate checkpoint is in .candidate_cache/. No "
                "off-target screen, ranking, chemistry or export has run, and "
                "there is no results.csv. Run the screen job against this same "
                "output directory to finish it."
            ),
        }, fh, indent=2)
        fh.write("\n")

    _emit("prepared",
          f"Prepare-only: {len(prescreen_set)} candidates checkpointed; "
          f"stopping before the screen (wrote {PREPARED_MARKER})", 100,
          {"prepared_candidates": len(prescreen_set)})

    stage_timer.write_report(os.path.join(output_dir, "run_timings.json"))

    return {
        "pairs": [],
        "manifest": manifest,
        "stats": stats,
    }


def run_pipeline(
    config: RunConfig,
    cancel_event: Optional[threading.Event] = None,
    progress_callback: Optional[Callable[[ProgressEvent], None]] = None,
    target_ids: Optional[Set[str]] = None,
    screen_hits_by_probe: Optional[Dict[str, List[Dict[str, Any]]]] = None,
) -> Dict[str, Any]:
    """Execute the full probe design pipeline.

    Parameters
    ----------
    config : RunConfig
        Complete run configuration.
    cancel_event : threading.Event, optional
        Set this event to request cancellation.
    progress_callback : callable, optional
        Called with ProgressEvent at each stage boundary.
    target_ids : set of str, optional
        Target gene IDs for specificity. If None, uses all genes from FASTA.
    screen_hits_by_probe : dict, optional
        Pre-computed BLAST hits keyed by probe_id. If provided, BLAST step
        is skipped (the CLI passes these for --cached-hits).

    Returns
    -------
    dict with keys: pairs, manifest, stats

    On a ``config.prepare_only`` run the SHAPE is the same and the CONTENT is
    not: ``pairs`` is empty, ``stats['pairs_output']`` is 0, and
    ``stats['prepare_only']`` is True. The zero means "no probes were
    attempted", not "every candidate was rejected" — a caller that reads
    ``pairs_output == 0`` as a failed run (as ``cli.cmd_design`` does for a full
    run) must check ``prepare_only`` first. See :func:`_stage_prepare_only`.
    """
    def _check_cancel():
        if cancel_event and cancel_event.is_set():
            raise CancellationError("Pipeline cancelled by user")

    # Times every stage off the progress-event boundaries below; the totals go
    # into the manifest so a multi-day HPC run's per-stage cost is recoverable
    # from the results alone, without the job log.
    stage_timer = _StageTimer()

    def _emit(stage: str, msg: str, pct: Optional[float] = 0.0,
              data: Optional[Dict] = None):
        stage_timer.enter(stage, msg, pct)
        ev = ProgressEvent(stage, msg, pct, data)
        logger.info("%s", ev)
        if progress_callback:
            progress_callback(ev)

    _prepare_config(config)

    output_dir = config.output_dir
    os.makedirs(output_dir, exist_ok=True)
    # Append-as-you-go, so a walltime kill still says where it died. A resumed
    # chain link appends to the same file: the links read as one timeline.
    stage_timer.open_sink(os.path.join(output_dir, "run_timeline.jsonl"))

    # Progress events, the log and the timing records all name the engine that
    # is ACTUALLY screening. `_scr_stage` is also the elapsed_seconds key, so a
    # manifest says which screen the run paid for rather than always "blast".
    _scr = screen_labels(config.screen.engine)
    _scr_stage = _scr["stage"]

    # The pre-screen cap, the checkpoint save and the pre-screen artefacts are
    # all announced under the SCREEN's stage name, because on a normal run they
    # are the screen stage's opening moves. On a --prepare-only run no screen
    # happens at all, so filing them there would put a `kmer` entry with real
    # seconds against it into a run_timings.json for a process that never built
    # a k-mer index. Identical to `_scr_stage` on every run that is not
    # prepare-only.
    _prep_stage = STAGE_PREPARE if config.prepare_only else _scr_stage

    # ---- Stage 1: Ingest ----
    _ingest = _stage_ingest(config, target_ids, _emit, _check_cancel)
    gene_seqs = _ingest.gene_seqs
    target_ids = _ingest.target_ids
    gene_aware = _ingest.gene_aware
    gene_of = _ingest.gene_of
    group_members = _ingest.group_members

    # ---- Resume: reload the post-cap candidate set, skipping tile+filter ----
    # A whole-transcriptome run generates hundreds of millions of tiles before
    # BLAST even starts. On an HPC resubmission with an unchanged input + config,
    # reload the checkpointed candidates and jump straight to BLAST (which then
    # resumes from its own chunk cache). Gated to the non-multi-strain path
    # (multi-strain rebuilds a consensus that must run) and only when we'd
    # actually BLAST (otherwise there's no long stage to resume).
    _ckpt = _load_candidate_checkpoint(
        config, gene_seqs, output_dir,
        screen_hits_by_probe=screen_hits_by_probe, _emit=_emit,
    )
    _use_cand_cache = _ckpt.use_cand_cache
    cand_fp = _ckpt.fingerprint
    cand_loaded = _ckpt.loaded
    passing = _ckpt.passing
    tiles_total = _ckpt.tiles_total
    filter_stats = _ckpt.filter_stats

    # ---- Stage 2: Tile (lazy generators; filtered streaming in Stage 3) ----
    _emit("tile", "Generating tiles", 15)
    _check_cancel()

    alignment: Optional[StrainAlignment] = None
    # Specificity strain rule; overridden to PAN_STRAIN by constitutive
    # (isoform) mode inside the multi-strain branch below.
    consensus_strain_rule = StrainRule.EXACT_SET
    # Which input records every candidate is REQUIRED to bind — the scope of
    # the filter funnel's junction-mismatch rule. Resolved from the consensus
    # strategy once the alignment exists (below); stays None on a single-
    # sequence run, which leaves that rule inert.
    on_target_records: Optional[FrozenSet[str]] = None
    if config.multi_strain:
        _ms = _stage_multi_strain_tiling(config, gene_seqs, output_dir, _emit)
        alignment = _ms.alignment
        gene_seqs = _ms.gene_seqs
        target_ids = _ms.target_ids
        on_target_records = _ms.on_target_records
        consensus_strain_rule = _ms.consensus_strain_rule
        candidate_iter = _ms.candidate_iter
    else:
        candidate_iter = _plain_tile_iter(
            gene_seqs, length=config.tile_length, step=config.tile_step,
        )

    # ---- Pre-screen cap + artefact plan ----
    _plan = _plan_prescreen(
        config, gene_seqs, output_dir,
        gene_aware=gene_aware, gene_of=gene_of,
        screen_hits_by_probe=screen_hits_by_probe, cand_loaded=cand_loaded,
    )
    # Only the two fields the MANIFEST reads are unpacked here; the cap and the
    # artefact writer are handed the plan whole.
    _prescreen_capture = _plan.capture
    _tally = _plan.tally
    # Set by the streamed filter below, not part of the plan: it records what
    # HAPPENED, not what was planned.
    _stream_capped = False

    # ---- Stage 3: Tile + Filter (streamed) ----
    # Tiles are generated lazily and filtered as they go, so a whole-
    # transcriptome run only ever holds the passing set + a capped reject
    # sample in memory rather than every tile.
    if not cand_loaded:
        _flt = _stage_filter(
            config, gene_seqs, candidate_iter, output_dir,
            plan=_plan, on_target_records=on_target_records,
            _scr_stage=_prep_stage, _emit=_emit, _check_cancel=_check_cancel,
        )
        passing = _flt.passing
        filter_stats = _flt.filter_stats
        tiles_total = _flt.tiles_total
        _stream_capped = _flt.stream_capped

    if not passing:
        _raise_no_candidates(gene_seqs, config.tile_length)

    # ---- Stage 3b: pre-screen per-gene cap (batch fallback) ----
    # The normal path already streamed the cap above (``_stream_capped``). This
    # batch form covers the cases that path skips — chiefly multi-strain runs,
    # which keep heavy per-strain state and aren't a whole-transcriptome memory
    # concern. Same result as the streamed cap, just materialised first.
    if not cand_loaded and not _stream_capped:
        passing = _apply_batch_cap(
            passing, plan=_plan, _scr_stage=_prep_stage, _emit=_emit,
        )

    # ---- The candidate set is now SETTLED ----
    # Every path above has converged here: resume-from-checkpoint, the streamed
    # filter, the parallel filter, the streamed cap and the batch cap. Nothing
    # below narrows it again (proved mechanically: `passing` has no further
    # assignment in this function).
    #
    # It gets its own name because `passing` is a funnel accumulator — it means
    # "whatever is still standing", which is honest while the funnel is running
    # and useless 400 lines later, where the reader needs to know WHICH
    # narrowings have already happened. From here down the answer is: all of
    # them. (`passing` is rebound several times above, meaning something
    # different at each.)
    #
    # NOTE the one place this name is deliberately NOT used: the manifest and
    # summary still report `len(...)` of this set as "filter_passed", which is
    # the post-cap count, not the filter's own. That gap is real, intentional
    # and already documented — see _write_filter_report and the comment above
    # the first _reconcile_degenerate_families call. Renaming here does not
    # change it; it just stops the name pretending the gap isn't there.
    prescreen_set = passing

    _save_candidate_checkpoint(
        output_dir, cand_fp, prescreen_set,
        use_cand_cache=_use_cand_cache, cand_loaded=cand_loaded,
        tiles_total=tiles_total, filter_stats=filter_stats,
        _scr_stage=_prep_stage, _emit=_emit,
    )

    # ---- What BLAST will actually screen ----
    # A consensus parent is a mixed-base string: blastn has no sensible
    # interpretation of it, and screening a single "representative" resolution
    # would be worse than not screening at all, because a variant specific to
    # strain B can hit a human transcript where strain A's variant does not —
    # exactly the risk the screen exists to catch. So the concrete variants are
    # screened SEPARATELY and their verdicts roll up onto the parent afterwards
    # (all children must pass).
    #
    # Defined HERE, above the artefacts below, because those document themselves
    # as the set BLAST screened and that has to stay literally true. Identical to
    # ``passing`` on every run with no degenerate parents.
    _screen = [c for c in prescreen_set if c.probe_role != PROBE_ROLE_CONSENSUS]
    _n_parents = len(prescreen_set) - len(_screen)
    if _n_parents:
        _emit(_prep_stage,
              f"Screening {len(_screen)} concrete variants; {_n_parents} "
              f"degenerate parent(s) inherit their variants' verdicts", 44,
              {"degenerate_parents": _n_parents})

    # ---- Pre-screen candidate artefacts (opt-in) ----
    _write_prescreen_artefacts(
        config, output_dir, _screen,
        plan=_plan, _scr_stage=_prep_stage, _emit=_emit,
    )

    # ---- STOP HERE on a prepare-only run ----
    # This is the split point, and it is here rather than a few lines either
    # side for two reasons. Above it, the candidate checkpoint is saved and the
    # pre-screen artefacts are written, so the screen job has everything it
    # needs. Below it, `_stage_screen` builds the k-mer index and forks
    # `screen_workers` processes — the work whose memory the split exists to
    # keep away from the tiling high-water mark.
    #
    # Nothing between the checkpoint save and this return narrows `passing` /
    # `prescreen_set`, so the set a resumed screen job loads is the set this
    # run would have screened. That is what makes the split result-neutral, and
    # it is asserted end-to-end in tests/test_prepare_only.py
    # (test_resume_after_prepare_matches_single_process).
    if config.prepare_only:
        return _stage_prepare_only(
            config, output_dir,
            prescreen_set=prescreen_set,
            screen_set=_screen,
            gene_seqs=gene_seqs,
            gene_aware=gene_aware,
            group_members=group_members,
            tiles_total=tiles_total,
            filter_stats=filter_stats,
            stage_timer=stage_timer,
            cand_loaded=cand_loaded,
            _prescreen_capture=_prescreen_capture,
            _tally=_tally,
            _emit=_emit,
            alignment=alignment,
        )

    # ---- Stage 4: off-target screen ----
    _sc = _stage_screen(
        config, _screen, output_dir,
        target_ids=target_ids,
        gene_of=gene_of,
        group_members=group_members,
        consensus_strain_rule=consensus_strain_rule,
        screen_hits_by_probe=screen_hits_by_probe,
        n_prescreen=len(prescreen_set),
        _scr=_scr,
        _emit=_emit,
        _check_cancel=_check_cancel,
    )
    # Only the four fields the LATER stages read are unpacked here; the rest of
    # _ScreenResult is handed to Stage 5 whole (`sc=_sc`), which is the only
    # reader of the screen's own bookkeeping.
    screen_ran = _sc.screen_ran
    _spec_kwargs = _sc.spec_kwargs
    _half_len = _sc.half_len
    _subjects = _sc.subjects

    # ---- Stage 5: Specificity ----
    _spec = _stage_specificity(
        config, prescreen_set, _screen,
        target_ids=target_ids, sc=_sc, n_parents=_n_parents,
        gene_aware=gene_aware, _emit=_emit, _check_cancel=_check_cancel,
    )
    specific = _spec.specific
    gated_low_iso = _spec.gated_low_iso


    # ---- Stage 5b: null-model calibration of the off-target gate ----
    # What fraction of composition-matched *noise* would this gate have
    # rejected? Without that number the run's own rejection count is
    # uninterpretable — a gate can reject 60% of windows because the probes are
    # unspecific, or because the database is large enough that random 50-mers
    # hit it too. Only meaningful when BLAST actually ran.
    null_model_stats: Optional[Dict[str, Any]] = None
    if config.null_model_shuffles > 0 and screen_ran and not config.screen.skip_screen:
        try:
            null_model_stats = _run_null_model(
                config, _screen, output_dir,
                spec_kwargs=_spec_kwargs, target_ids=target_ids,
                half_len=_half_len, emit=_emit, cancel=_check_cancel,
                subjects=_subjects,
            )
        except FlexprobeError:
            raise
        except Exception:
            # Calibration is a reported diagnostic, not a gate. A failure here
            # must never lose an otherwise-good design run.
            logger.warning("Null-model calibration failed", exc_info=True)

    # ---- Stage 6: Rank ----
    final_ranked = _stage_rank(
        specific, config, gene_seqs,
        gene_aware=gene_aware, _emit=_emit, _check_cancel=_check_cancel,
    )

    # ---- Stage 7: Chemistry & Pair ----
    _chem_result = _stage_chemistry_and_pair(
        final_ranked, config, _emit=_emit, _check_cancel=_check_cancel,
    )
    pairs = _chem_result.pairs
    bc_labels = _chem_result.bc_labels
    multi_barcode = _chem_result.multi_barcode
    probe_bc_label = _chem_result.probe_bc_label
    probe_bc_seq = _chem_result.probe_bc_seq
    barcode_pairs = _chem_result.barcode_pairs

    # ---- Stage 8: Export ----
    export, overwrite = _stage_export(
        final_ranked, pairs, config, output_dir,
        barcode_pairs=barcode_pairs, multi_barcode=multi_barcode,
        _emit=_emit, _check_cancel=_check_cancel,
    )

    # ---- Provenance ----
    return _stage_provenance_and_summary(
        config, output_dir,
        pairs=pairs,
        final_ranked=final_ranked,
        passing=prescreen_set,
        specific=specific,
        gene_seqs=gene_seqs,
        gene_aware=gene_aware,
        group_members=group_members,
        alignment=alignment,
        bc_labels=bc_labels,
        probe_bc_label=probe_bc_label,
        probe_bc_seq=probe_bc_seq,
        export=export,
        overwrite=overwrite,
        stage_timer=stage_timer,
        tiles_total=tiles_total,
        filter_stats=filter_stats,
        gated_low_iso=gated_low_iso,
        screen_ran=screen_ran,
        null_model_stats=null_model_stats,
        _prescreen_capture=_prescreen_capture,
        _subjects=_subjects,
        _tally=_tally,
        _emit=_emit,
    )
