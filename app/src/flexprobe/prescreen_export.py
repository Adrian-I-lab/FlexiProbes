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

"""Persist the pre-screen candidate set as a first-class analysis artefact.

Why this exists
---------------
``results.csv`` contains only the candidates that SURVIVED the off-target
screen. Any statistic computed over it — rank enrichment above all — is
conditioned on the outcome, so a null model built from it is a null model of the
survivors, not of the candidates. Correcting that needs the denominator: the set
of windows that were actually *considered*.

The pipeline already materialises that set (it is what goes into BLAST) and
already checkpoints it, but the checkpoint is an internal resume format —
``.candidate_cache/`` is dot-prefixed, schema-versioned, deleted on success, and
free to change shape whenever the resume logic does. Analysis must not depend on
it. So this module writes a stable, self-describing, gzipped CSV instead, and
the cache stays what it is: a cache.

What gets written
-----------------
``candidates_prescreen.csv.gz`` — one row per candidate that entered BLAST, i.e.
the set AFTER the ``screen_top_per_gene`` cap. This is exactly the denominator
for any statement about the specificity screen.

On a degenerate-consensus run that set is the VARIANTS, not the consensus
parents. A parent carries IUPAC ambiguity, is never orderable and is never
queried, so it is excluded upstream (see ``pipeline``'s ``_screen``) and does
not appear here — the file stays literally "what BLAST saw". But note what that
does to the unit: one qualifying window contributes as many rows as it has
distinct per-record resolutions, so a null model built over these rows is a null
over screened *oligos*, not over windows. Group by ``parent_probe_id`` first if
the window is the unit you mean.

``prescreen_cap_summary.csv.gz`` — one row per gene: how many candidates passed
the filters BEFORE the cap, and how many the cap kept. The cap is itself a
selection step applied before BLAST, and it ranks on ``prescreen_score``, which
correlates with the outcome — so a null that ignores it is still biased. These
counts, together with the per-row ``prescreen_score`` and ``rank_in_gene`` in the
file above, give the cap's truncation point per gene and let it be modelled
rather than assumed away. Genes are keyed the same way the cap keys them (the
Ensembl gene for gene-aware runs, else the record id).

``candidates_precap.csv.gz`` + ``prescreen_transcripts.csv.gz`` — OPTIONAL, off
unless asked for. EVERY window that passed the chemistry filters, before the
cap. This is the better denominator, because the cap is itself a selection step
and it ranks on a score that correlates with the outcome — but it is 7-11x the
row count of the capped set on the old engine (74M rows across the three
mammalian transcriptomes) and more than that under v0.4.0's wider funnel.

That size is why this file does NOT use the 12-column schema below. It carries
four columns, and everything else is reconstructed:

  t                the transcript's integer index into
                   ``prescreen_transcripts.csv.gz``. The transcript id is a
                   ~20-character string repeated on every row and would be, by
                   a wide margin, the largest single cost in the file; the
                   dictionary encoding removes it. (Rejected alternative:
                   Parquet, which does this natively and better — but pyarrow
                   is not a flexiprobes dependency and so is not in the HPC
                   container, and adding a ~100 MB wheel to the desktop app for
                   one optional artefact is not a trade worth making. The
                   measured bytes per row made it unnecessary.)
  start            window start on the transcript
  prescreen_score   the score the cap ranks on, 4 dp
  junction_base    the mRNA base at the ligation junction — the variable the
                   hard/soft TN comparison turns on, one character

  end              NOT stored: it is always ``start + tile_length``, and
                   tile_length is a scalar in the config.
  probe_id         NOT stored: it is ``f"{transcript_id}_t{start:06d}"``, the
                   tiler's own rule, so it joins to results.csv.
  sequence         NOT stored: slice the input FASTA. Storing it would dominate
                   the file by an order of magnitude.
  filter outcome   NOT stored: every row in this file passed every filter, by
                   construction. The rejected windows' outcomes are already
                   summarised, per reason, in ``filter_report.json``.

``prescreen_index.json`` states which of the two candidate sets a given run
captured, plus the schemas and the reconstruction rules above. A gzipped CSV
cannot carry its own metadata and "pre-cap or post-cap?" must never be a
guess.

Columns of the post-cap file, and why each one is here
-----------------------------------------------------
Every column is either an identifier needed to join back to ``results.csv`` /
the input FASTA, or an input to ``prescreen_score`` — the quantity the cap ranks
on. Nothing else is carried, because at 4M+ rows per species every column has a
real cost.

  probe_id        join key to results.csv; unique per candidate
  transcript_id   the FASTA record the window came from
  gene            the cap's grouping key (Ensembl gene, or the record id)
  start, end      window coordinates on the transcript
  gc_content      reported alongside the per-half columns below; NOT itself a
                  prescreen_score input (GC proximity is scored per half — see
                  ``ranking._gc_proximity_per_half_score``)
  gc_left         prescreen_score term (GC proximity AND GC balance across the
                  two halves)
  gc_right        prescreen_score term (GC proximity AND GC balance across the
                  two halves)
  homopolymer_max prescreen_score term (homopolymer penalty)
  junction_base   the mRNA base at the ligation junction; drives the junction
                  multiplier in prescreen_score, and is the variable the
                  hard/soft TN comparison turns on
  prescreen_score  the score the cap ranked on — recomputed here rather than
                  read off the candidate, because it is not persisted on the
                  dataclass
  rank_in_gene    1-based rank by prescreen_score within the gene, so the cap's
                  cut point is explicit in the file

``tn_junction_ok`` is deliberately NOT written: it is derivable
(``junction_base == 'A'``) and, unlike ``junction_base``, it is not persisted in
the resume checkpoint — so a resumed run would emit it as False for every row
and the file would silently disagree with a clean run's.
"""
from __future__ import annotations

import csv
import gzip
import os
from typing import (
    Any, Callable, Dict, IO, Iterable, List, Optional, Tuple, Union,
)

from .models import ProbeCandidate

CANDIDATES_FILE = "candidates_prescreen.csv.gz"
PRECAP_FILE = "candidates_precap.csv.gz"
TRANSCRIPTS_FILE = "prescreen_transcripts.csv.gz"
SUMMARY_FILE = "prescreen_cap_summary.csv.gz"
INDEX_FILE = "prescreen_index.json"

# Which candidate set a run captured. Recorded in prescreen_index.json and in the
# run stats, so it is never inferred from which files happen to exist.
CAPTURE_NONE = "none"
CAPTURE_POSTCAP = "post_cap"
CAPTURE_PRECAP = "pre_cap"

# The pre-cap dump is written from inside the filter stream, so its compression
# cost is paid on the pipeline's critical path at ~100M rows. Level 5 measured
# within a few percent of level 9's ratio on this data (integer columns, sorted
# by transcript then start, so the redundancy is local) at a fraction of the CPU.
PRECAP_COMPRESSLEVEL = 5

_COLUMNS = (
    "probe_id", "transcript_id", "gene", "start", "end",
    "gc_content", "gc_left", "gc_right", "homopolymer_max",
    "junction_base", "prescreen_score", "rank_in_gene",
)

# Compact schema — see the module docstring for what is reconstructed and why.
_PRECAP_COLUMNS = ("t", "start", "prescreen_score", "junction_base")

_TRANSCRIPT_COLUMNS = ("t", "transcript_id", "gene", "length")

_SUMMARY_COLUMNS = (
    "gene", "candidates_passed_filters", "candidates_screened", "kept_fraction",
)


def _seq_length_unset(_transcript_id: str) -> int:
    """Stand-in for a missing ``seq_length_of`` on a tally that isn't writing.

    Only reachable if the "write_rows implies seq_length_of" invariant enforced
    in ``PrecapTally.__init__`` is ever broken; failing loudly beats scoring
    every candidate against a silent length of 0.
    """
    raise RuntimeError(
        "PrecapTally has no seq_length_of; row writing should be disabled")


def _row(c: ProbeCandidate, gene: str, score: float,
         rank: Union[int, str]) -> List:
    # `rank` is Union[int, str] because the unranked (streamed) path writes an
    # empty cell rather than a placeholder integer — see write_candidates.
    return [
        c.probe_id, c.gene_id, gene, c.start, c.end,
        f"{c.gc_content:.4f}", f"{c.gc_left:.4f}", f"{c.gc_right:.4f}",
        c.homopolymer_max, c.junction_base, f"{score:.6f}", rank,
    ]


def write_candidates(
    output_dir: str,
    candidates: Iterable[ProbeCandidate],
    *,
    gene_key: Callable[[ProbeCandidate], str],
    seq_length_of: Callable[[str], int],
    weights=None,
    filename: str = CANDIDATES_FILE,
    ranked: bool = True,
) -> str:
    """Write the pre-screen candidate set to ``output_dir/<filename>``.

    ``ranked`` computes each candidate's 1-based rank within its gene, which
    needs the whole set grouped in memory — fine for the capped set (a few
    million rows, already resident) and wrong for the pre-cap stream, which is
    written incrementally with ``rank_in_gene`` left empty.
    """
    from .ranking import DEFAULT_WEIGHTS, prescreen_score

    w = weights or DEFAULT_WEIGHTS
    path = os.path.join(output_dir, filename)
    tmp = path + ".tmp"
    by_gene: Dict[str, List] = {}
    for c in candidates:
        gene = gene_key(c)
        score = prescreen_score(c, seq_length_of(c.gene_id), w)
        by_gene.setdefault(gene, []).append((score, c))

    with gzip.open(tmp, "wt", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(_COLUMNS)
        for gene in sorted(by_gene):
            rows = by_gene[gene]
            if ranked:
                # Same ordering the cap uses, so rank_in_gene reproduces the
                # cap's own ordering rather than an approximation of it.
                rows.sort(key=lambda sc: (-sc[0], sc[1].start, sc[1].probe_id))
            for i, (score, c) in enumerate(rows, start=1):
                wr.writerow(_row(c, gene, score, i if ranked else ""))
    os.replace(tmp, path)
    return path


class PrecapTally:
    """Streaming per-gene tally — and optional compact dump — of the pre-cap set.

    Counted BEFORE the ``screen_top_per_gene`` cap, which is the whole point: the
    cap runs inside the filter stream (a bounded per-gene accumulator), so the
    pre-cap set is never materialised anywhere and this counter is the only place
    its size is ever observed.

    Counting is O(1) per candidate and always on. Writing rows is the expensive
    part — 74M+ rows across the three mammalian transcriptomes on the old engine,
    and the v0.4.0 funnel is wider — so the row dump is a separate opt-in and
    uses the compact schema described in the module docstring rather than the
    12-column one used for the (much smaller) post-cap set.
    """

    def __init__(
        self,
        gene_key: Callable[[ProbeCandidate], str],
        *,
        output_dir: Optional[str] = None,
        seq_length_of: Optional[Callable[[str], int]] = None,
        weights=None,
        write_rows: bool = False,
        compresslevel: int = PRECAP_COMPRESSLEVEL,
    ):
        self._gene_key = gene_key
        self.counts: Dict[str, int] = {}
        self.rows_written = 0
        self._fh: Optional[IO[str]] = None
        self._writer: Optional[Any] = None   # csv writer has no public type
        self._path: Optional[str] = None
        self._output_dir = output_dir
        self._seq_length_of = seq_length_of
        # Non-Optional views of the two dependencies the row-writing path needs.
        # __init__ refuses to enable write_rows unless both were supplied, so on
        # that path these always hold the real values. They exist so the write
        # path's types are non-Optional without a narrowing check inside add(),
        # which runs once per candidate — ~10^8 times on a whole-transcriptome
        # run. The sentinels are unreachable while write_rows is off, and raise
        # rather than return a silent 0 if that ever stops being true.
        self._row_dir: str = output_dir or ""
        self._row_len_of: Callable[[str], int] = (
            seq_length_of if seq_length_of is not None else _seq_length_unset
        )
        self._weights = weights
        # transcript_id -> (index, gene, length), assigned in first-seen order.
        # This is the dictionary encoding: the transcript id is a ~15-25 char
        # string repeated on every one of ~100M rows, and by far the dominant
        # cost if written inline. Emitted once, at close, as its own table.
        self._transcripts: Dict[str, Tuple[int, str, int]] = {}
        if write_rows:
            if output_dir is None or seq_length_of is None:
                raise ValueError(
                    "write_rows needs output_dir and seq_length_of")
            self._path = os.path.join(output_dir, PRECAP_FILE)
            self._fh = gzip.open(self._path + ".tmp", "wt", newline="",
                                 encoding="utf-8", compresslevel=compresslevel)
            self._writer = csv.writer(self._fh)
            self._writer.writerow(_PRECAP_COLUMNS)

    def _transcript_index(self, c: ProbeCandidate, gene: str) -> int:
        rec = self._transcripts.get(c.gene_id)
        if rec is None:
            # Only ever reached from add()'s writing branch, so _row_len_of is
            # the caller-supplied function rather than the sentinel.
            rec = (len(self._transcripts), gene, self._row_len_of(c.gene_id))
            self._transcripts[c.gene_id] = rec
        return rec[0]

    def add(self, c: ProbeCandidate) -> None:
        gene = self._gene_key(c)
        self.counts[gene] = self.counts.get(gene, 0) + 1
        if self._writer is not None:
            from .ranking import DEFAULT_WEIGHTS, prescreen_score
            w = self._weights or DEFAULT_WEIGHTS
            score = prescreen_score(c, self._row_len_of(c.gene_id), w)
            self._writer.writerow((
                self._transcript_index(c, gene),
                c.start,
                f"{score:.4f}",
                c.junction_base,
            ))
            self.rows_written += 1

    def add_many(self, cands: Iterable[ProbeCandidate]) -> None:
        for c in cands:
            self.add(c)

    def close(self) -> Optional[str]:
        """Finalise the row dump and write its transcript dictionary."""
        path = self._path
        # `path` is set together with `_fh` and never independently; binding it
        # locally is what lets the rename below be typed.
        if self._fh is not None and path is not None:
            self._fh.close()
            os.replace(path + ".tmp", path)
            self._fh = None
            self._writer = None
            self._write_transcripts()
        return self._path

    def _write_transcripts(self) -> str:
        path = os.path.join(self._row_dir, TRANSCRIPTS_FILE)
        tmp = path + ".tmp"
        with gzip.open(tmp, "wt", newline="", encoding="utf-8") as fh:
            wr = csv.writer(fh)
            wr.writerow(_TRANSCRIPT_COLUMNS)
            for tid, (idx, gene, length) in sorted(
                    self._transcripts.items(), key=lambda kv: kv[1][0]):
                wr.writerow((idx, tid, gene, length))
        os.replace(tmp, path)
        return path


def write_index(
    output_dir: str,
    *,
    capture: str,
    tile_length: int,
    screen_top_per_gene: int,
    screened_rows: int,
    precap_rows: Optional[int] = None,
) -> str:
    """Write ``prescreen_index.json`` — the artefacts' self-description.

    A gzipped CSV cannot carry its own metadata, and "which candidate set did
    this run actually capture?" is the single question an analyst must never
    have to infer from filenames. So it is stated here, in words, next to the
    files: the capture mode, the schema of each file, and the two reconstruction
    rules (probe_id, and end) that let the compact pre-cap schema be joined back
    to results.csv.
    """
    import json

    files = {
        CANDIDATES_FILE: {
            "rows": screened_rows,
            "content": "every candidate the off-target screen actually saw "
                       "(post-cap). The denominator for any statement about "
                       "the specificity screen.",
            "columns": list(_COLUMNS),
        },
        SUMMARY_FILE: {
            "rows": None,
            "content": "per gene: candidates that passed the filters BEFORE "
                       "the screen_top_per_gene cap vs how many the cap kept.",
            "columns": list(_SUMMARY_COLUMNS),
        },
    }
    if capture == CAPTURE_PRECAP:
        files[PRECAP_FILE] = {
            "rows": precap_rows,
            "content": "EVERY window that passed the chemistry filters, before "
                       "the cap. Compact schema: the transcript id is "
                       "dictionary-encoded (see "
                       f"{TRANSCRIPTS_FILE}) and `end` / `probe_id` are "
                       "reconstructed, not stored.",
            "columns": list(_PRECAP_COLUMNS),
        }
        files[TRANSCRIPTS_FILE] = {
            "rows": None,
            "content": "transcript dictionary for the compact schema.",
            "columns": list(_TRANSCRIPT_COLUMNS),
        }
    payload = {
        "capture": capture,
        "capture_meaning": {
            CAPTURE_PRECAP: "pre-cap: every filter-passing window was recorded, "
                            "so the null can account for the "
                            "screen_top_per_gene cap as well as the screen.",
            CAPTURE_POSTCAP: "post-cap: only the windows the cap let through to "
                             "BLAST were recorded. The cap is itself an "
                             "outcome-correlated selection step; use "
                             f"{SUMMARY_FILE} to model it.",
        }[capture],
        "screen_top_per_gene": screen_top_per_gene,
        "tile_length": tile_length,
        "reconstruction": {
            "end": "start + tile_length",
            "probe_id": "f'{transcript_id}_t{start:06d}' — the same rule the "
                        "tiler uses, so it joins to results.csv.probe_id",
            "sequence": "not stored; slice the input FASTA at [start, end)",
        },
        "files": files,
    }
    path = os.path.join(output_dir, INDEX_FILE)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
    os.replace(tmp, path)
    return path


def write_cap_summary(
    output_dir: str,
    precap_counts: Dict[str, int],
    kept: Iterable[ProbeCandidate],
    *,
    gene_key: Callable[[ProbeCandidate], str],
) -> str:
    """Write the per-gene pre-cap vs screened counts.

    ``candidates_passed_filters`` is the pre-cap count; ``candidates_screened``
    is how many of them the cap let through to BLAST. Their ratio is the
    per-gene sampling fraction a null model has to account for — without it,
    the cap looks like it never happened.
    """
    kept_counts: Dict[str, int] = {}
    for c in kept:
        g = gene_key(c)
        kept_counts[g] = kept_counts.get(g, 0) + 1

    path = os.path.join(output_dir, SUMMARY_FILE)
    tmp = path + ".tmp"
    with gzip.open(tmp, "wt", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(_SUMMARY_COLUMNS)
        for gene in sorted(set(precap_counts) | set(kept_counts)):
            passed = precap_counts.get(gene, 0)
            screened = kept_counts.get(gene, 0)
            frac = f"{screened / passed:.6f}" if passed else ""
            wr.writerow([gene, passed, screened, frac])
    os.replace(tmp, path)
    return path
