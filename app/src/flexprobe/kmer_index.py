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

"""Exhaustive k-mismatch search over a reference, by pigeonhole seeding.

This is the search engine behind the k-mer off-target screen
(``screen.KmerScreen``). It answers exactly one question, exhaustively:

    where in the reference does this query sit within ``max_mm`` mismatches?

It answers it by indexing every k-mer of the reference, enumerating the
mismatch neighbourhood of each of the query's seeds, and verifying **every**
candidate offset by a FULL-LENGTH comparison. No local alignment, no trimming,
no score threshold, nothing reconstructed.

WHY THIS EXISTS
---------------
``blastn`` reports the maximum-scoring *local* segment, but the 10x off-target
rule is stated over the whole 25 nt probe half. Everything expensive in the
BLAST path is compensation for that gap: the 9 nt pushdown floor and the large
hit volume it lets through, ``blast.SubjectWindows`` and its ``blastdbcmd`` retrieval, and the
``conservative`` recount basis that over-rejects when retrieval fails. A search
that is end-to-end by construction has none of those, because the count it
returns IS the number of mismatches over the full query.

THE GUARANTEE, AND ITS LIMIT
----------------------------
Pigeonhole: if a query of length L carries at most ``max_mm`` mismatches and is
covered by ``s`` non-overlapping seeds, then some seed carries at most
``max_mm // s``. Seeding at that per-seed budget therefore CANNOT miss a hit.
Callers pick ``seed_offsets`` and ``seed_mm`` to satisfy this; :func:`search`
does not verify it for you, and getting it wrong silently loses hits. The two
configurations the screen uses, both checked against brute force in
``tests/test_kmer_oracle.py``:

    25 nt half, max_mm=4, seeds [0, 13], seed_mm=2      (two 12 nt seeds)
    50 nt window, max_mm=8, seeds [0,12,25,37], seed_mm=2

The limit: this is **Hamming distance only**. A 1 nt indel shifts the register
and reads as many mismatches. That is a deliberate design choice, not an
oversight: the 10x off-target rule is stated as a mismatch count, and the
BLAST engine remains available for a gapped search.

RECORD CONTAINMENT
------------------
The reference is ONE array with :data:`PAD` sentinels between records, so a
reported offset has to be certified as a site before it can be named. Two
things do that, and both are needed:

1. ``query_len <= PAD`` (:func:`assert_pad_covers`) keeps a window from
   reaching across a whole separator into the record after next.
2. :class:`RecordBounds`, passed as ``bounds`` to :func:`verify`,
   :func:`search` and :func:`brute`, requires the window to lie ENTIRELY inside
   one record.

(2) is not implied by (1). Sentinels only charge a mismatch, and a mismatch is
affordable: a window starting in record i's pad with 4 sentinels and 21 bases of
record i+1 costs 4 of a budget of 4, passes, and is then attributed to record
**i** by :meth:`Reference.record_of`. Such a window is a chimaera no probe can
hybridise to, and it is discarded rather than re-attributed — it is not record
i+1's site either, since part of it is not sequence at all.
"""
from __future__ import annotations

import functools
import hashlib
import json
import logging
import os
from typing import (
    Any, Dict, Iterable, List, NamedTuple, Optional, Sequence, Tuple,
)

import numpy as np

from . import reference_qc
from .reference_qc import ReferenceComposition
from .engine_identity import current_stamp
from .errors import KmerIndexStampMismatchError

logger = logging.getLogger("flexprobe.kmer_index")

# Base -> 2-bit code. 4 means "not a base": an N, an IUPAC ambiguity code, or a
# record separator. It can never equal a query code (0-3), so a window covering
# one always counts a mismatch there — which is what makes PAD work as a
# boundary and what makes an N conservative rather than invisible.
INVALID = 4
CODE = np.full(256, INVALID, dtype=np.uint8)
for _i, _ch in enumerate("ACGT"):
    CODE[ord(_ch)] = _i
    CODE[ord(_ch.lower())] = _i

# Sentinel bases inserted between records. MUST be >= the longest query length,
# or a window could span two records and report an alignment that does not
# exist in any single sequence. 64 covers the 50 nt tile with room to spare;
# :func:`load_reference` asserts against the queries it is asked about via
# :func:`assert_pad_covers`.
#
# THE PAD ALONE IS NOT A BOUNDARY, ONLY A COST. A sentinel can never equal a
# query code, so it charges a MISMATCH — but a mismatch is a price, not a veto.
# A window that STARTS inside record i's pad and runs on into record i+1 pays
# only as many mismatches as it covers sentinels, and with a budget of 4 it can
# buy four of them and still be reported. `Reference.record_of` then maps its
# start back to record **i**, so a site made almost entirely of record i+1's
# sequence is attributed to record i. Rare, but real, and always in the
# same direction.
#
# So the pad is backed by an explicit containment filter: see
# :class:`RecordBounds` and the ``bounds`` argument of :func:`verify` /
# :func:`search`. `assert_pad_covers` remains necessary (it is what keeps a
# window from reaching PAST the next record's start) but it is not
# sufficient on its own.
PAD = 64

# Bumped whenever the on-disk index layout or the encoding changes, so an index
# written by an older build is rebuilt rather than misread. Same role as
# ``blast._DB_BUILD_RECIPE``.
INDEX_FORMAT = 1

# HARD limit, and it is a correctness limit rather than a memory one. Reference
# positions are stored as int32 — 4 bytes per indexed k-mer, which is what keeps
# a transcriptome index near 1.5 GB instead of 3 GB. A reference longer than the
# int32 range would WRAP TO NEGATIVE POSITIONS and the search would return
# silent garbage: no exception, no warning, just wrong answers. So it is refused.
#
# The limit happens to fall exactly where the intended scope does. Every
# mammalian GENOME exceeds it (human ~3.1 Gb, mouse ~2.7 Gb) and every
# transcriptome is far below it (human cDNA ~358 Mb of sequence, ~372 Mb of
# indexed positions once inter-record padding is counted), so "the reference
# is too big to index" and "the reference is a genome, which this engine is not for"
# are the same condition.
MAX_INDEXABLE_BASES = 2 ** 31 - 1

# SOFT limit: above this the index still works, but its footprint is worth
# stating before a user waits for it. Cost is roughly 5 bytes per base (4 for
# the position array, 1 for the sequence) plus a fixed 8*4**k offsets array.
LARGE_REFERENCE_WARN_BASES = 1_000_000_000

_COMPLEMENT_CODE = np.array([3, 2, 1, 0, INVALID], dtype=np.uint8)


def encode_seq(seq: str) -> np.ndarray:
    """ACGT string -> uint8 code array (non-ACGT becomes :data:`INVALID`)."""
    return CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]


def revcomp_codes(codes: np.ndarray) -> np.ndarray:
    """Reverse complement of a code array, preserving :data:`INVALID`."""
    return _COMPLEMENT_CODE[codes[::-1]]


def check_indexable(n_bases: int, *, name: str = "") -> None:
    """Refuse a reference this index cannot address correctly.

    See :data:`MAX_INDEXABLE_BASES`. Raised BEFORE any work so a user learns at
    once rather than after a long build, and worded so the fix is obvious: this
    engine screens against a transcriptome, and a genome is the wrong input for
    it, not merely a large one.
    """
    what = f" ({name})" if name else ""
    if n_bases > MAX_INDEXABLE_BASES:
        raise ValueError(
            f"reference{what} is {n_bases / 1e9:.2f} Gbases, which exceeds the "
            f"{MAX_INDEXABLE_BASES / 1e9:.2f} Gbase limit of the k-mer screen's "
            "32-bit position index.\n"
            "Every mammalian genome exceeds this and every transcriptome is far "
            "below it, so this almost always means a GENOME was supplied where "
            "a transcriptome was intended. The off-target rule is stated over "
            "transcripts, and the pipeline is transcript-oriented throughout "
            "(gene-aware mode, isoform coverage).\n"
            "Fix: screen against a transcriptome FASTA (e.g. Ensembl cdna.all "
            "or a GENCODE transcripts file), or use --screen-engine blast, "
            "which has no such limit."
        )
    if n_bases > LARGE_REFERENCE_WARN_BASES:
        logger.warning(
            "Reference%s is %.2f Gbases; the k-mer index will need roughly "
            "%.1f GB. Build once and it is reused from disk thereafter.",
            what, n_bases / 1e9, (n_bases * 5) / 1e9,
        )


def assert_pad_covers(query_len: int) -> None:
    """Guard the one assumption that would silently fabricate alignments.

    NECESSARY, NOT SUFFICIENT. This bounds how far a window can reach — with
    ``query_len <= PAD`` no window can touch two records' real sequence at once
    — but it does not stop a window from being made of sentinels plus ONE
    record's sequence, which is still not a site. :class:`RecordBounds` is what
    stops that; see the comment on :data:`PAD`.
    """
    if query_len > PAD:
        raise ValueError(
            f"query length {query_len} exceeds the record separator PAD={PAD}; "
            "a window could span two records and report an alignment present in "
            "neither. Raise PAD and rebuild the index."
        )


class RecordBounds(NamedTuple):
    """Half-open ``[start, end)`` of every record's REAL sequence in ``seq``.

    ``starts`` is :attr:`Reference.starts`; ``ends[i]`` is where record ``i``'s
    sequence stops and its pad begins. Passed to :func:`verify` / :func:`search`
    to enforce the one invariant a padded concatenation cannot enforce by
    itself: a reported hit's whole window lies inside a single record.

    Two arrays rather than a per-base validity mask deliberately. A mask costs
    one byte per base — ~370 MB over GRCh38 cdna.all's padded 372 Mbase array,
    on top of the ~1.55 GB index and the reference itself — where this costs
    16 bytes per RECORD (3.3 MB for that file's 207,249 records) and is built
    once, lazily, from arrays the parse already produced.
    """

    starts: np.ndarray
    ends: np.ndarray

    def contains_window(self, pos: np.ndarray, length: int) -> np.ndarray:
        """Boolean mask: does ``seq[p:p+length]`` lie within one record?

        Vectorised, and it must stay that way — a Python loop here would dwarf
        the search it guards. ``searchsorted`` is O(n log R) and the rest is one
        gather and one comparison, so on a transcriptome-sized boundary array it
        costs on the order of 100 ns per candidate. That is a small fraction of
        a ``search`` call, and it buys the containment invariant.

        ``pos >= starts[rec]`` holds by construction of ``searchsorted``, so only
        the end needs testing. A start inside a pad is caught by the same
        comparison: ``ends[rec]`` is then already behind ``pos``.
        """
        rec = np.searchsorted(self.starts, pos, side="right") - 1
        return np.asarray(pos + length <= self.ends[rec])


class Reference:
    """Concatenated reference sequence plus its record boundaries.

    Records are joined with :data:`PAD` sentinel bases so that a window can
    never span two of them (see :func:`assert_pad_covers`). Positions are
    offsets into :attr:`seq`; :meth:`record_of` maps them back to record indices.

    :attr:`bounds` carries the record ends as well as the starts, because
    ``record_of`` alone cannot tell a position inside record ``i``'s sequence
    from one inside record ``i``'s trailing pad — both answer ``i``. Every
    search that reports positions from THIS reference must pass
    :attr:`bounds`; see :data:`PAD`.
    """

    def __init__(self, seq: np.ndarray, names: List[str],
                 starts: np.ndarray,
                 composition: Optional[ReferenceComposition] = None,
                 pad: int = PAD) -> None:
        self.seq = seq
        self.names = names
        self.starts = starts
        # The separator actually used when `seq` was built, not the module
        # default: `load_reference` takes `pad` as a parameter, and the record
        # ends are wrong by exactly that difference if the two disagree.
        self.pad = int(pad)
        # What the FASTA contained, and what was left out of `seq` — see
        # :mod:`flexprobe.reference_qc`. Optional because a Reference can be
        # built directly from arrays (tests, ad hoc scripts) with no file
        # behind it, and an absent tally is honest where a zeroed one would be
        # a claim.
        self.composition = composition

    def __len__(self) -> int:
        return int(len(self.seq))

    @property
    def bounds(self) -> RecordBounds:
        """``[start, end)`` of every record's real sequence. Computed once.

        ``ends[i] = starts[i+1] - pad`` for every record but the last, whose
        sequence ends ``pad`` bases before the array does — `load_reference`
        appends a pad block after EVERY record, the final one included.
        """
        cached = getattr(self, "_bounds", None)
        if cached is None:
            starts = np.asarray(self.starts, dtype=np.int64)
            ends = np.empty(len(starts), dtype=np.int64)
            if len(starts):
                ends[:-1] = starts[1:] - self.pad
                ends[-1] = len(self.seq) - self.pad
            cached = RecordBounds(starts, ends)
            self._bounds = cached
        return cached

    def record_of(self, pos: np.ndarray) -> np.ndarray:
        """Index into :attr:`names` of the record containing each position.

        CAUTION: a position inside record ``i``'s trailing pad also answers
        ``i``, because the pad is on record ``i``'s side of ``starts[i+1]``.
        That is only safe for a position a search has already certified as
        lying inside real sequence — which is what :attr:`bounds` is for.
        """
        return np.searchsorted(self.starts, pos, side="right") - 1

    def name_of(self, pos: int) -> str:
        return self.names[int(self.record_of(np.asarray([pos]))[0])]


def load_reference(
    path: str, pad: int = PAD, *,
    exclude_non_primary: bool = reference_qc.EXCLUDE_ALT_SCAFFOLDS_DEFAULT,
) -> Reference:
    """Parse a FASTA into a :class:`Reference`.

    Sequence lines are concatenated per record; record ids are the full header
    line minus ``>``. Non-ACGT characters become :data:`INVALID`.

    ⚠️ **A REFERENCE N IS NOT ALWAYS CONSERVATIVE.** It is tempting to say that
    INVALID "counts as a mismatch everywhere", so an N in the reference can only
    make the screen stricter. That is true of :func:`verify`, which charges a
    mismatch. It is NOT true of the seeded :func:`search`, because
    :func:`_kmer_codes` returns ``~bad`` and :meth:`KmerIndex.__init__` indexes
    only ``np.nonzero(valid)``, so a reference k-mer containing any non-ACGT
    base **is never indexed at all**. The site is not penalised, it is
    invisible: no seed can propose that offset, whatever its mismatch count.

    With the 25 nt plan's two 12 nt seeds at offsets 0 and 13, TWO reference Ns,
    one falling in each seed window, remove every candidate offset and hide the
    site. That is 12 x 12 = 144 of the 300 possible two-N placements, while no
    single-N placement hides a site. :func:`brute` finds the sites ``search``
    misses. Position 12 lies in neither seed.

    The mirror of this on the QUERY side is also real and is handled the same
    way: :func:`search` skips any seed containing INVALID (see its ``continue``
    on an ``(seed == INVALID).any()``), and ``max_n_per_half`` permits one
    degenerate base by default. Query ambiguity that is EXPANDED before
    searching is safe; ambiguity that reaches a seed is not.

    Prevalence on a real transcriptome is **UNVERIFIED** — the mechanism is
    demonstrated and bounded, not counted in GRCh38.

    COMPOSITION. Every header is classified as it goes past — primary assembly
    vs ALT/PATCH/scaffold, and gene biotype — and the tally lands on
    :attr:`Reference.composition` and in the run manifest. This is in the SAME
    pass and only in the header branch, so per-base cost is untouched and the
    added parse time is small beside an index build. See
    :mod:`flexprobe.reference_qc` for the breakdown.

    ``exclude_non_primary`` (default ON) drops ALT haplotype and PATCH scaffold
    records from the searchable reference. They are near-verbatim copies of
    primary regions, so on GRCh38 they produce spurious perfect-50bp
    off-target calls, i.e. good probes rejected. It is a
    NO-OP unless the file proves it has a primary assembly to keep, so a viral
    panel or any other non-Ensembl FASTA is unaffected and unmentioned. See
    :mod:`flexprobe.reference_qc`.
    """
    names: List[str] = []
    # Per record, parallel to `names`: the encoded sequence, its assembly class
    # and its biotype. Held per record rather than concatenated on the fly
    # because the keep/drop decision is only sound once the WHOLE file has been
    # seen (a file with no primary records must not be emptied).
    records: List[np.ndarray] = []
    assemblies: List[int] = []
    biotypes: List[str] = []
    cur: List[str] = []
    comp = reference_qc.ReferenceComposition(path)

    def flush() -> None:
        if not names:
            return
        arr = encode_seq("".join(cur))
        records.append(arr)
        asm = assemblies[-1]
        # Re-parsed only for the small minority of records that are not on a
        # primary chromosome (~7.6% of GRCh38 cdna.all), purely to report
        # distinct scaffolds per class and which region types were seen.
        region_type, region = (
            reference_qc.location_of(names[-1])
            if asm != reference_qc.PRIMARY else ("", ""))
        comp.add(asm, biotypes[-1], len(arr), region, region_type)

    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                cur.clear()
                header = line[1:].strip()
                names.append(header)
                asm, biotype = reference_qc.classify_header(header)
                assemblies.append(asm)
                biotypes.append(biotype)
            else:
                cur.append(line.strip())
    flush()
    if not names:
        raise ValueError(f"no FASTA records found in {path}")

    dropping = reference_qc.should_exclude(comp, exclude_non_primary)
    comp.finalise(excluded_non_primary=dropping)
    reference_qc.register(comp)
    # Before anything is dropped: a whole reference landing outside PRIMARY is
    # a parsing failure, not a reference without chromosomes.
    reference_qc.warn_if_no_primary(comp)

    kept: List[str] = []
    starts: List[int] = []
    chunks: List[np.ndarray] = []
    total = 0
    padblock = np.full(pad, INVALID, dtype=np.uint8)
    for name, asm, arr in zip(names, assemblies, records):
        if dropping and asm in reference_qc.EXCLUDABLE:
            continue
        kept.append(name)
        starts.append(total)
        chunks.append(arr)
        chunks.append(padblock)
        total += len(arr) + pad
    # Release the dropped records BEFORE the concatenate allocates the final
    # array, so excluding never costs peak memory (the kept ones stay alive
    # through `chunks`, and none of this copies).
    records.clear()

    seq = np.concatenate(chunks) if chunks else np.empty(0, dtype=np.uint8)
    logger.info("Loaded reference %s: %.1f Mbases, %d records",
                os.path.basename(path), len(seq) / 1e6, len(kept))
    if dropping:
        # One line, and only when something was actually dropped: it changes
        # which off-targets exist, so it belongs in the log next to the load.
        logger.info(
            "Excluded %d ALT/PATCH records (%.1f Mbases, %.1f%% of %s) from "
            "the off-target reference; %d records searched, of which %d are "
            "unplaced/unlocalised scaffolds that are RETAINED because they "
            "duplicate no primary sequence.",
            comp.records[reference_qc.ALT_PATCH],
            comp.bases[reference_qc.ALT_PATCH] / 1e6,
            100.0 * comp.alt_patch_base_fraction,
            os.path.basename(path), comp.records_searched,
            comp.records[reference_qc.UNPLACED],
        )
    else:
        reference_qc.warn_if_searching_alt(comp)
    return Reference(seq, kept, np.asarray(starts, dtype=np.int64), comp,
                     pad=pad)


def _kmer_codes(seq: np.ndarray, k: int) -> Tuple[np.ndarray, np.ndarray]:
    """``(codes, valid)`` for every k-mer start position in ``seq``."""
    n = len(seq) - k + 1
    if n <= 0:
        return (np.empty(0, dtype=np.int64), np.empty(0, dtype=bool))
    dt: Any = np.int32 if 4 ** k < 2 ** 31 else np.int64
    codes = np.zeros(n, dtype=dt)
    bad = np.zeros(n, dtype=bool)
    for i in range(k):
        col = seq[i:i + n]
        codes *= 4
        codes += np.where(col < INVALID, col, 0)
        bad |= (col == INVALID)
    return codes, ~bad


def code_of(seed: np.ndarray) -> int:
    """Integer code of a k-mer given as base codes. Assumes all are valid."""
    c = 0
    for b in seed:
        c = c * 4 + int(b)
    return c


def neighbourhood(seed: np.ndarray, max_mm: int) -> np.ndarray:
    """Every k-mer code within ``max_mm`` substitutions of ``seed``.

    Sizes are small and fixed: for k=12, 1 / 37 / 631 codes at max_mm 0 / 1 / 2.
    Computed arithmetically from the base code (changing position j from b to b'
    shifts the code by ``(b'-b) * 4**(k-1-j)``) rather than by materialising
    strings.
    """
    if max_mm < 0:
        raise ValueError("max_mm must be >= 0")
    if max_mm > 2:
        # Deliberate: the enumeration below is written out to two levels because
        # that is what the screen needs, and a general recursive version would be
        # both slower and an invitation to call it at a budget whose candidate
        # volume is not viable (a 3-mismatch 12-mer neighbourhood is 6,571 codes,
        # roughly ten times the candidates of the 631-code 2-mismatch one).
        raise NotImplementedError(
            f"neighbourhood() supports max_mm <= 2, got {max_mm}"
        )
    k = len(seed)
    base = code_of(seed)
    place = [4 ** (k - 1 - j) for j in range(k)]
    out: List[int] = [base]
    if max_mm >= 1:
        for j in range(k):
            for b in range(4):
                if b != seed[j]:
                    out.append(base + (b - int(seed[j])) * place[j])
    if max_mm >= 2:
        for j1 in range(k):
            for b1 in range(4):
                if b1 == seed[j1]:
                    continue
                d1 = base + (b1 - int(seed[j1])) * place[j1]
                for j2 in range(j1 + 1, k):
                    for b2 in range(4):
                        if b2 != seed[j2]:
                            out.append(d1 + (b2 - int(seed[j2])) * place[j2])
    return np.asarray(out, dtype=np.int64)


def _sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@functools.lru_cache(maxsize=1)
def builder_sha256() -> str:
    """SHA-256 of this module's source: the code that BUILDS an index.

    This file only, not the whole engine, so a ranking or filter change does
    not force a rebuild. ``INDEX_FORMAT`` tracks only the on-disk layout, so
    without this an index built by different builder code would be reused
    silently. A build that ships no ``.py`` source returns ``"unavailable"``.
    """
    try:
        with open(__file__, "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()
    except OSError:
        return "unavailable"


def _miss_reason(meta: Dict[str, Any], sha: str, k: int,
                 primary_only: bool) -> str:
    """Which cache-key field disqualified an index that IS on disk.

    Reported to the operator, not consumed by anything, so it names the first
    mismatch rather than enumerating them all. ``primary_only`` is checked last
    because the other three are the ordinary "the reference changed" cases,
    while a primary_only mismatch means the wrong VARIANT was staged and wants
    the louder message the caller then emits.
    """
    if meta.get("source_sha256") != sha:
        return "the reference has changed since the index was built"
    if meta.get("format") != INDEX_FORMAT:
        return (f"index layout v{meta.get('format')!r}, "
                f"this build reads v{INDEX_FORMAT}")
    if int(meta.get("k", -1)) != k:
        return f"index is k={meta.get('k')!r}, this run needs k={k}"
    if bool(meta.get("primary_only", False)) != primary_only:
        return ("the index on disk is "
                f"{'primary-only' if not primary_only else 'full-file'}, "
                f"this run needs "
                f"{'primary-only' if primary_only else 'full-file'}")
    return "unknown"


class KmerIndex:
    """Exact k-mer index over a reference: bucket code -> sorted positions.

    Stored CSR-style — one sorted ``pos`` array plus an ``offsets`` array of
    ``4**k + 1`` entries — so a lookup is a slice and costs no allocation.

    MEMORY. For k=12 over a 357.6 Mbase transcriptome -- 372 M indexed
    positions, since padding is addressed too -- the resident cost is
    ~1.55 GB (1.42 GB positions + 134 MB offsets), plus the reference itself.
    The BUILD peak is higher (~7 GB, dominated by the argsort), which is why
    :meth:`build_or_load` persists the result and :meth:`load` memory-maps it:
    an index is built once per reference and then costs page cache, not heap.
    """

    def __init__(self, seq: np.ndarray, k: int) -> None:
        check_indexable(len(seq))
        self.k = k
        codes, valid = _kmer_codes(seq, k)
        pos = np.nonzero(valid)[0]
        c = codes[pos]
        order = np.argsort(c, kind="stable")
        self.pos: np.ndarray = pos[order].astype(np.int32)
        counts = np.bincount(c, minlength=4 ** k)
        self.offsets: np.ndarray = np.zeros(4 ** k + 1, dtype=np.int64)
        np.cumsum(counts, out=self.offsets[1:])
        del codes, valid, pos, c, order, counts

    # -- lookup ------------------------------------------------------------

    def lookup(self, code: int) -> np.ndarray:
        """Reference positions of one exact k-mer code."""
        return self.pos[self.offsets[code]:self.offsets[code + 1]]

    def lookup_many(self, codes: np.ndarray) -> np.ndarray:
        """Positions for many bucket codes at once, concatenated.

        Vectorised rather than looped: a 2-mismatch 12-mer neighbourhood is 631
        buckets, and the screen issues four of those per probe half, so per-call
        Python overhead is the thing that matters here.
        """
        if len(codes) == 0:
            return np.empty(0, dtype=np.int32)
        lo = self.offsets[codes]
        n = self.offsets[codes + 1] - lo
        keep = n > 0
        if not keep.any():
            return np.empty(0, dtype=np.int32)
        lo, n = lo[keep], n[keep]
        total = int(n.sum())
        # Build a flat gather index: for each bucket, its start offset repeated
        # `n` times, plus a running 0..n-1 within the bucket.
        starts = np.repeat(lo - np.concatenate([[0], np.cumsum(n)[:-1]]), n)
        return self.pos[starts + np.arange(total, dtype=np.int64)]

    # -- persistence -------------------------------------------------------

    def save(self, prefix: str, *, source_sha: str = "",
             primary_only: bool = False) -> None:
        """Write the index beside ``prefix`` for later :meth:`load`."""
        os.makedirs(os.path.dirname(os.path.abspath(prefix)) or ".",
                    exist_ok=True)
        np.save(prefix + ".pos.npy", self.pos)
        np.save(prefix + ".off.npy", self.offsets)
        with open(prefix + ".meta.json", "w") as fh:
            json.dump({"k": self.k, "format": INDEX_FORMAT,
                       "n_positions": int(len(self.pos)),
                       # Belt-and-braces against the prefix: an index built
                       # over primary records only describes a DIFFERENT
                       # sequence than the same FASTA's full index, and the
                       # difference is invisible in the position array itself.
                       "primary_only": bool(primary_only),
                       "source_sha256": source_sha,
                       # Compared on reuse; see build_or_load.
                       "builder_sha256": builder_sha256(),
                       # For the record only. Not compared, so an engine
                       # change outside this file does not force a rebuild.
                       "engine_stamp": current_stamp()}, fh)

    @classmethod
    def load(cls, prefix: str, *, mmap: bool = True) -> "KmerIndex":
        """Load a saved index, memory-mapping the position array by default.

        ``mmap`` is the reason persistence exists: the positions array is by far
        the largest component, and mapping it keeps the process's resident set
        to the offsets array while the OS pages in only the buckets actually
        touched.
        """
        with open(prefix + ".meta.json") as fh:
            meta = json.load(fh)
        if meta.get("format") != INDEX_FORMAT:
            raise ValueError(
                f"index at {prefix} is format {meta.get('format')}, this build "
                f"reads {INDEX_FORMAT}; rebuild it"
            )
        obj = cls.__new__(cls)
        obj.k = int(meta["k"])
        obj.pos = np.load(prefix + ".pos.npy",
                          mmap_mode="r" if mmap else None)
        obj.offsets = np.load(prefix + ".off.npy")
        return obj

    @staticmethod
    def cache_prefix(fasta_path: str, k: int, *,
                     cache_dir: Optional[str] = None,
                     fallback_dir: Optional[str] = None,
                     primary_only: bool = False) -> str:
        """Where this reference's index lives.

        BESIDE THE FASTA by default, and that is the whole point of keying it on
        the source SHA-256: an index is a property of the REFERENCE, not of the
        run that happened to need it, so every project screening against the
        same transcriptome should share one. It is also exactly where
        ``makeblastdb`` puts a BLAST database (``blast.ensure_db_for_path``), so
        the two engines keep their derived data in the same place.

        Scoping it to an output directory instead would silently defeat the
        cache: every run would rebuild from scratch and write its own ~1.5 GB
        copy, paying minutes of build time where a cache hit costs seconds.

        ``fallback_dir`` is used only when the reference's own directory is not
        writable (a read-only reference mount is the case worth handling), so
        the run still works without silently losing reuse everywhere else.

        ``primary_only`` is part of the NAME, not merely of the metadata,
        because an index built after ALT/PATCH records were removed has
        DIFFERENT positions for the same FASTA: reusing one for the other would
        not fail, it would return hits at wrong coordinates in wrong records.
        The same argument applies to a change in WHICH records are removed,
        which is why the suffix is ``.noaltpatch`` rather than ``.primary``.
        Two names means the two variants coexist and each is reused, and — since
        the un-suffixed name is unchanged — every index already on disk stays
        valid for the runs that ask for the same thing it holds.
        """
        base = os.path.splitext(os.path.basename(fasta_path))[0]
        if primary_only:
            # ".noaltpatch", NOT ".primary". The suffix names WHAT WAS
            # REMOVED. An index named ".primary" was built with
            # unplaced/unlocalised scaffolds removed as well, so reusing one
            # here would search a reference missing records this build keeps,
            # at coordinates that no longer match. The distinct name means such
            # files are never looked up, and can be deleted.
            base += ".noaltpatch"
        if cache_dir:
            out_dir = cache_dir
        else:
            beside = os.path.dirname(os.path.abspath(fasta_path))
            if os.access(beside, os.W_OK):
                out_dir = beside
            elif fallback_dir:
                logger.info(
                    "Reference directory %s is not writable; caching the k-mer "
                    "index under %s instead. It will not be reused by other "
                    "runs — make the reference directory writable to share it.",
                    beside, fallback_dir,
                )
                out_dir = fallback_dir
            else:
                out_dir = beside          # will fail loudly on save
        os.makedirs(out_dir, exist_ok=True)
        return os.path.join(out_dir, f"{base}.k{k}")

    @classmethod
    def build_or_load(
        cls, fasta_path: str, k: int, *,
        cache_dir: Optional[str] = None, mmap: bool = True,
        fallback_dir: Optional[str] = None,
        exclude_non_primary: bool = (
            reference_qc.EXCLUDE_ALT_SCAFFOLDS_DEFAULT),
        rebuild: bool = False,
    ) -> Tuple["KmerIndex", Reference]:
        """Return ``(index, reference)``, reusing a cached index when valid.

        The cache is keyed on the source FASTA's SHA-256 plus
        :data:`INDEX_FORMAT`, so editing the reference or changing the layout
        rebuilds rather than silently misreads — the same contract as
        ``blast.build_db``'s ``.sha256`` sidecar. See :meth:`cache_prefix` for
        WHERE it goes and why that matters.

        ``exclude_non_primary`` participates in the cache IDENTITY, not just in
        the load: the FASTA's hash is unchanged by dropping records from it, so
        a cache keyed on the hash alone would hand a full-file index to a run
        that asked for primary-only and silently report hits at positions that
        belong to other records. It is keyed on what the exclusion ACTUALLY did
        (``ref.composition``), not on what was requested, so the no-op case —
        every reference without Ensembl location fields — keeps using and
        producing exactly the index it always did.

        An index that matches on all of the above must ALSO carry this
        build's :func:`builder_sha256`. A missing or different one raises
        :class:`~flexprobe.errors.KmerIndexStampMismatchError` and deletes
        nothing, unless ``rebuild`` (CLI ``--rebuild-kmer-index``) asks for a
        rebuild. A changed reference always rebuilds.
        """
        ref = load_reference(fasta_path,
                             exclude_non_primary=exclude_non_primary)
        primary_only = bool(ref.composition
                            and ref.composition.excluded_non_primary)
        prefix = cls.cache_prefix(fasta_path, k, cache_dir=cache_dir,
                                  fallback_dir=fallback_dir,
                                  primary_only=primary_only)
        sha = _sha256_file(fasta_path)
        if ref.composition is not None:
            # Already paid for by the cache key, so the manifest gets the
            # reference's hash for free — and a k-mer run has no BLAST database
            # sidecars for `_blast_db_sha256s` to digest, so without this its
            # reference is recorded by PATH alone.
            ref.composition.source_sha256 = sha
        meta_path = prefix + ".meta.json"
        reason = "no index beside the reference"
        if os.path.exists(meta_path):
            try:
                with open(meta_path) as fh:
                    meta = json.load(fh)
                if (meta.get("source_sha256") == sha
                        and meta.get("format") == INDEX_FORMAT
                        and int(meta.get("k", -1)) == k
                        # Absent from older index metadata, which always
                        # described a full-file index, so the default is False
                        # and those stay reusable by the runs they are valid for.
                        and bool(meta.get("primary_only", False)) == primary_only):
                    cached_builder = meta.get("builder_sha256")
                    if cached_builder == builder_sha256():
                        logger.info("k-mer index cache hit for %s", prefix)
                        return cls.load(prefix, mmap=mmap), ref
                    if not rebuild:
                        raise KmerIndexStampMismatchError(
                            prefix, cached_builder, builder_sha256())
                    reason = ("rebuild requested for an index built by "
                              "different k-mer index code")
                else:
                    reason = _miss_reason(meta, sha, k, primary_only)
            except (OSError, ValueError, json.JSONDecodeError):
                logger.debug("Unusable k-mer index at %s; rebuilding",
                             prefix, exc_info=True)
                reason = "index metadata unreadable"
        # Staging a prebuilt index next to a reference is a supported workflow
        # (shipping one on a USB stick, seeding a new machine, an HPC bundle),
        # and the commonest way to get it wrong is to stage the OTHER variant:
        # `--exclude-alt-scaffolds` is ON by default, so a default run against
        # GRCh38 cdna.all wants `.noaltpatch.k12.*` and walks straight past a
        # plain `.k12.*` sitting right beside it. If that were silent, the only
        # symptom would be a run taking minutes instead of seconds, which reads
        # as slowness rather than as a misconfiguration. Say which prefix was
        # wanted, and say it at WARNING when the wrong variant is what is there.
        sibling = cls.cache_prefix(fasta_path, k, cache_dir=cache_dir,
                                   fallback_dir=fallback_dir,
                                   primary_only=not primary_only)
        if os.path.exists(sibling + ".meta.json"):
            logger.warning(
                "k-mer index MISS (%s). This run needs '%s.*' but the only "
                "index present is '%s.*', which it cannot use. Rebuilding, "
                "which costs roughly 5 min and ~7 GB peak RAM. To reuse a "
                "staged index, ship the '%s' variant, or re-run with %s.",
                reason,
                os.path.basename(prefix), os.path.basename(sibling),
                "primary-only" if primary_only else "full-file",
                "--no-exclude-alt-scaffolds" if primary_only
                else "--exclude-alt-scaffolds")
        else:
            logger.info("k-mer index miss (%s); wanted '%s.*'",
                        reason, os.path.basename(prefix))
        logger.info("Building %d-mer index for %s", k,
                    os.path.basename(fasta_path))
        idx = cls(ref.seq, k)
        idx.save(prefix, source_sha=sha, primary_only=primary_only)
        return idx, ref


def verify(
    seq: np.ndarray, starts: np.ndarray, query: np.ndarray, max_mm: int,
    *, bounds: Optional[RecordBounds] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """``(starts, mismatches)`` for candidates within ``max_mm``, full-length.

    This is where the answer is actually decided. Everything upstream only
    proposes offsets; every one of them is compared across all ``len(query)``
    positions here, so a hit's mismatch count is a count and not an inference.

    ``bounds`` — pass it whenever ``seq`` is a padded concatenation of records,
    which means every :class:`Reference`. It is the ONLY thing that stops a
    window made of separator sentinels plus a neighbour's sequence from being
    reported as a site and attributed to the preceding record; a mismatch
    budget cannot do it, because sentinels only cost mismatches (see
    :data:`PAD`). Omitting it says "this array has no record structure", which
    is true of the bare arrays the oracle tests search and of nothing else.

    THE FILTER RUNS BEFORE THE WINDOW GATHER, so a discarded candidate costs a
    comparison rather than ``len(query)`` byte loads.
    """
    L = len(query)
    empty = (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int16))
    if len(starts) == 0:
        return empty
    cand = np.unique(starts.astype(np.int64))
    cand = cand[(cand >= 0) & (cand + L <= len(seq))]
    if bounds is not None and len(cand):
        cand = cand[bounds.contains_window(cand, L)]
    if len(cand) == 0:
        return empty
    win = seq[cand[:, None] + np.arange(L)]
    mm = (win != query[None, :]).sum(axis=1).astype(np.int16)
    keep = mm <= max_mm
    return cand[keep], mm[keep]


def search(
    seq: np.ndarray, index: KmerIndex, query: np.ndarray, max_mm: int,
    seed_offsets: Sequence[int], seed_mm: int,
    *, bounds: Optional[RecordBounds] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """All ``(start, mismatches)`` where ``query`` sits within ``max_mm``.

    ``seed_offsets`` must be non-overlapping k-length windows of the query whose
    count ``s`` satisfies ``max_mm // s <= seed_mm`` — otherwise the pigeonhole
    guarantee does not hold and hits are silently lost. See the module docstring
    for the two configurations the screen uses.

    ``bounds`` is forwarded to :func:`verify` and is REQUIRED for correctness on
    any :class:`Reference`. THIS IS THE LAYER THE CONTAINMENT RULE BELONGS IN:
    the index proposes a seed position that is always inside a record (a k-mer
    containing a sentinel is never indexed at all), and it is this function's
    ``hits - off`` shift that can walk a candidate backwards into the preceding
    pad. The defect is created here, so it is refused here, once, for every
    caller — rather than in each place a position is turned into a record name.
    """
    k = index.k
    cands: List[np.ndarray] = []
    for off in seed_offsets:
        if off < 0 or off + k > len(query):
            raise ValueError(
                f"seed offset {off} with k={k} does not fit a query of length "
                f"{len(query)}"
            )
        seed = query[off:off + k]
        if (seed == INVALID).any():
            continue          # an N in the seed: this seed proposes nothing
        hits = index.lookup_many(neighbourhood(seed, seed_mm))
        if len(hits):
            cands.append(hits.astype(np.int64) - off)
    if not cands:
        return (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int16))
    return verify(seq, np.concatenate(cands), query, max_mm, bounds=bounds)


def brute(
    seq: np.ndarray, query: np.ndarray, max_mm: int,
    chunk: int = 8_000_000, *, bounds: Optional[RecordBounds] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Ground truth by definition: every offset, compared at full length.

    ``bounds`` carries the SAME meaning as in :func:`verify`, and for the same
    reason: over a padded concatenation the naive definition — "every offset" —
    includes offsets that are not sites at all, so brute force without it is
    not the ground truth for a :class:`Reference`, only for a bare array. An
    oracle comparison against a padded reference must pass it to BOTH sides.

    Slow and unusable in production. It exists so :func:`search` — which is a
    METHOD — can be checked against the DEFINITION (``tests/test_kmer_oracle.py``),
    and so engine disagreements can be adjudicated by something neither engine
    had a hand in.

    CHUNKED, because the obvious one-liner is not merely slow, it is fatal: the
    comparison window is an ``(n, L)`` array, and against a 357.6 Mbase
    transcriptome (372 M indexed positions, which is what the 9.3 GB below is
    computed over) at L=25 that is 9.3 GB for a single query. Adjudicating a
    handful of disputed probes would OOM the machine before answering anything.
    ``chunk`` bounds it to a few hundred MB at no cost to the result.
    """
    L = len(query)
    n = len(seq) - L + 1
    if n <= 0:
        return (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int16))
    offs = np.arange(L)
    pos_out: List[np.ndarray] = []
    mm_out: List[np.ndarray] = []
    for lo in range(0, n, chunk):
        hi = min(lo + chunk, n)
        win = seq[np.arange(lo, hi)[:, None] + offs]
        mm = (win != query[None, :]).sum(axis=1)
        keep = np.nonzero(mm <= max_mm)[0]
        if len(keep) and bounds is not None:
            keep = keep[bounds.contains_window(keep + lo, L)]
        if len(keep):
            pos_out.append(keep + lo)
            mm_out.append(mm[keep].astype(np.int16))
    if not pos_out:
        return (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int16))
    return np.concatenate(pos_out), np.concatenate(mm_out)


# Seeding configurations used by the off-target screen. Named here so callers
# cannot get the pigeonhole arithmetic wrong by hand; both are pinned against
# brute force in tests/test_kmer_oracle.py.
#   name -> (query_len, max_mm, seed_offsets, seed_mm)
SEED_PLANS: Dict[str, Tuple[int, int, Tuple[int, ...], int]] = {
    "half25_mm4": (25, 4, (0, 13), 2),
    # The near-miss budget. 5 // 2 == 2, so the same two seeds cover it for
    # free — see KmerScreen.MAX_MISMATCH_REPORTED for why one band above the
    # disqualifying threshold is worth having.
    "half25_mm5": (25, 5, (0, 13), 2),
    "window50_mm8": (50, 8, (0, 12, 25, 37), 2),
}


def plan_for(query_len: int, max_mm: int) -> Tuple[Tuple[int, ...], int]:
    """``(seed_offsets, seed_mm)`` for a known configuration."""
    for qlen, mm, offs, smm in SEED_PLANS.values():
        if qlen == query_len and mm == max_mm:
            return offs, smm
    raise KeyError(
        f"no seeding plan for query_len={query_len} max_mm={max_mm}; add one to "
        "SEED_PLANS and pin it against brute force in tests/test_kmer_oracle.py"
    )


def iter_seed_plans() -> Iterable[Tuple[str, int, int, Tuple[int, ...], int]]:
    """``(name, query_len, max_mm, seed_offsets, seed_mm)`` for every plan."""
    for name, (qlen, mm, offs, smm) in SEED_PLANS.items():
        yield name, qlen, mm, offs, smm
