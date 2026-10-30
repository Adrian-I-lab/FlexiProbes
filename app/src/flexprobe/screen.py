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

"""Off-target search engines behind a common interface.

A :class:`Screen` answers one question for a chunk of candidates: *where does
each probe half match the reference, closely enough to matter?* It owns the
search and nothing else — chunking, resume caching, streaming, progress and
cancellation all live in ``pipeline._run_screen_search`` and are shared by every
engine.

WHY THE SEAM IS HERE
--------------------
``blastn`` reports a maximum-scoring *local* segment, but the 10x off-target
rule is stated over the whole 25 nt probe half. Everything expensive in the
BLAST path is compensation for that mismatch — the 9 nt pushdown floor and the
large hit volume it lets through, ``blast.SubjectWindows`` and its ``blastdbcmd`` retrieval, and
the ``conservative`` recount basis that over-rejects when retrieval fails. An
engine that searches end-to-end has none of those. Putting the seam at "produce
hits for a chunk" lets a second engine replace all of it without touching the
specificity rules, the cache, or the pipeline's control flow.

THE HIT CONTRACT
----------------
Every engine emits the same dicts ``specificity`` already consumes, keyed by
probe id: ``qseqid`` (``<probe_id>__L`` / ``__R``), ``sseqid``, ``pident``,
``length``, ``mismatch``, ``qstart``/``qend``, ``sstart``/``send``, ``evalue``,
``bitscore``, ``qseq``, ``sseq``. An engine that aligns end-to-end sets
``qstart=1``/``qend=half_len`` and supplies the full half as ``qseq``, which is
what makes ``specificity.full_half_mismatch_detail`` return ``RECOUNT_EXACT``
with no subject retrieval at all.

A probe's list may also carry ONE junction-duplex record, which is not an
alignment: see ``KmerScreen._junction_scan`` and
``specificity.split_junction_records``. Engines are not required to produce one,
and ``specificity`` reads its absence as "not measured" rather than as zero.
"""
from __future__ import annotations

import inspect
import itertools
import logging
import os
from typing import Any, Dict, List, Optional, Protocol, Sequence, Tuple

import numpy as np

from . import blast
from . import reference_qc
from .constants import (
    BLAST_MAX_TARGET_SEQS, JUNCTION_DUPLEX_MAX_SITES,
    REFERENCE_STRAND_BOTH, REFERENCE_STRAND_PLUS,
)
from .models import ProbeCandidate
# The junction-duplex measurement lives with the rest of the specificity rules,
# not here: this module owns the SEARCH, and what counts as a ligatable duplex
# is a rule. Imported rather than reimplemented so there is exactly one
# definition of it (specificity does not import screen, so this is not a cycle).
from .specificity import (
    JUNCTION_DUPLEX_KEY, UNSCREENABLE_HALF_KEY,
    gene_scoping_degraded,
    junction_clean_duplex as _junction_clean_duplex,
)

logger = logging.getLogger("flexprobe.screen")

ENGINE_BLAST = "blast"
ENGINE_KMER = "kmer"
ENGINES = (ENGINE_BLAST, ENGINE_KMER)

#: What each engine is CALLED in progress events, logs and the timing records.
#
# The names follow the engine actually in use. A log that says "BLAST chunk"
# for a screen that never invoked blastn is a label a reader will repeat in a
# methods section, so the k-mer engine gets its own wording.
#
# ``stage`` doubles as the key in ``elapsed_seconds`` / ``run_timings.json``, so
# a run says which screen it paid for. On-disk names (``.screen_cache/``,
# ``candidates_prescreen.csv.gz``, ``config.screen``) are deliberately NOT renamed:
# they are format contracts that resume and import depend on, and a cosmetic
# rename there would strand every existing run.
SCREEN_LABELS: Dict[str, Dict[str, str]] = {
    ENGINE_BLAST: {
        "stage": ENGINE_BLAST,
        "chunk": "BLAST chunk",
        "noun": "BLAST",
        "target": "DB(s)",
        "search": "blastn-short against {n} DB(s) — may be slow",
    },
    ENGINE_KMER: {
        "stage": ENGINE_KMER,
        "chunk": "k-mer chunk",
        "noun": "k-mer screen",
        "target": "reference(s)",
        "search": "k-mer screen against {n} reference(s)",
    },
}


def screen_labels(engine: str) -> Dict[str, str]:
    """Progress/log vocabulary for ``engine``; falls back to BLAST's wording.

    An unknown engine is a caller error that `make_screen` reports properly, so
    this only has to avoid raising a second, less useful exception first.
    """
    return SCREEN_LABELS.get(engine, SCREEN_LABELS[ENGINE_BLAST])

Hit = Dict[str, Any]
HitsByProbe = Dict[str, List[Hit]]

#: Query-id suffix carried by a junction-duplex record. Deliberately NOT
#: ``__L``/``__R``: ``specificity.hit_half`` must not attribute it to a half,
#: because an unattributable hit is charged to BOTH halves and this is not a hit
#: at all. The record is identified by :data:`specificity.JUNCTION_DUPLEX_KEY`,
#: never by parsing this suffix — the suffix is there so a record is legible in
#: a dumped cache file.
JUNCTION_QUERY_SUFFIX = "__J"

#: Query-id suffix carried by an unscreenable-half record. Same reasoning as
#: :data:`JUNCTION_QUERY_SUFFIX`: deliberately NOT ``__L``/``__R``, because
#: ``specificity.hit_half`` must not attribute it to a half — it is not a hit.
UNSCREENABLE_QUERY_SUFFIX = "__U"

#: Most concrete oligos one degenerate probe half may be expanded into before
#: the screen refuses to search it at all.
#:
#: 16 is chosen against the filter that decides which halves get here.
#: ``FilterThresholds.max_n_per_half`` defaults to 1, so a shipped-defaults run
#: never exceeds 4 (one 'N'), and 16 still covers a half carrying two 'N's or
#: four two-fold codes. Past that the enumeration stops being a search of the
#: probe the caller ordered — it is a search of up to a thousand oligos none of
#: which will be synthesised — and it costs 4x more per extra ambiguous base.
#: Above the cap the half is reported as unscreenable and the pair is REJECTED;
#: see ``specificity.STATUS_UNSCREENABLE_HALF``.
MAX_HALF_RESOLUTIONS = 16


def expand_degenerate_half(
    seq: str, cap: int = MAX_HALF_RESOLUTIONS,
) -> Optional[List[str]]:
    """Concrete ACGT oligos a probe half stands for, or None if it is too many.

    WHY THIS EXISTS. A k-mer seed containing a non-ACGT base proposes nothing,
    so a screen that searched the half as written would effectively SKIP any
    half carrying one, and the pair would keep whatever verdict its other half
    earned. With ``max_n_per_half`` defaulting to 1 that is reachable on a
    default run: a probe with one 'N' in its LHS could be certified specific
    while a perfect off-target sat under that LHS, unsearched. Expanding the
    code to its resolutions is what keeps the guarantee that a hit is never
    lost by seeding: every resolution is searched exhaustively, so no site is
    lost.

    THE UNION IS THE SAFE DIRECTION. The half's hits are the union over
    resolutions, keeping the LOWEST mismatch count found at each site. An 'N'
    means "this position is unknown", and the strictest honest reading of an
    unknown base opposite a reference base is that it could pair — so the site
    is scored at its best-case duplex and the screen can only reject more, never
    less. That is the same direction every other estimate in the screen errs in
    (see ``specificity.off_target_subject_key``, ``full_half_mismatches``).

    PER-POSITION, NOT PER-RECORD. ``consensus.resolve_degenerate_window`` fans a
    consensus window out into one oligo per input STRAIN, because that is what
    will be synthesised. This is a different question: here we do not know which
    resolution the ordered oligo will carry (the code may not even come from a
    consensus), so the cartesian product over the codes is the right enumeration
    — it is a superset of the per-record resolutions, and a superset is again
    the safe direction for a screen.

    Returns None when the product would exceed ``cap``, or when a character is
    not an IUPAC code at all (a 'U', a stray letter) — in both cases the caller
    must fail closed rather than search a subset. A pure-ACGT half returns
    ``[seq]``, so callers may use this unconditionally if they prefer.
    """
    from .consensus import IUPAC_BASES

    choices: List[str] = []
    total = 1
    for ch in seq.upper():
        bases = IUPAC_BASES.get(ch)
        if bases is None:
            return None
        opts = "".join(sorted(bases))
        choices.append(opts)
        total *= len(opts)
        if total > cap:
            return None
    return ["".join(t) for t in itertools.product(*choices)]


def _unscreenable_record(probe_id: str, qid: str, qseq: str) -> Hit:
    """The record that says "this half was not searched", for the hit list.

    Not an alignment and not evidence: ``specificity.split_screen_records``
    takes it back out, and ``evaluate_specificity`` turns it into a rejection.
    """
    return {
        UNSCREENABLE_HALF_KEY: True,
        "qseqid": f"{probe_id}{UNSCREENABLE_QUERY_SUFFIX}",
        # Which half, for a human reading a dumped cache file. Never parsed as
        # an attribution — the pair is rejected outright either way.
        "half": qid.rsplit("__", 1)[-1],
        "qseq": qseq,
    }


def _best_hit_per_site(hits: List[Hit]) -> List[Hit]:
    """Collapse hits from several resolutions of one half onto their sites.

    Two resolutions of the same degenerate half land on the same reference
    position with different mismatch counts. That is one off-target site, not
    two, and the count that matters is the smallest — the best duplex the
    ordered oligo could form there. Insertion order is preserved, so the output
    is deterministic given a deterministic resolution order.
    """
    best: Dict[Tuple[str, int, int], Hit] = {}
    for h in hits:
        key = (h["sseqid"], int(h["sstart"]), int(h["send"]))
        cur = best.get(key)
        if cur is None or h["mismatch"] < cur["mismatch"]:
            best[key] = h
    return list(best.values())


def half_queries(cand: ProbeCandidate) -> List[tuple]:
    """``[(query_id, sequence), ...]`` — the two 25 nt halves of a candidate.

    The BLAST screen queries each half separately rather than the 50 bp tile:
    on-/off-target hits are judged per 25-mer,
    and the 10x rule is a disjunction over the two halves. The ``__L``/``__R``
    suffix is the ONLY record of which half a hit came from — see
    ``specificity.hit_half`` — so every engine must use it.
    """
    seq = cand.sequence
    half = len(seq) // 2
    return [(f"{cand.probe_id}__L", seq[:half]),
            (f"{cand.probe_id}__R", seq[half:])]


class Screen(Protocol):
    """Produce off-target hits for chunks of candidates.

    Lifecycle: ``prepare()`` once, ``search_chunk()`` per chunk, ``close()``
    once. ``fingerprint_parts()`` is folded into the resume-cache digest, so an
    engine MUST surface every parameter that changes which hits come back — a
    cache keyed without one is a cache that silently returns another run's
    answers.
    """

    name: str

    def prepare(self, db_prefixes: Sequence[str]) -> None: ...

    def fingerprint_parts(self) -> Dict[str, Any]: ...

    def search_chunk(
        self, chunk: Sequence[ProbeCandidate], chunk_index: int,
        output_dir: str,
    ) -> HitsByProbe: ...

    def close(self) -> None: ...


_REF_FASTA_EXTS = (".fa", ".fasta", ".fna", ".ffn", ".frn")


def _blast_reference_content_hash(configured_path: str) -> str:
    """Content signature of one configured BLAST reference, for the resume-
    cache fingerprint (``screen_cache.fingerprint``'s ``ref_content_hash``).

    Prefers hashing the SOURCE FASTA — the exact bytes ``blast.build_db``
    itself hashes to decide whether to rebuild — so a reference edited (or
    replaced) at an unchanged path is caught the moment this run starts,
    rather than depending on whether the on-disk ``.nin``/``.nhr``/``.nsq``
    happen to have been rebuilt already. Falls back to hashing the resolved
    DB's own index files (same recipe as ``pipeline._blast_db_sha256s``, kept
    separate here rather than shared: that helper runs once at manifest time
    over already-resolved prefixes, this runs before the screen even starts
    and needs its own resolution) when no FASTA sits beside ``configured_path``
    — a bare pre-built DB prefix with no source file of its own.

    Never raises: an unreadable or absent reference must not be what fails a
    fingerprint computation. Returns ``""`` (== "unknown") on any I/O problem,
    same convention as ``_blast_db_sha256s``'s "an absent key means unknown".
    """
    fasta = None
    if os.path.exists(configured_path) and not os.path.isdir(configured_path):
        fasta = configured_path
    else:
        for ext in _REF_FASTA_EXTS:
            if os.path.exists(configured_path + ext):
                fasta = configured_path + ext
                break
    if fasta is not None:
        try:
            return blast._sha256_file(fasta)
        except OSError:
            logger.debug("Could not hash reference FASTA %s", fasta,
                        exc_info=True)
            return ""
    # No source FASTA on disk -- hash whatever DB index files are there.
    prefix = blast.resolve_db_prefix(configured_path)
    if prefix is None:
        return ""
    import glob
    files = sorted(f for f in glob.glob(prefix + ".*") if os.path.isfile(f))
    if not files:
        return ""
    try:
        import hashlib
        h = hashlib.sha256()
        for f in files:
            h.update(f"{os.path.basename(f)}:{os.path.getsize(f)}\0".encode())
            with open(f, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
        return h.hexdigest()
    except OSError:
        logger.debug("Could not hash DB files at %s", prefix, exc_info=True)
        return ""


class BlastScreen:
    """``blastn-short`` against one or more BLAST databases.

    Owns the search parameters, the threshold pushdown and the truncation
    warning. The surrounding chunk/resume/stream machinery lives in the
    pipeline and is shared with the k-mer engine.

    ``blast`` module functions are called through the module (``blast.run_blastn``)
    rather than imported by name, so a test can patch ``flexprobe.blast.run_blastn``
    and have it take effect here.
    """

    name = ENGINE_BLAST

    def __init__(
        self, *, min_pident: float = 0.0, min_aln: int = 0,
        query_half_len: int = 25, num_threads: int = 1,
        max_target_seqs: int = BLAST_MAX_TARGET_SEQS,
        reference_strand: str = REFERENCE_STRAND_BOTH,
        derive_qcov: Optional[Any] = None,
        resolve_db: Optional[Any] = None,
        **_ignored: Any,
    ) -> None:
        _reject_unknown_kwargs("BlastScreen", _ignored)
        # Injected rather than imported: DB resolution builds databases with
        # makeblastdb and its result is also what the manifest describes, so it
        # stays owned by the pipeline. An engine that needs no BLAST database
        # simply never calls it.
        self._resolve_db = resolve_db
        self.query_half_len = query_half_len
        self.num_threads = num_threads
        self.max_target_seqs = max_target_seqs
        # See constants.REFERENCE_STRAND_DEFAULT: "both" (blastn's own default)
        # unless the caller has asserted the reference is sense-oriented cDNA.
        self.reference_strand = reference_strand
        # Push the post-filter thresholds into BLAST so it stops emitting the
        # huge volume of sub-threshold seeds we discard anyway. Kept no stricter
        # than the Python classifier: perc_identity is the exact identity floor,
        # qcov_hsp_perc is derived from the length floor with a safety margin so
        # a gapped boundary HSP is never dropped server-side.
        self.perc_identity: Optional[float] = (
            min_pident if min_pident > 0 else None
        )
        self.qcov_hsp_perc: Optional[float] = (
            derive_qcov(min_aln, query_half_len) if derive_qcov else None
        )
        self.db_arg = ""
        self._ref_content_hash = ""
        self._tools: Dict[str, str] = {}

    def prepare(self, db_paths: Sequence[str]) -> None:
        self._tools = blast.discover_blast()
        blast.check_version(self._tools["blastn"])
        db_prefixes = (self._resolve_db(list(db_paths)) if self._resolve_db
                       else list(db_paths))
        # blastn's -db value is a single space-separated list of prefixes, which
        # it splits on whitespace. A prefix containing spaces (e.g. an output dir
        # like ".../TRUBV L/db") would be mis-split into two bogus DBs; blastn
        # honours double quotes, so wrap those.
        self.db_arg = " ".join(
            f'"{p}"' if " " in p else p for p in db_prefixes
        )
        # Content signature per CONFIGURED path (not the resolved prefix: a
        # resolved prefix has already lost the FASTA extension that tells
        # `_blast_reference_content_hash` where to find the source file), so
        # `db_arg` staying byte-identical across a content edit at the same
        # path is caught here instead of silently reusing a stale cache — see
        # `screen_cache.fingerprint`'s `ref_content_hash` parameter.
        self._ref_content_hash = "|".join(
            _blast_reference_content_hash(p) for p in db_paths
        )

    def fingerprint_parts(self) -> Dict[str, Any]:
        return {"db_arg": self.db_arg,
                "perc_identity": self.perc_identity,
                "qcov_hsp_perc": self.qcov_hsp_perc,
                # Changes which strand(s) of the reference blastn searches, so
                # it changes which hits come back — a resumed run under a
                # different reference_strand must not replay this one's cache.
                "reference_strand": self.reference_strand,
                "ref_content_hash": self._ref_content_hash}

    def search_chunk(
        self, chunk: Sequence[ProbeCandidate], chunk_index: int,
        output_dir: str,
    ) -> HitsByProbe:
        query_fa = os.path.join(output_dir, f"_blast_query_{chunk_index:04d}.fa")
        with open(query_fa, "w") as fh:
            for c in chunk:
                for qid, qseq in half_queries(c):
                    fh.write(f">{qid}\n{qseq}\n")
        tsv_path = blast.run_blastn(
            self._tools["blastn"], query_fa, self.db_arg,
            # evalue and word_size are NOT passed. run_blastn owns the
            # complete-coverage configuration, and its defaults (a very
            # permissive evalue, a small word_size) are what keep short
            # disqualifying near-matches from being filtered out; a stricter
            # evalue or a larger word_size can miss them. An override here
            # would silently undo run_blastn's own defaults, so do not add one
            # at this call site.
            num_threads=self.num_threads,
            perc_identity=self.perc_identity, qcov_hsp_perc=self.qcov_hsp_perc,
            # max_hsps is NOT passed: -max_hsps 1 keeps only the best-SCORING
            # alignment per subject, and best-scoring is not most-disqualifying.
            # A (query, subject) pair can hide an alignment with FEWER
            # mismatches than the one kept, including an exact 25 nt match
            # masked behind a 2-mismatch hit. run_blastn defaults it to None
            # and omits the flag, which is the correct behaviour; for the same
            # reason as the evalue/word_size note above, do not override it
            # here.
            max_target_seqs=self.max_target_seqs,
            # None (omit -strand, blastn's own "both" default) unless the
            # caller has asserted the reference is sense-oriented cDNA — see
            # constants.REFERENCE_STRAND_DEFAULT.
            strand=(None if self.reference_strand == REFERENCE_STRAND_BOTH
                    else self.reference_strand),
        )
        valid_ids = {c.probe_id for c in chunk}
        hits: HitsByProbe = {}
        per_query_counts: Dict[str, int] = {}
        for h in blast.parse_hits_tsv(tsv_path):
            qid = h.get("qseqid", "")
            per_query_counts[qid] = per_query_counts.get(qid, 0) + 1
            base = qid.rsplit("__", 1)[0]   # strip the __L / __R half suffix
            if base in valid_ids:
                hits.setdefault(base, []).append(h)
        # If any query half returned exactly the cap, its off-target hit list may
        # be truncated — a genuinely non-specific probe could slip through.
        capped = [q for q, n in per_query_counts.items()
                  if n >= self.max_target_seqs]
        if capped:
            logger.warning(
                "%d query half/halves returned the max_target_seqs cap (%d); "
                "their off-target hit lists may be truncated and a non-specific "
                "probe could pass. Re-run with a larger max_target_seqs. "
                "First few: %s",
                len(capped), self.max_target_seqs, ", ".join(sorted(capped)[:5]),
            )
        for p in (query_fa, tsv_path):
            try:
                os.remove(p)
            except OSError:
                pass
        return hits

    def close(self) -> None:
        return None


class KmerScreen:
    """Exhaustive end-to-end search over a k-mer index.

    Every probe half is compared to every candidate subject position across all
    25 positions, so the mismatch count a hit carries IS the number of
    mismatches — not an inference from a trimmed local alignment. The
    consequences for ``specificity`` are the point of the whole engine:

    * ``qstart=1`` / ``qend=half_len`` means ``_flank_mismatches`` returns 0
      without touching a database, so ``full_half_mismatch_detail`` reports
      :data:`specificity.RECOUNT_EXACT` for every hit. No ``blastdbcmd``, no
      ``-parse_seqids``, no conservative fallback that over-rejects.
    * ``hit_query_coverage`` is 1.0 by construction, so ``offtarget_min_qcov``
      stops being a caveat and becomes a guarantee.

    SEARCH BUDGET. Hits are returned out to :data:`MAX_MISMATCH_REPORTED`
    mismatches, not merely to the disqualifying threshold. The extra band feeds
    the *near-miss* statistics (``specificity_margin``,
    ``offtarget_nearmiss_count``): a probe whose closest off-target sits at
    exactly ``min_mismatch`` is the tightest possible pass, and a screen that
    stopped at the threshold could not tell it from one with no off-target at
    all. 5 is not arbitrary — it is the widest budget the two-seed plan covers
    for free (``5 // 2 == 2``); going wider needs a 3-mismatch seed
    neighbourhood and roughly ten times the candidates (see
    ``kmer_index.neighbourhood``).

    THE JUNCTION SCAN is the second thing this engine produces, and it exists
    because the hit list above CANNOT answer the question it asks. Hits stop at
    :data:`MAX_MISMATCH_REPORTED` mismatches PER HALF, so the most dangerous
    shape of all — a non-target gene matching perfectly across the ligation nick
    while both halves carry 6+ mismatches out at the oligo ends — is invisible
    to every hit in the list. Post-processing the hits could never recover it.
    So the scan seeds separately, on the junction itself; see
    :meth:`_junction_scan`.

    BOTH STRANDS are searched, as ``blastn`` does: the reverse-complement query
    is searched separately and its hits are reported with ``sstart > send``,
    which is how ``specificity._flank_mismatches`` and the rest of the codebase
    already detect a minus-strand HSP.

    MULTI-DB. One index per configured reference rather than a single
    concatenated one. Indexes are cached on disk keyed by the FASTA's SHA-256,
    and one reference (a human transcriptome, say) is typically shared across
    many projects — per-DB indexes make that cache reusable, which a
    concatenation keyed on the *combination* would throw away.
    """

    name = ENGINE_KMER

    #: Report hits out to this many mismatches over the full half. See the
    #: class docstring — this is the near-miss band, not the decision threshold.
    MAX_MISMATCH_REPORTED = 5

    #: k for the index. 12 keeps bucket occupancy near 20 on a transcriptome,
    #: which is what makes the 631-code neighbourhood affordable.
    K = 12

    def __init__(
        self, *, query_half_len: int = 25,
        max_mismatch: int = MAX_MISMATCH_REPORTED,
        k: int = K, cache_dir: Optional[str] = None, mmap: bool = True,
        fallback_cache_dir: Optional[str] = None,
        max_hits_per_query: int = BLAST_MAX_TARGET_SEQS,
        exclude_non_primary: bool = (
            reference_qc.EXCLUDE_ALT_SCAFFOLDS_DEFAULT),
        reference_strand: str = REFERENCE_STRAND_BOTH,
        rebuild_kmer_index: bool = False,
        **_ignored: Any,
    ) -> None:
        _reject_unknown_kwargs("KmerScreen", _ignored)
        # Rebuild an index whose builder stamp does not match, instead of
        # stopping on KmerIndexStampMismatchError. Opt-in only.
        self.rebuild_kmer_index = rebuild_kmer_index
        self.query_half_len = query_half_len
        self.max_mismatch = max_mismatch
        self.k = k
        self.cache_dir = cache_dir
        # See constants.REFERENCE_STRAND_DEFAULT: "both" searches the half
        # itself AND its reverse complement against the reference (the
        # default, safe for a genome or negative-sense-viral reference);
        # "plus" is an explicit opt-in that skips the minus-strand arm of both
        # the half search and the junction scan below, for a reference the
        # caller asserts is sense-oriented cDNA.
        self.reference_strand = reference_strand
        # Drop ALT haplotype / PATCH scaffold records from the reference before
        # indexing it. Only this engine can honour it — BLAST searches a
        # prebuilt database, not the FASTA — see flexprobe.reference_qc.
        self.exclude_non_primary = exclude_non_primary
        # Used only when the reference's own directory is read-only; see
        # KmerIndex.cache_prefix. NOT the default location — an index belongs
        # beside its reference so every run shares it.
        self.fallback_cache_dir = fallback_cache_dir
        self.mmap = mmap
        self.max_hits_per_query = max_hits_per_query
        self._refs: List[Any] = []          # [(source_path, KmerIndex, Reference)]

    # -- lifecycle ---------------------------------------------------------

    def prepare(self, db_paths: Sequence[str]) -> None:
        from .kmer_index import KmerIndex, assert_pad_covers, plan_for

        assert_pad_covers(self.query_half_len)
        # Fail here rather than after building an index: an unsupported budget
        # would otherwise surface as silently-missing hits.
        plan_for(self.query_half_len, self.max_mismatch)
        for path in db_paths:
            fasta = self._as_fasta(path)
            idx, ref = KmerIndex.build_or_load(
                fasta, self.k, cache_dir=self.cache_dir, mmap=self.mmap,
                fallback_dir=self.fallback_cache_dir,
                exclude_non_primary=self.exclude_non_primary,
                rebuild=self.rebuild_kmer_index,
            )
            self._refs.append((fasta, idx, ref))
        if not self._refs:
            raise ValueError("the k-mer screen was given no reference to search")
        self._warn_gene_scoping()

    def _warn_gene_scoping(self) -> None:
        """Say so when per-gene off-target scoping cannot work on this input.

        [SAFE] — this changes nothing about the screen. It only names, once per
        run, a degradation the run would otherwise take silently: on a reference
        whose record ids cannot express a gene (Ensembl ``cdna.all`` being the
        common case),
        ``off_target_same_subject`` groups off-target hits per TRANSCRIPT, which
        is the permissive direction. See
        ``specificity.gene_scoping_degraded`` for the argument and the wording.

        Emitted here because ``prepare`` is where the reference headers first
        exist and runs exactly once. The BLAST engine cannot do the same check:
        it is handed a prebuilt database, never the FASTA, so it has no record
        ids to inspect until hits come back.
        """
        for path, _idx, ref in self._refs:
            msg = gene_scoping_degraded(ref.names)
            if msg:
                logger.warning("%s: %s", os.path.basename(path), msg)

    @staticmethod
    def _as_fasta(path: str) -> str:
        """The FASTA behind a configured reference path.

        A BLAST database is NOT usable here — the index is built from sequence,
        and a ``.nin``/``.nhr``/``.nsq`` set carries it only in a form that needs
        ``blastdbcmd`` to extract, which is exactly the dependency this engine
        exists to remove. Say so, with the fix, rather than failing obscurely
        later.
        """
        if os.path.exists(path) and not os.path.isdir(path):
            return path
        for ext in (".fa", ".fasta", ".fna", ".ffn", ".frn"):
            if os.path.exists(path + ext):
                return path + ext
        if blast.is_db_prefix(path):
            raise ValueError(
                f"{path!r} is a BLAST database, and the k-mer screen searches "
                "sequence directly. Point --blast-db at the FASTA the database "
                "was built from (the index is cached beside it after the first "
                "run), or select the blast engine."
            )
        raise ValueError(
            f"{path!r} is not a readable FASTA for the k-mer screen"
        )

    def fingerprint_parts(self) -> Dict[str, Any]:
        # The sources AND the search budget: widening max_mismatch changes which
        # hits come back, so a cache keyed without it would replay a narrower
        # search's answers.
        #
        # A reference whose ALT/PATCH records were dropped is a DIFFERENT set of
        # subjects under the same path, so the exclusion is folded into the
        # source term — otherwise a resumed run would replay verdicts made
        # against records this run cannot see. Marked per reference and only
        # when it actually applied, so every existing cache for a reference the
        # filter no-ops on (every viral panel) keeps its digest and stays warm.
        sources = "|".join(
            p + (":primary" if (ref.composition
                                and ref.composition.excluded_non_primary)
                 else "")
            for p, _, ref in self._refs
        )
        # Already paid for: `KmerIndex.build_or_load` hashes the source FASTA
        # to key its own on-disk index cache and stashes the result on
        # `ref.composition.source_sha256` (see kmer_index.py), so this is a
        # free reuse, not a re-hash. Without it a reference edited at an
        # unchanged path would keep `sources` byte-identical and the stale
        # resume cache would be silently read back — see `screen_cache.
        # fingerprint`'s `ref_content_hash` parameter.
        ref_content_hash = "|".join(
            (ref.composition.source_sha256 or "") if ref.composition else ""
            for _, _, ref in self._refs
        )
        return {"db_arg": sources,
                "perc_identity": float(self.max_mismatch),
                "qcov_hsp_perc": float(self.k),
                # Drops the minus-strand arm of both the half search and the
                # junction scan when "plus" — a resumed run under a different
                # reference_strand must not replay this one's cache.
                "reference_strand": self.reference_strand,
                "ref_content_hash": ref_content_hash}

    def close(self) -> None:
        self._refs.clear()

    # -- search ------------------------------------------------------------

    def search_chunk(
        self, chunk: Sequence[ProbeCandidate], chunk_index: int,
        output_dir: str,
    ) -> HitsByProbe:
        from .kmer_index import (
            INVALID, encode_seq, plan_for, revcomp_codes, search,
        )

        offs, seed_mm = plan_for(self.query_half_len, self.max_mismatch)
        hits: HitsByProbe = {}
        capped: List[str] = []
        for cand in chunk:
            for qid, qstr in half_queries(cand):
                q = encode_seq(qstr)
                if len(q) != self.query_half_len:
                    # Not an ambiguity problem: this candidate's half is not the
                    # length the engine was configured for, so there is no query
                    # to run at all. A caller/config error, reported by the
                    # length checks in `models.RunConfig.validate`.
                    continue
                # A DEGENERATE BASE DOES NOT REMOVE THE HALF FROM THE SEARCH.
                # Seeds containing one propose nothing, so skipping the half
                # would let the pair keep whatever verdict its other half earned
                # — a probe with the single 'N' the default `max_n_per_half`
                # admits could pass with a perfect off-target sitting under the
                # unsearched half. Instead the code is expanded to its concrete
                # resolutions and every one of them is searched exhaustively; see
                # `expand_degenerate_half` for the cap and why the union is the
                # safe direction.
                probes = [(qstr, q)]
                if (q == INVALID).any():
                    resolved = expand_degenerate_half(qstr)
                    if resolved is None:
                        # Too many resolutions (or not IUPAC at all). FAIL
                        # CLOSED: say the half was not searched, so
                        # `specificity` rejects the pair, rather than searching
                        # an arbitrary subset and calling it a screen.
                        hits.setdefault(cand.probe_id, []).append(
                            _unscreenable_record(cand.probe_id, qid, qstr)
                        )
                        continue
                    probes = [(r, encode_seq(r)) for r in resolved]
                out: List[Hit] = []
                for rstr, rq in probes:
                    # "plus" skips the minus-strand arm: a probe half is
                    # single-stranded DNA that hybridises to sense-strand
                    # reference sequence, so a hit that only matches the
                    # reference's reverse complement is not one it can physically
                    # form — see constants.REFERENCE_STRAND_DEFAULT for why that
                    # is only true for a sense-oriented cDNA reference, and so is
                    # an opt-in rather than the default.
                    strand_arms = (
                        ((False, rq),)
                        if self.reference_strand == REFERENCE_STRAND_PLUS
                        else ((False, rq), (True, revcomp_codes(rq)))
                    )
                    for _src, idx, ref in self._refs:
                        for minus, probe in strand_arms:
                            # `bounds=` is not optional in practice: without it
                            # a window starting in a record's trailing PAD can
                            # spend its mismatch budget on sentinels, pass, and
                            # be attributed by `_hit` to the record BEFORE the
                            # one its sequence actually came from. See
                            # kmer_index.RecordBounds.
                            pos, mm = search(
                                ref.seq, idx, probe, self.max_mismatch,
                                offs, seed_mm, bounds=ref.bounds,
                            )
                            for p, n in zip(pos.tolist(), mm.tolist()):
                                out.append(
                                    self._hit(qid, rstr, rq, ref, p, n, minus)
                                )
                if len(probes) > 1:
                    # Only a degenerate half can put two hits on one site, and
                    # the dedup costs a dict insert per hit — which at
                    # whole-transcriptome hit volumes is not something to pay on
                    # the 99.99% of halves that are pure ACGT.
                    out = _best_hit_per_site(out)
                if len(out) >= self.max_hits_per_query:
                    capped.append(qid)
                hits.setdefault(cand.probe_id, []).extend(out)
            jd = self._junction_scan(cand)
            if jd is not None:
                # Appended to the SAME list, because that list is the only
                # channel between an engine and `specificity` and the only thing
                # the resume cache persists. It is marked, not disguised:
                # `specificity.split_junction_records` takes it back out before
                # anything counts a mismatch.
                hits.setdefault(cand.probe_id, []).append(jd)
        if capped:
            logger.warning(
                "%d query half/halves hit the reporting cap (%d); their hit "
                "lists may be truncated. This usually means a probe half lies "
                "in a repeat. First few: %s",
                len(capped), self.max_hits_per_query,
                ", ".join(sorted(capped)[:5]),
            )
        return hits

    def _junction_scan(self, cand: ProbeCandidate) -> Optional[Hit]:
        """Where does this probe's LIGATION NICK sit inside unbroken homology?

        Returns one record — ``{"junction_duplex": True, "sites": [[sseqid,
        left, right], ...]}`` — or None when the question could not be asked at
        all, which ``specificity`` reports as -1 rather than as zero.

        WHY THIS CANNOT BE A POST-PROCESSING PASS OVER THE HIT LIST. The hit
        list stops at :data:`MAX_MISMATCH_REPORTED` (5) mismatches PER HALF. The
        probe this gate exists to catch — a non-target gene matching perfectly
        across the nick, with every mismatch pushed out to the oligo ends — has
        6 or more mismatches in each half and therefore produces NO HIT AT ALL.
        It is not under-described in the hit list, it is absent from it. Any
        amount of post-processing would return "clean" for exactly the probes
        that are not.

        SO IT SEEDS ON THE JUNCTION INSTEAD, which inverts the problem into the
        one thing a k-mer index is fastest at: an EXACT lookup. The seed is the
        ``k``-mer of the target window centred on the nick (k=12 -> 6 bases each
        side), looked up at zero mismatches — one bucket slice, no
        neighbourhood, ~1/600th of the work one half's own search costs.

        And the seeding is COMPLETE for every threshold this gate accepts. A
        site with at least ``k//2`` clean bases on each side of the nick
        necessarily contains that exact centred k-mer, so no such site can be
        missed; a site with fewer is one the metric cannot see, which is why
        ``models`` refuses a threshold below
        :data:`constants.OFF_TARGET_JUNCTION_DUPLEX_MIN_BP`. Every site this
        returns therefore has ``left >= 6`` and ``right >= 6`` by construction.
        (Pigeonhole contract of :func:`kmer_index.search`: one seed at
        ``max_mm=0`` satisfies ``0 // 1 <= 0``.)

        The run is then extended outward from the nick through CONTIGUOUS
        reference bases. That contiguity is the biology: the two runs live on one
        molecule in one gene, so a probe that clears the count rule on both
        halves separately can still be caught here — which is the entire point,
        since RTL needs hybridisation AND a sealed nick, not five mismatches
        somewhere. Records are separated by :data:`kmer_index.PAD` sentinels that
        can never pair, so a run stops at a transcript boundary instead of
        walking into the next transcript.

        BOTH STRANDS, like the half search, UNLESS ``reference_strand`` is
        "plus" — a minus-strand junction is exactly as phantom here as it is
        in the half search: the reference's sense sequence at that site is the
        reverse complement of the probe, not sequence the probe can hybridise
        to or ligate across. On a sense-oriented cDNA reference this is the
        larger of the two sources of minus-strand phantom off-targets: its
        record carries no per-site strand field, so no post-filter of
        the hit list can remove the contamination — the fix has to be here, in
        which frames get scanned at all. On the minus strand the probe's own
        5'/3' sense is the reverse-complement frame's, so the two counts are
        swapped back before they are reported: ``left`` always means "5' of the
        nick in the probe's frame".
        """
        from .kmer_index import INVALID, encode_seq, revcomp_codes, search

        seq = cand.sequence or ""
        n = len(seq)
        mid = n // 2                    # first base 3' of the nick
        k = self.k
        a = k // 2                      # seed bases 5' of the nick
        b = k - a                       # and 3' of it
        if n < 2 or mid < a or n - mid < b:
            return None                 # window too short to carry the seed
        w_plus = encode_seq(seq)
        # Best (left, right) per subject. Per subject rather than per site
        # because only the worst site of a gene can decide anything, and because
        # `specificity` drops target/isoform subjects afterwards — so the list
        # must stay long enough to still contain a non-target one.
        sites: Dict[str, Tuple[Tuple[int, int], Tuple[int, int]]] = {}
        flips = (
            (False,) if self.reference_strand == REFERENCE_STRAND_PLUS
            else (False, True)
        )
        for flip in flips:
            w = revcomp_codes(w_plus) if flip else w_plus
            m = n - mid if flip else mid
            seed = w[m - a:m + b]
            if bool((seed == INVALID).any()):
                # An N spanning the nick: the seed proposes nothing, so a scan
                # would report "measured, found nothing" for a window it never
                # searched. Unknown is the honest answer.
                return None
            w_str = _codes_to_str(w)
            for _src, idx, ref in self._refs:
                # `bounds=` is a no-op for THIS call and is passed anyway. At
                # max_mm=0 a window covering a PAD sentinel already fails (a
                # sentinel equals no query code, so it costs a mismatch the
                # budget cannot afford), and the seed's own 12-mer is never
                # indexed unless it is wholly inside a record. So the argument
                # removes nothing at max_mm=0, but the containment invariant should
                # not silently depend on a caller's mismatch budget staying 0,
                # and `record_of(p)` below is only sound while it holds.
                pos, _mm = search(ref.seq, idx, seed, 0, (0,), 0,
                                  bounds=ref.bounds)
                if len(pos) > JUNCTION_DUPLEX_MAX_SITES:
                    logger.warning(
                        "Probe %s: the 12-mer spanning its ligation junction "
                        "occurs %d times in %s (chance alone gives ~24), so the "
                        "junction sits in a repeat. Its clean-duplex columns are "
                        "reported as unmeasured (-1) rather than as the worst of "
                        "an arbitrary subset; the probe keeps its off-target "
                        "count-rule verdict.",
                        cand.probe_id, len(pos), os.path.basename(_src),
                    )
                    return None
                for p in pos.tolist():
                    lo = p + a - m               # window start in the reference
                    sub = _window_str(ref.seq, lo, n)
                    left, right = _junction_clean_duplex(w_str, sub, m)
                    if flip:
                        # The reverse-complement frame reads the probe backwards,
                        # so its "left" is the probe's right.
                        left, right = right, left
                    name = ref.names[int(ref.record_of(np.asarray([p]))[0])]
                    # Keyed on the MINIMUM, because the gate is per side; the
                    # total only breaks ties.
                    key = (min(left, right), left + right)
                    if key > sites.get(name, ((-1, -1), (0, 0)))[0]:
                        sites[name] = (key, (left, right))
        return {
            "qseqid": f"{cand.probe_id}{JUNCTION_QUERY_SUFFIX}",
            JUNCTION_DUPLEX_KEY: True,
            # Lists, not tuples: this round-trips through the JSON resume cache.
            "sites": [[name, lr[0], lr[1]] for name, (_key, lr) in
                      sorted(sites.items())],
        }

    def _hit(
        self, qid: str, qstr: str, q: Any, ref: Any, pos: int, mm: int,
        minus: bool,
    ) -> Hit:
        """One full-half hit, in the dict shape ``specificity`` consumes."""
        from .kmer_index import revcomp_codes

        L = self.query_half_len
        rec = int(ref.record_of(np.asarray([pos]))[0])
        local = pos - int(ref.starts[rec]) + 1        # 1-based within record
        window = ref.seq[pos:pos + L]
        if minus:
            # The reverse-complement query matched here, so in the PROBE's
            # reading frame the subject is the reverse complement of this
            # window. sstart > send is how the rest of the codebase detects
            # a minus-strand hit (see specificity._flank_mismatches).
            sseq = _codes_to_str(revcomp_codes(window))
            sstart, send = local + L - 1, local
        else:
            sseq = _codes_to_str(window)
            sstart, send = local, local + L - 1
        return {
            "qseqid": qid,
            "sseqid": ref.names[rec],
            "pident": 100.0 * (L - mm) / L,
            "length": L,
            "mismatch": mm,
            "gapopen": 0,
            # The HSP spans the WHOLE half. This is the field that makes
            # specificity._flank_mismatches return 0 and the recount exact.
            "qstart": 1, "qend": L,
            "sstart": sstart, "send": send,
            # Undefinable for a k-mer search, and NaN rather than 0 because 0 is
            # a meaningful bitscore and would read as a real, terrible hit.
            # results.csv records the engine so these are never misread.
            "evalue": float("nan"), "bitscore": float("nan"),
            "qseq": qstr, "sseq": sseq,
        }


_BASES = "ACGTN"
# Vectorised form of the same table, for the junction scan's per-site window
# conversion. That runs once per candidate reference site — tens of times per
# probe — so the per-character Python loop below is the wrong shape for it.
_BASE_BYTES = np.frombuffer(_BASES.encode("ascii"), dtype=np.uint8)


def _codes_to_str(codes: Any) -> str:
    return "".join(_BASES[int(c)] for c in codes)


def _window_str(seq: Any, lo: int, length: int) -> str:
    """``length`` reference bases from ``lo``, as a string, padded off the ends.

    Positions before the start or past the end of the reference are filled with
    ``N``. That is not a guess: there is no sequence there, so nothing can
    base-pair, and ``junction_clean_duplex`` ends a run at exactly such a
    position. Same reasoning as ``specificity._flank_mismatches`` giving a flank
    that runs off a transcript end — absent is known, not unknown.
    """
    hi = lo + length
    clipped = seq[max(0, lo):max(0, min(len(seq), hi))]
    body = _BASE_BYTES[clipped].tobytes().decode("ascii")
    pad_lo = "N" * max(0, -lo)
    return (pad_lo + body).ljust(length, "N")


def _reject_unknown_kwargs(cls_name: str, ignored: dict) -> None:
    """Accept the OTHER engine's parameters; refuse anything neither accepts.

    WHY THIS EXISTS AND WHY IT IS NOT JUST ``**_ignored``. ``make_screen``
    passes ONE kwarg set to whichever engine is selected, so each engine has to
    tolerate the other's parameters -- ``cache_dir`` means nothing to BLAST and
    ``min_pident`` means nothing to the k-mer engine. A bare ``**_ignored``
    buys that tolerance at the price of swallowing EVERYTHING, including
    typos and settings the caller believed were taking effect.

    For example, ``KmerScreen(off_target_min_mismatch=7)`` and a literal
    nonsense keyword would both be absorbed silently. The first is the
    dangerous one, because ``off_target_min_mismatch`` is a real setting
    elsewhere in the codebase, so a caller can reasonably believe the search
    budget moved when it did not.

    The accepted set is derived from the two signatures at call time rather
    than hardcoded, so it cannot drift as either engine gains a parameter.
    """
    if not ignored:
        return
    accepted = set()
    for cls in (BlastScreen, KmerScreen):
        accepted |= {
            name for name, param in inspect.signature(cls.__init__).parameters.items()
            if param.kind is not inspect.Parameter.VAR_KEYWORD and name != "self"
        }
    unknown = sorted(set(ignored) - accepted)
    if unknown:
        raise TypeError(
            f"{cls_name}: unknown keyword argument(s) {unknown}. "
            f"Neither screening engine accepts them, so nothing would have "
            f"used the value. Accepted across both engines: "
            f"{sorted(accepted)}"
        )

def make_screen(engine: str, **kwargs: Any) -> Screen:
    """Construct the screen for ``engine``.

    The single place engine names are resolved, so an unknown or not-yet-built
    engine fails HERE — before a DB is resolved, an index is built, or a chunk
    is searched — rather than partway through a long run.
    """
    if engine == ENGINE_BLAST:
        return BlastScreen(**kwargs)
    if engine == ENGINE_KMER:
        return KmerScreen(**kwargs)
    raise ValueError(
        f"unknown screen engine {engine!r}; expected one of {ENGINES}"
    )
