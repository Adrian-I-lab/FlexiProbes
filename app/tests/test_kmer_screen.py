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

"""The k-mer screen's hit contract, and its agreement with BLAST.

The point of this engine is that a hit's mismatch count IS the number of
mismatches over the whole 25 nt half — not an inference from a trimmed local
alignment. These tests pin that end of it: the hits it emits must drive
``specificity`` to ``RECOUNT_EXACT`` with no subject retrieval, and must reach
the same verdicts BLAST does on the cases BLAST gets right.
"""
from __future__ import annotations

import random
import shutil

import pytest

from flexprobe.models import ProbeCandidate
from flexprobe.screen import ENGINE_KMER, KmerScreen, make_screen
from flexprobe.specificity import (
    RECOUNT_EXACT, evaluate_specificity, full_half_mismatch_detail, hit_half,
    split_junction_records,
)

HALF = 25


def _rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _cand(seq, pid="p0", gene="g"):
    return ProbeCandidate(gene_id=gene, probe_id=pid, sequence=seq,
                          start=0, end=len(seq))


def _aln(hits):
    """The ALIGNMENTS in a probe's hit list, without the junction-duplex record.

    ``search_chunk`` also appends one non-alignment record per probe (the
    junction scan — see ``KmerScreen._junction_scan``), which carries no
    qseq/sseq, no half suffix and no mismatch count. Every test below is about
    the alignment contract, so it filters exactly the way production code does.
    """
    return split_junction_records(hits)[0]


def _screen(tmp_path, fasta_text, **kw):
    fa = tmp_path / "ref.fa"
    fa.write_text(fasta_text)
    s = KmerScreen(cache_dir=str(tmp_path / "idx"), mmap=False, **kw)
    s.prepare([str(fa)])
    return s


# ---------------------------------------------------------------------------
# The hit contract — what makes the recount exact
# ---------------------------------------------------------------------------

def test_hits_span_the_full_half_so_the_recount_is_exact(tmp_path):
    """The whole reason the engine exists.

    ``qstart=1``/``qend=25`` makes ``_flank_mismatches`` return 0 without
    touching a database, so every hit is RECOUNT_EXACT with ``subjects=None``.
    Under BLAST this needed blastdbcmd, and degraded to a lower-bound count
    (which over-rejects) whenever retrieval failed.
    """
    rng = random.Random(11)
    ref = _rand_seq(3000, rng)
    probe = ref[500:550]
    s = _screen(tmp_path, f">REC1\n{ref}\n")
    hits = _aln(s.search_chunk([_cand(probe)], 0, str(tmp_path))["p0"])
    assert hits, "the probe's own locus must be found"
    for h in hits:
        assert h["qstart"] == 1 and h["qend"] == HALF
        assert h["length"] == HALF
        assert len(h["qseq"]) == HALF and len(h["sseq"]) == HALF
        n, basis = full_half_mismatch_detail(
            h, HALF, query_seq=h["qseq"], subjects=None,
        )
        assert basis == RECOUNT_EXACT
        assert n == h["mismatch"]


def test_mismatch_count_equals_column_disagreements(tmp_path):
    """The reported count must equal what qseq/sseq actually disagree on."""
    rng = random.Random(3)
    ref = _rand_seq(2000, rng)
    probe = list(ref[100:150])
    for p in (2, 20, 30, 47):
        probe[p] = "ACGT"[("ACGT".index(probe[p]) + 1) % 4]
    s = _screen(tmp_path, f">REC1\n{ref}\n")
    hits = _aln(s.search_chunk([_cand("".join(probe))], 0, str(tmp_path))["p0"])
    assert hits
    for h in hits:
        observed = sum(1 for a, b in zip(h["qseq"], h["sseq"]) if a != b)
        assert h["mismatch"] == observed
        assert h["pident"] == pytest.approx(100.0 * (HALF - observed) / HALF)


def test_bitscore_and_evalue_are_nan_not_zero(tmp_path):
    """Undefinable for a k-mer search. 0 is a MEANINGFUL bitscore and would
    read as a real, terrible hit; NaN cannot be mistaken for a measurement."""
    rng = random.Random(5)
    ref = _rand_seq(1500, rng)
    s = _screen(tmp_path, f">REC1\n{ref}\n")
    hits = _aln(s.search_chunk([_cand(ref[200:250])], 0, str(tmp_path))["p0"])
    assert hits
    for h in hits:
        assert h["bitscore"] != h["bitscore"]      # NaN
        assert h["evalue"] != h["evalue"]


def test_both_halves_are_queried_and_labelled(tmp_path):
    rng = random.Random(7)
    ref = _rand_seq(2000, rng)
    s = _screen(tmp_path, f">REC1\n{ref}\n")
    hits = _aln(s.search_chunk([_cand(ref[300:350])], 0, str(tmp_path))["p0"])
    assert {hit_half(h) for h in hits} == {"L", "R"}


# ---------------------------------------------------------------------------
# Strandedness
# ---------------------------------------------------------------------------

def test_minus_strand_hits_are_found_and_flagged(tmp_path):
    """At the DEFAULT setting (``reference_strand="both"``), blastn's own
    default searches both
    strands; so must this, since a genome or negative-sense-viral reference
    can carry a genuine off-target there. ``sstart > send`` is how the rest of
    the codebase detects it (specificity._flank_mismatches).

    Whether such a hit is real or a phantom depends on what kind of reference
    it came from (see constants.REFERENCE_STRAND_DEFAULT); this test only
    pins that the search finds and correctly flags it when asked to look —
    ``test_reference_strand.py`` covers the "should this count at all"
    question and its ``reference_strand="plus"`` opt-out.
    """
    rng = random.Random(13)
    ref = _rand_seq(2000, rng)
    probe = _revcomp(ref[800:850])          # matches only the minus strand
    s = _screen(tmp_path, f">REC1\n{ref}\n")
    hits = _aln(s.search_chunk([_cand(probe)], 0, str(tmp_path))["p0"])
    assert hits
    minus = [h for h in hits if h["sstart"] > h["send"]]
    assert minus, "a reverse-complement match must be reported as minus strand"
    for h in minus:
        assert h["mismatch"] == 0
        assert h["qseq"] == h["sseq"]       # subject put in the probe's frame


def test_plus_strand_hits_report_ascending_coordinates(tmp_path):
    rng = random.Random(17)
    ref = _rand_seq(2000, rng)
    s = _screen(tmp_path, f">REC1\n{ref}\n")
    hits = _aln(s.search_chunk([_cand(ref[400:450])], 0, str(tmp_path))["p0"])
    plus = [h for h in hits if h["mismatch"] == 0 and h["sstart"] < h["send"]]
    assert plus
    for h in plus:
        assert h["send"] - h["sstart"] == HALF - 1


def test_subject_coordinates_are_record_local(tmp_path):
    """Coordinates must be 1-based within their own record, not offsets into
    the concatenated reference — otherwise every record after the first reports
    positions that do not exist in it."""
    rng = random.Random(19)
    a, b = _rand_seq(300, rng), _rand_seq(300, rng)
    s = _screen(tmp_path, f">A\n{a}\n>B\n{b}\n")
    hits = _aln(s.search_chunk([_cand(b[100:150])], 0, str(tmp_path))["p0"])
    onb = [h for h in hits if h["sseqid"] == "B"]
    assert onb
    for h in onb:
        assert 1 <= min(h["sstart"], h["send"]) <= len(b)
        assert 1 <= max(h["sstart"], h["send"]) <= len(b)


def test_no_hit_is_assembled_out_of_a_records_trailing_pad(tmp_path):
    """REGRESSION, at the layer that names the gene.

    The engine-level case is in ``test_kmer_oracle.py``; this one asserts the
    SCREEN cannot emit it, because that is where a position becomes an
    ``sseqid`` and therefore where the error does its damage. ``KmerScreen``
    turns positions into hits in two places — the per-half search and
    ``_junction_scan`` — and either one forgetting to bound its search would
    show up here.

    ``lhs`` is four filler bases then recB's opening 21, so its best window
    starts four bases inside recA's PAD: four sentinels bought with the whole
    4-mismatch budget, 21 bases of recB, and ``record_of`` answering **recA**.
    Before the fix this test reports a hit on recA at ``sstart`` past the end of
    recA — a coordinate that does not exist in the sequence it names.
    """
    rec_a = "A" * 300
    rec_b = ("GTTACCGATCAGTTACCGATCAGTTACCGATCAGGGCCTTAACCGGTTAA"
             + "C" * 250)
    lengths = {"recA": len(rec_a), "recB": len(rec_b)}
    probe = "TTTT" + rec_b[:46]
    assert len(probe) == 50
    s = _screen(tmp_path, f">recA\n{rec_a}\n>recB\n{rec_b}\n")
    hits = _aln(s.search_chunk([_cand(probe)], 0, str(tmp_path))["p0"])
    assert hits, "the genuine recB hit must still be found"
    assert any(h["sseqid"] == "recB" for h in hits)
    for h in hits:
        n = lengths[h["sseqid"]]
        assert 1 <= min(h["sstart"], h["send"]) <= n, h
        assert 1 <= max(h["sstart"], h["send"]) <= n, h
        # The sequence a hit reports must be sequence, not separator.
        assert "N" not in h["sseq"], h


# ---------------------------------------------------------------------------
# Multi-DB, N handling, and input validation
# ---------------------------------------------------------------------------

def test_searches_every_configured_reference(tmp_path):
    rng = random.Random(23)
    a, b = _rand_seq(500, rng), _rand_seq(500, rng)
    (tmp_path / "a.fa").write_text(f">A\n{a}\n")
    (tmp_path / "b.fa").write_text(f">B\n{b}\n")
    s = KmerScreen(cache_dir=str(tmp_path / "idx"), mmap=False)
    s.prepare([str(tmp_path / "a.fa"), str(tmp_path / "b.fa")])
    # A probe whose halves come from different references: both must be seen.
    probe = a[100:125] + b[200:225]
    hits = _aln(s.search_chunk([_cand(probe)], 0, str(tmp_path))["p0"])
    assert {h["sseqid"] for h in hits} >= {"A", "B"}


def test_a_blast_database_is_refused_with_the_fix_named(tmp_path, monkeypatch):
    """A BLAST DB carries sequence only in a form needing blastdbcmd — the very
    dependency this engine removes. Fail with the fix, not obscurely later."""
    from flexprobe import blast as fp_blast
    monkeypatch.setattr(fp_blast, "is_db_prefix", lambda _p: True)
    s = KmerScreen(cache_dir=str(tmp_path / "idx"))
    with pytest.raises(ValueError, match="BLAST database"):
        s.prepare([str(tmp_path / "nonexistent_db")])


def test_no_reference_is_an_error(tmp_path):
    s = KmerScreen(cache_dir=str(tmp_path / "idx"))
    with pytest.raises(ValueError, match="no reference"):
        s.prepare([])


def test_probe_half_containing_n_is_expanded_not_skipped(tmp_path):
    """A half whose seeds cannot be built is RESOLVED, never dropped.

    Skipping such a half (so it contributes no hits) is worse than the
    problem it avoids: the pair then keeps whatever verdict its OTHER half
    earns, so a probe carrying the single 'N' that ``max_n_per_half`` admits by
    default can be certified specific with a perfect off-target sitting under
    the half that was never searched. So seeds that propose nothing are
    replaced by seeds that propose everything: the code is expanded to its
    concrete resolutions and each is searched exhaustively. See
    ``screen.expand_degenerate_half``.
    """
    rng = random.Random(29)
    ref = _rand_seq(2000, rng)
    probe = "N" + ref[600:624] + ref[625:650]
    s = _screen(tmp_path, f">REC1\n{ref}\n")
    hits = _aln(s.search_chunk([_cand(probe)], 0, str(tmp_path)).get("p0", []))
    left = [h for h in hits if hit_half(h) == "L"]
    assert left, "the N-carrying half must still be searched"
    # ref[599] is one of the four resolutions of the 'N', so the best duplex
    # this half can form at 599 is perfect — and the union over resolutions
    # keeps that best case, which is the direction that over-rejects.
    assert min(h["mismatch"] for h in left) == 0
    # One site, one hit: the four resolutions are collapsed back onto it.
    starts = [h["sstart"] for h in left if h["mismatch"] == 0]
    assert len(starts) == len(set(starts))


def test_near_miss_band_is_reported_above_the_threshold(tmp_path):
    """Hits are returned to MAX_MISMATCH_REPORTED, not to the disqualifying
    threshold: specificity_margin needs the band just past it to tell a tight
    pass from a probe with no off-target at all."""
    rng = random.Random(31)
    ref = _rand_seq(3000, rng)
    half = list(ref[1000:1025])
    for p in (1, 6, 12, 18, 23):                     # exactly 5 mismatches
        half[p] = "ACGT"[("ACGT".index(half[p]) + 1) % 4]
    probe = "".join(half) + ref[1025:1050]
    s = _screen(tmp_path, f">REC1\n{ref}\n")
    hits = _aln(s.search_chunk([_cand(probe)], 0, str(tmp_path))["p0"])
    assert 5 in {h["mismatch"] for h in hits if hit_half(h) == "L"}


# ---------------------------------------------------------------------------
# Agreement with BLAST
# ---------------------------------------------------------------------------

def _brute_min_mismatch(ref: str, query: str) -> int:
    """Fewest mismatches between ``query`` and any locus, either strand.

    Ground truth by definition. Quadratic, hence the small reference and the
    ``slow`` marker.
    """
    rc = query.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    best = len(query)
    for probe in (query, rc):
        for i in range(len(ref) - len(probe) + 1):
            best = min(best, sum(1 for a, b in zip(ref[i:i + len(probe)], probe)
                                 if a != b))
    return best


@pytest.mark.slow
@pytest.mark.skipif(
    shutil.which("blastn") is None or shutil.which("makeblastdb") is None,
    reason="ncbi-blast+ not installed",
)
def test_engine_disagreements_are_adjudicated_by_brute_force(tmp_path):
    """Where the engines differ, BRUTE FORCE decides — not the other engine.

    Cross-engine *agreement* is the wrong contract, and asserting it hid a real
    result the first time this ran. BLAST is heuristic: it reports the
    maximum-scoring LOCAL segment, so it can report a short exact core at one
    locus and never surface the genuinely closest one. Observed here on a probe
    whose right half sits 4 mismatches from the reference: BLAST reported an
    8 nt / 100% core elsewhere, recounted it at 13 mismatches, and passed the
    probe. The k-mer screen found the true 4-mismatch locus and rejected it.

    So the real contract is directional:

    * k-mer rejects, BLAST passes  -> EXPECTED. The exhaustive screen caught a
      non-specific probe BLAST missed. Verified against brute force below.
    * BLAST rejects, k-mer passes  -> the exhaustive screen cannot legitimately
      MISS a hit, so this is either BLAST's conservative recount over-rejecting
      or a k-mer bug. Adjudicated, never assumed.
    """
    from flexprobe.blast import SubjectWindows, build_db, discover_blast

    rng = random.Random(20260811)
    ref = _rand_seq(6000, rng)
    fa = tmp_path / "ref.fa"
    fa.write_text(f">OFFTGT\n{ref}\n")
    db = str(tmp_path / "db")
    build_db(discover_blast()["makeblastdb"], str(fa), db, title="off")

    cands = []
    for i in range(12):
        start = rng.randrange(0, len(ref) - 50)
        seq = list(ref[start:start + 50])
        for p in rng.sample(range(50), rng.randint(0, 10)):
            seq[p] = "ACGT"[("ACGT".index(seq[p]) + 1) % 4]
        cands.append(_cand("".join(seq), pid=f"p{i}"))

    km = make_screen(ENGINE_KMER, cache_dir=str(tmp_path / "idx"), mmap=False)
    km.prepare([str(fa)])
    k_hits = km.search_chunk(cands, 0, str(tmp_path))

    bl = make_screen("blast", min_pident=75.0, min_aln=9, query_half_len=HALF,
                     derive_qcov=lambda a, b: 0.0)
    bl.prepare([db])
    b_hits = bl.search_chunk(cands, 0, str(tmp_path))
    subjects = SubjectWindows([db])

    kmer_caught_extra = 0
    for c in cands:
        ck = _cand(c.sequence, pid=c.probe_id)
        cb = _cand(c.sequence, pid=c.probe_id)
        # The junction-duplex gate is switched OFF for this comparison on
        # purpose. It is OUR rule, not 10x's, and only the k-mer engine can
        # measure it — BLAST emits no junction scan — so leaving it on would
        # compare the two engines on a criterion only one of them has, and
        # every disagreement it produced would be a property of the engines'
        # capabilities rather than of their search. What this test adjudicates
        # is the COUNT rule, against brute-force truth.
        _no_jd = dict(off_target_junction_duplex_bp=0)
        evaluate_specificity(ck, k_hits.get(c.probe_id, []), set(), **_no_jd)
        evaluate_specificity(cb, b_hits.get(c.probe_id, []), set(),
                             subjects=subjects, **_no_jd)
        if ck.specificity_pass == cb.specificity_pass:
            continue
        truth = {w: _brute_min_mismatch(ref, s) for w, s in
                 (("L", c.sequence[:HALF]), ("R", c.sequence[HALF:]))}
        both_close = all(v <= 4 for v in truth.values())
        if not ck.specificity_pass and cb.specificity_pass:
            # k-mer rejected. Truth must agree that both halves are close.
            assert both_close, (
                f"{c.probe_id}: k-mer rejected but truth says {truth}"
            )
            kmer_caught_extra += 1
        else:
            pytest.fail(
                f"{c.probe_id}: BLAST rejected and the EXHAUSTIVE screen "
                f"passed; truth={truth}. The k-mer search cannot miss a hit, "
                f"so this is a bug in it."
            )
    # Not asserted as a fixed number — it depends on the seed — but the case
    # above did occur, and this records that the direction is the expected one.
    assert kmer_caught_extra >= 0


@pytest.mark.slow
@pytest.mark.skipif(
    shutil.which("blastn") is None or shutil.which("makeblastdb") is None,
    reason="ncbi-blast+ not installed",
)
def test_kmer_finds_every_close_off_target_blast_does(tmp_path):
    """The k-mer screen is exhaustive, so its <=4-mismatch set must CONTAIN
    BLAST's. A probe BLAST rejects and the k-mer screen passes would mean the
    exhaustive search missed something — the one failure mode it cannot have.
    """
    from flexprobe.blast import SubjectWindows, build_db, discover_blast

    rng = random.Random(4578)
    ref = _rand_seq(8000, rng)
    fa = tmp_path / "ref.fa"
    fa.write_text(f">OFFTGT\n{ref}\n")
    db = str(tmp_path / "db")
    build_db(discover_blast()["makeblastdb"], str(fa), db, title="off")

    cands = []
    for i in range(15):
        start = rng.randrange(0, len(ref) - 50)
        seq = list(ref[start:start + 50])
        for p in rng.sample(range(50), rng.randint(0, 8)):
            seq[p] = "ACGT"[("ACGT".index(seq[p]) + 1) % 4]
        cands.append(_cand("".join(seq), pid=f"p{i}"))

    km = make_screen(ENGINE_KMER, cache_dir=str(tmp_path / "idx"), mmap=False)
    km.prepare([str(fa)])
    k_hits = km.search_chunk(cands, 0, str(tmp_path))
    bl = make_screen("blast", min_pident=75.0, min_aln=9, query_half_len=HALF,
                     derive_qcov=lambda a, b: 0.0)
    bl.prepare([db])
    b_hits = bl.search_chunk(cands, 0, str(tmp_path))
    subjects = SubjectWindows([db])

    for c in cands:
        halves = {"L": c.sequence[:HALF], "R": c.sequence[HALF:]}
        # (half, subject, mismatches) that BLAST's exact recount puts at <= 4
        blast_close = set()
        for h in b_hits.get(c.probe_id, []):
            w = hit_half(h)
            n, basis = full_half_mismatch_detail(
                h, HALF, query_seq=halves.get(w, ""), subjects=subjects,
            )
            if basis == RECOUNT_EXACT and n <= 4:
                blast_close.add((w, h["sseqid"], n))
        kmer_close = {
            (hit_half(h), h["sseqid"], h["mismatch"])
            for h in _aln(k_hits.get(c.probe_id, [])) if h["mismatch"] <= 4
        }
        missing = blast_close - kmer_close
        assert not missing, f"{c.probe_id}: exhaustive search missed {missing}"
