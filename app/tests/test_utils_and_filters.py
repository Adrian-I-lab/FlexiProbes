"""Shared utils (safe_tm) + the filter funnel."""
import importlib.util

import pytest

from flexprobe._utils import fmt_tm, is_nan, safe_tm
from flexprobe.filters import apply_filters
from flexprobe.io_fasta import parse_fasta
from flexprobe.models import FilterThresholds, ProbeCandidate

# safe_tm delegates to primer3.calc_tm and returns NaN when primer3 cannot be
# imported. primer3-py is a required dependency, so that only happens in a
# broken environment. Skip the "real Tm value" assertions there rather than
# fail the whole module; the NaN-fallback behaviour is still covered by
# test_safe_tm_ambiguous_is_nan.
_HAS_PRIMER3 = importlib.util.find_spec("primer3") is not None
requires_primer3 = pytest.mark.skipif(
    not _HAS_PRIMER3, reason="primer3 not importable (required dependency)"
)


@requires_primer3
def test_safe_tm_clean_sequence():
    # 61.538 °C with primer3-py 2.3.1, the version constraints.txt pins. The
    # 0.5 °C tolerance absorbs rounding differences between primer3-py 2.x
    # releases (the unpinned CI leg) but not a change to the Tm model or its
    # salt and oligo concentrations, which moves this sequence by degrees.
    v = safe_tm("CACTGCCACTGTCGGTTAGTTTAAT")
    assert not is_nan(v)
    assert v == pytest.approx(61.538, abs=0.5)


@requires_primer3
def test_safe_tm_lowercase_ok():
    assert abs(safe_tm("cactgccactgtcggttagtttaat")
               - safe_tm("CACTGCCACTGTCGGTTAGTTTAAT")) < 1e-6


def test_safe_tm_ambiguous_is_nan():
    assert is_nan(safe_tm("GTGATGGTGGGCANCCATGCTTTAT"))
    assert is_nan(safe_tm(""))


def test_fmt_tm():
    assert fmt_tm(float("nan")) == "n/a"
    assert fmt_tm(None) == "n/a"
    assert fmt_tm(62.345) == "62.34" or fmt_tm(62.345) == "62.35"


def _candidate(seq):
    return ProbeCandidate(probe_id="t", sequence=seq, start=0, end=len(seq))


def test_filter_gc_once_matches_expected():
    # 50-bp, perfectly balanced halves at 0.48 GC.
    seq = ("ATGCATGCATGCATGCATGCATGCA" "TGCATGCATGCATGCATGCATGCAT")
    c = _candidate(seq)
    thr = FilterThresholds(gc_min=0.30, gc_max=0.70, max_n_per_half=0,
                           max_homopolymer_run=6, require_tn_junction=False)
    apply_filters([c], thr)
    assert abs(c.gc_left - 0.48) < 0.05
    assert abs(c.gc_content - (c.gc_left + c.gc_right) / 2) < 1e-6


def test_filter_tn_junction():
    seq = "A" * 50
    c = _candidate(seq[:25] + "C" + seq[26:])   # pos 25 == C -> fails TN
    thr = FilterThresholds(gc_min=0.0, gc_max=1.0, max_n_per_half=50,
                           max_homopolymer_run=99, require_tn_junction=True)
    counts = apply_filters([c], thr)
    assert "tn_junction" in c.filter_fail_reasons
    assert counts["tn_junction"] == 1


def test_filter_homopolymer_and_ncount():
    seq = "G" * 50
    c = _candidate(seq)
    thr = FilterThresholds(gc_min=0.0, gc_max=1.0, max_n_per_half=0,
                           max_homopolymer_run=4, require_tn_junction=False)
    apply_filters([c], thr)
    assert "homopolymer" in c.filter_fail_reasons
    assert c.homopolymer_max == 50


def test_filter_rejects_gap_tile():
    # A tile spanning an alignment gap ('-') from an aligned/consensus FASTA
    # must be dropped up front (and not crash GC/chemistry scoring), so gaps
    # never leak into a designed probe. Lock-in for the gap-handling fix.
    seq = "ACGT" * 6 + "-" + "ACGTACGTACGTACGTACGTACGTA"   # 50 bp incl. one gap
    assert len(seq) == 50
    c = _candidate(seq)
    thr = FilterThresholds(gc_min=0.0, gc_max=1.0, max_n_per_half=50,
                           max_homopolymer_run=99, require_tn_junction=False)
    counts = apply_filters([c], thr)
    assert "gap_or_invalid" in c.filter_fail_reasons
    assert counts["gap_or_invalid"] == 1
    assert not c.passed_filters


def test_iter_fasta_streams_records(tmp_path):
    # Audit #42: iter_fasta yields records lazily (one at a time) and parse_fasta
    # is just the list-materialising wrapper over it.
    from flexprobe.io_fasta import iter_fasta
    import types
    p = tmp_path / "multi.fasta"
    p.write_text(">g1 desc\nACGTACGT\n>g2\nTTTTGGGG\n>g3\nCCCCAAAA\n")
    it = iter_fasta(str(p))
    assert isinstance(it, types.GeneratorType)
    first = next(it)
    assert first[0] == "g1" and first[2] == "ACGTACGT"
    rest = list(it)
    assert [r[0] for r in rest] == ["g2", "g3"]
    # parse_fasta returns the same records as a materialised list.
    assert [r[0] for r in parse_fasta(str(p))] == ["g1", "g2", "g3"]


def test_parse_fasta_tolerates_gaps(tmp_path):
    # Aligned / consensus FASTAs carry '-' gaps; ingest must not abort on them
    # (gaps are kept so tile coordinates stay faithful and the filter funnel
    # drops gap-containing windows). Lock-in for the gap-handling fix.
    p = tmp_path / "gapped.fasta"
    p.write_text(">seg|Large|ACC1\nACGTACGT--ACGTACGT\nACGT----ACGT\n")
    recs = parse_fasta(str(p))
    assert len(recs) == 1
    gid, _desc, seq = recs[0]
    # The id is the full header token — NOT truncated at '|', which would
    # merge distinct pipe-delimited records (see test_fasta_identity.py).
    assert gid == "seg|Large|ACC1"
    assert "-" in seq            # gaps preserved (coordinates intact)
