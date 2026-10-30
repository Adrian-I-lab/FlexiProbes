"""Streaming filter equivalence (P1)."""
from flexprobe.filters import apply_filters, evaluate_candidate, new_filter_counters
from flexprobe.models import FilterThresholds, ProbeCandidate
from flexprobe.pipeline import _stream_filter
from flexprobe.tiling import tile_sequence

from conftest import SYNTH_200


def _tiles():
    for start, end, sub in tile_sequence(SYNTH_200, length=50, step=1):
        yield ProbeCandidate(gene_id="G", probe_id=f"G_t{start}",
                             sequence=sub, start=start, end=end)


THR = FilterThresholds(gc_min=0.40, gc_max=0.75, max_n_per_half=1,
                       max_homopolymer_run=5, require_tn_junction=True)


def test_stream_filter_matches_apply_filters():
    # Materialized reference.
    ref = list(_tiles())
    ref_counts = apply_filters(ref, THR)
    ref_pass = [c.probe_id for c in ref if c.passed_filters]

    # Streamed.
    passing, counters, rejects = _stream_filter(_tiles(), THR, reject_cap=25)

    assert counters["total"] == ref_counts["total"]
    assert counters["passed"] == ref_counts["passed"]
    assert counters["gc_range"] == ref_counts["gc_range"]
    assert counters["tn_junction"] == ref_counts["tn_junction"]
    assert [c.probe_id for c in passing] == ref_pass


def test_stream_reject_sample_is_capped():
    passing, counters, rejects = _stream_filter(_tiles(), THR, reject_cap=5)
    # No reason should appear more than the cap times in the sample.
    from collections import Counter
    reason_counts = Counter()
    for c in rejects:
        for r in c.filter_fail_reasons:
            reason_counts[r] += 1
    assert all(v <= 5 for v in reason_counts.values())
    # But total rejected (in counters) can far exceed the sample size.
    assert counters["total"] - counters["passed"] >= len(rejects)


def test_evaluate_candidate_single():
    c = ProbeCandidate(probe_id="x", sequence="G" * 50, start=0, end=50)
    counters = new_filter_counters()
    passed = evaluate_candidate(c, THR, counters)
    assert passed is False
    assert counters["total"] == 1
    assert "homopolymer" in c.filter_fail_reasons
