"""B7: ranking weights + specificity thresholds are configurable per run."""
import pytest

from conftest import FakeSubjects, aligned_hit
from flexprobe import constants
from flexprobe.errors import ConfigError
from flexprobe.models import ProbeCandidate, RunConfig, SpecificityThresholds
from flexprobe.specificity import batch_evaluate


def test_specificity_thresholds_default_to_constants():
    s = SpecificityThresholds()
    assert s.off_target_pident_min == constants.OFF_TARGET_PIDENT_MIN
    assert s.on_target_pident_min == constants.ON_TARGET_PIDENT_MIN
    RunConfig().validate()  # the new fields don't break a default config


def test_ranking_weights_unknown_key_rejected():
    cfg = RunConfig(ranking_weights={"bogus": 1.0})
    with pytest.raises(ConfigError):
        cfg.validate()
    RunConfig(ranking_weights={"gc_proximity": 0.5}).validate()  # valid subset


def _hit(sseqid, pident, length):
    return {"sseqid": sseqid, "pident": pident, "length": length}


def test_off_target_threshold_is_honoured():
    # A 88%-identity off-target hit in the LEGACY local-alignment mode:
    # rejected at the default 80% floor, but tolerated if the caller raises the
    # off-target floor to 90%. off_target_pident_min only decides anything in
    # that mode — the default full-half mode counts mismatches instead.
    hits = {"p": [_hit("TRUB|S", 99.0, 25), _hit("HUMAN|x", 88.0, 24)]}
    legacy = dict(off_target_full_half=False)

    c1 = ProbeCandidate(probe_id="p", sequence="A" * 50)
    assert batch_evaluate([c1], hits, {"TRUB"}, **legacy) == 0   # default floor
    assert c1.specificity_status == "failed_off_target"

    c2 = ProbeCandidate(probe_id="p", sequence="A" * 50)
    n = batch_evaluate([c2], hits, {"TRUB"},
                       off_target_pident_min=90.0, **legacy)
    assert n == 1 and c2.specificity_pass                     # relaxed: passes


def test_off_target_min_mismatch_is_honoured():
    # Full-half mode (the default). The count is a real count over the whole
    # half, so the hit is built from actual sequence: a 25 nt
    # off-target half with 4 mismatches, reported by BLAST over its full length.
    # 4 is one short of the 5 the 10x rule requires, so it disqualifies;
    # demanding only 4 lets it through.
    #
    # The off-target hit carries no __L/__R suffix, so it is charged to BOTH
    # halves — which is what makes it disqualifying under the both-halves rule.
    query = "ACGTACGTACGTACGTACGTACGTA"
    subject = "TCGTTCGTACGTACGAACGTACGTG"      # differs at 1, 5, 16, 25
    assert sum(1 for a, b in zip(query, subject) if a != b) == 4
    subjects = FakeSubjects({"HUMAN|x": subject})
    hits = {"p": [_hit("TRUB|S", 99.0, 25),
                  aligned_hit("HUMAN|x", query, subject, 1, 25, probe="p")]}

    c1 = ProbeCandidate(probe_id="p", sequence=query * 2)
    assert batch_evaluate([c1], hits, {"TRUB"}, subjects=subjects) == 0
    assert c1.specificity_status == "failed_off_target"
    assert c1.offtarget_recount_basis == "exact"

    c2 = ProbeCandidate(probe_id="p", sequence=query * 2)
    n = batch_evaluate([c2], hits, {"TRUB"}, off_target_min_mismatch=4,
                       subjects=subjects)
    assert n == 1 and c2.specificity_pass
