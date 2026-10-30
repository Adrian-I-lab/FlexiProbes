"""L3: a target's own isoforms must not be counted as off-target."""
from flexprobe.models import ProbeCandidate
from flexprobe.specificity import batch_evaluate, evaluate_specificity


def _hit(sseqid, pident, length):
    return {"sseqid": sseqid, "pident": pident, "length": length}


def test_isoform_hit_not_treated_as_off_target():
    # Probe is on-target for GENEX-001 but also strongly hits GENEX-002, an
    # isoform of the same gene. Without isoform awareness that's a rejection.
    hits = [_hit("GENEX-001|tx", 99.0, 25), _hit("GENEX-002|tx", 97.0, 25)]

    c1 = ProbeCandidate(probe_id="p", sequence="A" * 50)
    evaluate_specificity(c1, hits, {"GENEX-001"})
    assert not c1.specificity_pass
    assert c1.specificity_status == "failed_off_target"

    c2 = ProbeCandidate(probe_id="p", sequence="A" * 50)
    evaluate_specificity(c2, hits, {"GENEX-001"},
                         extra_on_target_ids={"GENEX-002"})
    assert c2.specificity_pass and c2.specificity_status == "passed"


def test_extra_on_target_does_not_satisfy_required_strains():
    # An isoform excused from off-target must NOT count toward the EXACT_SET
    # requirement of the actual target strains.
    hits = [_hit("OTHER_ISO|x", 99.0, 25)]
    c = ProbeCandidate(probe_id="p", sequence="A" * 50)
    evaluate_specificity(c, hits, {"REAL_TARGET"},
                         extra_on_target_ids={"OTHER_ISO"},
                         require_on_target=True)   # strict: REAL_TARGET must be hit
    assert not c.specificity_pass            # no on-target hit to REAL_TARGET


def test_batch_evaluate_forwards_extra_on_target():
    hits = {"p": [_hit("GX-1|a", 99.0, 25), _hit("GX-2|a", 96.0, 24)]}
    c = ProbeCandidate(probe_id="p", sequence="A" * 50)
    n = batch_evaluate([c], hits, {"GX-1"}, extra_on_target_ids={"GX-2"})
    assert n == 1
