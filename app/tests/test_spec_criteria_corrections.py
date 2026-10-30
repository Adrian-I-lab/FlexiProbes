"""Regression tests for the three 10x-spec corrections.

Each test pins one criterion to its verbatim source sentence, and each also
pins the opt-out that reproduces the pre-correction behaviour.

Sources (both on disk under Flex_probe_design_original_code/10X_documentation):
  • CG000621 Rev D — Visium spatial (p.2-3) + Chromium Flex v1 (p.12-13)
  • CG000839 Rev A — GEM-X Flex v2 (p.2)
"""
import pytest

from conftest import FakeSubjects, aligned_hit
from flexprobe._utils import safe_tm
from flexprobe.chemistry import get_chemistry, resolve_tn_mode
from flexprobe.constants import (
    TM_MIN_DEFAULT, TN_MODE_OFF, TN_MODE_SOFT, TN_MODE_STRICT_A,
)
from flexprobe.errors import ConfigError
from flexprobe.filters import evaluate_candidate, new_filter_counters
from flexprobe.models import (
    FilterThresholds, ModalityDefinition, ProbeCandidate,
)
from flexprobe.specificity import (
    evaluate_specificity, full_half_mismatches, hit_half,
    off_target_subject_key,
)


def _cand(seq="A" * 50, pid="p1"):
    # ``sequence`` is load-bearing for the full-half recount, which slices the
    # probe's two 25 nt halves out of it (specificity.full_half_mismatches).
    return ProbeCandidate(probe_id=pid, sequence=seq)


def _hit(sseqid, pident, length, *, half=None, probe="p1"):
    """A BLAST hit dict. ``half`` ('L'/'R') sets the __L/__R query-id suffix
    the pipeline uses to attribute a hit to one probe half."""
    qseqid = f"{probe}__{half}" if half else probe
    return {"qseqid": qseqid, "sseqid": sseqid,
            "pident": pident, "length": length,
            "bitscore": 50.0, "evalue": 1e-5}


def _run_filters(seq, thresholds):
    """Run the filter funnel on one 50-mer; return (candidate, counters)."""
    pc = ProbeCandidate(probe_id="p", sequence=seq)
    counters = new_filter_counters()
    evaluate_candidate(pc, thresholds, counters)
    return pc, counters


# ---------------------------------------------------------------------------
# Task 1 — the off-target rule is a disjunction over the two probe halves
#
#   CG000621 Rev D / CG000839 Rev A, verbatim:
#   "Matches to off-target genes should have at least five mismatches in at
#    least one of the LHS or RHS probes to prevent efficient hybridization."
# ---------------------------------------------------------------------------

def test_lhs_off_target_with_clean_rhs_now_passes():
    """THE headline fix: LHS hits an off-target, RHS does not -> the pair PASSES.

    RTL only yields a library molecule when both halves hybridise and the
    ligase seals the nick between them. A clean RHS means no duplex, no nick,
    no ligation product, no false signal — so the pair is specific and the old
    code rejecting it was inverting the published rule.
    """
    c = _cand()
    hits = [
        _hit("OTHERGENE", 100.0, 25, half="L"),   # LHS: perfect off-target
        _hit("TRUB|S", 100.0, 25, half="R"),      # RHS: on-target only
    ]
    evaluate_specificity(c, hits, {"TRUB"})
    assert c.specificity_pass is True
    assert c.specificity_status == "passed"


def test_rhs_off_target_with_clean_lhs_now_passes():
    """Symmetric to the above — the rule names either half."""
    c = _cand()
    hits = [
        _hit("TRUB|S", 100.0, 25, half="L"),
        _hit("OTHERGENE", 100.0, 25, half="R"),
    ]
    evaluate_specificity(c, hits, {"TRUB"})
    assert c.specificity_pass is True


def test_both_halves_off_target_still_rejected():
    """The rule only excuses a pair when at least ONE half is >=5 mismatches
    away. Both halves close to the same off-target => a real ligation product
    can form => reject."""
    c = _cand()
    hits = [
        _hit("OTHERGENE", 100.0, 25, half="L"),
        _hit("OTHERGENE", 100.0, 25, half="R"),
    ]
    evaluate_specificity(c, hits, {"TRUB"})
    assert c.specificity_pass is False
    assert c.specificity_status == "failed_off_target"


def test_either_half_mode_reproduces_previous_behaviour():
    """The opt-out restores the pre-correction verdict exactly."""
    c = _cand()
    hits = [_hit("OTHERGENE", 100.0, 25, half="L"),
            _hit("TRUB|S", 100.0, 25, half="R")]
    evaluate_specificity(c, hits, {"TRUB"},
                         off_target_require_both_halves=False)
    assert c.specificity_pass is False
    assert c.specificity_status == "failed_off_target"


def test_unattributable_hit_is_charged_to_both_halves():
    """A hit with no __L/__R suffix cannot be attributed. It must be treated
    conservatively (as failing both halves), never as a free pass — otherwise a
    legacy cached hit file would silently weaken the screen."""
    c = _cand()
    evaluate_specificity(c, [_hit("OTHERGENE", 100.0, 25)], {"TRUB"})
    assert c.specificity_pass is False


def test_hit_half_reads_the_query_suffix():
    assert hit_half(_hit("x", 100.0, 25, half="L")) == "L"
    assert hit_half(_hit("x", 100.0, 25, half="R")) == "R"
    assert hit_half(_hit("x", 100.0, 25)) == ""


# ---------------------------------------------------------------------------
# Task 1c — the per-subject reading (opt-in, for sensitivity analysis)
#
#   "Matches to off-target GENES should have at least five mismatches in at
#    least one of the LHS or RHS probes..."  — arguably scopes the disjunction
#   per off-target gene: an LHS on gene X and an RHS on gene Y are never
#   adjacent on one molecule, so the ligase has nothing to seal.
# ---------------------------------------------------------------------------

_SAME_SUBJ = dict(off_target_same_subject=True)
_WHOLE_PROBE = dict(off_target_same_subject=False)


def test_halves_hitting_different_genes_pass_by_default():
    """The whole point of the rule. LHS -> gene X, RHS -> gene Y: no single
    molecule carries both, so the ligase has nothing to seal and no library
    molecule forms on either gene. The default (per-gene) reading passes it;
    the pre-0.4.1 whole-probe reading rejects it.

    Direction matters and is asserted both ways: this pair is the 10.3% of the
    real 10x Flex v1 panel that the whole-probe reading discarded (17.7% ->
    7.4% rejected once scoped per gene), so a regression here is not a cosmetic
    default change — it silently removes a tenth of the designable sites.
    """
    hits = [_hit("GENE_X", 100.0, 25, half="L"),
            _hit("GENE_Y", 100.0, 25, half="R")]

    default = _cand()
    evaluate_specificity(default, hits, {"TRUB"})
    assert default.specificity_pass is True

    whole_probe = _cand()
    evaluate_specificity(whole_probe, hits, {"TRUB"}, **_WHOLE_PROBE)
    assert whole_probe.specificity_pass is False


def test_same_subject_reading_still_rejects_a_shared_off_target():
    """The flag must not weaken the real case: both halves on ONE gene can
    ligate and must still be rejected."""
    c = _cand()
    evaluate_specificity(
        c,
        [_hit("GENE_X", 100.0, 25, half="L"),
         _hit("GENE_X", 100.0, 25, half="R")],
        {"TRUB"}, **_SAME_SUBJ,
    )
    assert c.specificity_pass is False
    assert c.specificity_status == "failed_off_target"


@pytest.mark.parametrize("lhs_id,rhs_id", [
    ("KR013237.1", "KR013237.2"),                 # accession version differs
    ("ENSG001|ENST009", "ENSG001|ENST123"),       # different isoforms, one gene
    ("OROV_X|Large|ACC1", "OROV_X|Small|ACC2"),   # same strain, two segments
    ("GENE_X", "gene_x"),                         # case
])
def test_subject_grouping_is_robust_to_id_variants(lhs_id, rhs_id):
    """The same gene under different accessions/versions must land in ONE
    group, else a genuinely non-specific pair would slip through."""
    assert off_target_subject_key(lhs_id) == off_target_subject_key(rhs_id)
    c = _cand()
    evaluate_specificity(
        c,
        [_hit(lhs_id, 100.0, 25, half="L"), _hit(rhs_id, 100.0, 25, half="R")],
        {"TRUB"}, **_SAME_SUBJ,
    )
    assert c.specificity_pass is False, \
        "id variants of one subject must not split into two groups"


def test_same_subject_reading_is_failsafe_on_an_unknown_half():
    """Attribution gap 1: a hit with no __L/__R suffix is charged to both
    halves, so it alone rejects the pair even under the per-subject reading."""
    c = _cand()
    evaluate_specificity(c, [_hit("GENE_X", 100.0, 25)], {"TRUB"}, **_SAME_SUBJ)
    assert c.specificity_pass is False


@pytest.mark.parametrize("unknown_id", ["", "   "])
def test_same_subject_reading_is_failsafe_on_an_unknown_subject(unknown_id):
    """Attribution gap 2: a hit whose SUBJECT is unidentifiable is unioned into
    every group, so it can only make the screen stricter. Here the LHS hits an
    identified gene and the RHS an unidentifiable one — the pair must reject,
    because that unknown subject might be the very same gene."""
    c = _cand()
    evaluate_specificity(
        c,
        [_hit("GENE_X", 100.0, 25, half="L"),
         _hit(unknown_id, 100.0, 25, half="R")],
        {"TRUB"}, **_SAME_SUBJ,
    )
    assert c.specificity_pass is False


def test_same_subject_flag_is_inert_in_either_half_mode():
    """Documented interaction: the legacy mode never consults the second half,
    so the per-subject reading cannot apply."""
    hits = [_hit("GENE_X", 100.0, 25, half="L"),
            _hit("GENE_Y", 100.0, 25, half="R")]
    c = _cand()
    evaluate_specificity(c, hits, {"TRUB"},
                         off_target_require_both_halves=False, **_SAME_SUBJ)
    assert c.specificity_pass is False


def test_probe_surviving_only_on_the_per_subject_reading_reports_no_margin():
    """Such a probe has a disqualifying off-target on BOTH halves and no clean
    fallback half, so its reported specificity margin must collapse rather than
    look like a clean probe's."""
    c = _cand()
    evaluate_specificity(
        c,
        [_hit("GENE_X", 100.0, 25, half="L"),
         _hit("GENE_Y", 100.0, 25, half="R")],
        {"TRUB"}, **_SAME_SUBJ,
    )
    assert c.specificity_pass is True
    assert c.specificity_margin == 0.0


def test_same_subject_defaults_on():
    """Per-GENE scoping is the default.

    Pinned because the setting changes which probes pass: scoring the
    off-target disjunction over the whole probe rejects 17.7% of 10x's public
    Flex v1 panel against 7.4% per-gene, so a silent revert would quietly
    discard ~10% of designable sites. See
    constants.OFF_TARGET_SAME_SUBJECT_DEFAULT.
    """
    from flexprobe.models import SpecificityThresholds
    assert SpecificityThresholds().off_target_same_subject is True


# ---------------------------------------------------------------------------
# Task 1b — mismatches are counted over the FULL 25 nt half, not a sub-window
# ---------------------------------------------------------------------------

def test_clean_20nt_core_inside_a_5_mismatch_half_is_not_off_target():
    """The second mis-translation. A 25-mer with 5 mismatches clustered at its
    ends still contains a clean 20 nt core, which the legacy rule (>=80% over
    >=20 bp) scores as a disqualifying off-target. Counting over the full half
    gives 5 mismatches — exactly what the spec asks for — so it must be ignored.

    Six mismatches cannot coexist with a clean 20 nt core: 25 - 6 is 19. So
    this is the 5-mismatch boundary case. Because the count is made against
    real sequence rather than inferred from ``pident``, an inconsistent fixture
    cannot be written down at all: the two strings have to disagree in exactly
    five places.
    """
    # 5 mismatches, 3 at the head and 2 at the tail, leaving query positions
    # 4-23 identical. BLAST reports that 20 nt core and nothing else.
    #
    # The "the unaligned positions are mismatches" step is CHECKED against the
    # retrieved subject rather than assumed. Assuming it is catastrophically
    # wrong for a trimmed core (see
    # test_off_target_full_half_underestimates_mismatches); this case is where
    # the assumption happens to be true.
    query = "AAA" + "CGTACGTACGTACGTACGTA" + "GG"
    subject = "TTT" + "CGTACGTACGTACGTACGTA" + "CC"
    assert len(query) == 25
    assert sum(1 for a, b in zip(query, subject) if a != b) == 5
    subjects = FakeSubjects({"OTHERGENE": subject})

    core = aligned_hit("OTHERGENE", query, subject, 4, 23, half="L")
    assert core["pident"] == 100.0 and core["length"] == 20
    assert full_half_mismatches(core, 25, query_seq=query, subjects=subjects) == 5

    hits = [core, aligned_hit("OTHERGENE", query, subject, 4, 23, half="R")]
    both = _cand(seq=query * 2)
    evaluate_specificity(both, hits, {"TRUB"}, subjects=subjects)
    assert both.specificity_pass is True, \
        "6 mismatches over the full half cannot hybridise efficiently"
    assert both.offtarget_recount_basis == "exact"

    # ...and the legacy local-alignment mode rejects exactly this probe.
    legacy = _cand(seq=query * 2)
    evaluate_specificity(legacy, hits, {"TRUB"}, off_target_full_half=False,
                         subjects=subjects)
    assert legacy.specificity_pass is False
    assert legacy.offtarget_recount_basis == "legacy_local"


def test_four_mismatches_over_the_full_half_is_off_target():
    """The boundary: 4 mismatches is one short of the 5 the spec requires, so
    the half still hybridises and the match disqualifies."""
    h = _hit("OTHERGENE", 84.0, 25, half="L")     # 21/25 identities
    assert full_half_mismatches(h, 25) == 4
    c = _cand()
    evaluate_specificity(
        c,
        [h, _hit("OTHERGENE", 84.0, 25, half="R")],
        {"TRUB"},
    )
    assert c.specificity_pass is False


def test_off_target_pushdown_is_never_stricter_than_the_python_rule():
    """The BLAST pushdown exists to cut hit volume, not to decide. It must not
    demand more coverage/identity than the Python classifier would count —
    notably it must NOT be -qcov_hsp_perc 100, which would hide the end-trimmed
    HSPs that carry a genuine 4-mismatch off-target."""
    from flexprobe.models import SpecificityThresholds
    from flexprobe.pipeline import _derive_qcov_hsp_perc, _off_target_pushdown

    from flexprobe.constants import BLASTN_SHORT_MISMATCH_COST

    pident, aln = _off_target_pushdown(SpecificityThresholds(), 25)
    # 25 - COST*(5-1), DERIVED rather than typed. NOT 21: BLAST reports the
    # maximum-scoring LOCAL segment, so the floor has to be low enough to keep
    # a trimmed core. COST is reward + |penalty| and the screen sets both, so
    # this expectation tracks a scoring change instead of going stale against
    # it. At the current -penalty -1 that is 25 - 2*4 = 17; at blastn-short's
    # own -3 it was 25 - 4*4 = 9. See the exhaustive proof in
    # test_pushdown_keeps_every_hit_the_python_rule_would_reject.
    assert aln == 25 - BLASTN_SHORT_MISMATCH_COST * 4
    assert pident < 84.0                  # slack for gapped HSPs
    qcov = _derive_qcov_hsp_perc(aln, 25)
    assert qcov is not None and qcov < 100.0


def test_pushdown_keeps_every_hit_the_python_rule_would_reject():
    """EXHAUSTIVE: no disqualifying hit can be dropped server-side.

    This is the check that makes the recount worth having. If the BLAST floors
    discard the trimmed cores, ``full_half_mismatches`` never sees them and the
    whole fix is inert — which is exactly what the previous 21 nt floor did.

    For every arrangement of up to 4 mismatches in a 25-mer (all 15,275 of them,
    m = 0..4), work out the maximum-scoring local segment BLAST would report
    under the scoring the screen actually asks blastn for (imported below,
    not restated), then assert that segment clears
    both the length floor and the identity floor. Brute force over the real
    scoring, not a restatement of the formula the code uses.
    """
    from itertools import combinations

    from flexprobe.models import SpecificityThresholds
    from flexprobe.pipeline import _derive_qcov_hsp_perc, _off_target_pushdown

    # Imported, NOT restated. These three tests hardcoded penalty -3 and so
    # kept simulating blastn-short's defaults after the screen stopped using
    # them, which turned an exhaustive safety proof into a proof about a
    # configuration nothing runs.
    from flexprobe.constants import BLASTN_PENALTY, BLASTN_REWARD

    HALF, REWARD, PENALTY = 25, BLASTN_REWARD, BLASTN_PENALTY
    spec = SpecificityThresholds()
    pident_floor, aln_floor = _off_target_pushdown(spec, HALF)
    qcov_floor = _derive_qcov_hsp_perc(aln_floor, HALF)

    def best_local_segment(mismatch_at):
        """(length, identities) of the maximum-scoring local segment."""
        best = (0, 0, 0)          # (score, length, identities)
        for i in range(HALF):
            for j in range(i, HALF):
                n_mm = sum(1 for p in mismatch_at if i <= p <= j)
                length = j - i + 1
                score = (length - n_mm) * REWARD + n_mm * PENALTY
                if score > best[0]:
                    best = (score, length, length - n_mm)
        return best[1], best[2]

    checked = 0
    for m in range(spec.off_target_min_mismatch):      # 0..4 == disqualifying
        for positions in combinations(range(HALF), m):
            length, identities = best_local_segment(positions)
            checked += 1
            assert length >= aln_floor, (
                f"{m} mismatches at {positions} -> a {length} nt HSP, which the "
                f"{aln_floor} nt length floor would DROP server-side"
            )
            assert 100.0 * identities / length >= pident_floor, (
                f"{m} mismatches at {positions} -> {identities}/{length} "
                f"identity, below the {pident_floor:.1f}% floor"
            )
            assert 100.0 * length / HALF >= qcov_floor, (
                f"{m} mismatches at {positions} -> {length}/{HALF} coverage, "
                f"below the {qcov_floor:.1f}% qcov floor"
            )
    # C(25,0) + C(25,1) + ... + C(25,4) = 1 + 25 + 300 + 2300 + 12650.
    assert checked == 15276, "the sweep must actually be exhaustive"

    # Confirm the check COULD fail. Run the SAME sweep against the floor this
    # replaced. A guard that has not been watched fail is not a guard.
    #
    # The COUNT depends on the scoring, so it is recorded per penalty rather
    # than written as one magic number. On this exhaustive sweep:
    #   penalty -3, blastn-short's default : 8,259 of 15,276, being 54%
    #   penalty -1, what this screen asks  :   350 of 15,276, being 2.3%
    # The fall is the POINT of the change, not a weakening of the control. At
    # -1 the full-length alignment outscores short exact cores, so BLAST
    # reports something 21 nt or longer far more often and the old floor would
    # have discarded far less. The control still fires, which is its whole job.
    dropped = [
        p for m in range(spec.off_target_min_mismatch)
        for p in combinations(range(HALF), m)
        if best_local_segment(p)[0] < 21          # the old `keep` floor
    ]
    assert dropped, (
        "the old 21 nt floor must be shown to drop SOME disqualifying hit, "
        "otherwise this sweep cannot tell a good floor from a bad one"
    )
    expected = {-3: 8259, -1: 350}.get(PENALTY)
    if expected is not None:
        assert len(dropped) == expected, (
            f"at penalty {PENALTY} the old 21 nt floor dropped "
            f"{len(dropped)}/{checked}, expected {expected}. If the scoring "
            "changed deliberately, record the new number here."
        )


# ---------------------------------------------------------------------------
# Task 2 — the TN junction gate is chemistry-aware, hard by default
#
#   CG000621 Rev D (Visium p.3, Chromium Flex p.12), verbatim:
#   "The 25th nucleotide of the probe (3' most nucleotide of the LHS probe)
#    must be a T. The opposing nucleotide in the target RNA must be an A."
#   CG000839 Rev A (GEM-X Flex v2, p.2), verbatim:
#   "While preferred, the TN motif is not always strictly required, and other
#    nucleotide motifs can also function effectively. However, for a simplified
#    and robust design process, 10x Genomics recommends using the T at the
#    specified position."
# ---------------------------------------------------------------------------

# 50-mers differing only at index 25 — the base opposite the LHS 3' end.
_J = {b: ("ACGTACGTACGTACGTACGTACGTA" + b + "CGTACGTACGTACGTACGTACGTA")
      for b in "ACGT"}


@pytest.mark.parametrize("base", list("ACGT"))
def test_junction_base_is_recorded_for_ranking(base):
    """Both fields the ranking layer needs are populated whatever gate runs —
    that is what lets A-junction probes still sort to the top under `soft`."""
    pc, _ = _run_filters(_J[base], FilterThresholds(tn_junction_mode=TN_MODE_OFF))
    assert pc.junction_base == base
    assert pc.tn_junction_ok is (base == "A")


@pytest.mark.parametrize("base,expect_pass", [("A", True), ("C", False),
                                              ("G", False), ("T", False)])
def test_strict_gate_admits_only_adenine(base, expect_pass):
    thr = FilterThresholds(tn_junction_mode=TN_MODE_STRICT_A)
    pc, counters = _run_filters(_J[base], thr)
    assert ("tn_junction" in pc.filter_fail_reasons) is (not expect_pass)
    assert counters["tn_junction"] == (0 if expect_pass else 1)


@pytest.mark.parametrize("base,expect_pass", [("A", True), ("C", False),
                                              ("G", True), ("T", True)])
def test_soft_gate_rejects_only_cytosine(base, expect_pass):
    """The opt-in relaxation: A/G/T admitted, C rejected. Matches what 10x's own
    panels actually do (A ~51%, G ~36%, T ~12%, C ~1%)."""
    thr = FilterThresholds(tn_junction_mode=TN_MODE_SOFT)
    pc, _ = _run_filters(_J[base], thr)
    assert ("tn_junction" in pc.filter_fail_reasons) is (not expect_pass)


def test_default_gate_is_hard_not_soft():
    """Owner decision: hard is the default wherever the gate is engaged; the
    relaxation must always be asked for explicitly."""
    for version in ("v1", "v2"):
        chem = get_chemistry(ModalityDefinition(name="flex", version=version))
        thr = FilterThresholds(require_tn_junction=True)
        assert resolve_tn_mode(thr, chem) == TN_MODE_STRICT_A


def test_bare_thresholds_enforce_the_junction_adenine():
    """A bare FilterThresholds() now gates on the junction adenine.

    The default is True, matching every other entry point: the GUI modality
    presets set it True too. A False default would put the library and CLI at
    odds with the GUI, and a run from either would not reproduce a GUI run.
    """
    chem = get_chemistry(ModalityDefinition(name="flex", version="v2"))
    assert resolve_tn_mode(FilterThresholds(), chem) == TN_MODE_STRICT_A
    for base in "CGT":
        pc, _ = _run_filters(_J[base], FilterThresholds())
        assert "tn_junction" in pc.filter_fail_reasons
    pc, _ = _run_filters(_J["A"], FilterThresholds())
    assert "tn_junction" not in pc.filter_fail_reasons


def test_junction_gate_can_still_be_disabled():
    chem = get_chemistry(ModalityDefinition(name="flex", version="v2"))
    thr = FilterThresholds(require_tn_junction=False)
    assert resolve_tn_mode(thr, chem) == TN_MODE_OFF
    for base in "ACGT":
        pc, _ = _run_filters(_J[base], thr)
        assert "tn_junction" not in pc.filter_fail_reasons


def test_soft_gate_is_available_for_flex_only():
    for version in ("v1", "v2"):
        chem = get_chemistry(ModalityDefinition(name="flex", version=version))
        thr = FilterThresholds(tn_junction_mode=TN_MODE_SOFT)
        assert resolve_tn_mode(thr, chem) == TN_MODE_SOFT


def test_soft_gate_refused_for_visium():
    """CG000621's Visium section states the junction adenine flatly and neither
    document relaxes it for spatial, so asking for soft is an error rather than
    a silent downgrade."""
    chem = get_chemistry(ModalityDefinition(name="visium", version="v1",
                                            spatial=True))
    thr = FilterThresholds(tn_junction_mode=TN_MODE_SOFT)
    with pytest.raises(ConfigError, match="not available"):
        resolve_tn_mode(thr, chem)


def test_explicit_mode_wins_over_the_require_default():
    """require_tn_junction now DEFAULTS to True, so it is no longer evidence of
    intent. If resolve_tn_mode treated the combination as contradictory, the
    soft mode would be unreachable from a default-constructed config."""
    chem = get_chemistry(ModalityDefinition(name="flex", version="v2"))
    thr = FilterThresholds(tn_junction_mode=TN_MODE_SOFT)
    assert thr.require_tn_junction is True          # the default
    assert resolve_tn_mode(thr, chem) == TN_MODE_SOFT


def test_cli_refuses_an_explicitly_contradictory_junction_request():
    """The CLI *can* tell an explicitly-passed flag from a default, so that is
    where the genuine contradiction is caught."""
    from flexprobe.cli import main
    with pytest.raises(SystemExit) as exc:
        main(["design", "--fasta", "x.fa", "--out", "o", "--skip-screen",
              "--require-tn-junction", "--tn-junction-mode", "soft"])
    assert exc.value.code == 2                      # argparse usage error


def test_cli_default_enforces_the_junction_and_opt_out_flag_exists():
    from flexprobe.cli import build_parser
    p = build_parser()
    assert p.parse_args(
        ["design", "--fasta", "x.fa"]).require_tn_junction is None
    assert p.parse_args(
        ["design", "--fasta", "x.fa", "--require-tn-junction"]
    ).require_tn_junction is True
    assert p.parse_args(
        ["design", "--fasta", "x.fa", "--no-require-tn-junction"]
    ).require_tn_junction is False


# ---------------------------------------------------------------------------
# Task 3 — duplex stability is gated on per-half Tm, not the GC band
# ---------------------------------------------------------------------------

# A 50-mer whose halves are GC 0.36 each — below 10x's published 0.44 floor —
# but well above the empirical Tm floor. Exactly the window class the GC band
# discards on AT-rich (bunyavirus) targets.
_LOW_GC_HIGH_TM = "GGCATTAGCATTAGGCATTAGCAAT" "GGCATTAGCATTAGGCATTAGCAAT"

# The GC/Tm tests are about duplex stability only, so they switch the junction
# gate off rather than hand-tuning index 25 — otherwise the (now default) hard
# adenine gate would reject these windows for an unrelated reason.
_NO_TN = dict(require_tn_junction=False)


def test_low_gc_but_thermally_sound_window_now_passes():
    lo, hi = _LOW_GC_HIGH_TM[:25], _LOW_GC_HIGH_TM[25:]
    assert max(lo.count("G") + lo.count("C"),
               hi.count("G") + hi.count("C")) / 25 < 0.44
    assert min(safe_tm(lo), safe_tm(hi)) >= TM_MIN_DEFAULT

    pc, _ = _run_filters(_LOW_GC_HIGH_TM, FilterThresholds(**_NO_TN))
    assert "gc_range_left" not in pc.filter_fail_reasons
    assert "tm_range_left" not in pc.filter_fail_reasons
    assert pc.passed_filters is True
    # GC is still computed — the ranking layer and results.csv both read it.
    assert 0.0 < pc.gc_left < 0.44


def test_gc_band_opt_in_reproduces_previous_rejection():
    thr = FilterThresholds(enforce_gc_range=True, **_NO_TN)
    pc, counters = _run_filters(_LOW_GC_HIGH_TM, thr)
    assert "gc_range_left" in pc.filter_fail_reasons
    assert counters["gc_range"] == 1
    assert pc.passed_filters is False


def test_thermally_weak_window_is_rejected_by_the_tm_gate():
    """An AT-rich half well under the floor must be rejected — demoting GC must
    not leave duplex stability ungated."""
    weak = "ATATATATATATATATATATATATA" "ATATATATATATATATATATATATA"
    assert safe_tm(weak[:25]) < TM_MIN_DEFAULT
    pc, counters = _run_filters(weak, FilterThresholds(max_homopolymer_run=25, **_NO_TN))
    assert "tm_range_left" in pc.filter_fail_reasons
    assert counters["tm_range"] == 1


def test_tm_gate_can_be_disabled():
    weak = "ATATATATATATATATATATATATA" "ATATATATATATATATATATATATA"
    thr = FilterThresholds(tm_min=None, max_homopolymer_run=25, **_NO_TN)
    pc, _ = _run_filters(weak, thr)
    assert not any(r.startswith("tm_range") for r in pc.filter_fail_reasons)


def test_tm_floor_admits_every_probe_10x_ships():
    """The floor is calibrated to 10x's own panels: their observed per-half
    minimum is 56.32 C (human) / 56.35 C (mouse), so 56.0 rejects none of them.
    Guards against someone 'tidying' the constant upward."""
    assert TM_MIN_DEFAULT <= 56.32


def test_ambiguous_base_is_not_silently_rejected_by_the_tm_gate():
    """primer3 cannot score an ambiguity code, so safe_tm returns NaN. NaN means
    'cannot judge' and must not be read as 'too cold'."""
    amb = "GGCATTAGCATTAGGCATTAGCANT" "GGCATTAGCATTAGGCATTAGCAAT"
    thr = FilterThresholds(max_n_per_half=1, **_NO_TN)
    pc, _ = _run_filters(amb, thr)
    assert not any(r.startswith("tm_range") for r in pc.filter_fail_reasons)


# ---------------------------------------------------------------------------
# Task 4 — the strain-unique gate counts its mismatches per half, like the
#          off-target gate it borrows the number from
# ---------------------------------------------------------------------------

# ``RunConfig.strain_unique_min_mismatches`` defaults to 5 *because* of the same
# CG000621/CG000839 sentence pinned at the top of this file. That sentence scopes
# the 5 to a probe half, so scoping it to the whole 50-mer instead put the two
# gates in contradiction: a window with 3 LHS + 2 RHS mismatches was called
# strain-unique, while evaluate_specificity — correctly — treats a half with
# only 3 or 4 mismatches as still able to hybridise. Both halves bind, the
# ligase seals the nick, and the probe reports on the strain it was meant to
# exclude.


def test_strain_unique_default_counts_per_half():
    from flexprobe.constants import STRAIN_UNIQUE_PER_HALF_DEFAULT
    from flexprobe.models import RunConfig
    assert STRAIN_UNIQUE_PER_HALF_DEFAULT is True
    assert RunConfig().strain_unique_per_half is True


def test_strain_unique_gate_agrees_with_the_off_target_gate():
    """The two gates must give the same verdict on the same half.

    3 LHS + 2 RHS is 5 mismatches over the window. The off-target classifier
    calls each half a live hit (< 5 mismatches each), so the pair still ligates
    on that sequence — therefore the strain-unique gate must not accept it.
    """
    from flexprobe.constants import OFF_TARGET_MIN_MISMATCH
    from flexprobe.consensus import (
        build_consensus, score_probe_against_all_strains,
        strain_unique_start_positions,
    )
    ref = ("ACGT" * 13)[:50]
    other = list(ref)
    for p in (2, 7, 13, 30, 41):                 # 3 in the LHS, 2 in the RHS
        other[p] = "T" if other[p] != "T" else "A"
    aln = build_consensus([("A", ref), ("B", "".join(other))])

    # Per-half mismatch counts, from the same scorer that fills the
    # lhs_mm__* / rhs_mm__* columns of results.csv.
    score = {s.strain_id: s for s in score_probe_against_all_strains(
        aln, 0, tile_length=50, probe_seq=ref)}["B"]
    assert (score.lhs_mismatches, score.rhs_mismatches) == (3, 2)
    assert score.total_mismatches == 5           # 5 over the window as a whole

    # Both halves sit under the off-target floor, so the off-target side calls
    # each of them a live hit: the pair hybridises and ligates on strain B.
    assert score.lhs_mismatches < OFF_TARGET_MIN_MISMATCH
    assert score.rhs_mismatches < OFF_TARGET_MIN_MISMATCH

    # So the strain-unique side must agree it is not discriminating.
    assert strain_unique_start_positions(
        aln, target_strain_idx=0, tile_length=50, tile_step=1,
        min_mismatches_vs_others=OFF_TARGET_MIN_MISMATCH,
    ) == []


def test_whole_window_mode_reproduces_previous_behaviour():
    from flexprobe.consensus import (
        build_consensus, strain_unique_start_positions,
    )
    ref = ("ACGT" * 13)[:50]
    other = list(ref)
    for p in (2, 7, 13, 30, 41):
        other[p] = "T" if other[p] != "T" else "A"
    aln = build_consensus([("A", ref), ("B", "".join(other))])
    assert strain_unique_start_positions(
        aln, target_strain_idx=0, tile_length=50, tile_step=1,
        min_mismatches_vs_others=5, per_half=False,
    ) == [(0, 5)]


def test_cli_defaults_to_per_half_and_opt_out_flag_exists():
    from flexprobe.cli import _build_run_config, build_parser
    p = build_parser()
    assert _build_run_config(
        p.parse_args(["design", "--fasta", "x.fa"])
    ).strain_unique_per_half is True
    assert _build_run_config(
        p.parse_args(["design", "--fasta", "x.fa",
                      "--strain-unique-whole-window"])
    ).strain_unique_per_half is False


def test_per_half_setting_round_trips_through_a_saved_config(tmp_path):
    """Provenance: run_config.json must record which reading a run used, so a
    panel designed either way can be told apart after the fact."""
    import json
    from flexprobe.config_io import load_config, save_config
    from flexprobe.models import RunConfig
    path = tmp_path / "run_config.json"
    save_config(RunConfig(fasta_path="x.fa", strain_unique_per_half=False),
                str(path))
    assert json.loads(path.read_text())["strain_unique_per_half"] is False
    assert load_config(str(path)).strain_unique_per_half is False


def test_threshold_above_the_half_length_is_rejected_not_silently_empty():
    """In per-half mode the count is capped by the half, so a threshold of 30 on
    a 50-mer can never be met — that must be a config error, not a run that
    quietly yields zero probes."""
    from flexprobe.models import RunConfig
    cfg = RunConfig(fasta_path="x.fa", multi_strain=True,
                    consensus_strategy="strain_unique",
                    strain_unique_min_mismatches=30)
    with pytest.raises(ConfigError) as exc:
        cfg.validate()
    assert "exceeds the 25 bp probe half" in str(exc.value)

    # Same threshold is legitimate under the whole-window reading.
    cfg.strain_unique_per_half = False
    cfg.validate()                               # must not raise


# ---------------------------------------------------------------------------
# The recount is the only screen — there is no cheaper mode to fall back to
# ---------------------------------------------------------------------------

def test_pushdown_floor_is_always_the_trimmed_core_floor():
    """One floor, and it is the one that can see a trimmed core.

    A 21 nt floor is not offered. Checked against exhaustive ground truth, it
    discards a large share of usable probes, because demanding >=72% query
    coverage makes BLAST report a stretched, shifted alignment whose implied
    offset misplaces the reconstructed flanks. If this ever reads 21 again,
    that floor is back.
    """
    from flexprobe.constants import BLASTN_SHORT_MISMATCH_COST
    from flexprobe.models import SpecificityThresholds
    from flexprobe.pipeline import _off_target_pushdown
    _, aln = _off_target_pushdown(SpecificityThresholds(), 25)
    # Derived from the screen's own scoring, so a penalty change moves this
    # with it. The regression this guards against is the floor reading 21.
    assert aln == 25 - BLASTN_SHORT_MISMATCH_COST * 4
    assert aln != 21


def test_there_is_no_recount_mode_knob_left_to_set():
    """A knob whose every setting but one is known wrong is a trap."""
    from flexprobe.models import SpecificityThresholds
    assert not hasattr(SpecificityThresholds(), "off_target_recount_mode")


def test_large_screen_still_warns_about_its_cost():
    from flexprobe.models import SpecificityThresholds
    from flexprobe.pipeline import _warn_recount_cost
    from flexprobe.constants import RECOUNT_EXACT_WARN_HALVES
    spec = SpecificityThresholds()
    assert _warn_recount_cost(spec, RECOUNT_EXACT_WARN_HALVES - 1) is None
    msg = _warn_recount_cost(spec, RECOUNT_EXACT_WARN_HALVES)
    assert msg is not None and "hit volume" in msg
