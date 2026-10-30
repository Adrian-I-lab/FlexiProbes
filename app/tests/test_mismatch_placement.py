"""Graded positional scoring of ON-target mismatches.

Two things are under test here:

* ``constants.mismatch_position_weight`` — the linear ramp that replaced the
  binary ``DTM_JUNCTION_BONUS`` step in ``consensus._binding_score_for_window``.
  Its whole claim is that it grades the interior of the probe WITHOUT moving
  either endpoint, so the endpoint checks below are the load-bearing ones.
* ``ranking``'s seventh sub-score, which finally lets that placement reach the
  composite score instead of only results.csv.

The hard requirement running through all of it: a single-sequence run cannot
have on-target mismatches, so none of this may move a single-sequence score.
"""
import pytest

from flexprobe.constants import (
    DTM_PER_MISMATCH, JUNCTION_WINDOW,
    MISMATCH_WEIGHT_AT_END, MISMATCH_WEIGHT_AT_NICK,
    mismatch_distance_from_nick, mismatch_position_weight,
)
from flexprobe.consensus import _binding_score_for_window
from flexprobe.models import ProbeCandidate, RunConfig
from flexprobe.ranking import (
    DEFAULT_WEIGHTS, MISMATCH_PLACEMENT_DTM_FLOOR, RankingWeights,
    _central_placement_score, _gc_balance_score, _gc_proximity_per_half_score,
    _homopolymer_penalty_score, _junction_base_multiplier,
    _on_target_mismatch_placement_score, _specificity_margin_score,
    compute_composite_score,
)

HALF = 25
LENGTH = 2 * HALF


def _window(mismatch_at=()):
    """A 50-mer with 'A' at index 25 (so the TN-junction gate is inert)."""
    seq = list(("ACGTG" * 10)[:LENGTH])
    seq[HALF] = "A"
    return "".join(seq)


def _mutate(seq, positions):
    out = list(seq)
    for p in positions:
        out[p] = "T" if out[p] != "T" else "G"
    return "".join(out)


def _candidate(mismatch_positions=None, **kw):
    """A candidate carrying one strain's worth of on-target mismatch data."""
    c = ProbeCandidate(
        probe_id=kw.pop("probe_id", "p"), sequence=_window(),
        gc_content=0.52, gc_left=0.52, gc_right=0.52,
        homopolymer_max=2, specificity_margin=40.0, start=100, end=150, **kw,
    )
    if mismatch_positions is not None:
        c.strain_binding = [{
            "strain_id": "sB",
            "total_mismatches": len(mismatch_positions),
            "mismatch_positions": list(mismatch_positions),
        }]
    return c


# ---------------------------------------------------------------------------
# (iii) The ramp preserves today's endpoints exactly.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("position", [HALF - 1, HALF])       # 24 and 25
def test_mismatch_at_the_nick_still_costs_six_degrees(position):
    """The old cost at the junction was DTM_PER_MISMATCH + DTM_JUNCTION_BONUS."""
    seq = _window()
    score = _binding_score_for_window(
        seq, _mutate(seq, [position]), "sB", half_length=HALF)
    assert score.total_mismatches == 1
    assert score.estimated_dTm == pytest.approx(-6.0)


@pytest.mark.parametrize("position", [0, LENGTH - 1])        # 0 and 49
def test_mismatch_at_the_probe_end_still_costs_two_degrees(position):
    """The old cost outside the junction window was DTM_PER_MISMATCH alone."""
    seq = _window()
    score = _binding_score_for_window(
        seq, _mutate(seq, [position]), "sB", half_length=HALF)
    assert score.total_mismatches == 1
    assert score.estimated_dTm == pytest.approx(-2.0)


def test_weight_endpoints_are_the_two_old_constants():
    assert mismatch_position_weight(HALF - 1, HALF) == MISMATCH_WEIGHT_AT_NICK
    assert mismatch_position_weight(HALF, HALF) == MISMATCH_WEIGHT_AT_NICK
    assert mismatch_position_weight(0, HALF) == MISMATCH_WEIGHT_AT_END
    assert mismatch_position_weight(LENGTH - 1, HALF) == MISMATCH_WEIGHT_AT_END
    assert DTM_PER_MISMATCH * MISMATCH_WEIGHT_AT_NICK == -6.0
    assert DTM_PER_MISMATCH * MISMATCH_WEIGHT_AT_END == -2.0


def test_distance_from_nick_is_symmetric_about_the_junction():
    assert mismatch_distance_from_nick(HALF - 1, HALF) == 0
    assert mismatch_distance_from_nick(HALF, HALF) == 0
    assert mismatch_distance_from_nick(0, HALF) == HALF - 1
    assert mismatch_distance_from_nick(LENGTH - 1, HALF) == HALF - 1
    for d in range(HALF):
        assert mismatch_distance_from_nick(HALF - 1 - d, HALF) == d
        assert mismatch_distance_from_nick(HALF + d, HALF) == d
    # Out-of-window offsets clamp rather than run off the ramp.
    assert mismatch_distance_from_nick(-5, HALF) == HALF - 1
    assert mismatch_distance_from_nick(LENGTH + 5, HALF) == HALF - 1


def test_ramp_is_monotone_and_has_no_cliff_at_the_junction_window_edge():
    """The point of the ramp: a gradient everywhere, a step nowhere."""
    weights = [mismatch_position_weight(p, HALF) for p in range(LENGTH)]
    assert all(0.99 <= w <= 3.01 for w in weights)
    # Rises monotonically to the nick across the LHS half, falls across the RHS.
    assert weights[:HALF] == sorted(weights[:HALF])
    assert weights[HALF:] == sorted(weights[HALF:], reverse=True)
    # The old code stepped 3x -> 1x between the last in-window base and the
    # first out-of-window one. Now consecutive bases differ by one ramp step.
    step = (MISMATCH_WEIGHT_AT_NICK - MISMATCH_WEIGHT_AT_END) / (HALF - 1)
    edge = HALF + JUNCTION_WINDOW           # 29: first base the old step spared
    assert weights[edge - 1] - weights[edge] == pytest.approx(step)


def test_grading_does_not_change_any_mismatch_count():
    """Counts are reported unweighted; only the PRICE of a mismatch moved."""
    seq = _window()
    positions = [2, HALF - 2, HALF + 1, LENGTH - 3]
    score = _binding_score_for_window(
        seq, _mutate(seq, positions), "sB", half_length=HALF)
    assert score.total_mismatches == 4
    assert score.mismatch_positions == positions
    assert score.lhs_mismatches == 2 and score.rhs_mismatches == 2
    # bp 21..28 on a 50-mer, exactly as before.
    assert score.junction_mismatches == 2
    assert score.likely_binds is False


# ---------------------------------------------------------------------------
# (ii) Same mismatch COUNT, different placement -> different score.
# ---------------------------------------------------------------------------

def test_end_proximal_mismatch_outscores_a_nick_proximal_one():
    at_end = _candidate([0], probe_id="end")
    at_nick = _candidate([HALF], probe_id="nick")
    s_end = compute_composite_score(at_end, seq_length=1000)
    s_nick = compute_composite_score(at_nick, seq_length=1000)
    # Identical on every other axis, and the same single mismatch.
    assert at_end.strain_binding[0]["total_mismatches"] == \
        at_nick.strain_binding[0]["total_mismatches"]
    assert s_end > s_nick
    # 1 - 2/12 = 0.8333 against 1 - 6/12 = 0.5.
    assert at_end.score_components["mismatch_placement"] == pytest.approx(5 / 6)
    assert at_nick.score_components["mismatch_placement"] == pytest.approx(0.5)


def test_score_falls_monotonically_as_one_mismatch_walks_toward_the_nick():
    scores = []
    for p in range(HALF):                       # 0 (probe end) .. 24 (nick)
        c = _candidate([p], probe_id=f"p{p:02d}")
        scores.append(compute_composite_score(c, seq_length=1000))
    assert scores == sorted(scores, reverse=True)
    assert scores[0] > scores[-1]


def test_effect_is_visible_but_bounded_by_its_weight():
    """Worst case cannot move a score by more than its share of the total."""
    clean = _candidate([])
    ruined = _candidate([HALF - 1, HALF])       # 2 nick-adjacent -> the floor
    s_clean = compute_composite_score(clean, seq_length=1000)
    s_ruined = compute_composite_score(ruined, seq_length=1000)
    assert ruined.score_components["mismatch_placement"] == 0.0
    gap = s_clean - s_ruined
    assert gap == pytest.approx(
        DEFAULT_WEIGHTS.mismatch_placement / DEFAULT_WEIGHTS.total())
    # ~9% of the range: visible next to specificity_margin's 0.30, not dominant.
    assert 0.05 < gap < 0.15


def test_worst_strain_decides():
    """A probe is only as good as the strain it binds least well."""
    c = _candidate()
    c.strain_binding = [
        {"strain_id": "sA", "total_mismatches": 0, "mismatch_positions": []},
        {"strain_id": "sB", "total_mismatches": 1,
         "mismatch_positions": [HALF]},
        {"strain_id": "sC", "total_mismatches": 1, "mismatch_positions": [0]},
    ]
    assert _on_target_mismatch_placement_score(c) == pytest.approx(0.5)


def test_no_mismatches_scores_one():
    c = _candidate([])
    assert _on_target_mismatch_placement_score(c) == 1.0


def test_score_bottoms_out_at_the_floor_and_never_goes_negative():
    many = _candidate([HALF - 3, HALF - 2, HALF - 1, HALF, HALF + 1, HALF + 2])
    assert _on_target_mismatch_placement_score(many) == 0.0
    # Sanity: that really is past the documented floor.
    worst = DTM_PER_MISMATCH * sum(
        mismatch_position_weight(p, HALF)
        for p in many.strain_binding[0]["mismatch_positions"]
    )
    assert worst < MISMATCH_PLACEMENT_DTM_FLOOR


# ---------------------------------------------------------------------------
# (i) A single-sequence run is untouched. HARD REQUIREMENT.
# ---------------------------------------------------------------------------

def _legacy_score(c, seq_length, w):
    """The composite score exactly as it was computed before this term existed.

    Six weighted sub-scores over the sum of those six weights, then the
    TN-junction multiplier. Written out longhand on purpose: importing the new
    implementation to check the new implementation would prove nothing.

    ``gc_proximity`` is scored per half (``_gc_proximity_per_half_score`, the
    worst of the two halves) rather than on whole-window ``gc_content`` — that
    is a separate, independently-tested fix (audit M7 / ``test_ranking.py``'s
    per-half tests), not something this "before the seventh term" fixture is
    meant to hold at its pre-fix value. This mirrors exactly what
    ``compute_composite_score`` does for that sub-score.
    """
    six = (
        w.gc_proximity * _gc_proximity_per_half_score(c.gc_left, c.gc_right)
        + w.gc_balance * _gc_balance_score(c.gc_left, c.gc_right)
        + w.specificity_margin * _specificity_margin_score(c.specificity_margin)
        + w.homopolymer_penalty * _homopolymer_penalty_score(c.homopolymer_max)
        + w.central_placement * _central_placement_score(
            c.start, c.end, seq_length)
        + w.isoform_coverage * max(0.0, min(1.0, c.isoform_coverage))
    )
    denom = (w.gc_proximity + w.gc_balance + w.specificity_margin
             + w.homopolymer_penalty + w.central_placement + w.isoform_coverage)
    return _junction_base_multiplier(c.sequence) * six / denom


@pytest.mark.parametrize("gc,gc_l,gc_r,hp,margin,iso,start,end,seq_len", [
    (0.52, 0.52, 0.52, 2, 40.0, 0.0, 100, 150, 1000),
    (0.44, 0.30, 0.70, 7, 0.0, 1.0, 0, 50, 1000),
    (0.72, 0.72, 0.72, 4, 95.0, 0.5, 950, 1000, 1000),
    (0.10, 1.00, 0.00, 40, -5.0, 0.0, 10, 60, 100),
])
def test_single_sequence_score_is_bit_identical(
    gc, gc_l, gc_r, hp, margin, iso, start, end, seq_len,
):
    """No strain_binding == no on-target mismatch data == term does not exist.

    Bit-identical, not approximately equal: the term is dropped from the
    numerator AND the denominator, so the arithmetic is the same arithmetic.
    """
    for weights in (DEFAULT_WEIGHTS, RankingWeights(isoform_coverage=0.20)):
        c = ProbeCandidate(
            probe_id="p", sequence=_window(), gc_content=gc, gc_left=gc_l,
            gc_right=gc_r, homopolymer_max=hp, specificity_margin=margin,
            isoform_coverage=iso, start=start, end=end,
        )
        assert not c.strain_binding            # a single-sequence run
        s = compute_composite_score(c, seq_len, weights)
        assert s == _legacy_score(c, seq_len, weights)
        # Reported as inert, and reported at all, so the CSV column is never
        # blank for a run that simply has nothing to penalise.
        assert c.score_components["mismatch_placement"] == 1.0


def test_single_sequence_score_ignores_the_new_weight_entirely():
    """Turning the weight up to absurd values must not move the score."""
    baseline = None
    for w in (0.0, 0.10, 0.75, 5.0):
        c = ProbeCandidate(
            probe_id="p", sequence=_window(), gc_content=0.6, gc_left=0.5,
            gc_right=0.7, homopolymer_max=3, specificity_margin=30.0,
            start=100, end=150,
        )
        s = compute_composite_score(
            c, 1000, RankingWeights(mismatch_placement=w))
        if baseline is None:
            baseline = s
        assert s == baseline


def test_a_conserved_probe_with_zero_mismatches_is_not_the_same_case():
    """Multi-strain, but nothing to penalise: the term applies and is 1.0.

    The score is renormalised by the extra weight relative to a single-sequence
    run. That is correct and harmless — the term applies uniformly to every
    candidate in such a run, so it cannot re-order them — and it is what keeps
    "no mismatch data" distinguishable from "mismatch data, no mismatches".
    """
    c = _candidate([])
    s = compute_composite_score(c, 1000)
    assert c.score_components["mismatch_placement"] == 1.0
    # The seventh term enters the sum at its full weight (sub-score 1.0) and
    # the denominator grows by the same weight: an order-preserving affine map
    # of the old score, not a re-ranking.
    w_mm = DEFAULT_WEIGHTS.mismatch_placement
    assert s == pytest.approx(
        (_legacy_score(c, 1000, DEFAULT_WEIGHTS) + w_mm)
        / DEFAULT_WEIGHTS.total()
    )


# ---------------------------------------------------------------------------
# (iv) The [0, 1] invariant survives the seventh term.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("weights", [
    DEFAULT_WEIGHTS,
    RankingWeights(mismatch_placement=0.0),
    RankingWeights(mismatch_placement=1.0),
    RankingWeights(mismatch_placement=10.0),
    RankingWeights(gc_proximity=0.0, gc_balance=0.0, specificity_margin=0.0,
                   homopolymer_penalty=0.0, central_placement=0.0,
                   mismatch_placement=1.0),
    RankingWeights(gc_proximity=0.0, gc_balance=0.0, specificity_margin=0.0,
                   homopolymer_penalty=0.0, central_placement=0.0,
                   mismatch_placement=0.0),
])
def test_composite_stays_in_unit_range_with_the_new_weight(weights):
    positions_sets = [
        [], [0], [HALF], [HALF - 1, HALF], list(range(LENGTH)),
        [LENGTH - 1], [3, 24, 25, 46],
    ]
    for positions in positions_sets:
        for gc, gc_l, gc_r in ((0.0, 0.0, 1.0), (0.52, 0.5, 0.5),
                               (1.0, 1.0, 0.0)):
            for hp, margin, iso in ((0, -10.0, -1.0), (40, 1e6, 5.0)):
                for start, end, seq_len in ((0, 50, 100), (25, 75, 100),
                                            (10, 10, 100), (0, 50, 0)):
                    for jb in ("A", "G", "T", "C", "N", ""):
                        seq = (("ACGTA" * 5 + jb + "ACGTA" * 5)[:LENGTH]
                               if jb else "")
                        c = ProbeCandidate(
                            probe_id="p", sequence=seq, gc_content=gc,
                            gc_left=gc_l, gc_right=gc_r, homopolymer_max=hp,
                            specificity_margin=margin, isoform_coverage=iso,
                            start=start, end=end,
                        )
                        c.strain_binding = [{
                            "strain_id": "sB",
                            "total_mismatches": len(positions),
                            "mismatch_positions": list(positions),
                        }]
                        s = compute_composite_score(c, seq_len, weights)
                        assert 0.0 <= s <= 1.0, (s, weights, positions, jb)
                        sub = c.score_components["mismatch_placement"]
                        assert 0.0 <= sub <= 1.0, sub


def test_degenerate_half_lengths_do_not_blow_up():
    """A window too short to have an interior still returns a finite weight."""
    for half in (0, 1):
        assert mismatch_position_weight(0, half) == MISMATCH_WEIGHT_AT_NICK
        assert mismatch_distance_from_nick(0, half) == 0
    c = ProbeCandidate(probe_id="p", sequence="", gc_content=0.5,
                       gc_left=0.5, gc_right=0.5)
    c.strain_binding = [{"strain_id": "s", "total_mismatches": 1,
                         "mismatch_positions": [0]}]
    assert 0.0 <= compute_composite_score(c, 100) <= 1.0


# ---------------------------------------------------------------------------
# The weight is user-overridable.
# ---------------------------------------------------------------------------

def test_mismatch_placement_is_an_accepted_ranking_weight_key(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">G\n" + "ACGT" * 50 + "\n")
    cfg = RunConfig(fasta_path=str(fa), output_dir=str(tmp_path / "out"),
                    ranking_weights={"mismatch_placement": 0.05})
    cfg.validate()                              # must not raise
    with pytest.raises(Exception) as exc:
        RunConfig(fasta_path=str(fa), output_dir=str(tmp_path / "out"),
                  ranking_weights={"mismatch_plaecment": 0.05}).validate()
    assert "mismatch_placement" in str(exc.value)
