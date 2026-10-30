"""Unit tests for the composite ranking engine."""
from flexprobe.models import ProbeCandidate
from flexprobe.ranking import (
    DEFAULT_WEIGHTS, RankingWeights, _central_placement_score,
    _gc_balance_score, _gc_proximity_per_half_score, _gc_proximity_score,
    _homopolymer_penalty_score, compute_composite_score, rank_candidates,
)


def test_primer3_scoring_populates_thermo_fields():
    # Audit #21/#23: with primer3_scoring on, per-half Tm balance and the self-
    # and hetero-dimer ΔG are computed and recorded on the candidate.
    import importlib.util
    import pytest
    if importlib.util.find_spec("primer3") is None:
        pytest.skip("primer3 not installed")
    c = ProbeCandidate(probe_id="p", sequence="ACGTACGTACGTACGTACGTACGTA" * 2,
                       gc_content=0.5, gc_left=0.5, gc_right=0.5, start=25, end=75)
    rank_candidates([c], seq_length=100, primer3_scoring=True)
    assert not c.reduced_confidence          # clean ACGT sequence
    assert c.tm_balance >= 0.0
    # ΔG fields are populated (0.0 means "no dimer", still a real computed value)
    assert isinstance(c.self_dimer_dg, float)
    assert isinstance(c.hetero_dimer_dg, float)


def test_weights_sum_to_one():
    # The six terms that apply to EVERY run still sum to exactly 1.0, so a
    # single-sequence composite score is the plain weighted sum with no
    # renormalisation — unchanged by the addition of mismatch_placement.
    # total() is 1.10 because it also carries that seventh weight, which only a
    # run with on-target mismatch data can use; compute_composite_score leaves
    # it out of the denominator again on any candidate without such data.
    always_on = DEFAULT_WEIGHTS.total() - DEFAULT_WEIGHTS.mismatch_placement
    assert abs(always_on - 1.0) < 1e-9
    assert abs(DEFAULT_WEIGHTS.mismatch_placement - 0.10) < 1e-9


def test_gc_proximity_peaks_at_measured_tenx_ideal():
    # The ideal GC is the EMPIRICAL centre of 10x's shipped probe sets (median
    # per-half GC = 0.52 over 879,312 halves in the eight Chromium probe-set
    # CSVs), not the midpoint of the 44-72% acceptance band (0.58). 0.58 sat
    # above the 75th percentile of 10x's own probes.
    from flexprobe.ranking import _GC_IDEAL
    assert abs(_GC_IDEAL - 0.52) < 1e-9
    assert _gc_proximity_score(_GC_IDEAL) == 1.0
    assert _gc_proximity_score(0.58) < _gc_proximity_score(_GC_IDEAL)
    assert _gc_proximity_score(0.25) < _gc_proximity_score(0.55)
    assert _gc_proximity_score(0.0) == 0.0
    # Every GC in [0, 1] stays inside the unit range.
    assert all(0.0 <= _gc_proximity_score(g / 100.0) <= 1.0
               for g in range(0, 101))


def test_gc_balance():
    assert _gc_balance_score(0.5, 0.5) == 1.0
    assert _gc_balance_score(0.4, 0.6) < 1.0


# ---------------------------------------------------------------------------
# GC proximity is scored PER HALF, not on whole-window GC (audit finding M7).
#
# ``candidate.gc_content`` is the mean of the two halves (with equal 25/25
# halves it is *exactly* that mean — see ``filters.py`` ~line 483), so scoring
# it directly let a GC-rich half and a GC-poor half cancel out: a window whose
# halves sit at opposite edges of the acceptance band could still average to
# the empirical ideal and score a perfect 1.0. These tests pin the fixed
# behaviour and, where useful, contrast it explicitly with what the retired
# whole-window computation would have produced.
# ---------------------------------------------------------------------------

def test_gc_proximity_per_half_matches_whole_window_when_balanced():
    # Equal halves: gc_content == gc_left == gc_right, so per-half scoring must
    # reproduce the whole-window score exactly. This is the case the documented
    # rule and the filter both agree on, and it must not change.
    for half in (0.30, 0.44, 0.52, 0.58, 0.72, 0.90):
        assert _gc_proximity_per_half_score(half, half) == _gc_proximity_score(half)


def test_gc_proximity_per_half_uses_the_worse_half_not_the_average():
    import pytest

    # The exact pair demonstrated in the audit: halves at the GC-rich edge
    # (0.72) and below the band (0.32) average to whole-window GC == 0.52,
    # the empirical ideal, which the old whole-window computation scored 1.0.
    gc_left, gc_right = 0.72, 0.32
    whole_window_gc = (gc_left + gc_right) / 2.0
    assert whole_window_gc == pytest.approx(0.52)
    old_behaviour = _gc_proximity_score(whole_window_gc)
    assert old_behaviour == 1.0                      # what the bug scored

    new_behaviour = _gc_proximity_per_half_score(gc_left, gc_right)
    assert new_behaviour == pytest.approx(
        min(_gc_proximity_score(gc_left), _gc_proximity_score(gc_right)))
    assert new_behaviour == 0.0                       # both halves at/past the edge
    assert new_behaviour < old_behaviour


def test_gc_proximity_per_half_is_the_min_not_the_mean_of_the_two_halves():
    # One half exactly ideal, the other clearly not. The mean would still
    # reward the good half; the worst-of-two rule (matching the filter's
    # per-half AND gate and the mismatch_placement "worst strain decides"
    # precedent) must score exactly as if only the bad half existed.
    ideal, bad = 0.52, 0.30
    per_half = _gc_proximity_per_half_score(ideal, bad)
    assert per_half == _gc_proximity_score(bad)
    assert per_half < _gc_proximity_score(ideal)
    mean_gc = (ideal + bad) / 2.0
    assert per_half < _gc_proximity_score(mean_gc)     # strictly worse than averaging
    # Symmetric in the two halves.
    assert _gc_proximity_per_half_score(bad, ideal) == per_half


def test_gc_proximity_per_half_window_near_ideal_but_one_half_far_scores_low():
    # Whole-window GC is a perfect 0.52 (equal to _GC_IDEAL); one half alone
    # (0.80) is nowhere near it. A window this lopsided must score far below a
    # window whose halves are both genuinely near the ideal, even though the
    # two windows share the same whole-window GC.
    lopsided = ProbeCandidate(
        probe_id="lopsided", sequence="G" * 25 + "A" * 25,
        gc_content=0.52, gc_left=0.80, gc_right=0.24,
        homopolymer_max=0, specificity_margin=50.0, start=100, end=150,
    )
    balanced = ProbeCandidate(
        probe_id="balanced", sequence="ACGTA" * 10,
        gc_content=0.52, gc_left=0.52, gc_right=0.52,
        homopolymer_max=0, specificity_margin=50.0, start=100, end=150,
    )
    compute_composite_score(lopsided, seq_length=1000)
    compute_composite_score(balanced, seq_length=1000)
    assert lopsided.score_components["gc_proximity"] == 0.0
    assert balanced.score_components["gc_proximity"] == 1.0
    assert lopsided.composite_score < balanced.composite_score


def test_gc_proximity_per_half_reorders_ranking_vs_whole_window_average():
    # Same scenario driven through rank_candidates end-to-end (not just the
    # sub-score), because the M7 bug was specifically about ranking/selection
    # order: the lopsided window must not tie with, let alone beat, a window
    # that is actually GC-ideal on both halves.
    lopsided = ProbeCandidate(
        probe_id="lopsided", gc_content=0.52, gc_left=0.72, gc_right=0.32,
        homopolymer_max=0, specificity_margin=50.0, start=100, end=150,
    )
    ideal = ProbeCandidate(
        probe_id="ideal", gc_content=0.52, gc_left=0.52, gc_right=0.52,
        homopolymer_max=0, specificity_margin=50.0, start=100, end=150,
    )
    ranked = rank_candidates([lopsided, ideal], seq_length=1000)
    assert [c.probe_id for c in ranked] == ["ideal", "lopsided"]


def test_prescreen_score_also_scores_gc_proximity_per_half():
    # prescreen_score computes gc_proximity independently of compute_composite_
    # score (it runs before BLAST); it must not have kept the whole-window bug.
    lopsided = ProbeCandidate(
        gene_id="g", probe_id="lopsided", sequence="A" * 50,
        gc_content=0.52, gc_left=0.80, gc_right=0.24, homopolymer_max=0,
        start=0, end=50,
    )
    ideal = ProbeCandidate(
        gene_id="g", probe_id="ideal", sequence="A" * 50,
        gc_content=0.52, gc_left=0.52, gc_right=0.52, homopolymer_max=0,
        start=0, end=50,
    )
    assert prescreen_score(lopsided, seq_length=1000) < prescreen_score(
        ideal, seq_length=1000)


def test_homopolymer_penalty():
    # Graded INSIDE the passing range. The chemistry filter rejects runs above
    # max_homopolymer_run (default 4), so a term that only grades above 4 is
    # constant on every candidate the ranker ever sees, while still holding
    # its share of the default weight.
    assert _homopolymer_penalty_score(2) == 1.0
    assert _homopolymer_penalty_score(3) == 0.75
    assert _homopolymer_penalty_score(4) == 0.5
    assert _homopolymer_penalty_score(5) < _homopolymer_penalty_score(4)
    assert _homopolymer_penalty_score(100) == 0.0
    # Monotone non-increasing across the whole domain, including the old
    # threshold, so relaxing the filter cannot make a longer run score better.
    runs = list(range(1, 12))
    scores = [_homopolymer_penalty_score(r) for r in runs]
    assert all(a >= b for a, b in zip(scores, scores[1:])), scores


def test_homopolymer_penalty_is_not_constant_within_the_passing_range():
    """A scoring term must vary over the inputs the filter lets through.

    A scoring term that returns one value
    for every candidate that survives filtering contributes no ordering and
    silently dilutes the other terms by its weight.
    """
    from flexprobe.models import FilterThresholds

    ceiling = FilterThresholds().max_homopolymer_run
    reachable = {_homopolymer_penalty_score(r) for r in range(1, ceiling + 1)}
    assert len(reachable) > 1, (
        f"score is constant at {reachable} for every run length the filter "
        f"admits (1..{ceiling}); it can never change a rank")


def test_central_placement_edges():
    assert _central_placement_score(0, 50, 0) == 0.0     # unknown gene length
    assert _central_placement_score(10, 10, 100) == 0.0  # zero-width
    # Terminal 15% is strongly discouraged; everything inboard of it is not.
    L = 10_000
    edge = _central_placement_score(0, 50, L)            # midpoint at 0.25%
    inner = _central_placement_score(1000, 1050, L)      # midpoint at 10.25%
    assert edge < inner < 1.0
    assert edge < 0.01                                   # the very end ~ zero


def test_central_placement_is_flat_across_the_middle():
    """No preference between one interior position and another.

    This is the whole point of the plateau: a peak at the exact midpoint made
    every top-ranked window on a long transcript come from a narrow central
    band, so the chosen probes clustered there. Positions from 15% to 85% must
    score identically for spacing — not ranking — to decide the spread.
    """
    L = 10_000
    scores = {
        pos: _central_placement_score(pos, pos + 50, L)
        for pos in (1500, 2500, 5000, 7000, 8400)
    }
    assert set(scores.values()) == {1.0}, scores
    # ...and the plateau edges are exactly END_EXCLUSION_FRAC in from each end.
    assert _central_placement_score(1475, 1525, L) == 1.0        # mid = 0.150
    assert _central_placement_score(1375, 1425, L) < 1.0         # mid = 0.140
    assert _central_placement_score(8475, 8525, L) == 1.0        # mid = 0.850
    assert _central_placement_score(8575, 8625, L) < 1.0         # mid = 0.860


def test_central_placement_is_symmetric():
    import pytest
    L = 4_000
    for off in (10, 100, 300, 599):
        assert _central_placement_score(off, off + 50, L) == pytest.approx(
            _central_placement_score(L - off - 50, L - off, L))


def test_composite_score_in_unit_range():
    c = ProbeCandidate(probe_id="p", gc_content=0.5, gc_left=0.5, gc_right=0.5,
                       homopolymer_max=3, specificity_margin=50.0,
                       start=25, end=75)
    s = compute_composite_score(c, seq_length=100)
    assert 0.0 <= s <= 1.0
    assert c.composite_score == s


# ---------------------------------------------------------------------------
# composite_score ∈ [0, 1] — the contract the primer3 scaling bug broke.
#
# The shipped viral panel (n=177) ran with composite scores from -1.58 to +0.25,
# 89% of them negative, because `sd_penalty = max(0, dg / -10000)` treated
# primer3's cal/mol ΔG as if it were kcal/mol: a routine -13,560 cal/mol
# self-dimer became a 1.356 penalty against a score capped at 1.0. Nothing
# asserted the range, so nothing caught it. These tests do.
# ---------------------------------------------------------------------------

def test_dg_penalty_is_kcal_scaled_and_capped():
    from flexprobe.ranking import (
        DG_FLOOR_KCAL, DG_PENALTY_MAX, DG_SCALE_KCAL, _dg_penalty,
    )
    # primer3 returns cal/mol. Weak/absent structure -> no penalty.
    assert _dg_penalty(0.0) == 0.0
    assert _dg_penalty(500.0) == 0.0                     # positive ΔG
    assert _dg_penalty(DG_FLOOR_KCAL * 1000) == 0.0      # exactly at the floor
    # Saturates at the documented per-term maximum, and never exceeds it.
    assert _dg_penalty(DG_SCALE_KCAL * 1000) == DG_PENALTY_MAX
    # A real extreme like this one, unclamped, would give a penalty of 1.356,
    # larger than the entire score being nudged.
    assert _dg_penalty(-13559.9) == DG_PENALTY_MAX
    assert _dg_penalty(-1e9) == DG_PENALTY_MAX
    # Monotone and bounded across the whole plausible range.
    prev = 0.0
    for dg_cal in range(0, -30001, -100):
        p = _dg_penalty(float(dg_cal))
        assert 0.0 <= p <= DG_PENALTY_MAX
        assert p >= prev
        prev = p
    assert _dg_penalty(float("nan")) == 0.0


def test_composite_score_unit_range_over_input_grid():
    """0 <= composite_score <= 1 for every combination of feature extremes."""
    weight_sets = [
        DEFAULT_WEIGHTS,
        RankingWeights(gc_proximity=0.0, gc_balance=0.0, specificity_margin=1.0,
                       homopolymer_penalty=0.0, central_placement=0.0),
        RankingWeights(gc_proximity=3.0, gc_balance=2.0, specificity_margin=5.0,
                       homopolymer_penalty=1.0, central_placement=1.0,
                       isoform_coverage=4.0),
    ]
    for w in weight_sets:
        for gc in (0.0, 0.2, 0.44, 0.52, 0.58, 0.72, 1.0):
            for gc_l, gc_r in ((0.0, 1.0), (0.5, 0.5), (1.0, 0.0)):
                for hp in (0, 4, 5, 40):
                    for margin in (-10.0, 0.0, 50.0, 1e6):
                        for iso in (-1.0, 0.0, 0.5, 1.0, 5.0):
                            for start, end, seq_len in (
                                (0, 50, 100), (25, 75, 100), (900, 950, 1000),
                                (0, 50, 0), (10, 10, 100),
                            ):
                                for jb in ("A", "G", "T", "C", "N", ""):
                                    # 50-mer whose index-25 base is `jb`
                                    seq = (("ACGTA" * 5 + jb + "ACGTA" * 5)[:50]
                                           if jb else "")
                                    c = ProbeCandidate(
                                        probe_id="p", sequence=seq,
                                        gc_content=gc, gc_left=gc_l,
                                        gc_right=gc_r, homopolymer_max=hp,
                                        specificity_margin=margin,
                                        isoform_coverage=iso, start=start,
                                        end=end,
                                    )
                                    s = compute_composite_score(c, seq_len, w)
                                    assert 0.0 <= s <= 1.0, (s, w, gc, jb)


# ---------------------------------------------------------------------------
# TN ligation junction (mRNA base at index 25).
#
# Default = hard filter, so every survivor has 'A' and the gate must be exactly
# inert. Opt-in soft mode relaxes the floor to "not C", and then the ranker has
# to put A-junction probes decisively above G/T ones.
# ---------------------------------------------------------------------------

def _seq_with_junction(base, filler="ACGTA"):
    """50-mer whose mRNA index-25 base is *base*."""
    seq = (filler * 5 + base + filler * 5)[:50]
    assert seq[len(seq) // 2] == base
    return seq


def test_junction_multiplier_ranks_A_above_G_above_T_above_C():
    from flexprobe.ranking import (
        JUNCTION_BASE_MULTIPLIERS, _junction_base_multiplier,
    )
    m = {b: _junction_base_multiplier(_seq_with_junction(b)) for b in "ACGT"}
    assert m["A"] == 1.0
    assert m["A"] > m["G"] > m["T"] > m["C"] == 0.0
    # G is preferred over T because 10x's own shipped panels use G at the
    # junction ~3x as often as T (35.72% vs 12.05% of 439,656 probes), not
    # because G and T are chemically interchangeable.
    assert 2.0 < m["G"] / m["T"] < 3.5
    assert all(0.0 <= v <= 1.0 for v in JUNCTION_BASE_MULTIPLIERS.values())
    # Unverifiable base -> no credit; absent sequence -> inert, not zeroed.
    assert _junction_base_multiplier(_seq_with_junction("N")) == 0.0
    assert _junction_base_multiplier("") == 1.0


def test_junction_gate_is_exactly_inert_in_hard_mode():
    """Hard mode: every survivor has A, so the gate must change nothing."""
    common = dict(gc_content=0.52, gc_left=0.40, gc_right=0.60,
                  homopolymer_max=5, specificity_margin=30.0, start=10, end=60)
    probes = [
        ProbeCandidate(probe_id="a", sequence=_seq_with_junction("A"), **common),
        ProbeCandidate(probe_id="b", sequence=_seq_with_junction("A", "GGCTA"),
                       gc_content=0.60, gc_left=0.60, gc_right=0.60,
                       homopolymer_max=2, specificity_margin=90.0,
                       start=25, end=75),
        ProbeCandidate(probe_id="c", sequence=_seq_with_junction("A", "TTATA"),
                       gc_content=0.44, gc_left=0.20, gc_right=0.70,
                       homopolymer_max=9, specificity_margin=1.0,
                       start=0, end=50),
    ]
    for p in probes:
        compute_composite_score(p, seq_length=100)
        assert p.score_components["junction_base"] == 1.0
    # Identical candidates without any sequence must score bit-for-bit the same.
    for p in probes:
        bare = ProbeCandidate(
            probe_id=p.probe_id, gc_content=p.gc_content, gc_left=p.gc_left,
            gc_right=p.gc_right, homopolymer_max=p.homopolymer_max,
            specificity_margin=p.specificity_margin, start=p.start, end=p.end,
        )
        assert compute_composite_score(bare, seq_length=100) == p.composite_score
    ranked = rank_candidates(probes, seq_length=100)
    assert [c.probe_id for c in ranked] == ["b", "a", "c"]


def test_soft_mode_worst_A_junction_probe_beats_best_G_or_T_probe():
    """The gate is a floor-vs-ceiling separation, not a tie-breaker."""
    def score(base, good):
        feats = (dict(gc_content=0.52, gc_left=0.52, gc_right=0.52,
                      homopolymer_max=0, specificity_margin=1e6,
                      start=25, end=75)
                 if good else
                 dict(gc_content=0.72, gc_left=0.0, gc_right=1.0,
                      homopolymer_max=40, specificity_margin=0.0,
                      start=0, end=50))
        c = ProbeCandidate(probe_id=base, sequence=_seq_with_junction(base),
                           **feats)
        compute_composite_score(c, seq_length=100)
        return c.composite_score

    best_a, worst_a = score("A", True), score("A", False)
    assert best_a == 1.0
    for other in ("G", "T", "C"):
        assert score(other, True) < best_a
        # ...and the *best possible* non-A probe stays below a realistic
        # A-junction probe. The re-scored viral panel spans 0.468-0.956, so
        # anything at or above that floor is safe.
        assert score(other, True) < 0.468
    # Only a genuinely worthless A probe (every sourced axis at its floor) can
    # be caught, which is the intended and acceptable failure mode.
    assert worst_a < score("G", True)


def test_soft_mode_ranking_puts_every_A_probe_above_every_non_A_probe():
    cands = []
    for i, (base, gc, hp, margin) in enumerate([
        ("A", 0.44, 9, 1.0), ("A", 0.52, 0, 90.0), ("A", 0.70, 6, 20.0),
        ("G", 0.52, 0, 100.0), ("G", 0.52, 0, 95.0),
        ("T", 0.52, 0, 100.0), ("C", 0.52, 0, 100.0),
    ]):
        cands.append(ProbeCandidate(
            probe_id=f"{base}{i}", sequence=_seq_with_junction(base),
            gc_content=gc, gc_left=gc, gc_right=gc, homopolymer_max=hp,
            specificity_margin=margin, start=25, end=75))
    ranked = rank_candidates(cands, seq_length=100)
    bases = [c.probe_id[0] for c in ranked]
    assert bases[:3] == ["A", "A", "A"], bases
    assert bases[3:5] == ["G", "G"], bases
    assert bases[5:] == ["T", "C"], bases
    assert all(0.0 <= c.composite_score <= 1.0 for c in ranked)


def test_soft_mode_prescreen_cap_prefers_A_junction_windows():
    """The pre-screen per-gene cap must not fill up with G/T windows."""
    from flexprobe.ranking import prescreen_score, select_top_per_gene
    cands = []
    for i, base in enumerate("GGGGAATT"):
        cands.append(ProbeCandidate(
            gene_id="g", probe_id=f"p{i}_{base}",
            sequence=_seq_with_junction(base), start=i * 100, end=i * 100 + 50,
            gc_content=0.52, gc_left=0.52, gc_right=0.52, homopolymer_max=0))
    assert all(0.0 <= prescreen_score(c, 1000) <= 1.0 for c in cands)
    kept = select_top_per_gene(cands, gene_key=lambda c: c.gene_id,
                               seq_length_of=lambda _g: 1000, n=2)
    assert sorted(c.probe_id[-1] for c in kept) == ["A", "A"]


def test_junction_gate_survives_the_primer3_adjustment(monkeypatch):
    """Gate + bounded nudge together still land inside [0, 1]."""
    for base in "ACGT":
        for dg in (-1e6, 0.0, 1e6):
            _with_fake_primer3(monkeypatch, _FakePrimer3(dg, 60.0, 0.0, 200.0))
            c = ProbeCandidate(probe_id="p", sequence=_seq_with_junction(base),
                               gc_content=0.52, gc_left=0.52, gc_right=0.52,
                               homopolymer_max=0, specificity_margin=1e6,
                               start=25, end=75)
            rank_candidates([c], seq_length=100, primer3_scoring=True)
            assert 0.0 <= c.composite_score <= 1.0, (base, dg)
            # A C-junction probe is gated to zero and the nudge cannot revive it.
            if base == "C":
                assert c.composite_score == 0.0


class _FakeThermo:
    def __init__(self, dg):
        self.dg = dg


class _FakePrimer3:
    """primer3 stand-in that returns caller-chosen ΔG / Tm values."""

    def __init__(self, dg, tm=60.0, tm_lhs=60.0, tm_rhs=60.0):
        self._dg, self._tm = dg, tm
        self._halves = (tm_lhs, tm_rhs)
        self._n = 0

    def calc_tm(self, seq, **_kw):
        if len(seq) > 30:
            return self._tm
        self._n += 1
        return self._halves[self._n % 2]

    def calc_hairpin(self, seq, **_kw):
        return _FakeThermo(self._dg)

    def calc_homodimer(self, seq, **_kw):
        return _FakeThermo(self._dg)

    def calc_heterodimer(self, a, b, **_kw):
        return _FakeThermo(self._dg)


def _with_fake_primer3(monkeypatch, fake):
    import importlib.machinery
    import sys
    import types
    mod = types.ModuleType("primer3")
    mod.__spec__ = importlib.machinery.ModuleSpec("primer3", loader=None)
    for name in ("calc_tm", "calc_hairpin", "calc_homodimer", "calc_heterodimer"):
        setattr(mod, name, getattr(fake, name))
    monkeypatch.setitem(sys.modules, "primer3", mod)


def test_primer3_adjustment_cannot_push_score_out_of_unit_range(monkeypatch):
    """Extreme ΔG must not be able to drive composite_score out of [0, 1].

    Driven through a primer3 stand-in so the invariant is asserted against
    absurd inputs (-1e6 cal/mol, +1e6 cal/mol, wild Tm) rather than only the
    values primer3 happens to return today.
    """
    from flexprobe.ranking import (
        PRIMER3_ADJUSTMENT_MAX, PRIMER3_ADJUSTMENT_MIN,
    )
    seq = "ACGTACGTACGTACGTACGTACGTA" * 2
    cases = [
        (-1e6, 20.0, 0.0, 200.0),     # absurd structure + huge Tm imbalance
        (-13559.9, 60.0, 60.0, 60.0),  # the panel's real worst self-dimer
        (-1000.0, 60.0, 60.0, 60.0),
        (0.0, 60.0, 60.0, 60.0),
        (1e6, 60.0, 60.0, 60.0),       # positive ΔG (no structure)
    ]
    # Include feature-poor candidates, whose base score is near 0 and so are
    # the ones a large penalty would push negative.
    profiles = [
        dict(gc_content=0.52, gc_left=0.52, gc_right=0.52, homopolymer_max=0,
             specificity_margin=1e6),
        dict(gc_content=1.0, gc_left=0.0, gc_right=1.0, homopolymer_max=50,
             specificity_margin=0.0),
    ]
    for dg, tm, tm_l, tm_r in cases:
        for prof in profiles:
            _with_fake_primer3(monkeypatch, _FakePrimer3(dg, tm, tm_l, tm_r))
            c = ProbeCandidate(probe_id="p", sequence=seq, start=25, end=75,
                               **prof)
            rank_candidates([c], seq_length=100, primer3_scoring=True)
            assert 0.0 <= c.composite_score <= 1.0, (dg, prof, c.composite_score)
            adj = c.score_components["primer3_adjustment"]
            assert PRIMER3_ADJUSTMENT_MIN <= adj <= PRIMER3_ADJUSTMENT_MAX


def test_primer3_adjustment_stays_in_unit_range_on_real_sequences():
    import importlib.util
    import pytest
    if importlib.util.find_spec("primer3") is None:
        pytest.skip("primer3 not installed")
    seqs = [
        "ACGTACGTACGTACGTACGTACGTA" * 2,
        "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC",   # GC-rich
        "ATATATATATATATATATATATATATATATATATATATATATATATAT",   # AT-rich
        "GGGGGGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCC",  # self-complementary
        "CGGACTACATCTACATTTGATCCGGAGGCAGCATATGTGGCATTTGAAGC",  # real panel probe
        "TTCCAGTCGAATCCAGTGCCGGACAACGGTCTTACCCTGCACCGTCTGTC",  # real panel probe
    ]
    for i, s in enumerate(seqs):
        for gc, gl, gr, hp, margin in ((0.52, 0.52, 0.52, 0, 100.0),
                                       (1.0, 0.0, 1.0, 50, 0.0)):
            c = ProbeCandidate(probe_id=f"p{i}", sequence=s, gc_content=gc,
                               gc_left=gl, gc_right=gr, homopolymer_max=hp,
                               specificity_margin=margin, start=25, end=75)
            rank_candidates([c], seq_length=100, primer3_scoring=True)
            assert 0.0 <= c.composite_score <= 1.0, (s, c.composite_score)


def test_primer3_terms_use_the_synthesised_oligos_not_the_mrna_window():
    """The 50-mer mRNA window is never synthesised; the 25-mer halves are."""
    import importlib.util
    import pytest
    if importlib.util.find_spec("primer3") is None:
        pytest.skip("primer3 not installed")
    import primer3 as p3

    from flexprobe.chemistry import CHEMISTRIES
    from flexprobe.constants import (
        TM_DNA_CONC, TM_DNTP_CONC, TM_DV_CONC, TM_MV_CONC,
    )
    from flexprobe.ranking import _target_oligo_halves

    seq = "CGGACTACATCTACATTTGATCCGGAGGCAGCATATGTGGCATTTGAAGC"
    # The halves ranking uses must be exactly what the chemistry module builds.
    _m5, _m3, target_lhs, target_rhs = CHEMISTRIES["flex_v1"].split_tile(seq)
    assert _target_oligo_halves(seq) == (target_lhs, target_rhs)

    salt = dict(mv_conc=TM_MV_CONC, dv_conc=TM_DV_CONC, dntp_conc=TM_DNTP_CONC,
                dna_conc=TM_DNA_CONC)
    c = ProbeCandidate(probe_id="p", sequence=seq, gc_content=0.48,
                       gc_left=0.48, gc_right=0.48, start=25, end=75)
    rank_candidates([c], seq_length=100, primer3_scoring=True)
    expected_self = min(p3.calc_homodimer(target_lhs, **salt).dg,
                        p3.calc_homodimer(target_rhs, **salt).dg)
    assert abs(c.self_dimer_dg - expected_self) < 1e-6
    assert abs(c.hetero_dimer_dg
               - p3.calc_heterodimer(target_lhs, target_rhs, **salt).dg) < 1e-6
    # ...and NOT the mRNA-window value the pre-fix code recorded.
    assert abs(c.self_dimer_dg - p3.calc_homodimer(seq, **salt).dg) > 1.0


def test_composite_score_exposes_component_subscores():
    # Audit #6/#7: the per-component 0-1 sub-scores are recorded so a rank is
    # auditable, and the weighted+normalised sum reproduces composite_score.
    from flexprobe.ranking import DEFAULT_WEIGHTS
    c = ProbeCandidate(probe_id="p", gc_content=0.5, gc_left=0.5, gc_right=0.6,
                       homopolymer_max=5, specificity_margin=40.0,
                       start=25, end=75)
    s = compute_composite_score(c, seq_length=100, weights=DEFAULT_WEIGHTS)
    comps = c.score_components
    assert set(comps) == {
        "gc_proximity", "gc_balance", "specificity_margin",
        "homopolymer_penalty", "central_placement", "isoform_coverage",
        "mismatch_placement", "junction_base",
    }
    assert all(0.0 <= v <= 1.0 for v in comps.values())
    w = DEFAULT_WEIGHTS
    # This candidate carries no strain_binding, so mismatch_placement is
    # reported (1.0, "nothing to penalise") but weighted at 0 and left out of
    # the denominator — which is why the recombination below is the same
    # six-term expression it always was. See test_mismatch_placement.py.
    recombined = comps["junction_base"] * (
        w.gc_proximity * comps["gc_proximity"]
        + w.gc_balance * comps["gc_balance"]
        + w.specificity_margin * comps["specificity_margin"]
        + w.homopolymer_penalty * comps["homopolymer_penalty"]
        + w.central_placement * comps["central_placement"]
        + w.isoform_coverage * comps["isoform_coverage"]
    ) / (w.total() - w.mismatch_placement)
    assert abs(recombined - s) < 1e-9


def test_rank_candidates_deterministic_and_assigns_rank():
    a = ProbeCandidate(probe_id="a", gc_content=0.5, gc_left=0.5, gc_right=0.5,
                       homopolymer_max=2, start=40, end=90)
    b = ProbeCandidate(probe_id="b", gc_content=0.2, gc_left=0.2, gc_right=0.9,
                       homopolymer_max=8, start=0, end=50)
    out = rank_candidates([b, a], seq_length=100)
    assert [c.probe_id for c in out] == ["a", "b"]   # higher score first
    assert out[0].rank == 1 and out[1].rank == 2


def test_rank_tie_break_by_probe_id():
    # identical features -> equal score -> tie broken by probe_id ascending
    common = dict(gc_content=0.5, gc_left=0.5, gc_right=0.5, homopolymer_max=1,
                  start=10, end=60)
    z = ProbeCandidate(probe_id="z", **common)
    a = ProbeCandidate(probe_id="a", **common)
    out = rank_candidates([z, a], seq_length=100)
    assert [c.probe_id for c in out] == ["a", "z"]


def test_custom_weights_change_score():
    c = ProbeCandidate(probe_id="p", gc_content=0.5, gc_left=0.5, gc_right=0.5,
                       homopolymer_max=3, specificity_margin=0.0, start=25, end=75)
    s_default = compute_composite_score(c, 100, DEFAULT_WEIGHTS)
    s_spec = compute_composite_score(
        c, 100, RankingWeights(gc_proximity=0.0, gc_balance=0.0,
                               specificity_margin=1.0, homopolymer_penalty=0.0,
                               central_placement=0.0))
    assert s_default != s_spec


# ---- pre-screen per-gene cap (whole-transcriptome speedup) ----

from flexprobe.ranking import prescreen_score, select_top_per_gene  # noqa: E402


def _cand(pid, gene, start, gc=0.5, gc_l=0.5, gc_r=0.5, hp=0):
    return ProbeCandidate(
        gene_id=gene, probe_id=pid, sequence="A" * 50, start=start, end=start + 50,
        gc_content=gc, gc_left=gc_l, gc_right=gc_r, homopolymer_max=hp,
    )


def test_prescreen_score_in_unit_range_and_blast_free():
    c = _cand("p", "g", 100, gc=0.5, gc_l=0.5, gc_r=0.5, hp=0)
    s = prescreen_score(c, seq_length=200)
    assert 0.0 <= s <= 1.0
    # specificity_margin / isoform_coverage must NOT influence the pre-score
    c.specificity_margin = 99.0
    c.isoform_coverage = 0.1
    assert prescreen_score(c, seq_length=200) == s


def test_select_top_per_gene_caps_each_gene():
    # gene A: 5 candidates, gene B: 2 candidates; cap at 3
    cands = [_cand(f"A_{i}", "A", i * 50) for i in range(5)]
    cands += [_cand(f"B_{i}", "B", i * 50) for i in range(2)]
    kept = select_top_per_gene(
        cands, gene_key=lambda c: c.gene_id,
        seq_length_of=lambda _g: 1000, n=3,
    )
    from collections import Counter
    per_gene = Counter(c.gene_id for c in kept)
    assert per_gene["A"] == 3      # capped
    assert per_gene["B"] == 2      # <= n, kept whole


def test_select_top_per_gene_keeps_best_and_is_deterministic():
    # one strong candidate (GC right at 0.5, balanced) + weaker ones
    strong = _cand("g_best", "g", 500, gc=0.50, gc_l=0.50, gc_r=0.50, hp=0)
    weak = [_cand(f"g_{i}", "g", i * 10, gc=0.72, gc_l=0.72, gc_r=0.30, hp=6)
            for i in range(10)]
    cands = weak[:5] + [strong] + weak[5:]
    kept = select_top_per_gene(
        cands, gene_key=lambda c: c.gene_id,
        seq_length_of=lambda _g: 1000, n=1,
    )
    assert [c.probe_id for c in kept] == ["g_best"]
    # deterministic across calls
    kept2 = select_top_per_gene(
        cands, gene_key=lambda c: c.gene_id,
        seq_length_of=lambda _g: 1000, n=1,
    )
    assert [c.probe_id for c in kept] == [c.probe_id for c in kept2]


def test_select_top_per_gene_preserves_input_order():
    cands = [_cand(f"g_{i}", "g", i * 50) for i in range(6)]
    kept = select_top_per_gene(
        cands, gene_key=lambda c: c.gene_id,
        seq_length_of=lambda _g: 1000, n=10,   # n > count -> all kept, same order
    )
    assert [c.probe_id for c in kept] == [c.probe_id for c in cands]


def test_select_top_per_gene_n_zero_is_passthrough():
    cands = [_cand(f"g_{i}", "g", i) for i in range(3)]
    assert select_top_per_gene(
        cands, gene_key=lambda c: c.gene_id,
        seq_length_of=lambda _g: 100, n=0,
    ) == cands
