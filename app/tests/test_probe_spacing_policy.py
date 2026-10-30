# flexiprobes — probe design for 10x Genomics Flex v1 / Flex v2 / Visium.
# Copyright (C) 2026 Adrian Ilich and contributors
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""The 150 bp target / 40 bp fallback spacing policy, and the spread search.

Covers the three claims the policy makes:
  1. a roomy gene gets probes spread far wider than the 150 bp target;
  2. a cramped gene still gets its three probes, by dropping to 40 bp;
  3. whenever (2) happens the gene is FLAGGED, never silently emitted.
"""
from __future__ import annotations

import pytest

from flexprobe.auto_select import (
    SPACING_LADDER, auto_pick_per_gene, auto_pick_with_report,
    ladder_from, spread_gap_for,
)
from flexprobe.constants import (
    GAP_TIER_GREEN, GAP_TIER_ORANGE, GAP_TIER_RED, MIN_PROBE_GAP_BP,
    MIN_PROBE_GAP_FLOOR_BP, PROBE_GAP_RELAX_MIN_PROBES, RELAXED_PROBE_GAP_BP,
    gap_tier,
)
from flexprobe.isoform_select import coverage_pick, gene_index_from_manifest


def _row(pid, start, score, gene="ENST1"):
    return {"probe_id": pid, "gene_id": gene, "start": str(start),
            "end": str(start + 50), "composite_score": f"{score:.4f}"}


def _min_gap(rows_by_id, ids):
    ivs = sorted((int(rows_by_id[i]["start"]), int(rows_by_id[i]["end"]))
                 for i in ids)
    return min(ivs[k + 1][0] - ivs[k][1] for k in range(len(ivs) - 1))


# --------------------------------------------------------------------------- #
# Defaults
# --------------------------------------------------------------------------- #
def test_defaults_are_the_documented_policy():
    assert MIN_PROBE_GAP_BP == 150
    assert RELAXED_PROBE_GAP_BP == 40
    assert PROBE_GAP_RELAX_MIN_PROBES == 3
    assert RELAXED_PROBE_GAP_BP < MIN_PROBE_GAP_BP


# --------------------------------------------------------------------------- #
# 1. Spread
# --------------------------------------------------------------------------- #
def test_spread_gap_scales_with_the_gene():
    # Short transcript: even spread is tighter than the target, target wins.
    assert spread_gap_for(600, 3) >= MIN_PROBE_GAP_BP
    # Long transcript: three probes over the middle 70% of 12 kb is kilobases.
    wide = spread_gap_for(12_000, 3)
    assert wide > 3_000
    # More probes on the same transcript => each pair sits closer.
    assert spread_gap_for(12_000, 6) < wide
    # Degenerate inputs fall back to the floor rather than exploding.
    assert spread_gap_for(0, 3) >= MIN_PROBE_GAP_BP
    assert spread_gap_for(12_000, 1) >= MIN_PROBE_GAP_BP


def test_long_gene_probes_do_not_cluster_at_the_target():
    """A 12 kb transcript tiled everywhere: picks must spread, not bunch.

    Scores here *decrease* away from the centre, which is precisely the shape
    the old central-placement tent produced. Under a bare 150 bp floor the
    greedy walk would take the three highest-scoring windows — all within a few
    hundred bases of the midpoint. The spread search must beat that.
    """
    L = 12_000
    rows = []
    for i, start in enumerate(range(0, L - 50, 100)):
        mid = start + 25
        rows.append(_row(f"p{i}", start, 1.0 - abs(mid - L / 2) / L))
    by_id = {r["probe_id"]: r for r in rows}

    pick = auto_pick_with_report(rows, n_per_gene=3,
                                gene_length_of=lambda g: L)
    ids = pick.selections["ENST1"]
    assert len(ids) == 3
    gap = _min_gap(by_id, ids)
    assert gap > 2_000, f"probes still clustered: closest pair {gap} bp"
    assert not pick.spacing["ENST1"].tight
    assert pick.relaxed_genes == []


def test_spread_never_costs_a_probe():
    """Widening is an aspiration; the probe count is the constraint.

    Only four windows exist and they are packed into the first 800 bp of a long
    transcript, so the spread target is unreachable. Selection must still return
    three probes at the plain target rather than two well-spread ones.
    """
    rows = [_row(f"p{i}", i * 200, 1.0 - i * 0.01) for i in range(4)]
    pick = auto_pick_with_report(rows, n_per_gene=3,
                                gene_length_of=lambda g: 20_000)
    assert len(pick.selections["ENST1"]) == 3
    assert pick.spacing["ENST1"].gap_used == MIN_PROBE_GAP_BP


# --------------------------------------------------------------------------- #
# 2 + 3. Relaxation, and the flag that must accompany it
# --------------------------------------------------------------------------- #
def test_relaxes_to_forty_when_three_probes_will_not_fit_at_the_target():
    # Six 50-mers on a 60 bp grid: the widest three-probe layout has 70 bp gaps,
    # so the gene can only ever offer two probes at the 150 bp target.
    rows = [_row(f"p{i}", i * 60, 1.0 - i * 0.01) for i in range(6)]
    by_id = {r["probe_id"]: r for r in rows}

    # Disabling relaxation means removing BOTH lower rungs — raising only
    # relaxed_gap still leaves the 15 bp floor rung below it.
    strict = auto_pick_with_report(rows, n_per_gene=3,
                                   relaxed_gap=MIN_PROBE_GAP_BP,
                                   floor_gap=MIN_PROBE_GAP_BP)
    assert len(strict.selections["ENST1"]) == 2

    pick = auto_pick_with_report(rows, n_per_gene=3)
    ids = pick.selections["ENST1"]
    sp = pick.spacing["ENST1"]
    assert len(ids) == 3
    # It stopped on the orange rung, but the spread search inside that rung
    # still pushed the floor above 40 — gap_used is what was ACHIEVED, not the
    # rung's nominal value, so the report never understates the spacing.
    assert RELAXED_PROBE_GAP_BP <= sp.gap_used < MIN_PROBE_GAP_BP
    assert sp.tier == GAP_TIER_ORANGE
    assert _min_gap(by_id, ids) >= RELAXED_PROBE_GAP_BP


def test_relaxed_genes_are_flagged():
    rows = [_row(f"p{i}", i * 60, 1.0 - i * 0.01) for i in range(6)]
    pick = auto_pick_with_report(rows, n_per_gene=3)
    sp = pick.spacing["ENST1"]
    assert sp.relaxed is True
    assert sp.tight is True
    assert sp.min_gap is not None and sp.min_gap < MIN_PROBE_GAP_BP
    assert pick.relaxed_genes == ["ENST1"]
    assert pick.summary()["genes_below_min_gap"] == 1


def test_comfortable_gene_is_not_flagged():
    rows = [_row(f"p{i}", i * 1_000, 1.0 - i * 0.01) for i in range(4)]
    pick = auto_pick_with_report(rows, n_per_gene=3, gene_length_of=lambda g: 4_000)
    sp = pick.spacing["ENST1"]
    assert sp.relaxed is False and sp.tight is False
    assert pick.relaxed_genes == []
    assert pick.summary()["genes_below_min_gap"] == 0


def test_no_relaxation_when_the_target_already_yields_three():
    """A gene that clears the target keeps its spacing even if packing tighter
    would have fitted more probes. Relaxation is a rescue, not an optimiser."""
    rows = [_row(f"p{i}", i * 200, 1.0 - i * 0.01) for i in range(10)]
    pick = auto_pick_with_report(rows, n_per_gene=3,
                                gene_length_of=lambda g: 2_000)
    assert pick.spacing["ENST1"].gap_used >= MIN_PROBE_GAP_BP
    assert not pick.spacing["ENST1"].tight


def test_relax_below_bounds_the_rescue():
    """With relax_below=2 a gene that already produced two probes is not
    rescued, so the tighter third probe is never emitted."""
    rows = [_row(f"p{i}", i * 60, 1.0 - i * 0.01) for i in range(6)]
    pick = auto_pick_with_report(rows, n_per_gene=3, relax_below=2)
    assert len(pick.selections["ENST1"]) == 2
    assert pick.spacing["ENST1"].relaxed is False


def test_thin_wrapper_returns_just_the_selections():
    rows = [_row(f"p{i}", i * 1_000, 1.0 - i * 0.01) for i in range(4)]
    plain = auto_pick_per_gene(rows, n_per_gene=3)
    report = auto_pick_with_report(rows, n_per_gene=3)
    assert plain == report.selections


def test_gene_length_falls_back_to_the_furthest_probe():
    """Without a length map the spread target still works, just conservatively."""
    rows = [_row(f"p{i}", i * 500, 1.0 - i * 0.01) for i in range(8)]
    pick = auto_pick_with_report(rows, n_per_gene=3)
    assert pick.spacing["ENST1"].gene_length == 7 * 500 + 50
    assert len(pick.selections["ENST1"]) == 3


# --------------------------------------------------------------------------- #
# The coverage selector runs the same policy
# --------------------------------------------------------------------------- #
def _cov_row(pid, start, score, tid="ENST1"):
    return {"probe_id": pid, "gene_id": tid, "start": str(start),
            "end": str(start + 50), "composite_score": f"{score:.4f}",
            "specificity_pass": "TRUE", "isoform_hit_ids": tid,
            "mrna_target_50bp": "A" * 50}


def test_coverage_pick_relaxes_and_flags():
    rows = [_cov_row(f"p{i}", i * 60, 1.0 - i * 0.01) for i in range(6)]
    index = gene_index_from_manifest({"G1": ["ENST1"]})

    strict = coverage_pick(rows, index, n_per_gene=3,
                           relaxed_gap=MIN_PROBE_GAP_BP,
                           floor_gap=MIN_PROBE_GAP_BP)
    assert strict.summary["probes"] == 2
    assert strict.summary["genes_below_min_gap"] == 0

    pick = coverage_pick(rows, index, n_per_gene=3)
    assert pick.summary["probes"] == 3
    assert pick.summary["genes_below_min_gap"] == 1
    flag = pick.flags[0]
    assert "spacing_relaxed" in flag.reasons
    assert flag.gap_used_bp == RELAXED_PROBE_GAP_BP
    assert flag.min_gap_bp is not None and flag.min_gap_bp < MIN_PROBE_GAP_BP


def test_coverage_pick_records_spacing_on_comfortable_genes_too():
    rows = [_cov_row(f"p{i}", i * 1_000, 1.0 - i * 0.01) for i in range(4)]
    index = gene_index_from_manifest({"G1": ["ENST1"]})
    pick = coverage_pick(rows, index, n_per_gene=3)
    assert pick.summary["probes"] == 3
    assert pick.summary["genes_below_min_gap"] == 0
    assert all("spacing_relaxed" not in f.reasons for f in pick.flags)


@pytest.mark.parametrize("gap", [0, 25, 150, 400])
def test_greedy_walk_is_monotone_in_the_gap(gap):
    """Lowering the floor can never lose a pick — the property the binary
    search and the relaxation fallback both rest on."""
    rows = [_row(f"p{i}", i * 120, 1.0 - i * 0.001) for i in range(30)]
    counts = [
        len(auto_pick_with_report(rows, n_per_gene=8, min_gap=g,
                                  relaxed_gap=g).selections.get("ENST1", []))
        for g in (gap, max(0, gap - 20))
    ]
    assert counts[1] >= counts[0]


# --------------------------------------------------------------------------- #
# Interior preference — spread must not be bought by walking to the ends
# --------------------------------------------------------------------------- #
def test_spread_stays_inside_the_middle_70_percent():
    """Maximising a minimum gap, unchecked, pushes the outer probes to the very
    ends of the transcript. The interior-first pass is what stops it."""
    L = 6_000
    rows = [_row(f"p{i}", s, 0.5) for i, s in enumerate(range(0, L - 50, 50))]
    by_id = {r["probe_id"]: r for r in rows}

    pick = auto_pick_with_report(rows, n_per_gene=3, gene_length_of=lambda g: L)
    ids = pick.selections["ENST1"]
    assert len(ids) == 3
    for pid in ids:
        mid = (int(by_id[pid]["start"]) + int(by_id[pid]["end"])) / 2 / L
        assert 0.15 <= mid <= 0.85, f"{pid} at {mid:.3f} is in the terminal 15%"
    # ...and it is still a wide spread, not a huddle in the exact middle.
    assert _min_gap(by_id, ids) > 1_500


def test_ends_are_used_when_the_interior_cannot_supply_the_quota():
    """Strongly discouraged, not forbidden. Three windows exist and two of them
    sit in the terminal 15%; the gene must still get all three."""
    L = 10_000
    rows = [_row("edge5", 100, 0.9),       # midpoint at 1.25%
            _row("mid", 5_000, 0.8),       # midpoint at 50.25%
            _row("edge3", 9_800, 0.7)]     # midpoint at 98.25%
    pick = auto_pick_with_report(rows, n_per_gene=3, gene_length_of=lambda g: L)
    assert sorted(pick.selections["ENST1"]) == ["edge3", "edge5", "mid"]


def test_interior_preference_does_not_cost_probes():
    """A gene whose only spaced windows are near the 5' end keeps all three."""
    L = 4_000
    rows = ([_row(f"e{i}", i * 60, 0.9) for i in range(5)]      # all 5'-terminal
            + [_row("m", 2_000, 0.4)])
    pick = auto_pick_with_report(rows, n_per_gene=3, gene_length_of=lambda g: L)
    assert len(pick.selections["ENST1"]) == 3


# --------------------------------------------------------------------------- #
# The third rung, the tiers, and the filter that removes rungs
# --------------------------------------------------------------------------- #
def test_ladder_is_the_documented_three_rungs():
    assert MIN_PROBE_GAP_FLOOR_BP == 15
    assert SPACING_LADDER == (
        (150, 0.15),
        (40, 0.15),
        (15, 0.10),
    )
    # strictly descending in gap, non-increasing in terminal band — a rung that
    # tightened either would be unreachable
    gaps = [g for g, _ in SPACING_LADDER]
    fracs = [f for _, f in SPACING_LADDER]
    assert gaps == sorted(gaps, reverse=True) and len(set(gaps)) == len(gaps)
    assert fracs == sorted(fracs, reverse=True)


def test_gap_tier_boundaries():
    assert gap_tier(150) == GAP_TIER_GREEN
    assert gap_tier(1_000) == GAP_TIER_GREEN
    assert gap_tier(149) == GAP_TIER_ORANGE
    assert gap_tier(40) == GAP_TIER_ORANGE
    assert gap_tier(39) == GAP_TIER_RED
    assert gap_tier(0) == GAP_TIER_RED
    assert gap_tier(-1) == GAP_TIER_RED
    # A single-probe gene has no pair, so no crowding — not a spacing problem.
    assert gap_tier(None) == GAP_TIER_GREEN


def test_third_rung_rescues_a_gene_the_second_cannot():
    """Windows on a 30 bp grid: gaps are 20 bp at best, under the orange floor.

    Two probes fit at 40 bp (0 and 90, gap 40); the third needs the 15 bp rung.
    """
    L = 4_000
    rows = [_row(f"p{i}", 1_000 + i * 80, 1.0 - i * 0.01) for i in range(3)]
    by_id = {r["probe_id"]: r for r in rows}

    two = auto_pick_with_report(rows, n_per_gene=3, floor_gap=RELAXED_PROBE_GAP_BP,
                                gene_length_of=lambda g: L)
    assert len(two.selections["ENST1"]) == 2

    three = auto_pick_with_report(rows, n_per_gene=3, gene_length_of=lambda g: L)
    ids = three.selections["ENST1"]
    assert len(ids) == 3
    gap = _min_gap(by_id, ids)
    assert MIN_PROBE_GAP_FLOOR_BP <= gap < RELAXED_PROBE_GAP_BP
    assert three.spacing["ENST1"].tier == GAP_TIER_RED
    assert three.spacing["ENST1"].end_frac_used == 0.10


def test_third_rung_also_loosens_the_terminal_band():
    """The last rung is the only one that will use the 10-15% band.

    All candidates sit at ~12% of the transcript, inside the terminal 15% the
    first two rungs avoid but outside the 10% the last one avoids.
    """
    L = 10_000
    rows = [_row(f"p{i}", 1_050 + i * 80, 1.0 - i * 0.01) for i in range(3)]
    mids = [(1_050 + i * 80 + 25) / L for i in range(3)]
    assert all(0.10 < m < 0.15 for m in mids), mids
    pick = auto_pick_with_report(rows, n_per_gene=3, gene_length_of=lambda g: L)
    sp = pick.spacing["ENST1"]
    assert len(pick.selections["ENST1"]) == 3
    assert sp.end_frac_used == 0.10


def test_tier_is_earned_not_inherited_from_the_rung():
    """A gene can stop on a lower rung and still come out green.

    Only two windows exist and they are 1 kb apart, so no rung reaches three
    probes and selection walks to the bottom — but the two it has are 1 kb
    apart, and calling that gene red would be a lie.
    """
    L = 6_000
    rows = [_row("a", 1_500, 0.9), _row("b", 3_500, 0.8)]
    pick = auto_pick_with_report(rows, n_per_gene=3, gene_length_of=lambda g: L)
    sp = pick.spacing["ENST1"]
    assert sp.n_probes == 2 and sp.min_gap == 1_950
    assert sp.tier == GAP_TIER_GREEN
    assert pick.tier_counts() == {"green": 2, "orange": 0, "red": 0}


def test_tier_counts_count_probes_and_always_sum():
    rows = ([_row(f"a{i}", i * 1_000, 0.9, gene="A") for i in range(3)]
            + [_row(f"b{i}", i * 60, 0.9, gene="B") for i in range(6)])
    pick = auto_pick_with_report(rows, n_per_gene=3,
                                 gene_length_of=lambda g: 6_000)
    counts = pick.tier_counts()
    assert set(counts) == {"green", "orange", "red"}
    assert sum(counts.values()) == sum(len(v) for v in pick.selections.values())
    assert counts["green"] == 3          # gene A, 1 kb apart
    assert counts["orange"] == 3         # gene B, 70 bp apart


def test_summary_reports_every_tier_even_at_zero():
    rows = [_row(f"a{i}", i * 1_000, 0.9) for i in range(3)]
    s = auto_pick_with_report(rows, n_per_gene=3,
                              gene_length_of=lambda g: 6_000).summary()
    for name in ("green", "orange", "red"):
        assert f"probes_{name}" in s and f"genes_{name}" in s
    assert s["probes_red"] == 0 and s["genes_red"] == 0
    assert s["floor_gap_bp"] == MIN_PROBE_GAP_FLOOR_BP


# --------------------------------------------------------------------------- #
# ladder_from — the mechanism behind the GUI spacing filter
# --------------------------------------------------------------------------- #
def test_ladder_from_drops_rungs_that_cannot_loosen():
    assert ladder_from(150, 40, 15) == SPACING_LADDER
    # "no red": the floor rung collapses into the orange one and is dropped
    assert ladder_from(150, 40, 40) == ((150, 0.15), (40, 0.15))
    # "green only": one rung
    assert ladder_from(150, 150, 150) == ((150, 0.15),)
    # a floor looser than asked for is still only reached after the rung above
    assert [g for g, _ in ladder_from(150, 100, 50)] == [150, 100, 50]


def test_filtering_red_costs_probes_not_spacing():
    """Excluding a tier removes probes; it never re-spaces the ones that remain.

    This is the promise the GUI filter makes — 'red probes are never candidates'
    rather than 'red probes get quietly nudged apart'.
    """
    L = 4_000
    rows = [_row(f"p{i}", 1_000 + i * 80, 1.0 - i * 0.01) for i in range(3)]

    allowed = auto_pick_with_report(rows, n_per_gene=3, gene_length_of=lambda g: L)
    no_red = auto_pick_with_report(
        rows, n_per_gene=3, ladder=ladder_from(150, 40, 40),
        gene_length_of=lambda g: L)

    assert allowed.tier_counts()["red"] == 3
    assert no_red.tier_counts()["red"] == 0
    assert len(no_red.selections["ENST1"]) < len(allowed.selections["ENST1"])
    # every surviving pick clears the orange floor
    assert (no_red.spacing["ENST1"].min_gap or 10 ** 9) >= RELAXED_PROBE_GAP_BP


def test_green_only_filter_yields_only_green():
    L = 4_000
    rows = ([_row(f"p{i}", 1_000 + i * 80, 0.9) for i in range(3)]
            + [_row("far", 3_200, 0.5)])
    pick = auto_pick_with_report(
        rows, n_per_gene=3, ladder=ladder_from(150, 150, 150),
        gene_length_of=lambda g: L)
    assert pick.tier_counts()["orange"] == 0
    assert pick.tier_counts()["red"] == 0


def test_coverage_pick_walks_the_same_three_rungs():
    rows = [_cov_row(f"p{i}", 1_000 + i * 80, 1.0 - i * 0.01) for i in range(3)]
    index = gene_index_from_manifest({"G1": ["ENST1"]})

    two = coverage_pick(rows, index, n_per_gene=3,
                        floor_gap=RELAXED_PROBE_GAP_BP)
    assert two.summary["probes"] == 2

    three = coverage_pick(rows, index, n_per_gene=3)
    assert three.summary["probes"] == 3
    assert three.summary["probes_red"] == 3
    assert three.summary["genes_red"] == 1
    assert three.flags[0].tier_name == "red"
    assert (three.summary["probes_green"] + three.summary["probes_orange"]
            + three.summary["probes_red"]) == three.summary["probes"]


def test_descending_a_rung_never_trades_a_probe_away():
    """A looser rung is kept only when it yields MORE probes.

    Otherwise a gene could swap a well-spaced pick for an equally-sized tighter
    one and be reported as relaxed for nothing.
    """
    L = 8_000
    rows = [_row("a", 2_000, 0.9), _row("b", 5_000, 0.8)]
    pick = auto_pick_with_report(rows, n_per_gene=3, gene_length_of=lambda g: L)
    assert sorted(pick.selections["ENST1"]) == ["a", "b"]
    assert pick.spacing["ENST1"].min_gap == 2_950
