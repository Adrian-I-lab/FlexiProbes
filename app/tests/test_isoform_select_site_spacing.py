"""Spacing in `isoform_select` is measured per SITE, not per emitted probe id.

`flexprobe.pipeline` builds every variant child of a degenerate consensus window
with the PARENT's own start/end, so a family's siblings are one probe at one
place. Measuring them against each other is a category error: the answer is
always -1 bp and a correctly spaced panel reports itself as entirely crowded.

The review dialog applies the same per-site rule. This file covers the non-GUI
callers.

EVERY assertion here has a control that makes the same check return the OTHER
answer. A check that says "not crowded" about every input would pass this file
while silencing real crowding, which is the failure mode worth guarding.
"""
from __future__ import annotations

from flexprobe.constants import PROBE_ROLE_CONSENSUS, PROBE_ROLE_VARIANT
from flexprobe.isoform_select import (
    _Cand, _closest_gap, _conflicts, coverage_pick, gene_index_from_manifest,
)

PROBE_LEN = 50
MIN_GAP = 150


def _cand(probe_id: str, start: int, site_key: str | None = None) -> _Cand:
    c = _Cand(probe_id=probe_id, transcript="T1", start=start,
              end=start + PROBE_LEN, target="", score=0.9,
              site_key=(site_key or probe_id))
    c.footprint["T1"] = start
    return c


# ── the defect ───────────────────────────────────────────────────────────────
def test_siblings_of_one_family_do_not_conflict():
    a, b = _cand("v00", 100, "pA"), _cand("v01", 100, "pA")
    assert _conflicts(a, b, PROBE_LEN, MIN_GAP) is False


def test_closest_gap_ignores_a_family_sharing_one_site():
    a, b = _cand("v00", 100, "pA"), _cand("v01", 100, "pA")
    # None, not -1: there is no pair whose spacing is defined, which is not the
    # same as a comfortable gap.
    assert _closest_gap([a, b], PROBE_LEN) is None


# ── CONTROLS: the same checks must still fire on real crowding ───────────────
def test_two_different_families_10bp_apart_still_conflict():
    """The control that stops this becoming 'never report crowding'."""
    a, b = _cand("v00", 100, "pA"), _cand("w00", 160, "pB")
    assert _conflicts(a, b, PROBE_LEN, MIN_GAP) is True
    assert _closest_gap([a, b], PROBE_LEN) == 10


def test_ordinary_probes_far_apart_do_not_conflict():
    a, b = _cand("x", 100), _cand("y", 400)
    assert _conflicts(a, b, PROBE_LEN, MIN_GAP) is False
    assert _closest_gap([a, b], PROBE_LEN) == 250


def test_ordinary_probes_close_together_still_conflict():
    a, b = _cand("x", 100), _cand("y", 160)
    assert _conflicts(a, b, PROBE_LEN, MIN_GAP) is True


# ── end to end, through the public entry point ───────────────────────────────
def _row(probe_id, start, role="", parent=""):
    return {
        "probe_id": probe_id, "gene_id": "G_A", "transcript_id": "T1",
        "start": str(start), "end": str(start + PROBE_LEN),
        "composite_score": "0.9", "screen_status": "pass",
        "isoform_hit_ids": "T1", "probe_role": role, "parent_probe_id": parent,
    }


def _index():
    return gene_index_from_manifest({"G_A": ["T1"]})


def test_choosing_one_child_takes_its_siblings_with_it():
    """`auto_select._walk` emits every child of a chosen parent and counts the
    family once. Coverage maximisation must agree, or the two selectors disagree
    about what a family is and one silently drops resolutions."""
    rows = [
        _row("pA", 100, PROBE_ROLE_CONSENSUS),
        _row("pA_v00", 100, PROBE_ROLE_VARIANT, "pA"),
        _row("pA_v01", 100, PROBE_ROLE_VARIANT, "pA"),
        _row("pA_v02", 100, PROBE_ROLE_VARIANT, "pA"),
    ]
    pick = coverage_pick(rows, _index(), n_per_gene=1, require_pass=False)
    chosen = {p for ids in pick.selections.values() for p in ids}
    assert "pA" not in chosen, "a consensus parent is never orderable"
    assert {"pA_v00", "pA_v01", "pA_v02"} <= chosen, (
        "all three resolutions of one degenerate window must survive, got %r" % chosen)


def test_a_family_spends_one_quota_slot_not_three():
    """n_per_gene=1 with one family plus one ordinary probe must still reach the
    ordinary probe only if a slot remains. The family is one site, so it spends
    one slot and no more."""
    rows = [
        _row("pA_v00", 100, PROBE_ROLE_VARIANT, "pA"),
        _row("pA_v01", 100, PROBE_ROLE_VARIANT, "pA"),
    ]
    pick = coverage_pick(rows, _index(), n_per_gene=1, require_pass=False)
    chosen = {p for ids in pick.selections.values() for p in ids}
    assert chosen == {"pA_v00", "pA_v01"}, (
        "one family, one slot, both children: got %r" % chosen)


def test_a_run_with_no_families_is_unchanged():
    """CONTROL for the whole change. With no degenerate rows every site_key is
    the probe's own id, so nothing above can alter the ordinary path."""
    rows = [_row("p1", 100), _row("p2", 400), _row("p3", 700)]
    pick = coverage_pick(rows, _index(), n_per_gene=2, require_pass=False)
    chosen = {p for ids in pick.selections.values() for p in ids}
    assert len(chosen) == 2, "two slots, two ordinary probes: got %r" % chosen
    assert chosen <= {"p1", "p2", "p3"}
