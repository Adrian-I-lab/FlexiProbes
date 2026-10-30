"""Multi-lane alignment view in the review dialog.

One lane per sequence — every isoform of the gene, or every input record of a
homology run — under the existing single summary track. These tests assert the
things a screenshot would show: how many lanes there are, which lanes a probe
is drawn on, that a run carrying none of that data still gets exactly the old
single track, and that dozens of isoforms scroll instead of growing the dialog.

Skipped automatically when PySide6 isn't installed. Runs headless via the
offscreen Qt platform.
"""
import csv
import importlib.util
import json
import os
import pathlib

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 not installed"
)


@pytest.fixture(scope="module")
def qapp():
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    yield app


# --- fixtures ---------------------------------------------------------------

def _write_run(out_dir, rows, params):
    """Write the minimal results.csv + manifest.json a ReviewDialog needs."""
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fields: list = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with open(out_dir / "results.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(out_dir / "manifest.json", "w") as fh:
        json.dump({"params": params}, fh)
    return str(out_dir)


def _probe_row(pid, gene, start, **extra):
    row = {
        "rank": str(len(pid)), "probe_id": pid, "gene_id": gene,
        "start": str(start), "end": str(start + 50),
        "gc_content_full": "0.50", "gc_lhs_half": "0.50", "gc_rhs_half": "0.50",
        "tm_target_lhs": "60", "tm_target_rhs": "60",
        "composite_score": "0.9", "specificity_pass": "True",
        "specificity_status": "passed", "isoform_hit_ids": "",
    }
    row.update({k: str(v) for k, v in extra.items()})
    return row


def _dialog(out_dir):
    from flexprobe_gui.review_dialog import ReviewDialog
    from flexprobe_gui.run_history import RunRecord
    return ReviewDialog(RunRecord(name="t", output_dir=out_dir,
                                  status="ready_for_review"))


def _select_all(dlg, gene, probe_ids):
    dlg.max_spin.setValue(max(len(probe_ids), 1))
    dlg._selections[gene] = list(probe_ids)
    dlg._refresh_viz()


# --- isoform lanes ----------------------------------------------------------

def _three_isoform_run(tmp_path):
    """One gene, three isoforms. p_all binds all three; p_two binds ENST1 and
    ENST3 only (ENST2 is where the user must SEE the gap)."""
    rows = [
        _probe_row("p_all", "ENST1", 100, isoform_hit_ids="ENST1;ENST2;ENST3"),
        _probe_row("p_two", "ENST1", 900, isoform_hit_ids="ENST1;ENST3"),
    ]
    return _write_run(tmp_path / "iso", rows, {
        "gene_lengths": {"ENST1": 2000, "ENST2": 1500, "ENST3": 2400},
        "gene_isoforms": {"ENSG1": ["ENST1", "ENST2", "ENST3"]},
        "strain_ids": [],
    })


def test_lane_count_matches_isoform_count(qapp, tmp_path):
    from flexprobe_gui.review_dialog import LaneSet
    dlg = _dialog(_three_isoform_run(tmp_path))
    try:
        _select_all(dlg, "ENST1", ["p_all", "p_two"])
        lanes = dlg.map_widget.lanes
        assert len(lanes) == 3, [ln.label for ln in lanes]
        # The source transcript leads, because its coordinates are the axis.
        assert lanes[0].label.startswith("ENST1")
        assert "[source]" in lanes[0].label
        assert [ln.label.split()[0] for ln in lanes] == ["ENST1", "ENST2", "ENST3"]
        # Each lane reports the isoform's OWN length, not the axis length.
        assert [ln.length for ln in lanes] == [2000, 1500, 2400]
        assert isinstance(dlg._build_lane_set("ENST1", []), LaneSet)
        # The coordinate system is stated in the UI, not left implied.
        assert dlg.coord_label.isVisibleTo(dlg)
        assert dlg.lane_legend.isVisibleTo(dlg)
        assert "source transcript ENST1" in dlg.coord_label.text()
    finally:
        dlg.close()


def test_probe_binding_two_of_three_isoforms_is_drawn_on_two_lanes(qapp, tmp_path):
    from flexprobe_gui.review_dialog import MARK_ABSENT, MARK_BINDS
    dlg = _dialog(_three_isoform_run(tmp_path))
    try:
        _select_all(dlg, "ENST1", ["p_all", "p_two"])
        lanes = dlg.map_widget.lanes
        state = {
            ln.label.split()[0]: {m.probe_id: m.state for m in ln.marks}
            for ln in lanes
        }
        # p_two is solid on ENST1 + ENST3 and explicitly absent on ENST2.
        drawn = [iso for iso, s in state.items() if s["p_two"] != MARK_ABSENT]
        assert drawn == ["ENST1", "ENST3"]
        assert state["ENST2"]["p_two"] == MARK_ABSENT
        # p_all binds everywhere — the contrast is the whole point of the view.
        assert all(s["p_all"] == MARK_BINDS for s in state.values())
        # Absence carries a reason, not just a missing rectangle.
        absent = next(m for m in lanes[1].marks if m.probe_id == "p_two")
        assert "isoform_hit_ids" in absent.detail
    finally:
        dlg.close()


def test_summary_tooltip_reports_binding_fraction(qapp, tmp_path):
    dlg = _dialog(_three_isoform_run(tmp_path))
    try:
        _select_all(dlg, "ENST1", ["p_all", "p_two"])
        assert dlg.map_widget._binding_summary("p_two") == " — binds 2 / 3 sequences"
        assert dlg.map_widget._binding_summary("p_all") == " — binds 3 / 3 sequences"
    finally:
        dlg.close()


# --- fallback: nothing extra to show ---------------------------------------

def test_run_without_isoform_data_falls_back_to_one_track(qapp, tmp_path):
    from flexprobe_gui.review_dialog import SUMMARY_BLOCK_H
    out = _write_run(tmp_path / "plain",
                     [_probe_row("p1", "GENE", 100)],
                     {"gene_lengths": {"GENE": 2000}, "strain_ids": []})
    dlg = _dialog(out)
    try:
        _select_all(dlg, "GENE", ["p1"])
        # No lane block at all: the map is the single track it has always been.
        assert dlg.map_widget.lanes == ()
        assert dlg.map_widget.content_height() == SUMMARY_BLOCK_H
        assert dlg.map_scroll.height() == SUMMARY_BLOCK_H
        # isVisibleTo, not isVisible: the dialog is never shown in a headless
        # test, so isVisible() is False for every widget and would pass here
        # whatever the code did.
        assert not dlg.coord_label.isVisibleTo(dlg)
        assert not dlg.lane_legend.isVisibleTo(dlg)
    finally:
        dlg.close()


def test_single_isoform_gene_gets_no_lanes(qapp, tmp_path):
    # A gene-aware run whose gene has exactly one transcript has nothing to
    # compare, so it must not sprout a one-lane block.
    out = _write_run(tmp_path / "solo",
                     [_probe_row("p1", "ENST1", 100, isoform_hit_ids="ENST1")],
                     {"gene_lengths": {"ENST1": 2000},
                      "gene_isoforms": {"ENSG1": ["ENST1"]}, "strain_ids": []})
    dlg = _dialog(out)
    try:
        _select_all(dlg, "ENST1", ["p1"])
        assert dlg.map_widget.lanes == ()
    finally:
        dlg.close()


def test_strain_ids_without_per_record_columns_fall_back(qapp, tmp_path):
    # An older homology run lists strain_ids in the manifest but carries no
    # mismatches__<id> columns. Lanes built from absent columns would be
    # fiction, so there must be none.
    out = _write_run(tmp_path / "oldstrain",
                     [_probe_row("p1", "CONS", 100)],
                     {"gene_lengths": {"CONS": 3000},
                      "strain_ids": ["A", "B"]})
    dlg = _dialog(out)
    try:
        _select_all(dlg, "CONS", ["p1"])
        assert dlg._has_strain_cols is False
        assert dlg.map_widget.lanes == ()
    finally:
        dlg.close()


# --- height cap -------------------------------------------------------------

def test_many_isoforms_scroll_instead_of_growing_the_widget(qapp, tmp_path):
    from flexprobe_gui.review_dialog import MAX_MAP_VIEWPORT_H
    isoforms = [f"ENST{i:02d}" for i in range(40)]
    out = _write_run(tmp_path / "many",
                     [_probe_row("p1", "ENST00", 100,
                                 isoform_hit_ids=";".join(isoforms))],
                     {"gene_lengths": {t: 2000 for t in isoforms},
                      "gene_isoforms": {"ENSG1": isoforms}, "strain_ids": []})
    dlg = _dialog(out)
    try:
        _select_all(dlg, "ENST00", ["p1"])
        assert len(dlg.map_widget.lanes) == 40
        # The content really is taller than the cap — otherwise this guard
        # would pass for the wrong reason.
        assert dlg.map_widget.content_height() > MAX_MAP_VIEWPORT_H
        # …and what the dialog actually gives the map stays capped.
        assert dlg.map_scroll.height() <= MAX_MAP_VIEWPORT_H
        assert dlg.map_scroll.maximumHeight() <= MAX_MAP_VIEWPORT_H
    finally:
        dlg.close()


def test_lane_count_is_capped_and_overflow_reported(qapp, tmp_path):
    from flexprobe_gui.review_dialog import MAX_LANES
    n = MAX_LANES + 5
    isoforms = [f"ENST{i:03d}" for i in range(n)]
    out = _write_run(tmp_path / "huge",
                     [_probe_row("p1", "ENST000", 100,
                                 isoform_hit_ids=";".join(isoforms))],
                     {"gene_lengths": {t: 2000 for t in isoforms},
                      "gene_isoforms": {"ENSG1": isoforms}, "strain_ids": []})
    dlg = _dialog(out)
    try:
        _select_all(dlg, "ENST000", ["p1"])
        assert len(dlg.map_widget.lanes) == MAX_LANES
        assert dlg.map_widget._overflow == 5
        assert "5 not drawn" in dlg.coord_label.text()
    finally:
        dlg.close()


# --- per-record lanes -------------------------------------------------------

def test_record_lanes_show_per_sequence_mismatch_state(qapp, tmp_path):
    from flexprobe_gui.review_dialog import MARK_ABSENT, MARK_BINDS, MARK_MISMATCH
    rows = [_probe_row(
        "p1", "CONS", 100,
        **{"mismatches__A": 0, "likely_binds__A": "True", "est_dTm__A": "-0.0",
           "mismatches__B": 2, "likely_binds__B": "True", "est_dTm__B": "-4.0",
           "mismatches__C": 9, "likely_binds__C": "False", "est_dTm__C": "-22.0"})]
    out = _write_run(tmp_path / "rec", rows,
                     {"gene_lengths": {"CONS": 3000},
                      "strain_ids": ["A", "B", "C"]})
    dlg = _dialog(out)
    try:
        _select_all(dlg, "CONS", ["p1"])
        lanes = dlg.map_widget.lanes
        assert [ln.label for ln in lanes] == ["A", "B", "C"]
        assert [ln.marks[0].state for ln in lanes] == \
            [MARK_BINDS, MARK_MISMATCH, MARK_ABSENT]
        # Tooltip detail carries mismatch count + ΔTm for that sequence.
        assert "2 mismatch(es)" in lanes[1].marks[0].detail
        assert "ΔTm -22.0" in lanes[2].marks[0].detail
        # Alignment coordinates, said out loud.
        assert "alignment (consensus) coordinates" in dlg.coord_label.text()
    finally:
        dlg.close()


def test_variant_probe_maps_onto_the_lane_it_resolves_to(qapp, tmp_path):
    # Degenerate-consensus contract: a variant is the concrete oligo for its own
    # record(s) only. Read defensively — the columns are absent in older runs.
    from flexprobe_gui.review_dialog import MARK_ABSENT
    rows = [_probe_row(
        "v_A", "CONS", 100,
        probe_role="variant", parent_probe_id="c1", resolves_to="A",
        **{"mismatches__A": 0, "likely_binds__A": "True",
           "mismatches__B": 0, "likely_binds__B": "True"})]
    out = _write_run(tmp_path / "var", rows,
                     {"gene_lengths": {"CONS": 3000}, "strain_ids": ["A", "B"]})
    dlg = _dialog(out)
    try:
        _select_all(dlg, "CONS", ["v_A"])
        a, b = dlg.map_widget.lanes
        # Lane B's own likely_binds says True, but the variant belongs to A —
        # resolves_to wins, otherwise a variant would look pan-strain.
        assert a.marks[0].state != MARK_ABSENT
        assert b.marks[0].state == MARK_ABSENT
        assert "resolves to A" in b.marks[0].detail
    finally:
        dlg.close()


def test_missing_probe_role_columns_use_likely_binds(qapp, tmp_path):
    from flexprobe_gui.review_dialog import MARK_ABSENT, MARK_BINDS
    rows = [_probe_row("p1", "CONS", 100,
                       **{"mismatches__A": 0, "likely_binds__A": "True",
                          "mismatches__B": 0, "likely_binds__B": "True"})]
    out = _write_run(tmp_path / "norole", rows,
                     {"gene_lengths": {"CONS": 3000}, "strain_ids": ["A", "B"]})
    dlg = _dialog(out)
    try:
        _select_all(dlg, "CONS", ["p1"])
        assert "probe_role" not in dlg._rows_by_gene["CONS"][0]
        assert [ln.marks[0].state for ln in dlg.map_widget.lanes] == \
            [MARK_BINDS, MARK_BINDS]
        assert all(m.state != MARK_ABSENT
                   for ln in dlg.map_widget.lanes for m in ln.marks)
    finally:
        dlg.close()


# --- the summary track must not move --------------------------------------

def test_lane_block_does_not_disturb_the_summary_track(qapp):
    """The top 110 px of the map are the single track this dialog has always
    drawn. Adding lanes below must leave it pixel-identical — that is the
    regression bar for every single-sequence run."""
    from PySide6.QtGui import QImage, QPixmap
    from flexprobe_gui.review_dialog import (
        LaneMark, LaneSet, MARK_ABSENT, MARK_BINDS, ProbeMapWidget,
        SUMMARY_BLOCK_H, SequenceLane,
    )
    probes = [{"probe_id": "a", "rank": 1, "start": 120, "end": 170},
              {"probe_id": "b", "rank": 2, "start": 900, "end": 950},
              {"probe_id": "c", "rank": 3, "start": 960, "end": 1010}]
    variants = [50, 300, 1400, 1900]
    marks = tuple(LaneMark(p["probe_id"], p["rank"], p["start"], p["end"],
                           MARK_BINDS if i else MARK_ABSENT)
                  for i, p in enumerate(probes))
    lane_set = LaneSet(lanes=(SequenceLane("iso1", marks),
                              SequenceLane("iso2", marks)),
                       coord_caption="x", kind="isoform")

    def shot(lanes):
        w = ProbeMapWidget()
        w.set_data("GENE", 2000, probes, variants, lanes)
        h = w.content_height()
        w.resize(900, h)
        pm = QPixmap(900, h)
        w.render(pm)
        return pm.toImage().convertToFormat(QImage.Format_RGBA8888)

    plain = shot(None)
    laned = shot(lane_set)
    assert plain.height() == SUMMARY_BLOCK_H
    assert laned.height() > SUMMARY_BLOCK_H          # lanes really were drawn
    assert laned.copy(0, 0, 900, SUMMARY_BLOCK_H) == plain
