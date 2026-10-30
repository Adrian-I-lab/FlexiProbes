"""Behavioural tests for the discovery-mode Review dialog.

These assert what the USER experiences — the numbers on the labels, which
probes survive a filter, what lands in the order CSV — not that a method was
called. Bugs in a dialog are caught by a suite that exercises its behaviour,
not by one that checks its wiring.

Covered here:
  * the run's OWN spacing configuration drives colours, labels, the filter and
    re-picks, including the legacy no-config fallback being SAID rather
    than silently substituted,
  * selection: the per-gene cap, the spacing filter refusing a manual tick, and
    the filter retro-actively dropping picks it no longer permits,
  * the candidate filter bar (GC / Tm / score / text),
  * the export path: Confirm & export writes order CSVs restricted to the
    selection,
  * the HPC import path.

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


# --- run-directory builders -------------------------------------------------

_ROW_DEFAULTS = {
    "gc_content_full": "0.50", "gc_lhs_half": "0.50", "gc_rhs_half": "0.50",
    "tm_target_lhs": "60", "tm_target_rhs": "60",
    "composite_score": "0.90", "specificity_pass": "True",
    "specificity_status": "passed", "isoform_hit_ids": "",
    "probe_role": "",
}


def _row(pid, gene, start, rank, **extra):
    row = dict(_ROW_DEFAULTS)
    row.update({
        "rank": str(rank), "probe_id": pid, "gene_id": gene,
        "start": str(start), "end": str(start + 50),
        # export_selected_probes needs real oligo halves to write an order CSV.
        "lhs_oligo": f"LHS_{pid}", "rhs_oligo": f"/5Phos/RHS_{pid}",
        "probe_rc_50bp": "A" * 50,
    })
    row.update({k: str(v) for k, v in extra.items()})
    return row


def _write_run(out_dir, rows, params=None, config=None):
    """results.csv + manifest.json, and optionally the run_config.json sidecar."""
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fields = []
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
        json.dump({"params": params or {}}, fh)
    if config is not None:
        from flexprobe.config_io import save_config
        save_config(config, str(out_dir / "run_config.json"))
    return str(out_dir)


def _dialog(out_dir, name="t"):
    from flexprobe_gui.review_dialog import ReviewDialog
    from flexprobe_gui.run_history import RunRecord
    return ReviewDialog(RunRecord(name=name, output_dir=str(out_dir),
                                  status="ready_for_review"))


def _wide_config(**over):
    """A RunConfig designed at 300 / 200 / 100 bp — nothing like the defaults."""
    from flexprobe.models import RunConfig
    kw = dict(fasta_path="x.fa", min_probe_gap_bp=300,
              relaxed_probe_gap_bp=200, floor_probe_gap_bp=100)
    kw.update(over)
    return RunConfig(**kw)


def _two_probes_250_apart(tmp_path, **run_kw):
    """One gene, two probes whose gap is exactly 250 bp.

    250 is the number that separates the two configurations: it is legal at the
    150 bp default minimum and illegal at a 300 bp one.
    """
    rows = [_row("pA", "G1", 0, 1), _row("pB", "G1", 300, 2)]
    assert 300 - (0 + 50) == 250
    return _write_run(tmp_path / "run", rows,
                      params={"gene_lengths": {"G1": 2000}}, **run_kw)


# ---------------------------------------------------------------------------
# The run's own spacing configuration
# ---------------------------------------------------------------------------

def test_run_configured_at_300bp_is_shown_as_300bp(qapp, tmp_path):
    """A 300 bp run must SAY 300 bp — in the legend, the filter and the
    tooltips. A label carrying the module default is the same bug as a filter
    carrying it."""
    out = _two_probes_250_apart(tmp_path, config=_wide_config())
    dlg = _dialog(out)

    assert (dlg._spacing.min_gap, dlg._spacing.relaxed_gap,
            dlg._spacing.floor_gap) == (300, 200, 100)
    assert dlg._spacing.from_run_config is True

    filter_labels = [dlg.spacing_filter.itemText(i)
                     for i in range(dlg.spacing_filter.count())]
    assert any("300 bp" in t for t in filter_labels), filter_labels
    assert any("200 bp" in t for t in filter_labels), filter_labels
    assert not any("150 bp" in t for t in filter_labels), filter_labels

    legend = _spacing_legend_text(dlg)
    assert "300 bp" in legend and "200" in legend and "100 bp" in legend
    assert "150" not in legend

    for tip in (dlg.spacing_filter.toolTip(), dlg.autopick_btn.toolTip()):
        assert "300 bp" in tip and "150 bp" not in tip

    dlg.deleteLater()


def _spacing_legend_text(dlg):
    from PySide6.QtWidgets import QLabel
    for lbl in dlg.findChildren(QLabel):
        if lbl.accessibleName() == "Probe spacing legend":
            return lbl.text()
    raise AssertionError("no spacing legend label found")


def test_filtering_at_300bp_hides_probes_250bp_apart(qapp, tmp_path):
    """The user's experience: with the run configured at 300 bp, choosing
    'Green only' must refuse a second probe 250 bp from the first — and at the
    150 bp default the very same pair must be allowed."""
    from flexprobe_gui.review_dialog import GAP_TIER_GREEN

    out = _two_probes_250_apart(tmp_path, config=_wide_config())
    dlg = _dialog(out)
    dlg.max_spin.setValue(2)
    _select_green_only(dlg, GAP_TIER_GREEN)
    assert dlg._spacing_floor() == 300

    dlg._selections["G1"] = ["pA"]
    assert dlg._blocked_by_spacing_filter("G1", "pB") == 250

    # ...and the same pair already selected is pruned back to one.
    dlg._selections["G1"] = ["pA", "pB"]
    dlg._on_spacing_filter_changed(dlg.spacing_filter.currentIndex())
    assert dlg._selections.get("G1") == ["pA"]
    dlg.deleteLater()


def test_same_pair_is_legal_at_the_150bp_default(qapp, tmp_path):
    """Control for the test above: nothing about the 250 bp pair itself is
    illegal — it is the run's configuration that makes it so."""
    from flexprobe_gui.review_dialog import GAP_TIER_GREEN

    out = _two_probes_250_apart(tmp_path / "default")
    dlg = _dialog(out)
    dlg.max_spin.setValue(2)
    _select_green_only(dlg, GAP_TIER_GREEN)
    assert dlg._spacing_floor() == 150
    dlg._selections["G1"] = ["pA"]
    assert dlg._blocked_by_spacing_filter("G1", "pB") is None
    dlg.deleteLater()


def _select_green_only(dlg, green_tier):
    for i in range(dlg.spacing_filter.count()):
        if dlg.spacing_filter.itemData(i)[1] == green_tier:
            dlg.spacing_filter.setCurrentIndex(i)
            return
    raise AssertionError("no green-only spacing filter entry")


def test_300bp_run_colours_a_250bp_pair_orange_not_green(qapp, tmp_path):
    """The map's colours come from the run's rungs too. 250 bp is green at the
    default minimum and orange at a 300 bp one — so this is the tier the user
    sees on the bars, and the tally under them."""
    out = _two_probes_250_apart(tmp_path, config=_wide_config())
    dlg = _dialog(out)
    dlg.max_spin.setValue(2)
    dlg._selections["G1"] = ["pA", "pB"]
    dlg._refresh_viz()
    assert dlg.map_widget.has_tight is True
    assert dlg.map_widget.has_overlap is False
    assert dlg.selection_tier_counts() == {"green": 0, "orange": 2, "red": 0}
    dlg.deleteLater()

    dflt = _dialog(_two_probes_250_apart(tmp_path / "default"))
    dflt.max_spin.setValue(2)
    dflt._selections["G1"] = ["pA", "pB"]
    dflt._refresh_viz()
    assert dflt.map_widget.has_tight is False
    assert dflt.selection_tier_counts() == {"green": 2, "orange": 0, "red": 0}
    dflt.deleteLater()


def test_repick_uses_the_runs_ladder_not_the_default_one(qapp, tmp_path):
    """'Pick best N (spaced)' must walk the run's rungs."""
    out = _two_probes_250_apart(tmp_path, config=_wide_config())
    dlg = _dialog(out)
    rungs = [gap for gap, _tol in dlg._pick_ladder()]
    assert rungs == [300, 200, 100]
    dlg.deleteLater()


def test_missing_run_config_falls_back_and_says_so(qapp, tmp_path):
    """A legacy or hand-made output folder has no run_config.json. The dialog
    may fall back to the module defaults, but it must NOT pretend they are the
    run's own settings."""
    out = _two_probes_250_apart(tmp_path / "legacy")
    dlg = _dialog(out)
    assert dlg._spacing.from_run_config is False
    assert dlg._spacing.min_gap == 150

    note = dlg.spacing_source_label.text()
    assert dlg.spacing_source_label.isVisibleTo(dlg)
    assert "No run_config.json" in note
    assert "application defaults" in note
    assert "150" in note
    # and the same caveat reaches the tooltips
    assert "application defaults" in dlg.spacing_filter.toolTip()
    dlg.deleteLater()


def test_present_run_config_is_named_as_the_source(qapp, tmp_path):
    out = _two_probes_250_apart(tmp_path, config=_wide_config())
    dlg = _dialog(out)
    note = dlg.spacing_source_label.text()
    assert "run_config.json" in note and "No run_config.json" not in note
    assert "300" in note
    dlg.deleteLater()


# ---------------------------------------------------------------------------
# Selection behaviour
# ---------------------------------------------------------------------------

def _tick(dlg, probe_id, checked=True):
    """Tick/untick a probe's Select box the way a user would — through the
    table item, so the dialog's own itemChanged handling really runs."""
    from PySide6.QtCore import Qt
    from flexprobe_gui.review_dialog import SELECT_COL
    for row in range(dlg.table.rowCount()):
        item = dlg.table.item(row, SELECT_COL)
        if item.data(Qt.UserRole) == probe_id:
            item.setCheckState(Qt.Checked if checked else Qt.Unchecked)
            return item
    raise AssertionError(f"{probe_id} is not in the visible table")


def _spread_run(tmp_path, n=5, step=1000):
    """One gene, ``n`` probes 950 bp apart — legal at every rung."""
    rows = [_row(f"p{i}", "G1", i * step, i + 1) for i in range(n)]
    return _write_run(tmp_path / "spread", rows,
                      params={"gene_lengths": {"G1": n * step + 100}})


def test_ticking_a_box_selects_the_probe_and_updates_the_counter(qapp, tmp_path):
    dlg = _dialog(_spread_run(tmp_path))
    dlg.max_spin.setValue(3)
    _tick(dlg, "p0")
    _tick(dlg, "p2")
    assert dlg._selections["G1"] == ["p0", "p2"]
    assert "Selected for G1: 2 / 3" in dlg.counter_label.text()
    _tick(dlg, "p0", checked=False)
    assert dlg._selections["G1"] == ["p2"]
    dlg.deleteLater()


def test_the_per_gene_cap_refuses_the_extra_tick(qapp, tmp_path, monkeypatch):
    """Not "warns about" — refuses: the box goes back to unchecked and the
    probe never enters the selection."""
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QMessageBox
    from flexprobe_gui.review_dialog import SELECT_COL
    monkeypatch.setattr(QMessageBox, "information",
                        staticmethod(lambda *a, **k: QMessageBox.Ok))

    dlg = _dialog(_spread_run(tmp_path))
    dlg.max_spin.setValue(2)
    _tick(dlg, "p0")
    _tick(dlg, "p1")
    item = _tick(dlg, "p2")
    assert dlg._selections["G1"] == ["p0", "p1"]
    assert item.checkState() == Qt.Unchecked
    assert dlg.table.item(item.row(), SELECT_COL) is item
    dlg.deleteLater()


def test_lowering_the_cap_trims_existing_selections(qapp, tmp_path):
    dlg = _dialog(_spread_run(tmp_path))
    dlg.max_spin.setValue(4)
    for pid in ("p0", "p1", "p2", "p3"):
        _tick(dlg, pid)
    dlg.max_spin.setValue(2)
    assert dlg._selections["G1"] == ["p0", "p1"]
    dlg.deleteLater()


def test_selections_survive_switching_genes_and_back(qapp, tmp_path):
    rows = [_row("a1", "GA", 0, 1), _row("a2", "GA", 1000, 2),
            _row("b1", "GB", 0, 1)]
    out = _write_run(tmp_path / "two", rows,
                     params={"gene_lengths": {"GA": 2000, "GB": 2000}})
    dlg = _dialog(out)
    dlg.max_spin.setValue(2)
    assert dlg.gene_combo.count() == 2
    _tick(dlg, "a1")
    dlg.gene_combo.setCurrentIndex(1)
    _tick(dlg, "b1")
    dlg.gene_combo.setCurrentIndex(0)
    assert dlg._selections == {"GA": ["a1"], "GB": ["b1"]}
    assert "Total across genes: 2" in dlg.counter_label.text()
    dlg.deleteLater()


def test_autopick_respects_the_per_gene_maximum(qapp, tmp_path):
    dlg = _dialog(_spread_run(tmp_path, n=6))
    dlg.max_spin.setValue(3)
    dlg._on_autopick()
    assert len(dlg._selections["G1"]) == 3
    dlg.deleteLater()


# ---------------------------------------------------------------------------
# The candidate filter bar
# ---------------------------------------------------------------------------

def _visible_probe_ids(dlg):
    from PySide6.QtCore import Qt
    from flexprobe_gui.review_dialog import SELECT_COL
    return {dlg.table.item(r, SELECT_COL).data(Qt.UserRole)
            for r in range(dlg.table.rowCount())}


def _mixed_run(tmp_path):
    rows = [
        _row("good", "G1", 0, 1, gc_content_full="0.50",
             tm_target_lhs="65", tm_target_rhs="66", composite_score="0.90"),
        _row("gc_low", "G1", 1000, 2, gc_content_full="0.20",
             tm_target_lhs="65", tm_target_rhs="66", composite_score="0.90"),
        _row("tm_low", "G1", 2000, 3, gc_content_full="0.50",
             tm_target_lhs="65", tm_target_rhs="48", composite_score="0.90"),
        _row("score_low", "G1", 3000, 4, gc_content_full="0.50",
             tm_target_lhs="65", tm_target_rhs="66", composite_score="0.10"),
    ]
    return _write_run(tmp_path / "mixed", rows,
                      params={"gene_lengths": {"G1": 5000}})


def test_gc_filter_hides_rows_outside_the_band(qapp, tmp_path):
    dlg = _dialog(_mixed_run(tmp_path))
    assert _visible_probe_ids(dlg) == {"good", "gc_low", "tm_low", "score_low"}
    dlg.f_gc_min.setValue(0.40)
    assert "gc_low" not in _visible_probe_ids(dlg)
    dlg.deleteLater()


def test_tm_floor_needs_BOTH_halves_above_it(qapp, tmp_path):
    """The weaker half limits hybridisation, so a 65/48 probe must fail a 50 C
    floor even though one half clears it."""
    dlg = _dialog(_mixed_run(tmp_path))
    dlg.f_tm_min.setValue(50.0)
    assert "tm_low" not in _visible_probe_ids(dlg)
    assert "good" in _visible_probe_ids(dlg)
    dlg.deleteLater()


def test_score_floor_and_text_filter(qapp, tmp_path):
    dlg = _dialog(_mixed_run(tmp_path))
    dlg.f_score_min.setValue(0.50)
    assert "score_low" not in _visible_probe_ids(dlg)
    dlg.f_score_min.setValue(0.0)
    dlg.f_text.setText("gc")
    assert _visible_probe_ids(dlg) == {"gc_low"}
    dlg.deleteLater()


def test_clear_restores_every_row_and_the_counter_says_so(qapp, tmp_path):
    dlg = _dialog(_mixed_run(tmp_path))
    dlg.f_gc_min.setValue(0.40)
    dlg.f_score_min.setValue(0.50)
    assert "showing 2/4 candidates" in dlg.counter_label.text()
    dlg.f_clear.click()
    assert len(_visible_probe_ids(dlg)) == 4
    assert "showing 4/4 candidates" in dlg.counter_label.text()
    dlg.deleteLater()


def test_a_filtered_out_probe_stays_selected(qapp, tmp_path):
    """Hiding a row is a view change, not a de-selection — losing the pick
    silently would be a data-loss bug."""
    dlg = _dialog(_mixed_run(tmp_path))
    dlg.max_spin.setValue(2)
    _tick(dlg, "gc_low")
    dlg.f_gc_min.setValue(0.40)
    assert "gc_low" not in _visible_probe_ids(dlg)
    assert dlg._selections["G1"] == ["gc_low"]
    dlg.deleteLater()


# ---------------------------------------------------------------------------
# The export path
# ---------------------------------------------------------------------------

def test_confirm_exports_only_the_selected_probes(qapp, tmp_path, monkeypatch):
    """What the user gets to order. The order CSV must contain the probes they
    ticked and nothing else, and the run must be reported as confirmed."""
    from PySide6.QtWidgets import QMessageBox
    monkeypatch.setattr(QMessageBox, "information",
                        staticmethod(lambda *a, **k: QMessageBox.Ok))

    out = _spread_run(tmp_path, n=4)
    (pathlib.Path(out) / "manifest.json").write_text(
        json.dumps({"params": {"modality": "flex_v2",
                               "gene_lengths": {"G1": 4100}}}))
    dlg = _dialog(out)
    dlg.max_spin.setValue(2)
    _tick(dlg, "p1")
    _tick(dlg, "p3")

    emitted = []
    dlg.confirmed.connect(emitted.append)
    dlg._on_confirm()

    idt = pathlib.Path(out) / "idt_opools.csv"
    assert idt.exists(), sorted(p.name for p in pathlib.Path(out).iterdir())
    names = [r[1] for r in list(csv.reader(idt.open()))[1:]]
    assert sorted(names) == ["p1_LHS", "p1_RHS", "p3_LHS", "p3_RHS"]
    assert not any(n.startswith(("p0_", "p2_")) for n in names)

    saved = json.loads((pathlib.Path(out) / "selections.json").read_text())
    assert saved["selections"] == {"G1": ["p1", "p3"]}
    assert saved["note"] == "confirmed"
    assert emitted == [dlg._record.run_id]


def test_confirm_with_nothing_selected_exports_nothing(qapp, tmp_path,
                                                       monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    warned = []
    monkeypatch.setattr(QMessageBox, "warning",
                        staticmethod(lambda *a, **k: warned.append(a)))
    out = _spread_run(tmp_path)
    dlg = _dialog(out)
    dlg._on_confirm()
    assert warned
    assert not (pathlib.Path(out) / "idt_opools.csv").exists()
    dlg.deleteLater()


def test_save_draft_writes_selections_without_exporting(qapp, tmp_path,
                                                        monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    monkeypatch.setattr(QMessageBox, "information",
                        staticmethod(lambda *a, **k: QMessageBox.Ok))
    out = _spread_run(tmp_path)
    dlg = _dialog(out)
    dlg.max_spin.setValue(2)
    _tick(dlg, "p0")
    dlg._on_save_draft()
    saved = json.loads((pathlib.Path(out) / "selections.json").read_text())
    assert saved["selections"] == {"G1": ["p0"]}
    assert saved["note"] == "draft"
    assert not (pathlib.Path(out) / "idt_opools.csv").exists()
    dlg.deleteLater()


def test_a_saved_draft_is_reloaded_when_the_dialog_reopens(qapp, tmp_path,
                                                           monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    from flexprobe.select_export import clear_results_cache
    monkeypatch.setattr(QMessageBox, "information",
                        staticmethod(lambda *a, **k: QMessageBox.Ok))
    out = _spread_run(tmp_path)
    dlg = _dialog(out)
    dlg.max_spin.setValue(2)
    _tick(dlg, "p1")
    dlg._on_save_draft()
    dlg.deleteLater()

    clear_results_cache()
    again = _dialog(out)
    assert again._selections == {"G1": ["p1"]}
    again.deleteLater()


# ---------------------------------------------------------------------------
# The HPC import path
# ---------------------------------------------------------------------------

def test_importing_a_cluster_folder_registers_a_reviewable_run(
        qapp, tmp_path, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    from flexprobe_gui.hpc_dialog import ImportHpcDialog
    from flexprobe_gui.run_history import RunHistory
    monkeypatch.setattr(QMessageBox, "information",
                        staticmethod(lambda *a, **k: QMessageBox.Ok))

    out = _spread_run(tmp_path, n=3)
    history = RunHistory(str(tmp_path / "history.json"))
    dlg = ImportHpcDialog(history)
    dlg.src_edit.setText(out)
    dlg.name_edit.setText("cluster run")
    dlg._on_import()

    records = history.list_completed()
    assert len(records) == 1
    rec = records[0]
    assert rec.name == "cluster run"
    assert rec.status == "ready_for_review"
    assert (rec.stats or {}).get("imported") is True
    assert rec.output_dir == out

    # ...and the imported run really opens in the review dialog.
    from flexprobe_gui.review_dialog import ReviewDialog
    review = ReviewDialog(rec)
    assert review.gene_combo.count() == 1
    assert len(review._rows_by_gene["G1"]) == 3
    review.deleteLater()


def test_importing_an_archive_extracts_and_registers_it(qapp, tmp_path,
                                                        monkeypatch):
    import tarfile
    from PySide6.QtWidgets import QMessageBox
    from flexprobe_gui.hpc_dialog import ImportHpcDialog
    from flexprobe_gui.run_history import RunHistory
    monkeypatch.setattr(QMessageBox, "information",
                        staticmethod(lambda *a, **k: QMessageBox.Ok))

    out = _spread_run(tmp_path, n=2)
    archive = tmp_path / "flexiprobes_output.tar.gz"
    with tarfile.open(archive, "w:gz") as tf:
        tf.add(out, arcname="output")

    history = RunHistory(str(tmp_path / "history.json"))
    dlg = ImportHpcDialog(history)
    dlg.src_edit.setText(str(archive))
    dlg._on_import()

    records = history.list_completed()
    assert len(records) == 1
    assert os.path.exists(os.path.join(records[0].output_dir, "results.csv"))
    assert records[0].status == "ready_for_review"


def test_importing_a_folder_with_no_results_csv_is_refused(qapp, tmp_path,
                                                           monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    from flexprobe_gui.hpc_dialog import ImportHpcDialog
    from flexprobe_gui.run_history import RunHistory
    errors = []
    monkeypatch.setattr(QMessageBox, "critical",
                        staticmethod(lambda *a, **k: errors.append(a)))

    empty = tmp_path / "empty"
    empty.mkdir()
    history = RunHistory(str(tmp_path / "history.json"))
    dlg = ImportHpcDialog(history)
    dlg.src_edit.setText(str(empty))
    dlg._on_import()

    assert errors
    assert history.list_completed() == []


# ---------------------------------------------------------------------------
# Accessibility, and the spacing tier without colour
# ---------------------------------------------------------------------------

def test_every_interactive_control_has_an_accessible_name(qapp, tmp_path):
    """setAccessibleName appeared zero times in the whole GUI package: to a
    screen reader these controls were unlabelled boxes."""
    dlg = _dialog(_spread_run(tmp_path))
    controls = [
        dlg.gene_combo, dlg.max_spin, dlg.spacing_filter, dlg.autopick_btn,
        dlg.coverage_btn, dlg.f_gc_min, dlg.f_gc_max, dlg.f_tm_min,
        dlg.f_score_min, dlg.f_text, dlg.f_clear, dlg.table,
        dlg.save_btn, dlg.confirm_btn,
    ]
    unnamed = [c for c in controls if not c.accessibleName().strip()]
    assert unnamed == [], [type(c).__name__ for c in unnamed]
    dlg.deleteLater()


def test_label_field_pairs_are_joined_with_setBuddy(qapp, tmp_path):
    """A buddy is what gives a field its name from its visible label and makes
    the label's mnemonic move focus to it."""
    from PySide6.QtWidgets import QLabel
    dlg = _dialog(_spread_run(tmp_path))
    buddies = {lbl.buddy() for lbl in dlg.findChildren(QLabel)
               if lbl.buddy() is not None}
    for field in (dlg.gene_combo, dlg.max_spin, dlg.spacing_filter,
                  dlg.f_gc_min, dlg.f_tm_min, dlg.f_score_min):
        assert field in buddies, field.accessibleName()
    dlg.deleteLater()


def test_the_buttons_carry_mnemonics(qapp, tmp_path):
    dlg = _dialog(_spread_run(tmp_path))
    for btn in (dlg.autopick_btn, dlg.coverage_btn, dlg.save_btn,
                dlg.confirm_btn, dlg.f_clear):
        assert "&" in btn.text(), btn.text()
    # "Confirm && export": the doubled ampersand is a literal, the single one
    # is the mnemonic.
    assert dlg.confirm_btn.text() == "Confirm && &export"
    dlg.deleteLater()


def test_tab_order_follows_the_dialog_in_reading_order(qapp, tmp_path):
    """Qt's default tab chain follows widget CREATION order, which here
    interleaves the map's decorations with the controls."""
    dlg = _dialog(_spread_run(tmp_path))
    chain = [dlg.gene_combo, dlg.max_spin, dlg.spacing_filter,
             dlg.autopick_btn, dlg.coverage_btn, dlg.f_gc_min, dlg.f_gc_max,
             dlg.f_tm_min, dlg.f_score_min, dlg.f_text, dlg.f_clear,
             dlg.show_parents_check, dlg.table, dlg.save_btn, dlg.confirm_btn]
    # Neither a disabled nor a hidden control is a Tab stop. The coverage button
    # is disabled on a run with no isoform data, and the consensus-parent toggle
    # is hidden on a run with no degenerate families; both drop out here.
    chain = [w for w in chain if w.isEnabled() and w.isVisibleTo(dlg)]
    assert len(chain) >= 12
    for first, second in zip(chain, chain[1:]):
        assert _next_focusable(first) is second, (
            f"{first.accessibleName()} -> got "
            f"{_describe(_next_focusable(first))}, wanted "
            f"{second.accessibleName()}")
    dlg.deleteLater()


def _is_inside(widget, ancestor):
    w = widget
    while w is not None:
        if w is ancestor:
            return True
        w = w.parentWidget()
    return False


def _next_focusable(first, limit=80):
    """The next Tab stop after ``first``, skipping its own internals.

    A QSpinBox and a QComboBox each own a child QLineEdit that sits in the
    focus chain; those are part of the same Tab stop, not the next one.

    Hidden widgets are skipped too. Qt does that itself once the dialog is on
    screen, but these dialogs are never shown, so ``isVisible()`` is False for
    everything and the question has to be asked as ``isVisibleTo`` the window:
    would this widget be visible if its window were.
    """
    from PySide6.QtCore import Qt
    w = first
    for _ in range(limit):
        w = w.nextInFocusChain()
        if w is first or w is None:
            return None
        if _is_inside(w, first):
            continue
        if (w.focusPolicy() != Qt.NoFocus and w.isEnabled()
                and w.isVisibleTo(w.window())):
            return w
    return None


def _describe(w):
    if w is None:
        return "nothing"
    return w.accessibleName() or type(w).__name__


def test_the_spacing_tier_is_readable_without_colour(qapp, tmp_path):
    """Red/green fill alone excludes a colour-blind reader, a screen reader and
    any greyscale print. Each tier must also have a letter, a fill pattern and
    a name in words."""
    from flexprobe_gui.review_dialog import (
        GAP_TIER_GREEN, GAP_TIER_ORANGE, GAP_TIER_RED, _GAP_TIER_GLYPHS,
        _GAP_TIER_PATTERNS, tier_phrase,
    )
    tiers = (GAP_TIER_GREEN, GAP_TIER_ORANGE, GAP_TIER_RED)
    glyphs = [_GAP_TIER_GLYPHS[t] for t in tiers]
    assert glyphs == ["G", "O", "R"]
    assert len(set(glyphs)) == 3
    # The two tiers that need attention are also texture-distinct.
    assert _GAP_TIER_PATTERNS[GAP_TIER_ORANGE] != _GAP_TIER_PATTERNS[GAP_TIER_RED]
    assert _GAP_TIER_PATTERNS[GAP_TIER_GREEN] is None
    for t in tiers:
        phrase = tier_phrase(t)
        assert _GAP_TIER_GLYPHS[t] in phrase
        assert len(phrase.split(":", 1)[1].strip()) > 10   # actual words


def test_the_legend_names_the_letters_not_only_the_colours(qapp, tmp_path):
    dlg = _dialog(_spread_run(tmp_path))
    legend = _spacing_legend_text(dlg)
    for glyph, word in (("G", "green"), ("O", "orange"), ("R", "red")):
        assert f"<b>{glyph}</b>" in legend, legend
        assert word in legend
    dlg.deleteLater()


def test_the_probe_map_describes_its_selection_in_words(qapp, tmp_path):
    """The map is a QPainter drawing with no child widgets: without an
    accessible description it is an unlabelled rectangle."""
    dlg = _dialog(_spread_run(tmp_path))
    assert dlg.map_widget.accessibleName() == "Probe map"
    assert "no probes selected" in dlg.map_widget.accessibleDescription()
    dlg.max_spin.setValue(2)
    _tick(dlg, "p0")
    _tick(dlg, "p1")
    desc = dlg.map_widget.accessibleDescription()
    assert "G1" in desc and "2 selected probe(s)" in desc and "green" in desc
    dlg.deleteLater()


def test_the_map_tooltip_names_the_tier_in_words(qapp, tmp_path):
    """The hover text gives the spacing verdict in words, not only as a fill
    colour, alongside position and binding."""
    from PySide6.QtCore import QPointF
    from PySide6.QtGui import QMouseEvent
    from PySide6.QtCore import QEvent, Qt
    from PySide6.QtWidgets import QToolTip

    out = _two_probes_250_apart(tmp_path, config=_wide_config())
    dlg = _dialog(out)
    dlg.max_spin.setValue(2)
    dlg._selections["G1"] = ["pA", "pB"]
    dlg._refresh_viz()

    shown = []
    real = QToolTip.showText
    QToolTip.showText = staticmethod(lambda pos, text, *a: shown.append(text))
    try:
        w = dlg.map_widget
        w.resize(600, 200)
        track_x, track_w = w._track_span()
        # x for position 10 (inside probe pA, which spans 0-50)
        x = track_x + track_w * 10 / max(w._gene_length, 1)
        ev = QMouseEvent(QEvent.MouseMove, QPointF(x, 60), QPointF(x, 60),
                         Qt.NoButton, Qt.NoButton, Qt.NoModifier)
        w.mouseMoveEvent(ev)
    finally:
        QToolTip.showText = real

    assert shown, "no tooltip was shown"
    text = shown[0]
    assert "pA" in text
    assert "orange" in text and "tight" in text, text
    dlg.deleteLater()


# ---------------------------------------------------------------------------
# Degenerate consensus families in the review window
#
# A consensus parent and its variant children occupy the SAME [start, end) by
# construction: pipeline.py builds every child with the parent's start/end.
# They are one probe at one site: "one family = one pick" in auto_select.
#
# Two things follow, and the dialog got both wrong.
#   1. Spacing must be measured per SITE. Measuring it per emitted id compares
#      a family's siblings against each other, always reads -1 bp, and paints a
#      correctly spaced panel entirely red.
#   2. A parent is never orderable and never exported, so it does not belong in
#      the review table by default.
# ---------------------------------------------------------------------------

def _degenerate_run(tmp_path, n_variants=2, other_start=1000):
    """One gene: a degenerate family of ``n_variants`` at 0-50, plus one
    ordinary probe at ``other_start``.

    The family shares a single 50 bp footprint. The ordinary probe sits 950 bp
    clear of it, which is green at every rung this project ships. So the only
    way anything here can read as crowded is if siblings are measured against
    each other.
    """
    rows = [_row("pA", "G1", 0, 1, probe_role="consensus", n_variants=n_variants)]
    for i in range(n_variants):
        rows.append(_row(f"pA_v{i:02d}", "G1", 0, 1,
                         probe_role="variant", parent_probe_id="pA",
                         resolves_to=f"rec{i}", n_variants=n_variants))
    rows.append(_row("pB", "G1", other_start, 2))
    assert other_start - 50 == 950
    return _write_run(tmp_path / "degenerate", rows,
                      params={"gene_lengths": {"G1": 3000}})


def _table_probe_ids(dlg):
    """Every probe id the review table is currently showing, in row order."""
    from PySide6.QtCore import Qt
    from flexprobe_gui.review_dialog import SELECT_COL
    return [dlg.table.item(r, SELECT_COL).data(Qt.UserRole)
            for r in range(dlg.table.rowCount())]


def test_autopick_on_a_degenerate_gene_does_not_return_parents(qapp, tmp_path):
    """Guards the picker's contract, which the dialog relies on.

    ``auto_pick_per_gene`` ranks parents but emits their CHILDREN. If that ever
    breaks, the view-layer tests below would start passing for the wrong
    reason.
    """
    dlg = _dialog(_degenerate_run(tmp_path))
    dlg.max_spin.setValue(2)
    dlg._on_autopick()
    picked = dlg._selections["G1"]
    assert "pA" not in picked, picked
    assert set(picked) == {"pA_v00", "pA_v01", "pB"}, picked
    dlg.deleteLater()


def test_a_families_siblings_are_not_measured_against_each_other(qapp, tmp_path):
    """Sibling variants must not turn the whole set red.

    Two variants of one parent share a footprint by design, so the gap between
    them is -1 bp. Reading that as crowding is a category error. There is one
    site here, 950 bp from the next, which is green. A red bar on this panel is
    a false alarm on a correct design.
    """
    dlg = _dialog(_degenerate_run(tmp_path))
    dlg.max_spin.setValue(2)
    dlg._on_autopick()

    assert dlg.map_widget.has_overlap is False, (
        "sibling variants of one consensus parent were read as overlapping")
    assert dlg.map_widget.has_tight is False
    # 2 SITES (one family of 2 variants, one ordinary probe), counted per
    # site. Counting ids would give 3.
    assert dlg.selection_tier_counts() == {"green": 2, "orange": 0, "red": 0}
    assert dlg.overlap_warning.text() == ""
    dlg.deleteLater()


def test_real_crowding_between_two_families_is_still_red(qapp, tmp_path):
    """The negative control for the test above.

    Collapsing a family to one site must not collapse the check itself. Put a
    second family 10 bp from the first and the panel must still go red, or the
    family handling has simply switched the warning off.
    """
    rows = [
        _row("pA", "G1", 0, 1, probe_role="consensus", n_variants=2),
        _row("pA_v00", "G1", 0, 1, probe_role="variant",
             parent_probe_id="pA", resolves_to="rec0", n_variants=2),
        _row("pA_v01", "G1", 0, 1, probe_role="variant",
             parent_probe_id="pA", resolves_to="rec1", n_variants=2),
        _row("pC", "G1", 60, 2, probe_role="consensus", n_variants=2),
        _row("pC_v00", "G1", 60, 2, probe_role="variant",
             parent_probe_id="pC", resolves_to="rec0", n_variants=2),
        _row("pC_v01", "G1", 60, 2, probe_role="variant",
             parent_probe_id="pC", resolves_to="rec1", n_variants=2),
    ]
    assert 60 - 50 == 10           # 10 bp apart: under every rung this ships
    out = _write_run(tmp_path / "crowded", rows,
                     params={"gene_lengths": {"G1": 3000}})
    dlg = _dialog(out)
    dlg.max_spin.setValue(2)
    dlg._selections["G1"] = ["pA_v00", "pA_v01", "pC_v00", "pC_v01"]
    dlg._refresh_viz()

    assert dlg.map_widget.has_overlap is True, (
        "two DIFFERENT families 10 bp apart must still read as crowded")
    # Two families, so 2 red SITES (it read 4 when the tally counted ids).
    assert dlg.selection_tier_counts()["red"] == 2
    dlg.deleteLater()


def test_the_review_table_hides_consensus_parents_by_default(qapp, tmp_path):
    """Consensus parents must not appear in the review window.

    A parent is a mixed-base string. It is never ordered and never exported.
    Showing it next to its own children, on the same
    coordinates, reads as a duplicate the user must mentally discount.
    """
    dlg = _dialog(_degenerate_run(tmp_path))
    shown = _table_probe_ids(dlg)
    assert "pA" not in shown, shown
    assert shown == ["pA_v00", "pA_v01", "pB"], shown
    dlg.deleteLater()


def test_a_toggle_brings_the_parents_back(qapp, tmp_path):
    """Hiding them by default must not make them unreachable. The parent row
    carries the IUPAC sequence and the degenerate positions, which is the only
    place a user can see WHY a family has the variants it has."""
    dlg = _dialog(_degenerate_run(tmp_path))
    assert dlg.show_parents_check.isChecked() is False
    dlg.show_parents_check.setChecked(True)
    shown = _table_probe_ids(dlg)
    assert "pA" in shown, shown
    assert shown == ["pA", "pA_v00", "pA_v01", "pB"], shown
    dlg.show_parents_check.setChecked(False)
    assert "pA" not in _table_probe_ids(dlg)
    dlg.deleteLater()


def test_a_shown_parent_is_still_not_orderable(qapp, tmp_path):
    """With the toggle on, the parent is visible but must stay unselectable and
    must not enter the spacing arithmetic. Visibility is a view choice; being
    orderable is not."""
    from PySide6.QtCore import Qt
    from flexprobe_gui.review_dialog import SELECT_COL

    dlg = _dialog(_degenerate_run(tmp_path))
    dlg.max_spin.setValue(2)
    dlg.show_parents_check.setChecked(True)
    row = _table_probe_ids(dlg).index("pA")
    item = dlg.table.item(row, SELECT_COL)
    assert not (item.flags() & Qt.ItemIsUserCheckable), (
        "a consensus parent must not offer a tick box")
    dlg.deleteLater()


def test_a_run_with_no_degenerate_rows_is_unchanged(qapp, tmp_path):
    """The toggle must be hidden on an ordinary run. Every probe there is its
    own site, so there is nothing to hide and nothing to collapse."""
    dlg = _dialog(_spread_run(tmp_path, n=3))
    assert dlg.show_parents_check.isVisibleTo(dlg) is False
    assert _table_probe_ids(dlg) == ["p0", "p1", "p2"]
    dlg.max_spin.setValue(3)
    dlg._on_autopick()
    assert dlg.selection_tier_counts() == {"green": 3, "orange": 0, "red": 0}
    dlg.deleteLater()


def test_a_legacy_draft_holding_a_parent_id_is_still_read_as_one_site(
        qapp, tmp_path):
    """Before this change a parent carried a tick box, so a draft saved by an
    older build can hold a parent id alongside its own children.

    The exporter already copes, because ``select_export`` replaces a selected
    parent with its variants. But the MAP has to cope too, or such a draft reopens
    solid red: parent and children are three ids on one footprint.
    """
    dlg = _dialog(_degenerate_run(tmp_path))
    dlg.max_spin.setValue(2)
    dlg._selections["G1"] = ["pA", "pA_v00", "pA_v01", "pB"]
    dlg._refresh_viz()
    assert dlg.map_widget.has_overlap is False
    assert dlg.map_widget.has_tight is False
    assert dlg.selection_tier_counts()["red"] == 0
    dlg.deleteLater()


def test_the_parent_toggle_is_reachable_by_keyboard_when_it_is_there(
        qapp, tmp_path):
    """The control above is hidden on an ordinary run, so that test alone would
    also pass if it were hidden ALWAYS. On a degenerate run it must be visible,
    accessible and a real Tab stop between Clear and the table."""
    dlg = _dialog(_degenerate_run(tmp_path))
    assert dlg.show_parents_check.isVisibleTo(dlg) is True
    assert dlg.show_parents_check.accessibleName()
    assert _next_focusable(dlg.f_clear) is dlg.show_parents_check
    assert _next_focusable(dlg.show_parents_check) is dlg.table
    dlg.deleteLater()


# ---------------------------------------------------------------------------
# The per-gene limit is a limit on SITES, not on ids
#
# auto_select counts a degenerate family once ("one family = one pick") and
# emits its children, so a correct pick of 2 sites can be 3 ids. Counting ids
# instead would make the counter read "6 / 3" on a correct pick, make lowering
# the limit cut the id list and split a family, and refuse a manual tick while
# under the limit.
# ---------------------------------------------------------------------------

def test_the_counter_counts_sites_not_variant_ids(qapp, tmp_path):
    dlg = _dialog(_degenerate_run(tmp_path))
    dlg.max_spin.setValue(2)
    dlg._on_autopick()
    assert len(dlg._selections["G1"]) == 3          # 2 variants + 1 probe
    assert "Selected for G1: 2 / 2" in dlg.counter_label.text(), (
        dlg.counter_label.text())
    dlg.deleteLater()


def test_lowering_the_limit_never_splits_a_family_or_drops_a_site_in_range(
        qapp, tmp_path):
    dlg = _dialog(_degenerate_run(tmp_path))
    dlg.max_spin.setValue(3)
    dlg._selections["G1"] = ["pB", "pA_v00", "pA_v01"]   # 2 sites
    dlg.max_spin.setValue(2)                             # still 2 sites: keep all
    assert set(dlg._selections["G1"]) == {"pB", "pA_v00", "pA_v01"}
    dlg.max_spin.setValue(1)                             # keep the first SITE
    assert dlg._selections["G1"] == ["pB"], dlg._selections["G1"]
    dlg._selections["G1"] = ["pA_v00", "pA_v01", "pB"]
    dlg.max_spin.setValue(2)
    dlg.max_spin.setValue(1)
    assert dlg._selections["G1"] == ["pA_v00", "pA_v01"], (
        "a family must be kept or dropped whole")
    dlg.deleteLater()


def test_a_manual_tick_is_allowed_while_under_the_site_limit(
        qapp, tmp_path, monkeypatch):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QMessageBox
    from flexprobe_gui.review_dialog import SELECT_COL
    refused = []
    monkeypatch.setattr(QMessageBox, "information",
                        lambda *a, **k: refused.append(a[1] if len(a) > 1 else ""))
    dlg = _dialog(_degenerate_run(tmp_path))
    dlg.max_spin.setValue(2)
    dlg._selections["G1"] = ["pA_v00", "pA_v01"]         # 1 site of 2
    dlg._populate_table()
    row = _table_probe_ids(dlg).index("pB")
    dlg.table.item(row, SELECT_COL).setCheckState(Qt.Checked)
    assert refused == [], refused
    assert "pB" in dlg._selections["G1"]
    dlg.deleteLater()


# ---------------------------------------------------------------------------
# The spacing tally under the map counts SITES, each by its OWN nearest site
#
# Giving every selected id in a gene the tier of that gene's single closest
# pair, then counting ids, would let one crowded pair turn a gene's well-spaced
# probes red in the tally while the map drew them green, and a correct pick of
# 3 degenerate sites would read "green 6" beside "Selected: 3 / 3".
# ---------------------------------------------------------------------------

def test_the_tally_tiers_each_site_by_its_own_nearest_neighbour(qapp, tmp_path):
    rows = [_row("pA", "G1", 0, 1), _row("pB", "G1", 1000, 2),
            _row("pC", "G1", 1060, 3)]          # pB-pC 10 bp apart, pA far off
    dlg = _dialog(_write_run(tmp_path / "mixed", rows,
                             params={"gene_lengths": {"G1": 3000}}))
    dlg.max_spin.setValue(3)
    dlg._selections["G1"] = ["pA", "pB", "pC"]
    dlg._refresh_viz()
    assert dlg.selection_tier_counts() == {"green": 1, "orange": 0, "red": 2}
    dlg.deleteLater()


def test_the_tally_and_the_counter_count_the_same_sites(qapp, tmp_path):
    dlg = _dialog(_degenerate_run(tmp_path))
    dlg.max_spin.setValue(2)
    dlg._on_autopick()
    assert sum(dlg.selection_tier_counts().values()) == 2   # 2 sites, 3 ids
    assert "Selected for G1: 2 / 2" in dlg.counter_label.text()
    dlg.deleteLater()


# ---------------------------------------------------------------------------
# The Spacing filter must not measure a family's variants against each other
#
# _prune_to_floor and _blocked_by_spacing_filter compared raw footprints, and a
# family's variants share one, so siblings read as -1 bp apart. Even the
# default "Allow all" floor (15 bp) then refused the second variant of a site
# ("sits -1 bp from one you have already chosen") and any change of the
# Spacing box dropped variants from the order.
# ---------------------------------------------------------------------------

def test_ticking_the_second_variant_of_a_site_is_allowed(
        qapp, tmp_path, monkeypatch):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QMessageBox
    from flexprobe_gui.review_dialog import SELECT_COL
    refused = []
    monkeypatch.setattr(QMessageBox, "information",
                        lambda *a, **k: refused.append(a[2] if len(a) > 2 else ""))
    dlg = _dialog(_degenerate_run(tmp_path))
    dlg.max_spin.setValue(2)
    dlg._selections["G1"] = ["pA_v00"]
    dlg._populate_table()
    row = _table_probe_ids(dlg).index("pA_v01")
    dlg.table.item(row, SELECT_COL).setCheckState(Qt.Checked)
    assert refused == [], refused
    assert set(dlg._selections["G1"]) == {"pA_v00", "pA_v01"}
    dlg.deleteLater()


def test_changing_the_spacing_filter_keeps_every_variant_of_a_site(qapp, tmp_path):
    dlg = _dialog(_degenerate_run(tmp_path))
    dlg.max_spin.setValue(2)
    dlg._selections["G1"] = ["pA_v00", "pA_v01", "pB"]
    for i in range(dlg.spacing_filter.count()):
        dlg.spacing_filter.setCurrentIndex(i)
        assert set(dlg._selections["G1"]) == {"pA_v00", "pA_v01", "pB"}, (
            dlg.spacing_filter.currentText(), dlg._selections["G1"])
    dlg.deleteLater()


def test_the_spacing_filter_still_drops_a_genuinely_crowded_site(qapp, tmp_path):
    """Negative control: two DIFFERENT sites 10 bp apart must still be pruned."""
    rows = [_row("pA", "G1", 0, 1), _row("pB", "G1", 60, 2)]
    dlg = _dialog(_write_run(tmp_path / "crowd", rows,
                             params={"gene_lengths": {"G1": 3000}}))
    dlg.max_spin.setValue(2)
    dlg._selections["G1"] = ["pA", "pB"]
    for i in range(dlg.spacing_filter.count()):
        if dlg.spacing_filter.itemText(i).startswith("No red"):
            dlg.spacing_filter.setCurrentIndex(i)
    assert dlg._selections["G1"] == ["pA"], dlg._selections["G1"]
    dlg.deleteLater()
