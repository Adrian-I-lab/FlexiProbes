"""GUI runtime tests: column registry + the run_manager worker thread.

Skipped automatically when PySide6 isn't installed (e.g. the minimal CI/sandbox
env). Runs headless via the offscreen Qt platform.
"""
import importlib.util
import os
import time

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


# --- column registry --------------------------------------------------------

def test_column_registry_select_and_headers(qapp):
    from flexprobe_gui.review_dialog import _BASE_COLS, SELECT_COL
    labels = [c.label for c in _BASE_COLS]
    assert _BASE_COLS[SELECT_COL].kind == "select"
    # the three GC% columns are present and in order
    assert labels[labels.index("GC% 50bp"):labels.index("GC% 50bp") + 3] == \
        ["GC% 50bp", "GC% LHS", "GC% RHS"]


def test_make_cell_dispatch(qapp):
    from PySide6.QtCore import Qt
    from flexprobe_gui.review_dialog import ReviewDialog, _Col
    row = {"gc_lhs_half": "0.48", "composite_score": "0.91", "probe_id": "p7",
           "specificity_status": "not_screened"}
    gc = ReviewDialog._make_cell(_Col("GC% LHS", 68, "gc", "gc_lhs_half"),
                                 row, "p7", set())
    assert abs(gc.data(Qt.DisplayRole) - 48.0) < 1e-9      # fraction -> percent
    num = ReviewDialog._make_cell(_Col("Score", 80, "num", "composite_score"),
                                  row, "p7", set())
    assert abs(num.data(Qt.DisplayRole) - 0.91) < 1e-9
    txt = ReviewDialog._make_cell(_Col("Probe ID", 200, "text", "probe_id"),
                                  row, "p7", set())
    assert txt.text() == "p7"
    sel = ReviewDialog._make_cell(_Col("Select", 56, "select"), row, "p7", {"p7"})
    assert sel.checkState() == Qt.Checked


def test_tm_filter_uses_displayed_half_columns(qapp):
    # The "Tm >=" floor must act on the per-half Tms the table displays
    # (tm_target_lhs/rhs), requiring BOTH to clear it, not on tm_full_50,
    # which is never shown.
    from flexprobe_gui.review_dialog import ReviewDialog

    class _Spin:
        def __init__(self, v):
            self._v = v

        def value(self):
            return self._v

    class _Text:
        def text(self):
            return ""

    fake = type("F", (), {})()
    fake.f_gc_min = _Spin(0.0)
    fake.f_gc_max = _Spin(100.0)
    fake.f_tm_min = _Spin(50.0)
    fake.f_score_min = _Spin(0.0)
    fake.f_text = _Text()

    base = {"gc_content_full": "55", "composite_score": "0.5", "probe_id": "p1"}
    # A weak displayed half (45) must fail even though the invisible full-50
    # Tm (70) is comfortably above the floor.
    low_half = {**base, "tm_full_50": "70",
                "tm_target_lhs": "45", "tm_target_rhs": "60"}
    assert ReviewDialog._passes_filter(fake, low_half) is False
    # Both halves above the floor pass, even when full-50 is below it.
    ok = {**base, "tm_full_50": "40",
          "tm_target_lhs": "55", "tm_target_rhs": "58"}
    assert ReviewDialog._passes_filter(fake, ok) is True


# --- run_manager worker thread ---------------------------------------------

def test_run_manager_completes_one_shot(qapp, tmp_path):
    from flexprobe_gui.run_history import RunHistory, RunRecord, TERMINAL_STATUSES
    from flexprobe_gui.run_manager import RunManager
    from flexprobe.models import ScreenResources, ModalityDefinition, RunConfig

    fa = tmp_path / "toy.fasta"
    fa.write_text(
        ">toy|gene|ACC1\n"
        "GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGCATGCATGCATGCA"
        "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG\n"
    )
    out = tmp_path / "out"
    hist = RunHistory(path=str(tmp_path / "runs.json"))
    mgr = RunManager(hist)
    rec = RunRecord(name="smoke", output_dir=str(out), fasta_path=str(fa))
    cfg = RunConfig(fasta_path=str(fa), output_dir=str(out),
                    modality=ModalityDefinition(name="flex", version="v2"),
                    screen=ScreenResources(skip_screen=True))
    cfg.filters.require_tn_junction = True

    mgr.submit(rec, cfg)

    deadline = time.time() + 60
    while time.time() < deadline:
        qapp.processEvents()
        cur = hist.get(rec.run_id)
        if cur and cur.status in TERMINAL_STATUSES:
            break
        time.sleep(0.02)

    final = hist.get(rec.run_id)
    assert final is not None and final.status == "complete", \
        f"status={final.status if final else None} error={getattr(final,'error',None)}"
    assert (out / "results.csv").exists()


# --- bulk select / move / delete (home window) ------------------------------

def test_home_window_bulk_select_move(qapp, tmp_path, monkeypatch):
    # Redirect the on-disk stores into tmp so we never touch the real
    # ~/.flexprobe while the widget builds its RunHistory/ProjectHistory.
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path))
    from PySide6.QtWidgets import QInputDialog
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.run_history import RunRecord

    win = HomeWindow()
    try:
        proj = win.projects.create("Target")
        win._refresh_project_filter()
        for nm in ("r1", "r2"):
            win.manager.history.append(RunRecord(name=nm, status="complete"))
        win._refresh_completed()
        assert win.table.rowCount() == 2
        # Checkbox column present, nothing checked yet.
        assert win._checked_run_ids() == []
        assert not win.bulk_move_btn.isEnabled()
        assert not win.bulk_delete_btn.isEnabled()

        # Select-all checks every row and enables the actions.
        win._on_select_all_runs(True)
        assert len(win._checked_run_ids()) == 2
        assert win.bulk_move_btn.isEnabled() and win.bulk_delete_btn.isEnabled()
        assert win.select_all_chk.isChecked()

        # Move both checked runs into the project (stub the picker dialog).
        monkeypatch.setattr(
            QInputDialog, "getItem",
            staticmethod(lambda *a, **k: ("Target", True)),
        )
        win._on_bulk_move_runs()
        assert [r.project_id for r in win.manager.history.list_completed()] == \
            [proj.project_id, proj.project_id]
        # Rows rebuilt unchecked after the action.
        assert win._checked_run_ids() == []
    finally:
        win.close()


def test_home_window_bulk_delete(qapp, tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path))
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QMessageBox
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.run_history import RunRecord

    win = HomeWindow()
    try:
        for nm in ("a", "b", "c"):
            win.manager.history.append(RunRecord(name=nm, status="complete"))
        win._refresh_completed()
        win.table.item(0, 0).setCheckState(Qt.Checked)
        assert len(win._checked_run_ids()) == 1
        # Confirm the delete (no output folders → no checkbox in the dialog).
        monkeypatch.setattr(QMessageBox, "exec", lambda self: QMessageBox.Yes)
        win._on_bulk_delete_runs()
        assert len(win.manager.history.list_completed()) == 2
    finally:
        win.close()


# --- bulk select / delete (manage-projects dialog) --------------------------

def test_manage_projects_bulk_delete(qapp, tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path))
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QMessageBox
    from flexprobe_gui.manage_projects_dialog import ManageProjectsDialog
    from flexprobe_gui.project_history import ProjectHistory
    from flexprobe_gui.run_history import RunHistory

    ph = ProjectHistory(str(tmp_path / "projects.json"))
    for nm in ("Alpha", "Beta", "Gamma"):
        ph.create(nm)
    dlg = ManageProjectsDialog(ph, RunHistory(str(tmp_path / "runs.json")))
    try:
        assert dlg.table.rowCount() == 3
        assert not dlg.delete_btn.isEnabled()

        dlg._on_select_all(True)
        assert len(dlg._checked_ids()) == 3
        assert dlg.delete_btn.isEnabled()
        # Rename/describe need exactly one checked target.
        assert not dlg.rename_btn.isEnabled()
        dlg.table.item(1, 0).setCheckState(Qt.Unchecked)
        dlg.table.item(2, 0).setCheckState(Qt.Unchecked)
        assert dlg.rename_btn.isEnabled()

        # Re-check all, then confirm a bulk delete.
        dlg._on_select_all(True)
        monkeypatch.setattr(
            QMessageBox, "question",
            staticmethod(lambda *a, **k: QMessageBox.Yes),
        )
        dlg._on_delete()
        assert dlg.table.rowCount() == 0
        assert ph.list() == []
    finally:
        dlg.close()
