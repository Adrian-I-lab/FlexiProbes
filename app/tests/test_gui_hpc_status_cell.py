"""An imported-HPC run's Status cell must not paint its text twice.

``setCellWidget`` is an OVERLAY, not a replacement: the view keeps painting the
underlying item's ``DisplayRole`` and the widget is drawn on top. The HPC badge
widget is a plain ``QWidget`` with no background, so wherever the item also
carried the status text, both strings rendered — the item's at the cell's own
text origin, the widget's label a few pixels right of it, reading as bold,
doubled, smudged text next to the blue "HPC" pill.

Only imported runs take that branch, so locally-run designs never show the
problem.

The guard asserts the invariant that actually matters and is cheap to keep:
whenever a Status cell carries an overlay widget, the item beneath it must be
blank. It is deliberately not a screenshot test.
"""
import importlib.util
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 not installed"
)

STATUS_COL = 5


def _rows_by_name(win):
    return {win.table.item(r, 2).text(): r for r in range(win.table.rowCount())}


def test_imported_run_status_is_not_painted_twice(qapp, tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path))
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.run_history import RunRecord

    win = HomeWindow()
    try:
        win.manager.history.append(RunRecord(
            name="local", status="ready_for_review", stats={}))
        win.manager.history.append(RunRecord(
            name="imported", status="ready_for_review",
            stats={"imported": True}))
        win._refresh_completed()

        rows = _rows_by_name(win)
        local, imported = rows["local"], rows["imported"]

        # The local run is the control: plain item, no overlay, text visible.
        assert win.table.cellWidget(local, STATUS_COL) is None
        assert win.table.item(local, STATUS_COL).text() == "ready for review"

        # The imported run gets the badge widget, and the item under it must be
        # blank or the same words paint twice.
        assert win.table.cellWidget(imported, STATUS_COL) is not None
        assert win.table.item(imported, STATUS_COL).text() == ""

        # The status must still be legible — it moved into the widget, it did
        # not disappear. Guards against "fixing" this by dropping the text.
        from PySide6.QtWidgets import QLabel
        labels = win.table.cellWidget(imported, STATUS_COL).findChildren(QLabel)
        assert "ready for review" in [lb.text() for lb in labels]
        assert "HPC" in [lb.text() for lb in labels]
    finally:
        win.close()


def test_imported_failed_run_keeps_its_error_tooltip(qapp, tmp_path, monkeypatch):
    """Blanking the item must not take the failure reason with it.

    The tooltip lived on the item. An overlay widget swallows hover events, so
    a tooltip left only on the item beneath it is unreachable for exactly the
    runs that need it most.
    """
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path))
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.run_history import RunRecord

    win = HomeWindow()
    try:
        win.manager.history.append(RunRecord(
            name="boom", status="failed", error="Traceback: kmer engine died",
            stats={"imported": True}))
        win._refresh_completed()

        widget = win.table.cellWidget(0, STATUS_COL)
        assert widget is not None
        assert "kmer engine died" in widget.toolTip()
    finally:
        win.close()


def test_user_visible_branding_is_capitalised(qapp, tmp_path, monkeypatch):
    """The product name reads FlexiProbes wherever a user can see it.

    Asserted against the live widgets, not by grepping the source, because a
    source grep would also match the lowercase strings that MUST stay
    lowercase: the ``flexiprobes`` console-script name printed in the Tools
    dialog, ``flexiprobes.log``, ``.flexiprobes_DONE`` and
    ``flexiprobes_output.tar.gz``. Those are commands and filenames, and
    capitalising them would tell people to type something that does not exist.
    """
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path))
    from PySide6.QtWidgets import QLabel
    from flexprobe_gui.home_window import HomeWindow

    win = HomeWindow()
    try:
        assert win.windowTitle().startswith("FlexiProbes")
        wordmarks = [lb.text() for lb in win.findChildren(QLabel)
                     if lb.text().lower() == "flexiprobes"]
        assert wordmarks, "home-page wordmark not found"
        assert all(w == "FlexiProbes" for w in wordmarks), wordmarks
    finally:
        win.close()
