"""The home table's per-run buttons must be reachable without scrolling.

Review, Filters, Open and Re-run sit before the Output column, which holds a
full folder path. Placed after it, on a 1280 px window the path would push
every button off the right edge, so the app's main actions would need a
sideways scroll to find.
"""
import importlib.util
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 not installed"
)

LONG_PATH = "/path/to/" + "very_long_folder_name_" * 8 + "run"


@pytest.fixture()
def win(qapp, tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path))
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.run_history import RunRecord
    w = HomeWindow()
    w.manager.history.append(RunRecord(
        name="r1", status="complete", output_dir=LONG_PATH, stats={},
        started_at="2026-09-22T10:43:47.123456"))
    w._refresh_completed()
    w.resize(1280, 820)
    w.show()
    qapp.processEvents()
    yield w
    w.close()


def _header(win, col):
    return win.table.horizontalHeaderItem(col).text()


def test_actions_column_comes_before_the_output_path(win):
    from flexprobe_gui._home_completed import COL_ACTIONS, COL_OUTPUT
    assert COL_ACTIONS < COL_OUTPUT
    assert _header(win, COL_OUTPUT) == "Output"
    assert win.table.cellWidget(0, COL_ACTIONS) is not None


def test_action_buttons_are_visible_at_1280_px_without_scrolling(win):
    from flexprobe_gui._home_completed import COL_ACTIONS
    t = win.table
    right_edge = t.columnViewportPosition(COL_ACTIONS) + t.columnWidth(COL_ACTIONS)
    assert t.horizontalScrollBar().value() == 0
    assert right_edge <= t.viewport().width(), (right_edge, t.viewport().width())


def test_the_full_output_path_is_still_available_as_a_tooltip(win):
    from flexprobe_gui._home_completed import COL_OUTPUT
    # The folder does not exist here, so the "no longer exists" note is added
    # too. The path must survive it.
    assert LONG_PATH in win.table.item(0, COL_OUTPUT).toolTip()


def test_started_is_shown_as_a_readable_time_not_a_raw_timestamp(win):
    cell = win.table.item(0, 1)
    assert cell.text() == "2026-09-22 10:43"
    assert cell.toolTip() == "2026-09-22T10:43:47.123456"


def test_the_table_has_nothing_to_scroll_sideways_at_1280_px(
        qapp, tmp_path, monkeypatch):
    """Output takes exactly the width that is left, so a long path can never
    push the table past the window. The first fix capped Output at 360 px,
    which still overran on the screenshot pass.

    Built the way the screenshot was: the history is on disk BEFORE the window
    opens, the run folders EXIST (which adds the Review, Filters and Open
    buttons), and there are several rows of mixed status.
    """
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path / "home"))
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.run_history import RunHistory, RunRecord
    hist = RunHistory()
    for i, status in enumerate(["ready_for_review", "complete",
                                "ready_for_review", "failed"]):
        d = tmp_path / ("a_long_output_folder_name_for_run_" * 3 + str(i))
        d.mkdir(parents=True)
        hist.append(RunRecord(name=f"R{i:02d}_OROV_L_unique_BeAn", status=status,
                              output_dir=str(d) if status != "failed" else "",
                              modality_label="Flex v1 (Chromium)", stats={},
                              error="boom" if status == "failed" else None))
    hist.flush()
    w = HomeWindow()
    try:
        w.resize(1280, 820)
        w.show()
        qapp.processEvents()
        assert w.table.horizontalScrollBar().maximum() == 0, (
            w.table.horizontalScrollBar().maximum())
    finally:
        w.close()
