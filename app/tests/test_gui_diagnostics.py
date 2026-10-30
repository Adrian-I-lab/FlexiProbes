"""The GUI's log file and last-resort error dialog.

Every ``logger.warning`` in the engine and in ``run_manager`` must reach a log
file the user can open. In a packaged ``.app`` there is no console, so without
a handler and an error dialog an unhandled exception would take the window
with it and leave nothing behind to read or to quote.
"""
import logging
import os

import pytest

from flexprobe_gui import diagnostics
from flexprobe_gui.paths import flexprobe_home


@pytest.fixture(autouse=True)
def _no_leaked_handler():
    """Never leave a handler pointing at a tmp_path the next test deletes."""
    diagnostics.remove_file_logging()
    yield
    diagnostics.remove_file_logging()


# --- the log file -----------------------------------------------------------

def test_log_lands_under_flexprobe_home(tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path / "home"))

    path = diagnostics.install_file_logging()

    assert path == str(tmp_path / "home" / "logs" / "flexiprobes.log")
    assert os.path.isdir(str(tmp_path / "home" / "logs"))


@pytest.mark.parametrize("tree", ["flexprobe", "flexprobe_gui"])
def test_both_logger_trees_reach_the_file(tmp_path, monkeypatch, tree):
    """The engine runs inside the GUI's worker thread, so half the messages
    worth having come from ``flexprobe`` and half from ``flexprobe_gui``."""
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path / "home"))
    path = diagnostics.install_file_logging()

    logging.getLogger(f"{tree}.somewhere.deep").warning("mud %s", "flap")
    for h in logging.getLogger(tree).handlers:
        h.flush()

    assert "mud flap" in open(path).read()


def test_the_home_is_resolved_at_call_time_not_at_import(tmp_path, monkeypatch):
    """The suite redirects FLEXPROBE_HOME per test. A handler that bound the
    path at import would spend the whole run writing into a real ~/.flexprobe."""
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path / "first"))
    first = diagnostics.install_file_logging()
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path / "second"))
    second = diagnostics.install_file_logging()

    assert first != second
    assert second.startswith(str(tmp_path / "second"))

    logging.getLogger("flexprobe_gui.t").warning("after the move")
    for h in logging.getLogger("flexprobe_gui").handlers:
        h.flush()
    assert "after the move" in open(second).read()
    assert "after the move" not in open(first).read()


def test_installing_twice_does_not_double_the_handlers(tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path / "home"))
    diagnostics.install_file_logging()
    diagnostics.install_file_logging()

    for tree in diagnostics.LOGGER_TREES:
        ours = [h for h in logging.getLogger(tree).handlers
                if isinstance(h, logging.handlers.RotatingFileHandler)]
        assert len(ours) == 1


def test_the_log_rotates_rather_than_growing_without_bound(tmp_path,
                                                           monkeypatch):
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path / "home"))
    diagnostics.install_file_logging()
    handler = [h for h in logging.getLogger("flexprobe_gui").handlers
               if isinstance(h, logging.handlers.RotatingFileHandler)][0]
    assert handler.maxBytes > 0
    assert handler.backupCount > 0


def test_an_unwritable_home_does_not_stop_the_app(tmp_path, monkeypatch):
    # A file where the home directory should be: nothing under it can be made.
    blocked = tmp_path / "not-a-dir"
    blocked.write_text("")
    monkeypatch.setenv("FLEXPROBE_HOME", str(blocked))

    assert diagnostics.install_file_logging() is None
    logging.getLogger("flexprobe_gui.t").warning("still runs")   # must not raise


def test_flexprobe_home_helper_is_what_the_stores_use(tmp_path, monkeypatch):
    from flexprobe_gui.run_history import history_path
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path / "home"))
    assert history_path() == os.path.join(flexprobe_home(), "runs.json")
    assert flexprobe_home() == str(tmp_path / "home")


# --- the excepthook ---------------------------------------------------------

@pytest.fixture
def fresh_excepthook(monkeypatch):
    """Install the hook onto a scratch sys.excepthook and undo it after."""
    import sys
    monkeypatch.setattr(diagnostics, "_excepthook_installed", False)
    monkeypatch.setattr(sys, "excepthook", sys.__excepthook__)
    yield
    monkeypatch.setattr(diagnostics, "_excepthook_installed", False)


def _raise_through(hook):
    try:
        raise ValueError("a distinctive bang")
    except ValueError:
        import sys
        hook(*sys.exc_info())


def test_excepthook_logs_the_traceback(tmp_path, monkeypatch, caplog,
                                       fresh_excepthook):
    import sys
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path / "home"))
    path = diagnostics.install_file_logging()
    diagnostics.install_excepthook(report=lambda *a: True)

    with caplog.at_level(logging.CRITICAL, logger="flexprobe_gui.diagnostics"):
        _raise_through(sys.excepthook)
    for h in logging.getLogger("flexprobe_gui").handlers:
        h.flush()

    written = open(path).read()
    assert "a distinctive bang" in written
    assert "Traceback" in written


def test_excepthook_offers_the_user_something_to_copy(monkeypatch,
                                                      fresh_excepthook):
    import sys
    seen = []
    diagnostics.install_excepthook(report=lambda *a: seen.append(a) or True)

    _raise_through(sys.excepthook)

    assert len(seen) == 1
    title, text, detail = seen[0]
    assert "ValueError: a distinctive bang" in text
    # The detail is the whole traceback, which is the part a bug report needs.
    assert "Traceback" in detail and "a distinctive bang" in detail
    assert title


def test_excepthook_survives_a_reporter_that_itself_fails(monkeypatch,
                                                          fresh_excepthook):
    """It runs while Qt may already be tearing down. It reports a failure; it
    must never become one."""
    import sys

    def _boom(*_a):
        raise RuntimeError("QApplication is gone")

    diagnostics.install_excepthook(report=_boom)
    _raise_through(sys.excepthook)          # must not raise


def test_excepthook_leaves_keyboard_interrupt_alone(monkeypatch,
                                                    fresh_excepthook):
    import sys
    seen = []
    diagnostics.install_excepthook(report=lambda *a: seen.append(a) or True)
    try:
        raise KeyboardInterrupt
    except KeyboardInterrupt:
        sys.excepthook(*sys.exc_info())
    assert seen == []


def test_excepthook_defaults_to_the_copyable_dialog(fresh_excepthook):
    import sys
    diagnostics.install_excepthook()
    assert sys.excepthook is not sys.__excepthook__
    # Read the default out of the closure rather than calling it: a QApplication
    # left behind by another Qt test in the same process would turn a call into
    # a modal dialog with nobody to click it, and the suite would hang.
    cells = [c.cell_contents for c in (sys.excepthook.__closure__ or ())]
    assert diagnostics.show_error_dialog in cells


def test_installing_the_excepthook_twice_does_not_nest_it(fresh_excepthook):
    import sys
    diagnostics.install_excepthook(report=lambda *a: True)
    once = sys.excepthook
    diagnostics.install_excepthook(report=lambda *a: True)
    assert sys.excepthook is once
