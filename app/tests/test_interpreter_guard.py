"""The interpreter guard's own tests.

The guard exists to stop a Qt-less interpreter aborting the run. It is worth
testing because its failure mode is SILENCE: the first version checked
``sys.modules`` for PySide6, which is empty at collection time because the GUI
modules import Qt lazily, so it never fired and the crash it was written to
prevent still happened. A guard that does nothing looks exactly like a guard
that works.

These tests drive the decision logic with fakes rather than a real broken
interpreter, so they run identically everywhere.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
import _interpreter_guard as guard  # noqa: E402


class FakeItem:
    def __init__(self, path):
        self.fspath = path
        self.markers = []

    def add_marker(self, marker):
        self.markers.append(marker)


@pytest.fixture
def qt_test(tmp_path):
    p = tmp_path / "test_gui_thing.py"
    p.write_text("from PySide6.QtWidgets import QWidget\ndef test_x(): pass\n")
    return p


@pytest.fixture
def core_test(tmp_path):
    p = tmp_path / "test_core_thing.py"
    p.write_text("def test_y(): pass\n")
    return p


# --- source scan -----------------------------------------------------------

def test_qt_dependent_detects_pyside_import(qt_test):
    assert guard.qt_dependent(qt_test) is True


def test_qt_dependent_ignores_a_mere_mention(tmp_path):
    """This file is the reason: it discusses Qt without ever starting it."""
    p = tmp_path / "test_talks_about_it.py"
    p.write_text("# PySide6 is discussed here but never imported\n"
                 "MESSAGE = 'install PySide6 to run the GUI'\n")
    assert guard.qt_dependent(p) is False


def test_qt_dependent_detects_importorskip(tmp_path):
    p = tmp_path / "test_skipper.py"
    p.write_text("import pytest\npytest.importorskip(\"PySide6\")\n")
    assert guard.qt_dependent(p) is True


def test_qt_dependent_detects_gui_package(tmp_path):
    p = tmp_path / "test_z.py"
    p.write_text("from flexprobe_gui.app import main\n")
    assert guard.qt_dependent(p) is True


def test_qt_dependent_false_for_core_module(core_test):
    assert guard.qt_dependent(core_test) is False


def test_qt_dependent_survives_an_unreadable_path(tmp_path):
    assert guard.qt_dependent(tmp_path / "nope.py") is False


# --- venv discovery --------------------------------------------------------

def test_find_project_venv_locates_a_venv_above(tmp_path):
    venv = tmp_path / "app" / ".venv" / "bin"
    venv.mkdir(parents=True)
    (venv / "python").write_text("#!/bin/sh\n")
    found = guard.find_project_venv(tmp_path / "app" / "tests")
    assert found == venv / "python"


def test_find_project_venv_returns_none_when_absent(tmp_path):
    # The normal state in a fresh clone or a git worktree: .venv is gitignored.
    assert guard.find_project_venv(tmp_path) is None


# --- decision logic --------------------------------------------------------

def _check(items, tmp_path, monkeypatch, *, probes, ci=False, skip_env=False):
    """Run guard.check with qt_probe stubbed. ``probes`` maps executable -> error."""
    monkeypatch.setattr(guard, "qt_probe", lambda exe, **kw: probes.get(exe, None))
    monkeypatch.delenv("CI", raising=False)
    monkeypatch.delenv(guard.SKIP_ENV, raising=False)
    if ci:
        monkeypatch.setenv("CI", "true")
    if skip_env:
        monkeypatch.setenv(guard.SKIP_ENV, "1")

    class Boom(Exception):
        pass

    guard.check(items, tmp_path, Boom)
    return Boom


def test_core_only_run_is_never_blocked(core_test, tmp_path, monkeypatch):
    """No Qt tests collected: no probe, no error, no skip -- even if Qt is broken."""
    called = []
    monkeypatch.setattr(guard, "qt_probe",
                        lambda exe, **kw: called.append(exe) or "broken")
    item = FakeItem(core_test)
    guard.check([item], tmp_path, RuntimeError)
    assert called == [], "probed despite no Qt tests being collected"
    assert item.markers == []


def test_working_qt_leaves_everything_alone(qt_test, tmp_path, monkeypatch):
    item = FakeItem(qt_test)
    _check([item], tmp_path, monkeypatch, probes={})
    assert item.markers == []


def test_broken_qt_with_a_working_venv_raises(qt_test, tmp_path, monkeypatch):
    venv = tmp_path / ".venv" / "bin"
    venv.mkdir(parents=True)
    (venv / "python").write_text("#!/bin/sh\n")
    item = FakeItem(qt_test)
    with pytest.raises(Exception) as exc:
        _check([item], tmp_path, monkeypatch,
               probes={sys.executable: "no cocoa plugin"})
    assert "no cocoa plugin" in str(exc.value)
    assert "-m pytest" in str(exc.value), "message must give a runnable command"
    assert item.markers == [], "should fail fast, not skip"


def test_broken_qt_without_a_venv_skips_rather_than_aborting(
        qt_test, core_test, tmp_path, monkeypatch):
    qt_item, core_item = FakeItem(qt_test), FakeItem(core_test)
    _check([qt_item, core_item], tmp_path, monkeypatch,
           probes={sys.executable: "no cocoa plugin"})
    assert len(qt_item.markers) == 1, "Qt test should be skipped"
    assert "no cocoa plugin" in qt_item.markers[0].kwargs["reason"]
    assert core_item.markers == [], "core test must still run"


def test_broken_qt_in_ci_raises_instead_of_skipping(qt_test, tmp_path, monkeypatch):
    """A green CI build that quietly dropped the GUI tests is worse than a red one."""
    item = FakeItem(qt_test)
    with pytest.raises(Exception) as exc:
        _check([item], tmp_path, monkeypatch,
               probes={sys.executable: "libEGL.so.1 missing"}, ci=True)
    assert "libEGL.so.1 missing" in str(exc.value)
    assert "ci.yml" in str(exc.value), "message must point at the failing step"
    assert item.markers == []


def test_escape_hatch_skips_without_probing(qt_test, tmp_path, monkeypatch):
    called = []
    monkeypatch.setattr(guard, "qt_probe",
                        lambda exe, **kw: called.append(exe) or "broken")
    monkeypatch.setenv(guard.SKIP_ENV, "1")
    monkeypatch.delenv("CI", raising=False)
    item = FakeItem(qt_test)
    guard.check([item], tmp_path, RuntimeError)
    assert called == [], "escape hatch should not pay for the probe"
    assert len(item.markers) == 1


# --- the probe itself ------------------------------------------------------

def test_qt_probe_reports_a_missing_interpreter():
    err = guard.qt_probe("/definitely/not/a/python")
    assert err is not None and "could not run" in err
