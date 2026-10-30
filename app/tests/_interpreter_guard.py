"""Fail fast when the interpreter running the tests cannot start Qt.

WHY THIS EXISTS. ``python3`` on a developer machine is very often not the
interpreter this project is installed into. On the machine this was written for
it resolves to Anaconda, whose PySide6 cannot find its own ``cocoa`` platform
plugin; ``QApplication()`` then calls ``abort()``. pytest does not survive that:
the run dies with ``Fatal Python error: Aborted`` roughly 40 tests in, with no
failure report and no indication of which test or why.

That failure mode is actively misleading. It looks like the test suite is
broken, it lands in the middle of an unrelated file, and it reproduces on a
clean checkout, so the natural reading is "main is red". Anyone who runs a bare
``python3 -m pytest`` from that shell hits it.

The check is BEHAVIOURAL, not a heuristic. It does not sniff for Anaconda, a
Python version, or a path pattern, because the thing that matters is whether Qt
actually starts, and that can break for reasons that have nothing to do with
conda. It starts a real QApplication in a subprocess and looks at what happens.

It costs one subprocess (~0.3 s) and only when GUI tests were actually
collected. A core-only run never pays for it and is never blocked by it, since
the core suite passes under any interpreter.
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path
from typing import List, Optional

# Constructing a QWidget is the cheapest thing that exercises the full stack:
# the platform plugin has to be found and initialised before a widget exists.
# Importing PySide6 alone is NOT enough -- the import succeeds under Anaconda
# and the abort comes later, which is exactly why collection succeeds and the
# run dies mid-flight.
_QT_PROBE = (
    "from PySide6.QtWidgets import QApplication, QWidget\n"
    "QApplication([])\n"
    "QWidget()\n"
)

SKIP_ENV = "FLEXIPROBES_SKIP_QT_PRECHECK"


def qt_probe(executable: str, timeout: float = 120.0) -> Optional[str]:
    """Return None if ``executable`` can start Qt, else the reason it cannot.

    A crash (SIGABRT) is the expected failure here, not an exception, so the
    return code is what is checked -- ``subprocess`` reports a fatal signal as a
    negative code and raises nothing.
    """
    try:
        proc = subprocess.run(
            [executable, "-c", _QT_PROBE],
            capture_output=True, text=True, timeout=timeout,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return f"could not run {executable}: {exc}"
    if proc.returncode == 0:
        return None
    detail = (proc.stderr or proc.stdout or "").strip()
    if not detail:
        detail = f"exited with code {proc.returncode} and no output"
    # Qt's own diagnostic is the useful part and it is usually the first line;
    # the rest is a Python traceback or a signal dump that adds nothing.
    return detail.splitlines()[0].strip()


def find_project_venv(start: Path) -> Optional[Path]:
    """The project virtualenv interpreter, if one has been created.

    Looked for beside the package rather than anywhere on ``PATH``: the point is
    to name the interpreter this checkout is set up for, not whichever venv the
    caller happens to have active. Returns None in a fresh clone or a git
    worktree, where ``.venv`` is gitignored and so does not exist yet -- that is
    a normal state, not an error.
    """
    for base in (start, *start.parents):
        for rel in ("bin/python", "Scripts/python.exe"):
            cand = base / ".venv" / rel
            if cand.is_file():
                return cand
    return None


def _same_interpreter(a: Path, b: Path) -> bool:
    try:
        return a.resolve() == b.resolve()
    except OSError:
        return str(a) == str(b)


def wrong_interpreter_message(venv: Path, running: str, qt_error: str) -> str:
    """The message shown when a working alternative interpreter exists."""
    try:
        shown = venv.relative_to(Path.cwd())
    except ValueError:
        shown = venv
    return (
        "\n"
        "This interpreter cannot start Qt, so the GUI tests would abort the\n"
        "whole run instead of failing. Stopping before that happens.\n"
        "\n"
        f"  running   {running}\n"
        f"            (Python {sys.version.split()[0]})\n"
        f"  Qt says   {qt_error}\n"
        "\n"
        "The project virtualenv CAN start Qt. Verified just now. Use it:\n"
        "\n"
        f"    {shown} -m pytest\n"
        "\n"
        "The core suite does not need Qt and passes under any interpreter, so\n"
        f"to run it here anyway set {SKIP_ENV}=1 -- the GUI tests will be\n"
        "skipped rather than crashing the run.\n"
    )


def broken_qt_skip_reason(qt_error: str) -> str:
    """The message used when there is no better interpreter to point at."""
    return (
        f"Qt cannot start under this interpreter ({qt_error}); "
        "GUI tests skipped so they cannot abort the run. "
        "Create the project virtualenv to run them."
    )


def ci_qt_broken_message(qt_error: str) -> str:
    """Shown when Qt is broken in CI, where skipping would hide a real defect."""
    return (
        "\n"
        "Qt cannot start, and this is CI.\n"
        "\n"
        f"  Qt says   {qt_error}\n"
        "\n"
        "The GUI tests are NOT being skipped, because a green build that\n"
        "quietly ran 100+ fewer tests is worse than a red one. CI installs the\n"
        "Qt system libraries and sets QT_QPA_PLATFORM=offscreen on purpose, so\n"
        "this means that provisioning has broken. See the 'Install Qt system\n"
        "libraries' step in .github/workflows/ci.yml.\n"
    )


# An IMPORT of Qt, not a mention of it. A bare substring search also matched
# this guard's own tests and its docstrings, which discuss PySide6 at length
# without ever starting it -- and skipping the guard's tests on a broken
# interpreter would remove the check precisely where it earns its place.
_QT_IMPORT = re.compile(
    r"^\s*(?:from|import)\s+(?:PySide6|flexprobe_gui)\b"
    r"|importorskip\(\s*[\"\']PySide6",
    re.MULTILINE,
)


def qt_dependent(path: Path) -> bool:
    """Whether a test module drives Qt, by what it imports.

    Read from source rather than from the imported module, because the import
    is exactly what must not be trusted: PySide6 imports cleanly under a broken
    interpreter and only blows up later.
    """
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False
    return _QT_IMPORT.search(text) is not None


def check(items: List, conftest_dir: Path, raise_usage_error) -> None:
    """Guard entry point. Called from ``conftest.pytest_collection_modifyitems``.

    ``raise_usage_error`` is the EXCEPTION TYPE to raise (pytest.UsageError in
    practice), passed in rather than imported so this module stays importable,
    and unit-testable, without pytest.
    """
    # Which collected tests drive Qt. Determined by SOURCE SCAN, not by looking
    # for PySide6 in sys.modules: the GUI modules import it lazily, inside
    # fixtures and test bodies, so at collection time sys.modules says "no Qt
    # here" for a run that is about to abort in exactly those tests. That check
    # was the first version of this guard and it silently never fired.
    qt_items = _qt_items(items)
    if not qt_items:
        return                        # core-only run: costs nothing, blocks nothing

    if os.environ.get(SKIP_ENV):
        _mark_skipped(qt_items, f"{SKIP_ENV} is set; GUI tests skipped.")
        return

    qt_error = qt_probe(sys.executable)
    if qt_error is None:
        return                                  # Qt works here; nothing to do.

    venv = find_project_venv(conftest_dir)
    if (venv is not None
            and not _same_interpreter(venv, Path(sys.executable))
            and qt_probe(str(venv)) is None):
        raise raise_usage_error(
            wrong_interpreter_message(venv, sys.executable, qt_error))

    # In CI a silent skip would be worse than the crash: the build would go
    # green having run 106 fewer tests than it reports covering, and the Qt
    # provisioning could rot unnoticed. CI provisions Qt deliberately
    # (.github/workflows/ci.yml installs the system libraries and sets
    # QT_QPA_PLATFORM=offscreen), so Qt failing there is a real defect and has
    # to be loud.
    if os.environ.get("CI"):
        raise raise_usage_error(ci_qt_broken_message(qt_error))

    # Locally, Qt is broken and nothing better is available. Skipping beats
    # aborting: the ~740 core tests still run and still mean something.
    _mark_skipped(qt_items, broken_qt_skip_reason(qt_error))


def _qt_items(items: List) -> List:
    """The collected items whose module drives Qt."""
    cache: dict = {}
    out = []
    for item in items:
        path = Path(str(getattr(item, "fspath", "") or ""))
        if not path.name:
            continue
        if path not in cache:
            cache[path] = qt_dependent(path)
        if cache[path]:
            out.append(item)
    return out


def _mark_skipped(qt_items: List, reason: str) -> None:
    import pytest

    marker = pytest.mark.skip(reason=reason)
    for item in qt_items:
        item.add_marker(marker)
