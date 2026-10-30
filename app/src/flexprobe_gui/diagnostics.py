# flexiprobes — probe design for 10x Genomics Flex v1 / Flex v2 / Visium.
# Copyright (C) 2026 Adrian Ilich and contributors
#
# This file is part of flexiprobes.
#
# flexiprobes is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# flexiprobes is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""The GUI's diagnostics: a log file on disk, and a visible last-resort dialog.

Before this module existed the GUI had neither. Every ``logger.warning`` in the
engine and in ``run_manager`` — including the one that fires when a run's order
CSVs fail to write — went to a logger with no handler at all, which is to say
nowhere. In a packaged ``.app`` there is no console either, so an unhandled
exception took the window with it and left the user with nothing to report and
nothing to read. A comment in ``run_manager`` promised a "diagnostics log";
this is that file finally existing.

Two halves:

``install_file_logging``
    A size-capped rotating file under ``$FLEXPROBE_HOME/logs/``, attached to
    both the ``flexprobe`` (engine) and ``flexprobe_gui`` logger trees. Rotating
    rather than plain: this is a desktop app that may never be restarted for
    weeks, and an unbounded log in a user's home directory is its own bug.
``install_excepthook``
    Catches what nothing else caught, writes the traceback to that file, and
    puts it on screen in a form the user can SELECT AND COPY. A crash the user
    cannot quote is a crash that cannot be reported.
"""
from __future__ import annotations

import logging
import logging.handlers
import os
import sys
import traceback
from types import TracebackType
from typing import Callable, Optional, Type

from .paths import flexprobe_home

logger = logging.getLogger("flexprobe_gui.diagnostics")

#: Both trees, because half the interesting messages come from the engine
#: running inside the GUI's worker thread and half from the GUI itself.
LOGGER_TREES = ("flexprobe", "flexprobe_gui")

LOG_FILE_NAME = "flexiprobes.log"
#: 1 MB x 5 keeps roughly the last few thousand messages — enough to cover the
#: run that went wrong — while capping the app's footprint in a home directory.
MAX_BYTES = 1_000_000
BACKUP_COUNT = 5

_handler: Optional[logging.Handler] = None
_excepthook_installed = False

ErrorReporter = Callable[[str, str, str], bool]


def log_dir() -> str:
    """``$FLEXPROBE_HOME/logs``, resolved now — never cached at import."""
    return flexprobe_home("logs")


def log_path() -> str:
    return os.path.join(log_dir(), LOG_FILE_NAME)


def install_file_logging(level: int = logging.INFO) -> Optional[str]:
    """Attach the rotating log file to both logger trees. Idempotent.

    Returns the path being written to, or None if the location is not writable
    — an unwritable home is a reason to run without a log, never a reason to
    refuse to start.

    Re-resolves ``$FLEXPROBE_HOME`` on every call and reinstalls if it has
    moved, so a caller (or a test) that redirects the home directory gets a
    handler pointing at the new one rather than at wherever it pointed first.
    """
    global _handler
    try:
        path = log_path()
    except OSError:
        logger.warning("No writable FLEXPROBE_HOME; running without a log file",
                       exc_info=True)
        return None
    if _handler is not None:
        if getattr(_handler, "baseFilename", None) == os.path.abspath(path):
            return path
        remove_file_logging()
    try:
        handler: logging.Handler = logging.handlers.RotatingFileHandler(
            path, maxBytes=MAX_BYTES, backupCount=BACKUP_COUNT,
            encoding="utf-8",
        )
    except OSError:
        logger.warning("Could not open the log file %s; running without one",
                       path, exc_info=True)
        return None
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    ))
    for name in LOGGER_TREES:
        lg = logging.getLogger(name)
        lg.addHandler(handler)
        # Only ever lower the bar. Raising it would silently override a level
        # someone else set on purpose (the CLI's -v, say).
        if lg.level == logging.NOTSET or lg.level > level:
            lg.setLevel(level)
    _handler = handler
    return path


def remove_file_logging() -> None:
    """Detach and close the log handler. Mainly so tests do not leave one
    pointing at a temporary directory that is about to be deleted."""
    global _handler
    if _handler is None:
        return
    for name in LOGGER_TREES:
        logging.getLogger(name).removeHandler(_handler)
    try:
        _handler.close()
    except Exception:       # noqa: BLE001 — closing a log must not raise
        pass
    _handler = None


def show_error_dialog(title: str, text: str, detail: str) -> bool:
    """Show a copyable error box. Returns whether one was actually shown.

    Everything here is conditional and guarded on purpose: this is called from
    the excepthook, which runs at the worst possible moments — before Qt is up,
    after it has begun tearing down, and on the way out of a crash. It reports
    a failure; it must never become one.
    """
    try:
        from PySide6.QtCore import Qt
        from PySide6.QtWidgets import QApplication, QMessageBox
    except Exception:       # noqa: BLE001 — no Qt, no dialog; the log still has it
        return False
    try:
        app = QApplication.instance()
        if app is None or QApplication.closingDown():
            return False
        box = QMessageBox()
        box.setIcon(QMessageBox.Icon.Critical)
        box.setWindowTitle(title)
        box.setText(text)
        # The whole point: the user has to be able to select this and paste it
        # into a bug report. A screenshot of a traceback is not a bug report.
        box.setDetailedText(detail)
        box.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse
            | Qt.TextInteractionFlag.TextSelectableByKeyboard
        )
        box.setStandardButtons(QMessageBox.StandardButton.Ok)
        box.exec()
        return True
    except Exception:       # noqa: BLE001 — see the docstring
        logger.exception("Could not show the error dialog")
        return False


def install_excepthook(report: Optional[ErrorReporter] = None) -> None:
    """Log and display anything that reaches ``sys.excepthook``. Idempotent.

    ``report`` is injectable so the behaviour can be tested without putting a
    modal dialog on screen; it defaults to :func:`show_error_dialog`.
    """
    global _excepthook_installed
    if _excepthook_installed:
        return
    previous = sys.excepthook
    reporter: ErrorReporter = report or show_error_dialog

    def _hook(exc_type: Type[BaseException], exc: BaseException,
              tb: Optional[TracebackType]) -> None:
        if issubclass(exc_type, KeyboardInterrupt):
            previous(exc_type, exc, tb)     # Ctrl-C is a request, not a fault
            return
        detail = "".join(traceback.format_exception(exc_type, exc, tb))
        try:
            logger.critical("Unhandled exception\n%s", detail)
        except Exception:       # noqa: BLE001 — logging must not mask the crash
            pass
        try:
            previous(exc_type, exc, tb)     # keep stderr's copy for a dev run
        except Exception:       # noqa: BLE001
            pass
        try:
            reporter(
                "FlexiProbes hit an unexpected error",
                f"{exc_type.__name__}: {exc}\n\n"
                f"The full details are below and in the log file:\n"
                f"{log_path_or_unknown()}",
                detail,
            )
        except Exception:       # noqa: BLE001 — a failed report is not a crash
            logger.exception("Could not report an unhandled exception")

    sys.excepthook = _hook
    _excepthook_installed = True


def log_path_or_unknown() -> str:
    """The log's location for display, or a placeholder if it has none."""
    if _handler is not None:
        return str(getattr(_handler, "baseFilename", "(no log file)"))
    return "(no log file)"


def install(level: int = logging.INFO,
            report: Optional[ErrorReporter] = None) -> Optional[str]:
    """Install both halves. Call once, at GUI startup. Returns the log path."""
    path = install_file_logging(level)
    install_excepthook(report)
    return path
