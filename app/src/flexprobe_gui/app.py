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

"""Application entry-point and QApplication bootstrap."""
import logging
import sys


def _set_app_icon(app) -> None:
    """Set the window/dock icon from the packaged PNG, if present."""
    try:
        from importlib.resources import files
        from PySide6.QtGui import QIcon
        icon_path = files("flexprobe_gui") / "assets" / "icon.png"
        if icon_path.is_file():
            app.setWindowIcon(QIcon(str(icon_path)))
    except Exception:
        pass  # a missing icon is purely cosmetic; never block startup


def _install_history_error_dialog() -> None:
    """Tell the user when their run history could not be read.

    ``RunHistory`` reports a damaged ``runs.json`` through a module-level hook
    rather than raising a dialog itself — it is constructed in headless
    contexts too. The dialog is deferred onto the event loop with a zero-delay
    timer so it appears once the window is up, instead of blocking the
    construction of the very window it would have to sit on top of.
    """
    from PySide6.QtCore import QTimer

    from flexprobe_gui import diagnostics, run_history

    def _notify(err) -> None:
        detail = (f"file: {err.path}\n"
                  f"copy kept: {err.backup_path or '(none — see the log)'}")
        QTimer.singleShot(0, lambda: diagnostics.show_error_dialog(
            "FlexiProbes could not read your run history", err.message, detail))

    run_history.set_load_error_notifier(_notify)


def main(argv=None):
    """Launch the flexiprobes GUI."""
    # First thing, before anything can fail: a packaged .app has no console, so
    # until the log file exists every warning and every traceback goes nowhere.
    from flexprobe_gui import diagnostics
    log_file = diagnostics.install()

    try:
        from PySide6.QtCore import Qt
        from PySide6.QtWidgets import QApplication
    except ImportError:
        print("ERROR: PySide6 is required for the GUI.  pip install PySide6", file=sys.stderr)
        sys.exit(1)

    # Crisp rendering on fractional-scale displays (e.g. Windows 150%, Linux
    # fractional scaling). Must be set before the QApplication is constructed.
    QApplication.setHighDpiScaleFactorRoundingPolicy(
        Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
    )

    app = QApplication(argv or sys.argv)
    # Before any window exists, so every combo and spin box is polished under
    # it: the mouse wheel scrolls the page instead of changing their values.
    from flexprobe_gui import wheel_guard
    wheel_guard.install(app)
    app.setApplicationName("flexiprobes")
    app.setApplicationDisplayName("FlexiProbes")
    app.setOrganizationName("Adrian Ilich")
    _set_app_icon(app)
    if log_file:
        logging.getLogger("flexprobe_gui").info(
            "flexiprobes GUI starting; logging to %s", log_file)
    _install_history_error_dialog()

    from flexprobe_gui.home_window import HomeWindow
    win = HomeWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
