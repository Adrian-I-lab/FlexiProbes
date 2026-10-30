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

"""Palette-derived colors that adapt to the active light/dark system theme.

Qt 6 already follows the OS color scheme for standard widgets (this app sets no
global stylesheet that would override it). The only thing that breaks in dark
mode is *hardcoded* hex colors in per-widget stylesheets — e.g. a ``#555`` hint
that is unreadable on a dark background. These helpers return theme-appropriate
hex strings from the current QPalette so those stylesheets adapt automatically.
"""
from __future__ import annotations

from PySide6.QtGui import QPalette
from PySide6.QtWidgets import QApplication


def _role_hex(role: QPalette.ColorRole, fallback: str) -> str:
    app = QApplication.instance()
    if app is None:
        return fallback
    color = app.palette().color(role)
    return color.name() if color.isValid() else fallback


def muted_hex() -> str:
    """Hex for secondary / hint text — light-grey on light, dim on dark."""
    return _role_hex(QPalette.ColorRole.PlaceholderText, "#555555")


def link_hex() -> str:
    """Hex for hyperlink-style text — adapts to the active palette."""
    return _role_hex(QPalette.ColorRole.Link, "#1565c0")
