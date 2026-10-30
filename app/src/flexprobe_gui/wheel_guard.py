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

"""Stop the mouse wheel changing a value the user never clicked into.

On Windows (and under Fusion) a wheel event over a ``QComboBox`` or any
``QAbstractSpinBox`` steps its value, so scrolling the New query page quietly
rewrote whichever setting happened to pass under the cursor. The macOS style
ignores the wheel on those controls, which is why the defect never showed on
the machine the GUI was developed on.

One application-wide event filter fixes every dialog at once, on every
platform:

* a wheel event over an UNFOCUSED combo or spin box (or one of its children,
  such as a spin box's line edit) is taken away from the control and handed to
  the nearest enclosing scroll area, so the page scrolls instead;
* a focused control, and an open combo popup, keep the wheel as before;
* each such control is given ``Qt.StrongFocus`` when it is polished. Qt gives
  focus on a wheel event to a ``WheelFocus`` widget BEFORE the event reaches
  any filter, so without this the first wheel tick would focus the control and
  the filter would then let the value change.

The event is forwarded explicitly rather than ignored and left to propagate.
Qt propagates only spontaneous (window-system) wheel events to parents; a
synthesised one stops at its receiver, so an ignore-and-propagate filter
cannot be exercised by a test. Forwarding behaves the same for both. The cost
is that a scroll area already at its limit does not pass the wheel on to an
outer one.
"""

from __future__ import annotations

from PySide6.QtCore import QEvent, QObject, Qt
from PySide6.QtGui import QWheelEvent
from PySide6.QtWidgets import (
    QAbstractScrollArea, QAbstractSpinBox, QApplication, QComboBox, QWidget,
)

_GUARDED = (QComboBox, QAbstractSpinBox)


def _guarded_owner(widget: QWidget) -> QWidget | None:
    """The combo or spin box ``widget`` is, or belongs to, else ``None``.

    The walk stops at a window boundary, so a combo's popup list (a separate
    top-level window parented to the combo) is never treated as the combo.
    """
    w: QWidget | None = widget
    while w is not None:
        if isinstance(w, _GUARDED):
            return w
        if w.isWindow():
            return None
        w = w.parentWidget()
    return None


def _scroll_area_around(widget: QWidget) -> QAbstractScrollArea | None:
    """The nearest scroll area enclosing ``widget`` within its window."""
    w = widget.parentWidget()
    while w is not None:
        if isinstance(w, QAbstractScrollArea):
            return w
        if w.isWindow():
            return None
        w = w.parentWidget()
    return None


def _has_focus(owner: QWidget) -> bool:
    if owner.hasFocus():
        return True
    focused = QApplication.focusWidget()
    return focused is not None and owner.isAncestorOf(focused)


class WheelGuard(QObject):
    """Application event filter; see the module docstring."""

    def eventFilter(self, watched: QObject, event: QEvent) -> bool:
        etype = event.type()
        if etype == QEvent.Polish:
            if isinstance(watched, _GUARDED):
                _strong_focus(watched)
            return False
        if not (isinstance(event, QWheelEvent) and isinstance(watched, QWidget)):
            return False
        owner = _guarded_owner(watched)
        if owner is None or _has_focus(owner):
            return False
        if isinstance(owner, QComboBox) and owner.view().isVisible():
            return False
        area = _scroll_area_around(owner)
        if area is not None:
            viewport = area.viewport()
            forwarded = QWheelEvent(
                viewport.mapFromGlobal(event.globalPosition()),
                event.globalPosition(), event.pixelDelta(),
                event.angleDelta(), event.buttons(), event.modifiers(),
                event.phase(), event.inverted())
            QApplication.sendEvent(viewport, forwarded)
        # Accepted, so Qt does not also propagate a real wheel event to the
        # parent and scroll the page twice. With no scroll area it is simply
        # swallowed: there is nothing to scroll.
        event.accept()
        return True


def _strong_focus(widget: QWidget) -> None:
    """Click and Tab focus, never wheel focus. A deliberate ``NoFocus`` stays."""
    if widget.focusPolicy() != Qt.NoFocus:
        widget.setFocusPolicy(Qt.StrongFocus)


def install(app: QApplication) -> WheelGuard:
    """Install the guard on ``app`` once; return the installed filter.

    Parented to the application so it lives as long as the app does. Widgets
    that already exist were polished before the guard could see them, so they
    are given their focus policy here.
    """
    existing = app.findChild(WheelGuard)
    if existing is not None:
        return existing
    guard = WheelGuard(app)
    app.installEventFilter(guard)
    for widget in app.allWidgets():
        if isinstance(widget, _GUARDED):
            _strong_focus(widget)
    return guard
