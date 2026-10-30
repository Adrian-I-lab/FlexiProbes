"""The mouse wheel scrolls the page; it does not change an unfocused value.

On Windows, scrolling the New query page rewrote whichever combo box or spin
box passed under the cursor. The guard in ``flexprobe_gui.wheel_guard`` hands
those wheel events to the enclosing scroll area instead.

Every case asserts both halves of what the user sees: the value did not move
AND the page did (or the reverse, for a focused control). A guard that only
swallowed the wheel would pass the first half and leave the page frozen.

The combo is given the Fusion style explicitly. The macOS style ignores the
wheel on combos by itself, which would let the combo checks pass with no guard
at all; ``test_without_the_guard_the_wheel_changes_unfocused_values`` is the
control that shows these widgets really do react to the wheel.
"""
import importlib.util
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 not installed"
)


@pytest.fixture
def guard(qapp):
    from flexprobe_gui.wheel_guard import WheelGuard
    g = WheelGuard()
    qapp.installEventFilter(g)
    yield g
    qapp.removeEventFilter(g)


class _Page:
    """A short scroll area over a tall page holding a spin box and a combo."""

    def __init__(self):
        from PySide6.QtWidgets import (
            QComboBox, QDoubleSpinBox, QLineEdit, QScrollArea, QSpinBox,
            QStyleFactory, QVBoxLayout, QWidget,
        )
        self.area = QScrollArea()
        self.area.resize(300, 200)
        inner = QWidget()
        lay = QVBoxLayout(inner)
        # Somewhere for focus to sit that is not a guarded control.
        self.decoy = QLineEdit()
        self.spin = QSpinBox()
        self.spin.setRange(0, 100)
        self.spin.setValue(50)
        self.dspin = QDoubleSpinBox()
        self.combo = QComboBox()
        self.combo.addItems([str(i) for i in range(10)])
        self.combo.setCurrentIndex(5)
        # Held on the page: setStyle does not take ownership of the style.
        self._fusion = QStyleFactory.create("Fusion")
        self.combo.setStyle(self._fusion)
        for w in (self.decoy, self.spin, self.dspin, self.combo):
            lay.addWidget(w)
        lay.addSpacing(2000)
        self.area.setWidget(inner)
        self.area.show()
        self.area.activateWindow()
        self.focus(self.decoy)

    def focus(self, widget):
        from PySide6.QtWidgets import QApplication
        widget.setFocus()
        QApplication.processEvents()
        assert widget.hasFocus(), "offscreen focus failed"

    @property
    def scrolled(self):
        return self.area.verticalScrollBar().value()

    def close(self):
        self.area.close()
        self.area.deleteLater()


def _wheel(widget, dy=-120):
    """Send one wheel notch (down, by default) to ``widget``."""
    from PySide6.QtCore import QPoint, QPointF, Qt
    from PySide6.QtGui import QWheelEvent
    from PySide6.QtWidgets import QApplication
    local = QPoint(4, 4)
    event = QWheelEvent(
        QPointF(local), QPointF(widget.mapToGlobal(local)), QPoint(0, 0),
        QPoint(0, dy), Qt.NoButton, Qt.NoModifier, Qt.NoScrollPhase, False)
    QApplication.sendEvent(widget, event)
    return event


def test_without_the_guard_the_wheel_changes_unfocused_values(qapp):
    """Control: the widgets below are wheel-sensitive, so a pass elsewhere in
    this file means the guard stopped something, not that nothing happened."""
    page = _Page()
    _wheel(page.spin)
    _wheel(page.combo)
    assert not page.spin.hasFocus() and not page.combo.hasFocus()
    assert page.spin.value() == 49
    assert page.combo.currentIndex() == 6
    assert page.scrolled == 0
    page.close()


def test_unfocused_spin_box_scrolls_the_page_instead(qapp, guard):
    page = _Page()
    _wheel(page.spin)
    assert page.spin.value() == 50
    assert page.scrolled > 0
    page.close()


def test_unfocused_combo_box_scrolls_the_page_instead(qapp, guard):
    page = _Page()
    _wheel(page.combo)
    assert page.combo.currentIndex() == 5
    assert page.scrolled > 0
    page.close()


def test_wheel_over_a_spin_box_line_edit_scrolls_the_page(qapp, guard):
    """The cursor is usually over the spin box's text, i.e. its child
    QLineEdit, not the spin box itself."""
    page = _Page()
    before = page.dspin.value()
    _wheel(page.dspin.lineEdit())
    assert page.dspin.value() == before
    assert page.scrolled > 0
    page.close()


def test_scrolling_back_up_works_too(qapp, guard):
    page = _Page()
    for _ in range(3):
        _wheel(page.spin)
    low = page.scrolled
    assert low > 0
    _wheel(page.spin, dy=120)
    assert 0 <= page.scrolled < low
    assert page.spin.value() == 50
    page.close()


def test_a_focused_spin_box_still_takes_the_wheel(qapp, guard):
    page = _Page()
    page.focus(page.spin)
    _wheel(page.spin)
    assert page.spin.value() == 49
    assert page.scrolled == 0
    page.close()


def test_a_focused_combo_box_still_takes_the_wheel(qapp, guard):
    page = _Page()
    page.focus(page.combo)
    _wheel(page.combo)
    assert page.combo.currentIndex() == 6
    assert page.scrolled == 0
    page.close()


def test_guarded_widgets_get_strong_focus_when_polished(qapp, guard):
    """WheelFocus would let the first wheel tick focus the control before the
    filter sees the event. StrongFocus keeps click and Tab focus."""
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QSpinBox
    page = _Page()
    for w in (page.spin, page.dspin, page.combo):
        assert w.focusPolicy() == Qt.StrongFocus, type(w).__name__
    # A plain line edit is not the guard's business.
    assert page.decoy.focusPolicy() == Qt.StrongFocus
    # A control made unfocusable on purpose stays that way.
    off = QSpinBox()
    off.setFocusPolicy(Qt.NoFocus)
    off.ensurePolished()
    assert off.focusPolicy() == Qt.NoFocus
    page.close()
    off.deleteLater()


def test_an_open_combo_popup_keeps_the_wheel(qapp, guard):
    page = _Page()
    page.combo.showPopup()
    qapp.processEvents()
    assert page.combo.view().isVisible()
    view_port = page.combo.view().viewport()
    event = _wheel(view_port)
    assert guard.eventFilter(view_port, event) is False
    assert guard.eventFilter(page.combo, event) is False
    page.combo.hidePopup()
    page.close()


def test_a_control_outside_any_scroll_area_is_left_unchanged(qapp, guard):
    from PySide6.QtWidgets import QSpinBox, QVBoxLayout, QWidget
    host = QWidget()
    spin = QSpinBox(host)
    spin.setValue(5)
    QVBoxLayout(host).addWidget(spin)
    host.show()
    qapp.processEvents()
    spin.clearFocus()
    event = _wheel(spin)
    assert spin.value() == 5
    assert event.isAccepted()
    host.close()
    host.deleteLater()


def test_install_is_idempotent_and_fixes_existing_widgets(qapp):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QComboBox
    from flexprobe_gui import wheel_guard
    early = QComboBox()
    early.setFocusPolicy(Qt.WheelFocus)
    first = wheel_guard.install(qapp)
    try:
        assert wheel_guard.install(qapp) is first
        assert early.focusPolicy() == Qt.StrongFocus
    finally:
        qapp.removeEventFilter(first)
        first.setParent(None)
        first.deleteLater()
        early.deleteLater()
