"""New-query form layout checks.

* Hiding a homology-mode setting must hide the whole row, including the unit
  text, so "isoforms (0 = all)" and similar labels do not float in the form
  with nothing beside them.
* The filters group title "Filters & tile" must keep its "&" rather than lose
  it to Qt's shortcut marker and render as "Filters _tile".
* The dialog must open at least as wide as its own content, so it has no
  sideways scrollbar and no clipped Browse buttons.
"""
import importlib.util
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 not installed"
)

UNIT_TEXTS = ("isoforms (0 = all)", "per 25 nt half (0 = identical everywhere)",
              "bp (≥5 matches the 10x off-target threshold)")


@pytest.fixture(scope="module")
def qapp():
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    yield app


@pytest.fixture()
def dlg(qapp, tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXPROBE_HOME", str(tmp_path))
    from flexprobe_gui.new_query_dialog import NewQueryDialog
    d = NewQueryDialog()
    d.show()
    qapp.processEvents()
    yield d
    d.close()


def _visible_label_texts(d):
    from PySide6.QtWidgets import QLabel
    return [lb.text() for lb in d.findChildren(QLabel) if lb.isVisible()]


def test_single_sequence_mode_shows_no_orphan_unit_labels(dlg):
    shown = _visible_label_texts(dlg)
    for t in UNIT_TEXTS:
        assert t not in shown, f"orphan unit label visible: {t!r}"


def test_each_unit_label_appears_with_its_own_mode(dlg, qapp):
    from flexprobe_gui._nq_build import (
        STRAIN_MODE_HOMOLOGY, STRATEGY_CONSERVED, STRATEGY_CONSTITUTIVE,
        STRATEGY_UNIQUE)
    dlg.strain_mode_combo.setCurrentText(STRAIN_MODE_HOMOLOGY)
    expect = {STRATEGY_CONSERVED: UNIT_TEXTS[1],
              STRATEGY_CONSTITUTIVE: UNIT_TEXTS[0],
              STRATEGY_UNIQUE: UNIT_TEXTS[2]}
    for strategy, unit in expect.items():
        dlg.strategy_combo.setCurrentText(strategy)
        qapp.processEvents()
        shown = _visible_label_texts(dlg)
        assert unit in shown, (strategy, unit)
        for other in set(UNIT_TEXTS) - {unit}:
            assert other not in shown, (strategy, other)


def test_the_filters_group_title_keeps_its_ampersand(dlg):
    from PySide6.QtWidgets import QGroupBox
    titles = [g.title() for g in dlg.findChildren(QGroupBox)]
    filt = [t for t in titles if t.startswith("Filters")]
    assert filt, titles
    # A single "&" is a shortcut marker and is not drawn. The title must not
    # rely on one.
    assert "&" not in filt[0].replace("&&", ""), filt[0]


def test_the_dialog_opens_without_a_sideways_scrollbar(dlg, qapp):
    from PySide6.QtWidgets import QScrollArea
    area = dlg.findChildren(QScrollArea)[0]
    qapp.processEvents()
    assert area.horizontalScrollBar().maximum() == 0, (
        dlg.width(), area.widget().minimumSizeHint().width())
