"""The filter report must never show a negative count.

``filters.py`` checks every rule on every tile and counts a failing tile once
under EACH rule it fails. So the per-rule counts overlap, and their sum is
normally larger than the number of rejected tiles. Subtracting that sum from
``total - passed`` and printing the result as "rejected for a reason not
listed below" would give a negative number on most real runs.
"""
from __future__ import annotations

import json
import re

import pytest

pytest.importorskip("PySide6")


@pytest.fixture(scope="module")
def qapp():
    from PySide6.QtWidgets import QApplication
    app = QApplication.instance() or QApplication([])
    yield app


def _dialog(tmp_path, counts, thresholds=None):
    from flexprobe_gui.filter_report_dialog import FilterReportDialog
    (tmp_path / "filter_report.json").write_text(json.dumps({
        "counts": counts,
        "thresholds": thresholds or {},
        "examples": {},
    }))
    return FilterReportDialog(str(tmp_path), "run")


def _all_text(dlg) -> str:
    from PySide6.QtWidgets import QLabel, QTableWidget
    parts = [lbl.text() for lbl in dlg.findChildren(QLabel)]
    for t in dlg.findChildren(QTableWidget):
        for c in range(t.columnCount()):
            h = t.horizontalHeaderItem(c)
            parts.append(h.text() if h else "")
    return "\n".join(parts)


# Overlapping counts from a real run: 594 rejected tiles, 721 rule failures
# between them.
R01 = {"total": 767, "passed": 173, "tn_junction": 528, "tm_range": 125,
       "homopolymer": 68, "gc_range": 0}


def test_overlapping_rule_counts_never_show_a_negative_number(qapp, tmp_path):
    text = _all_text(_dialog(tmp_path, R01))
    assert not re.search(r"-\d", re.sub(r"<[^>]+>", "", text)), text
    assert "not listed below" not in text


def test_headline_states_rejections_as_total_minus_passed(qapp, tmp_path):
    text = re.sub(r"<[^>]+>", "", _all_text(_dialog(tmp_path, R01)))
    assert "767" in text and "173" in text
    assert "594 rejected" in text


def test_rule_column_says_it_counts_tiles_failing_that_rule(qapp, tmp_path):
    text = _all_text(_dialog(tmp_path, R01))
    assert "Tiles failing this rule" in text
    assert "can fail several rules" in text


def test_rejected_tiles_with_no_recorded_reason_are_still_flagged(qapp, tmp_path):
    # 6 rejected, but the rules account for only 2: an engine accounting gap,
    # which is the one case the old line was meant to catch.
    text = re.sub(r"<[^>]+>", "", _all_text(_dialog(
        tmp_path, {"total": 10, "passed": 4, "tm_range": 2})))
    assert "4 rejected tiles have no recorded reason" in text


def test_gc_gate_is_described_in_words_not_as_a_python_boolean(qapp, tmp_path):
    thr = {"gc_min": 0.44, "gc_max": 0.72, "enforce_gc_range": False,
           "tm_min": 56.0, "max_n_per_half": 1, "max_homopolymer_run": 4,
           "tn_junction_mode": "strict_a"}
    text = _all_text(_dialog(tmp_path, R01, thr))
    assert "gate: False" not in text and "gate: True" not in text
    assert "not enforced" in text


def test_a_switched_off_tm_gate_reads_off_not_none(qapp, tmp_path):
    thr = {"gc_min": 0.44, "gc_max": 0.72, "enforce_gc_range": False,
           "tm_min": None, "max_n_per_half": 1, "max_homopolymer_run": 4,
           "tn_junction_mode": "strict_a"}
    text = _all_text(_dialog(tmp_path, R01, thr))
    assert "None" not in text
    assert "Tm/half off" in text
