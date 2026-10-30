"""Regression tests for pct=None progress events.

On 13 Aug 2026 a whole-transcriptome cluster run died on its own warning: the
pre-screen "exact recount is expensive" notice is emitted with pct=None, and
every consumer formatted pct as a float unconditionally. The run had done no
work yet, so the crash was instant and perfectly reproducible — and the chain
resumed it 11 more times. These tests pin the contract at each sink.
"""
import pytest

from flexprobe.pipeline import ProgressEvent, _warn_recount_cost
from flexprobe.cli import _make_progress_printer


def test_progress_event_accepts_none_pct():
    ev = ProgressEvent("blast", "expensive", None)
    assert ev.pct is None


def test_repr_does_not_format_none():
    # logger.info("%s", ev) inside _emit hits __repr__ before anything else,
    # so this is the *first* thing that crashed on the cluster.
    assert "n/a" in repr(ProgressEvent("blast", "expensive", None))


def test_repr_still_shows_a_percentage_when_known():
    assert "45%" in repr(ProgressEvent("blast", "scoring", 45.0))


def test_cli_printer_handles_none_pct(capsys):
    _make_progress_printer()(ProgressEvent("blast", "expensive", None))
    out = capsys.readouterr().out
    assert "expensive" in out
    assert "--" in out


def test_cli_printer_column_width_is_stable(capsys):
    """The no-pct marker must occupy the same width as a percentage, or the
    job log stops being greppable by column."""
    p = _make_progress_printer()
    p(ProgressEvent("blast", "x", None))
    p(ProgressEvent("blast", "x", 100.0))
    a, b = capsys.readouterr().out.splitlines()
    assert a.index("]") == b.index("]")


def test_cli_printer_still_formats_floats(capsys):
    _make_progress_printer()(ProgressEvent("rank", "done", 100.0))
    assert "100.0%" in capsys.readouterr().out


class _Spec:
    off_target_recount = "exact"


def test_warn_recount_cost_emits_none_pct_without_raising():
    """The end-to-end shape of the original failure: the warning path calls
    emit() with pct=None and the printer must survive it."""
    seen = []

    def emit(stage, msg, pct=None, data=None):
        seen.append(pct)
        _make_progress_printer()(ProgressEvent(stage, msg, pct))

    msg = _warn_recount_cost(_Spec(), 10_000_000, emit)
    assert msg is not None
    assert seen == [None]


def test_pre_fix_printer_would_have_raised():
    """Guards the guard: if this stops raising, the test above proves nothing
    because the naive formatter it stands in for has changed meaning."""
    ev = ProgressEvent("blast", "expensive", None)
    with pytest.raises(TypeError):
        f"[{ev.pct:5.1f}%]"
