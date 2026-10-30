"""A corrupt runs.json must never cost the user their run history.

If ``RunHistory._load`` turned a parse failure into an empty list, the next
terminal run would write that empty list straight back over the file, and
every recorded run would be lost with nothing said and nothing kept. A
half-written runs.json is exactly what a power-off during a flush leaves
behind, so the corrupt file must be backed up and the failure logged.
"""
import json
import logging
import os

import pytest

from flexprobe_gui import run_history as rh
from flexprobe_gui.run_history import RunHistory, RunRecord


GOOD = [
    {"run_id": "aaa", "name": "run-a", "status": "complete"},
    {"run_id": "bbb", "name": "run-b", "status": "complete"},
    {"run_id": "ccc", "name": "run-c", "status": "failed"},
]


@pytest.fixture(autouse=True)
def _no_stray_notifier():
    rh.set_load_error_notifier(None)
    yield
    rh.set_load_error_notifier(None)


def _corrupt(tmp_path, text='[{"run_id": "aaa", "name": "run-a", "sta'):
    path = str(tmp_path / "runs.json")
    with open(path, "w") as fh:
        fh.write(text)
    return path


def test_corrupt_history_is_copied_aside_before_anything_overwrites_it(tmp_path):
    original = '[{"run_id": "aaa", "name": "run-a", "sta'
    path = _corrupt(tmp_path, original)

    h = RunHistory(path=path)
    # The app carries on with an empty list, which is correct.
    assert h.list_completed() == []

    bak = path + ".bak"
    assert os.path.exists(bak), "the corrupt file was discarded without a copy"
    with open(bak) as fh:
        assert fh.read() == original

    # ...and it survives the next write to runs.json.
    h.append(RunRecord(run_id="new", name="fresh", status="complete"))
    with open(bak) as fh:
        assert fh.read() == original


def test_a_second_corruption_does_not_clobber_the_first_evidence(tmp_path):
    first = '[{"run_id": "aaa"'
    path = _corrupt(tmp_path, first)
    RunHistory(path=path)

    second = '[{"run_id": "zzz"'
    _corrupt(tmp_path, second)
    RunHistory(path=path)

    with open(path + ".bak") as fh:
        assert fh.read() == first
    with open(path + ".bak.1") as fh:
        assert fh.read() == second


def test_valid_json_that_is_not_a_run_list_counts_as_corrupt(tmp_path):
    # A dict parses fine and then iterates as its KEYS, none of which are
    # records — so the old code quietly produced an empty history from a file
    # that was never a history at all.
    path = str(tmp_path / "runs.json")
    with open(path, "w") as fh:
        json.dump({"runs": GOOD}, fh)

    h = RunHistory(path=path)

    assert h.list_completed() == []
    assert os.path.exists(path + ".bak")
    assert h.load_error is not None


def test_a_healthy_history_is_left_alone(tmp_path):
    path = str(tmp_path / "runs.json")
    with open(path, "w") as fh:
        json.dump(GOOD, fh)

    h = RunHistory(path=path)

    assert sorted(r.name for r in h.list_completed()) == ["run-a", "run-b",
                                                          "run-c"]
    assert not os.path.exists(path + ".bak")
    assert h.load_error is None


def test_a_missing_history_is_not_a_corruption(tmp_path):
    h = RunHistory(path=str(tmp_path / "never-existed.json"))
    assert h.load_error is None
    assert not os.path.exists(str(tmp_path / "never-existed.json.bak"))


def test_the_corruption_is_logged_at_error(tmp_path, caplog):
    path = _corrupt(tmp_path)
    with caplog.at_level(logging.ERROR, logger="flexprobe_gui.run_history"):
        RunHistory(path=path)
    assert any(r.levelno >= logging.ERROR for r in caplog.records), (
        "a destroyed run history has to reach the diagnostics log"
    )


def test_the_corruption_is_surfaced_to_the_user(tmp_path):
    seen = []
    rh.set_load_error_notifier(seen.append)
    path = _corrupt(tmp_path)

    h = RunHistory(path=path)

    assert len(seen) == 1
    err = seen[0]
    assert err is h.load_error
    assert err.path == path
    assert err.backup_path == path + ".bak"
    assert "run" in err.message.lower()


def test_a_broken_notifier_cannot_stop_the_app_starting(tmp_path):
    def _boom(_err):
        raise RuntimeError("the UI is not up yet")

    rh.set_load_error_notifier(_boom)
    h = RunHistory(path=_corrupt(tmp_path))
    assert h.load_error is not None
