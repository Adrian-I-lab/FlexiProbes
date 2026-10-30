"""RunManager.resume() — re-queues an interrupted run into its own output dir.

Constructs a RunManager headless (no QApplication) and neuters the worker launch
so the queue can be inspected without starting a real pipeline thread.
"""
import os

from flexprobe import candidate_cache
from flexprobe.config_io import save_config
from flexprobe.models import (
    ScreenResources, ModalityDefinition, ProbeCandidate, RunConfig,
)
from flexprobe_gui.run_history import RunHistory, RunRecord
from flexprobe_gui.run_manager import RUN_CONFIG_SIDECAR, RunManager


def _interrupted_output(tmp_path):
    out = str(tmp_path / "out")
    os.makedirs(out)
    cfg = RunConfig(
        fasta_path="in.fa", output_dir=out, resume_screen=True,
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(db_paths=["/db"]),
    )
    cfg.export.overwrite = False           # resume must force this back on
    save_config(cfg, os.path.join(out, RUN_CONFIG_SIDECAR))
    candidate_cache.save(out, "fp", [ProbeCandidate(probe_id="g_t0", gene_id="g")])
    return out


def _manager(tmp_path):
    m = RunManager(RunHistory(path=str(tmp_path / "runs.json")))
    # Don't spin up a real worker thread; just record that the queue was kicked.
    m._started = []
    m._maybe_start_next = lambda: m._started.append(True)
    return m


def test_resume_requeues_into_same_dir_with_overwrite(tmp_path):
    out = _interrupted_output(tmp_path)
    m = _manager(tmp_path)
    rec = RunRecord(name="big", output_dir=out, status="interrupted")
    m.history.append(rec)

    assert m.resume(rec.run_id) is True
    assert m._started == [True]            # the queue was kicked
    assert len(m._queue) == 1
    item = m._queue[0]
    assert item["record"].run_id == rec.run_id
    assert item["record"].status == "queued"
    cfg = item["config"]
    assert cfg.output_dir == out           # same folder -> reuse checkpoints
    assert cfg.resume_screen is True
    assert cfg.export.overwrite is True    # replace any partial final output


def test_resume_false_without_config_sidecar(tmp_path):
    out = str(tmp_path / "out")
    os.makedirs(out)
    m = _manager(tmp_path)
    rec = RunRecord(name="x", output_dir=out, status="failed")
    m.history.append(rec)
    assert m.resume(rec.run_id) is False
    assert not m._queue


def test_resume_false_for_unknown_run(tmp_path):
    m = _manager(tmp_path)
    assert m.resume("does-not-exist") is False
