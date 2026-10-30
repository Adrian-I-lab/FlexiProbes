"""In-app resume: config sidecar, resumability detection, crash recovery.

Logic-only (no QApplication): exercises the persistence/recovery helpers that
back the GUI's Resume button, not the Qt widgets themselves.
"""
import os

from flexprobe import candidate_cache
from flexprobe.config_io import save_config
from flexprobe.models import (
    ScreenResources, ModalityDefinition, ProbeCandidate, RunConfig,
)

from flexprobe_gui.run_history import RunHistory, RunRecord
from flexprobe_gui.run_manager import (
    RUN_CONFIG_SIDECAR, is_resumable, load_run_config,
)


def _cfg(out):
    return RunConfig(
        fasta_path="in.fa", output_dir=out,
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(db_paths=["/db"]),
        tile_length=50, tile_step=3, screen_top_per_gene=25,
    )


# ---- config sidecar round-trip ----

def test_load_run_config_round_trip(tmp_path):
    out = str(tmp_path / "out")
    os.makedirs(out)
    save_config(_cfg(out), os.path.join(out, RUN_CONFIG_SIDECAR))
    cfg = load_run_config(out)
    assert cfg is not None
    assert cfg.tile_step == 3 and cfg.screen_top_per_gene == 25
    assert cfg.screen.db_paths == ["/db"]


def test_load_run_config_missing_returns_none(tmp_path):
    assert load_run_config(str(tmp_path / "nope")) is None


# ---- resumability detection ----

def _record(out):
    return RunRecord(name="r", output_dir=out, status="cancelled")


def test_is_resumable_requires_sidecar_and_checkpoint(tmp_path):
    out = str(tmp_path / "out")
    os.makedirs(out)
    rec = _record(out)
    # nothing yet
    assert not is_resumable(rec)
    # sidecar alone is not enough (no saved progress to resume)
    save_config(_cfg(out), os.path.join(out, RUN_CONFIG_SIDECAR))
    assert not is_resumable(rec)
    # sidecar + a checkpoint dir -> resumable
    candidate_cache.save(out, "fp", [ProbeCandidate(probe_id="g_t0", gene_id="g")])
    assert is_resumable(rec)


def test_is_resumable_false_when_checkpoints_cleared(tmp_path):
    out = str(tmp_path / "out")
    os.makedirs(out)
    save_config(_cfg(out), os.path.join(out, RUN_CONFIG_SIDECAR))
    candidate_cache.save(out, "fp", [ProbeCandidate(probe_id="g_t0", gene_id="g")])
    candidate_cache.clear(out)            # simulate a fully-successful finish
    assert not is_resumable(_record(out))


def test_is_resumable_false_without_output_dir():
    assert not is_resumable(RunRecord(name="r", output_dir="", status="failed"))


# ---- crash / power-off recovery ----

def test_running_record_recovers_as_interrupted(tmp_path):
    path = str(tmp_path / "runs.json")
    h = RunHistory(path=path)
    rec = RunRecord(name="big", output_dir=str(tmp_path / "out"))
    h.append(rec)                 # queued — not yet persisted
    rec.status = "running"
    h.update(rec)                 # running update doesn't force a save...
    h.flush()                     # ...but starting a run does (persist running)

    # Reload as if the app restarted after a crash mid-run.
    h2 = RunHistory(path=path)
    again = h2.get(rec.run_id)
    assert again is not None
    assert again.status == "interrupted"     # stale "running" -> resumable
    # interrupted is terminal, so it shows in the completed list
    assert any(r.run_id == rec.run_id for r in h2.list_completed())


def test_queued_only_run_not_persisted(tmp_path):
    # A run that never started (still queued at exit) writes no record.
    path = str(tmp_path / "runs.json")
    h = RunHistory(path=path)
    h.append(RunRecord(name="q", status="queued"))
    h2 = RunHistory(path=path)
    assert h2.list_completed() == [] and h2.list_active() == []
