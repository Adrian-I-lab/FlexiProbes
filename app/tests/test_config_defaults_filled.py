"""A config loaded from disk says which settings it did NOT carry.

``config_from_dict`` fills every missing key from today's dataclass defaults, so
an old config reruns under settings nobody chose. A config written by an older
release can lack ``prepare_only``, ``screen.reference_strand`` and
``strain_unique_junction_duplex_bp``, the last being the strain-unique
junction-duplex gate that did not exist when that config was written.
Reloading it switches that gate on with no trace.

Filling stays (it is what lets an old config load at all). What changes is
that it is never silent: the loader logs the filled keys, and a run started
from a loaded config writes them into ``manifest.json`` as
``config_filled_from_defaults``. A config built in memory (GUI, CLI flags) was
not loaded from anything, so its manifest records null there.
"""
import json
import logging
import os

from flexprobe.config_io import (
    config_from_dict, config_to_dict, filled_from_defaults, load_config,
    save_config,
)
from flexprobe.models import (
    ExportSettings, ModalityDefinition, RunConfig, ScreenResources,
)
from flexprobe.pipeline import run_pipeline

from conftest import SYNTH_200, write_fasta

# The three keys a config written by an older release lacks.
OLD_CONFIG_MISSING = ("prepare_only", "screen.reference_strand",
                      "strain_unique_junction_duplex_bp")


def _old_style_payload(tmp_path):
    fa = write_fasta(tmp_path / "g.fa", [("G1", SYNTH_200)])
    cfg = RunConfig(
        fasta_path=fa, output_dir=str(tmp_path / "out"),
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True),
    )
    data = config_to_dict(cfg)
    del data["prepare_only"]
    del data["strain_unique_junction_duplex_bp"]
    del data["screen"]["reference_strand"]
    return data


def test_a_complete_config_fills_nothing():
    cfg = RunConfig(fasta_path="x.fa")
    assert filled_from_defaults(config_from_dict(config_to_dict(cfg))) == []


def test_an_old_config_names_exactly_the_keys_it_lacked(tmp_path):
    cfg = config_from_dict(_old_style_payload(tmp_path))
    assert filled_from_defaults(cfg) == sorted(OLD_CONFIG_MISSING)


def test_a_missing_block_names_every_leaf_in_it():
    data = config_to_dict(RunConfig(fasta_path="x.fa"))
    n_export = len(data.pop("export"))
    filled = filled_from_defaults(config_from_dict(data))
    assert len([k for k in filled if k.startswith("export.")]) == n_export > 0


def test_a_legacy_blast_block_is_not_reported_as_missing():
    """0.6.0 renamed the block. The key moved, the values did not."""
    data = config_to_dict(RunConfig(fasta_path="x.fa"))
    data["blast"] = data.pop("screen")
    assert filled_from_defaults(config_from_dict(data)) == []


def test_filling_is_logged(tmp_path, caplog):
    path = tmp_path / "old.json"
    path.write_text(json.dumps(_old_style_payload(tmp_path)))
    with caplog.at_level(logging.WARNING, logger="flexprobe.config_io"):
        load_config(str(path))
    msg = " ".join(r.getMessage() for r in caplog.records)
    for key in OLD_CONFIG_MISSING:
        assert key in msg, msg


def test_a_config_built_in_memory_reports_none():
    assert filled_from_defaults(RunConfig(fasta_path="x.fa")) is None


def test_loaded_old_config_run_records_filled_keys_and_writes_full_config(tmp_path):
    path = tmp_path / "old.json"
    path.write_text(json.dumps(_old_style_payload(tmp_path)))
    cfg = load_config(str(path))
    run_pipeline(cfg)

    with open(os.path.join(cfg.output_dir, "manifest.json")) as fh:
        manifest = json.load(fh)
    assert manifest["config_filled_from_defaults"] == sorted(OLD_CONFIG_MISSING)

    # The run_config.json written beside the results is complete: every key
    # the old file lacked is now explicit, so a rerun of THIS file fills
    # nothing and is exactly what ran.
    written = os.path.join(cfg.output_dir, "run_config.json")
    with open(written) as fh:
        data = json.load(fh)
    assert "prepare_only" in data
    assert "strain_unique_junction_duplex_bp" in data
    assert "reference_strand" in data["screen"]
    assert filled_from_defaults(load_config(written)) == []
    assert load_config(written) == cfg


def test_in_memory_run_manifest_records_null(tmp_path):
    fa = write_fasta(tmp_path / "g.fa", [("G1", SYNTH_200)])
    cfg = RunConfig(
        fasta_path=fa, output_dir=str(tmp_path / "out"),
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True),
    )
    run_pipeline(cfg)
    with open(os.path.join(cfg.output_dir, "manifest.json")) as fh:
        assert json.load(fh)["config_filled_from_defaults"] is None


def test_save_config_of_a_loaded_config_round_trips_clean(tmp_path):
    """Saving never leaks the report into the file as if it were a setting."""
    cfg = config_from_dict(_old_style_payload(tmp_path))
    out = save_config(cfg, str(tmp_path / "again.json"))
    with open(out) as fh:
        assert "config_filled_from_defaults" not in fh.read()
