"""Every run record carries the ranking weights the run actually used.

``params.ranking_weights`` records only what the user overrode, and the
effective weights also depend on the run: gene-aware mode turns on
``isoform_coverage``, and a non-conserved multi-strain run turns off
``mismatch_placement``. So both files carry all seven under
``effective_ranking_weights``, from the same function
the ranking stage uses.
"""
import dataclasses
import json
import os

from flexprobe.config_io import load_config
from flexprobe.models import (
    ExportSettings, ModalityDefinition, RunConfig, ScreenResources,
)
from flexprobe.pipeline import effective_ranking_weights, run_pipeline
from flexprobe.ranking import DEFAULT_WEIGHTS, RankingWeights

from conftest import SYNTH_200, write_fasta

SEVEN = {f.name for f in dataclasses.fields(RankingWeights)}


def _run(tmp_path, **over):
    fa = write_fasta(tmp_path / "in.fa", [("G1", SYNTH_200), ("G2", SYNTH_200[::-1])])
    cfg = RunConfig(fasta_path=fa, output_dir=str(tmp_path / "out"),
                    modality=ModalityDefinition(name="flex", version="v2",
                                                spatial=False),
                    screen=ScreenResources(skip_screen=True),
                    export=ExportSettings(overwrite=True), workers=1, **over)
    run_pipeline(cfg)
    out = cfg.output_dir
    with open(os.path.join(out, "manifest.json")) as fh:
        manifest = json.load(fh)
    with open(os.path.join(out, "run_config.json")) as fh:
        run_config = json.load(fh)
    return cfg, manifest, run_config


def test_there_are_seven_weights():
    assert len(SEVEN) == 7


def test_overridden_weights_are_recorded_in_both_files(tmp_path):
    cfg, manifest, run_config = _run(
        tmp_path, group_by_gene=False,
        ranking_weights={"gc_proximity": 0.42, "central_placement": 0.11})
    expected = dataclasses.asdict(
        effective_ranking_weights(cfg, gene_aware=False))
    assert expected["gc_proximity"] == 0.42
    assert expected["central_placement"] == 0.11
    assert expected["gc_balance"] == DEFAULT_WEIGHTS.gc_balance

    for record in (manifest["params"], run_config):
        got = record["effective_ranking_weights"]
        assert set(got) == SEVEN
        assert got == expected


def test_the_gene_aware_default_is_recorded_as_applied(tmp_path):
    from flexprobe.ranking import DEFAULT_ISOFORM_COVERAGE_WEIGHT
    _cfg, manifest, run_config = _run(tmp_path, group_by_gene=True)
    applied = manifest["params"]["group_by_gene_applied"]
    want = DEFAULT_ISOFORM_COVERAGE_WEIGHT if applied else 0.0
    assert manifest["params"]["effective_ranking_weights"][
        "isoform_coverage"] == want
    assert run_config["effective_ranking_weights"] == \
        manifest["params"]["effective_ranking_weights"]


def test_run_config_still_reloads(tmp_path):
    cfg, _m, _rc = _run(tmp_path, ranking_weights={"gc_proximity": 0.42})
    again = load_config(os.path.join(cfg.output_dir, "run_config.json"))
    assert again.ranking_weights == {"gc_proximity": 0.42}
