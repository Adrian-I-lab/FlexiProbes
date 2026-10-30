"""group_by_gene defaults to True.

It must be a harmless no-op on input without Ensembl-style ``gene:`` header
fields (viral and custom inputs). That is not automatic: with grouping on and
no labels, each record becomes its own gene, so a probe that also hits ANOTHER
input record turns from on-target into off-target. So when no input record
carries a ``gene:`` field, the run switches grouping off and says so in one
log line, and its output is byte-identical to an explicit grouping-off run.
"""
import logging
import os

from conftest import SYNTH_200, write_fasta
from flexprobe.models import (
    ExportSettings, ModalityDefinition, RunConfig, ScreenResources,
)
from flexprobe.pipeline import run_pipeline

OTHER = SYNTH_200[::-1]


def _cfg(tmp_path, fa, name, **kw):
    return RunConfig(
        fasta_path=fa, output_dir=str(tmp_path / name),
        modality=ModalityDefinition(name="flex", version="v2"),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True), **kw)


def _manifest(cfg):
    import json
    with open(os.path.join(cfg.output_dir, "manifest.json")) as fh:
        return json.load(fh)


def test_the_default_is_true():
    assert RunConfig(fasta_path="x.fa").group_by_gene is True


def test_gene_labelled_input_is_grouped_by_default(tmp_path):
    fa = tmp_path / "tx.fa"
    fa.write_text(f">t1 gene:G1\n{SYNTH_200}\n>t2 gene:G1\n{SYNTH_200}\n"
                  f">t3 gene:G2\n{OTHER}\n")
    cfg = _cfg(tmp_path, str(fa), "lab")
    run_pipeline(cfg)
    params = _manifest(cfg)["params"]
    assert params["group_by_gene"] is True and params["group_by_gene_applied"] is True


def test_unlabelled_input_is_unchanged_and_says_so(tmp_path, caplog):
    fa = write_fasta(tmp_path / "v.fa", [("L", SYNTH_200), ("S", OTHER)])
    default = _cfg(tmp_path, fa, "default")
    explicit_off = _cfg(tmp_path, fa, "off", group_by_gene=False)
    with caplog.at_level(logging.INFO, logger="flexprobe"):
        run_pipeline(default)
    assert any("gene:" in r.getMessage() and "grouping" in r.getMessage()
               for r in caplog.records), [r.getMessage() for r in caplog.records]
    run_pipeline(explicit_off)
    params = _manifest(default)["params"]
    assert params["group_by_gene"] is True          # what was asked for
    assert params["group_by_gene_applied"] is False  # what ran
    with open(os.path.join(default.output_dir, "results.csv"), "rb") as a, \
            open(os.path.join(explicit_off.output_dir, "results.csv"), "rb") as b:
        assert a.read() == b.read()


def test_cli_default_and_opt_out():
    from flexprobe import cli
    p = cli.build_parser()
    on = cli._build_run_config(p.parse_args(["design", "--fasta", "a.fa"]))
    off = cli._build_run_config(p.parse_args(["design", "--fasta", "a.fa",
                                             "--no-group-by-gene"]))
    assert on.group_by_gene is True and off.group_by_gene is False
