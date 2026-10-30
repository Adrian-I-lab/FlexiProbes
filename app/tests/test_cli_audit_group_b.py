"""CLI flags whose behaviour must match their help.

Each check pins that a silent trap is loud, or that the help is true. The
section comments name the trap each group of tests closes.
"""
import pytest

from flexprobe.errors import ConfigError
from flexprobe.models import ModalityDefinition, RunConfig


@pytest.fixture(autouse=True)
def _restore_logging():
    """cli.main() calls setup_logging, which attaches a handler bound to the
    CURRENT stderr. Under capsys that stream is closed after the test, so a
    later test logging to the same tree would write into a closed file
    (seen as test_gui_diagnostics failing only after these tests). Put every
    touched logger back as it was."""
    import logging
    names = (None, "flexprobe", "flexprobe_gui")
    saved = {n: (list(logging.getLogger(n).handlers), logging.getLogger(n).level)
             for n in names}
    yield
    for n, (handlers, level) in saved.items():
        lg = logging.getLogger(n)
        lg.handlers[:] = handlers
        lg.setLevel(level)


def _cfg(**kw):
    kw.setdefault("modality", ModalityDefinition(name="flex", version="v2"))
    kw.setdefault("output_dir", "o")
    return RunConfig(fasta_path="x.fa", **kw)


# -- --tile-size other than the chemistry's window silently gave 0 probes --

def test_a_tile_length_the_chemistry_cannot_assemble_is_refused():
    with pytest.raises(ConfigError, match="tile_length 60"):
        _cfg(tile_length=60).validate()


def test_the_chemistry_window_still_validates():
    _cfg(tile_length=50).validate()


# -- a negative --conserved-max-degenerate-per-half was clamped to 0 --

def test_a_negative_degenerate_budget_is_refused():
    with pytest.raises(ConfigError, match="conserved_max_degenerate_per_half"):
        _cfg(conserved_max_degenerate_per_half=-1).validate()


def test_a_zero_degenerate_budget_is_allowed():
    _cfg(conserved_max_degenerate_per_half=0).validate()


# -- --config silently ignored every design flag it does not override --

def _config_file(tmp_path):
    from flexprobe.config_io import save_config
    # output_dir inside tmp_path: if the refusal ever regresses, the design run
    # that follows must not write into the working directory.
    return save_config(_cfg(output_dir=str(tmp_path / "out")),
                       str(tmp_path / "config.json"))


def test_a_design_flag_with_config_is_refused_by_name(tmp_path, capsys):
    from flexprobe import cli
    with pytest.raises(SystemExit) as exc:
        cli.main(["design", "--config", _config_file(tmp_path),
                  "--tile-step", "5", "--gc-min", "0.4"])
    assert exc.value.code == 2
    err = capsys.readouterr().err
    assert "ignored with --config" in err
    assert "--tile-step" in err and "--gc-min" in err


def test_the_honoured_overrides_are_not_flagged(tmp_path):
    from flexprobe import cli
    parser = cli.build_parser()
    args = parser.parse_args(
        ["design", "--config", _config_file(tmp_path), "--fasta", "a.fa",
         "--out", "o", "--workers", "2", "--force", "--no-resume",
         "--skip-screen", "--screen-top-per-gene", "0"])
    assert cli._flags_ignored_by_config(args) == []


def test_export_hpc_is_checked_too(tmp_path, capsys):
    from flexprobe import cli
    with pytest.raises(SystemExit):
        cli.main(["export-hpc", "--config", _config_file(tmp_path),
                  "--dest", str(tmp_path / "b"), "--image", "f.sif",
                  "--tile-step", "5"])
    err = capsys.readouterr().err
    assert "ignored with --config" in err and "--tile-step" in err


def test_override_list_matches_what_build_run_config_applies():
    """CONFIG_OVERRIDES must name exactly the args the --config branch reads."""
    import inspect
    import re
    from flexprobe import cli
    src = inspect.getsource(cli._build_run_config)
    branch = src.split('if getattr(args, "config", None):', 1)[1].split("return cfg", 1)[0]
    read = set(re.findall(r"args\.([a-z_]+)", branch))
    read |= set(re.findall(r'getattr\(args, "([a-z_]+)"', branch))
    assert read | {"config"} == set(cli.CONFIG_OVERRIDES)


def _help_of(dest, command="design"):
    from flexprobe import cli
    sub = cli.build_parser()._subparsers._group_actions[0].choices[command]
    return next(a.help for a in sub._actions if a.dest == dest)


def test_help_lists_the_real_overrides():
    text = _help_of("config")
    for flag in ("--off-target-ref", "--prepare-only", "--screen-workers",
                 "--no-resume", "--keep-candidate-cache"):
        assert flag in text, flag
    assert "--blast-db" not in text and "refused" in text


# -- --junction-window was stored and never read; the engine used 4 --

def _run_scenario_capturing(monkeypatch, tmp_path, scenario, module, func, window):
    import dataclasses
    import golden_scenarios as g
    from flexprobe.pipeline import run_pipeline
    seen = []
    real = getattr(module, func)

    def spy(*a, **k):
        seen.append(k.get("junction_window"))
        return real(*a, **k)

    monkeypatch.setattr(module, func, spy)
    cfg = g.SCENARIOS[scenario](str(tmp_path))
    cfg = dataclasses.replace(
        cfg, filters=dataclasses.replace(cfg.filters, junction_window=window))
    run_pipeline(cfg)
    return seen


def test_multi_strain_scoring_uses_the_configured_window(monkeypatch, tmp_path):
    import flexprobe.pipeline as pl
    seen = _run_scenario_capturing(monkeypatch, tmp_path, "multi_strain_capped",
                                   pl, "score_probe_against_all_strains", 6)
    assert seen and set(seen) == {6}


def test_specificity_uses_the_configured_window(monkeypatch, tmp_path):
    import flexprobe.specificity as sp
    # hit_junction_mismatches only runs on a qualifying off-target hit, which
    # this scenario has none of, so the spy sits on evaluate_specificity, which
    # every screened candidate passes through with the run's spec_kwargs.
    seen = _run_scenario_capturing(monkeypatch, tmp_path, "kmer_screened",
                                   sp, "evaluate_specificity", 6)
    assert seen and set(seen) == {6}


def test_evaluate_specificity_forwards_the_window(monkeypatch):
    import flexprobe.specificity as sp
    import inspect
    assert "junction_window" in inspect.signature(sp.evaluate_specificity).parameters


def test_an_impossible_window_is_refused():
    from flexprobe.models import FilterThresholds
    for bad in (0, 26):
        with pytest.raises(ConfigError, match="junction_window"):
            _cfg(filters=FilterThresholds(junction_window=bad)).validate()


# -- --no-idt / --no-cellranger were ignored by --auto-select-per-gene --

@pytest.mark.parametrize("flag,absent,present", [
    ("--no-idt", "idt_opools.csv", "cellranger_probe_set.csv"),
    ("--no-cellranger", "cellranger_probe_set.csv", "idt_opools.csv"),
])
def test_auto_select_honours_the_skip_flags(tmp_path, flag, absent, present):
    import os
    from conftest import SYNTH_200, write_fasta
    from flexprobe import cli
    fa = write_fasta(tmp_path / "g.fa", [("G1", SYNTH_200), ("G2", SYNTH_200[::-1])])
    out = tmp_path / "out"
    rc = cli.main(["design", "--fasta", str(fa), "--out", str(out),
                   "--skip-screen", "--auto-select-per-gene", "1", flag,
                   "--quiet", "--force"])
    assert rc == 0
    assert os.path.exists(out / "selections.json"), "auto-select did not run"
    assert not os.path.exists(out / absent)
    assert os.path.exists(out / present)


# -- help said BLAST "warns" under the junction gate; validate refuses --

def _blast(**kw):
    from flexprobe.models import ScreenResources
    from flexprobe.screen import ENGINE_BLAST
    return _cfg(screen=ScreenResources(engine=ENGINE_BLAST, db_paths=["db"], **kw))


def test_blast_with_the_gate_on_is_refused_with_the_way_out():
    with pytest.raises(ConfigError) as exc:
        _blast().validate()
    assert "--no-off-target-junction-duplex" in str(exc.value)


def test_blast_with_skip_screen_is_not_refused():
    """No screen runs, so no gate is claimed: nothing to refuse."""
    _blast(skip_screen=True).validate()


def test_help_says_refused_not_warned():
    for dest in ("screen_engine", "off_target_junction_duplex"):
        text = " ".join(_help_of(dest).split())
        assert "refused" in text, dest
        assert "a warning is printed" not in text, dest



# -- the spacing ladder, described as it runs --

def test_relaxed_gap_help_no_longer_claims_to_disable_fallback():
    for dest, cmd in (("relaxed_probe_gap", "design"), ("relaxed_gap", "select-coverage")):
        text = " ".join(_help_of(dest, cmd).split())
        assert "disable" not in text, (cmd, text)
        assert "this rung only" in text, (cmd, text)


def test_the_ladder_really_keeps_the_floor_rung():
    """The fact the help now states, checked against the code."""
    from flexprobe.auto_select import ladder_from
    assert [g for g, _ in ladder_from(150, 150, 15)] == [150, 15]
    assert [g for g, _ in ladder_from(150, 150, 150)] == [150]
    assert [g for g, _ in ladder_from(150, 40, 60)] == [150, 40]


def test_floor_help_no_longer_promises_an_absolute_floor():
    for dest, cmd in (("floor_probe_gap", "design"), ("floor_gap", "select-coverage")):
        text = " ".join(_help_of(dest, cmd).split())
        assert "never" not in text, (cmd, text)



def test_spacing_help():
    assert "--auto-select-per-gene" in _help_of("min_probe_gap")
    assert "isoform_hit_ids" in _help_of("max_mismatches", "select-coverage")
    assert "tie-break" in _help_of("preferred_gap", "select-coverage")
    assert "min(" in _help_of("relax_below", "select-coverage")


# -- export-hpc resource defaults the code never used --

def _export(tmp_path, capsys, *extra):
    from conftest import SYNTH_200, write_fasta
    from flexprobe import cli
    fa = write_fasta(tmp_path / "g.fa", [("G1", SYNTH_200)])
    rc = cli.main(["export-hpc", "--fasta", str(fa), "--skip-screen",
                   "--image", "f.sif", "--dest", str(tmp_path / "b"), *extra])
    assert rc == 0
    return capsys.readouterr().out


def test_split_screen_job_defaults_to_the_screen_memory(tmp_path, capsys):
    from flexprobe.hpc_export import DEFAULT_SCREEN_MEM_GB
    out = _export(tmp_path, capsys, "--split-prepare-screen")
    assert f"screen:  walltime=24:00:00 mem={DEFAULT_SCREEN_MEM_GB}gb" in out


def test_single_job_keeps_its_default_and_an_explicit_value_wins(tmp_path, capsys):
    from flexprobe.hpc_export import DEFAULT_MEM_GB
    out = _export(tmp_path, capsys)
    assert f"mem={DEFAULT_MEM_GB}gb" in out
    out = _export(tmp_path, capsys, "--split-prepare-screen", "--mem-gb", "300")
    assert "screen:  walltime=24:00:00 mem=300gb" in out


def test_prepare_defaults_to_eight_cores(tmp_path, capsys):
    out = _export(tmp_path, capsys, "--split-prepare-screen")
    assert "prepare: walltime=06:00:00 mem=96gb cores=8" in out


# -- validate-blast checked three sidecars and nothing the pipeline needs --

def _blast_dir():
    import os
    import platform
    import sys
    plat = {"darwin": "macos-arm64", "linux": "linux-x64"}.get(sys.platform)
    if plat is None or (plat == "macos-arm64" and platform.machine() != "arm64"):
        return None
    from pathlib import Path
    d = Path(__file__).resolve().parents[1] / "packaging" / "vendor" / "blast" / plat / "bin"
    return str(d) if (d / "blastdbcmd").exists() and os.access(d / "blastdbcmd", os.X_OK) else None


def _make_db(tmp_path, bdir, name, parse_seqids):
    import subprocess
    fa = tmp_path / f"{name}.fa"
    fa.write_text(">acc1 one\n" + "ACGT" * 40 + "\n>acc2 two\n" + "TTGCA" * 40 + "\n")
    cmd = [f"{bdir}/makeblastdb", "-in", str(fa), "-dbtype", "nucl",
           "-out", str(tmp_path / name)]
    if parse_seqids:
        cmd.append("-parse_seqids")
    subprocess.run(cmd, check=True, capture_output=True)
    return str(tmp_path / name)


def _validate(tmp_path, monkeypatch, bdir, db):
    import json
    from flexprobe import cli
    monkeypatch.setenv("FLEXPROBE_BLAST_DIR", bdir)
    out = tmp_path / "vb"
    rc = cli.main(["validate-blast", "--db", db, "--out", str(out)])
    return rc, json.loads((out / "blast_validation.json").read_text())


def test_a_db_without_parse_seqids_fails_validation(tmp_path, monkeypatch):
    bdir = _blast_dir()
    if bdir is None:
        pytest.skip("no vendored BLAST for this platform")
    rc, report = _validate(tmp_path, monkeypatch, bdir,
                           _make_db(tmp_path, bdir, "plain", False))
    assert rc == 1 and report["retrieval"] is False


def test_a_db_with_parse_seqids_passes(tmp_path, monkeypatch):
    bdir = _blast_dir()
    if bdir is None:
        pytest.skip("no vendored BLAST for this platform")
    rc, report = _validate(tmp_path, monkeypatch, bdir,
                           _make_db(tmp_path, bdir, "ok", True))
    assert rc == 0 and report["retrieval"] is True
    assert report["tools"].get("blastdbcmd")


def test_a_multi_volume_alias_counts_as_present(tmp_path):
    from flexprobe import cli
    (tmp_path / "big.nal").write_text("DBLIST big.00 big.01\n")
    assert cli._blast_db_layout(str(tmp_path / "big")) == "multi-volume"
    assert cli._blast_db_layout(str(tmp_path / "nothing")) is None


# -- --max-target-seqs does nothing under the (default) k-mer engine --

def test_max_target_seqs_help_says_blast_only():
    assert "--screen-engine blast only" in _help_of("max_target_seqs")


def test_max_target_seqs_under_kmer_warns():
    from flexprobe import cli
    p = cli.build_parser()
    a = p.parse_args(["design", "--fasta", "a.fa", "--max-target-seqs", "50"])
    assert any("--max-target-seqs" in w for w in cli._unused_flag_warnings(a))
    a = p.parse_args(["design", "--fasta", "a.fa", "--max-target-seqs", "50",
                      "--screen-engine", "blast",
                      "--no-off-target-junction-duplex"])
    assert cli._unused_flag_warnings(a) == []
    a = p.parse_args(["design", "--fasta", "a.fa"])
    assert cli._unused_flag_warnings(a) == []
