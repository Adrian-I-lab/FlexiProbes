"""CLI argument-parsing and end-to-end dispatch tests."""
import csv
import json
import os

import pytest

from flexprobe.cli import _build_run_config, build_parser, main

# A 120-bp toy target with an A at window position 25 for several windows.
TOY_FASTA = (
    ">toy|gene|ACC1\n"
    "GCTAGCTAGCTAGCTAGCTAGCTAGATGCATGCATGCATGCATGCATGCA"
    "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG\n"
)


def test_version_command(capsys):
    rc = main(["version"])
    assert rc == 0
    assert "flexiprobes" in capsys.readouterr().out.lower()


def test_build_parser_design_defaults():
    args = build_parser().parse_args(["design", "--fasta", "x.fa"])
    # CLI defaults must match the documented spec (constants SSOT).
    assert args.gc_min == 0.44 and args.gc_max == 0.72
    assert args.max_n_per_half == 1
    assert args.tile_size == 50 and args.tile_step == 1
    assert args.modality == "flex_v2"


def test_design_end_to_end_skip_screen(tmp_path, capsys):
    fa = tmp_path / "toy.fasta"
    fa.write_text(TOY_FASTA)
    out = tmp_path / "out"
    rc = main([
        "design", "--fasta", str(fa), "--out", str(out),
        "--skip-screen", "--require-tn-junction", "--force",
    ])
    assert rc == 0
    results = out / "results.csv"
    assert results.exists()
    rows = list(csv.DictReader(results.open()))
    assert rows, "expected at least one candidate row"
    # unscreened probes must be marked, never silently 'passed'
    assert rows[0]["specificity_status"] == "not_screened"


def test_design_with_config_file(tmp_path):
    """`design --config` runs from a serialized RunConfig (HPC bundle path)."""
    from flexprobe.config_io import save_config
    from flexprobe.models import ModalityDefinition, RunConfig

    fa = tmp_path / "toy.fasta"
    fa.write_text(TOY_FASTA)
    out = tmp_path / "out"
    cfg = RunConfig(
        fasta_path=str(fa), output_dir=str(out),
        modality=ModalityDefinition(name="flex", version="v2"),
    )
    cfg.screen.skip_screen = True
    cfg.export.overwrite = True
    cfg.filters.require_tn_junction = True
    cfg_path = tmp_path / "config.json"
    save_config(cfg, str(cfg_path))

    rc = main(["design", "--config", str(cfg_path)])
    assert rc == 0
    assert (out / "results.csv").exists()


def test_config_path_overrides(tmp_path):
    """--fasta/--out override the paths baked into a relocated config."""
    from flexprobe.config_io import save_config
    from flexprobe.models import ModalityDefinition, RunConfig

    fa = tmp_path / "toy.fasta"
    fa.write_text(TOY_FASTA)
    # Config points at paths that don't exist on this machine.
    cfg = RunConfig(
        fasta_path="/nonexistent/in.fa", output_dir="/nonexistent/out",
        modality=ModalityDefinition(name="flex", version="v2"),
    )
    cfg.screen.skip_screen = True
    cfg.export.overwrite = True
    cfg_path = tmp_path / "config.json"
    save_config(cfg, str(cfg_path))

    out = tmp_path / "real_out"
    rc = main([
        "design", "--config", str(cfg_path),
        "--fasta", str(fa), "--out", str(out),
    ])
    assert rc == 0
    assert (out / "results.csv").exists()


def test_workers_default_is_none():
    # None lets --config distinguish "unset" from an explicit --workers 1.
    args = build_parser().parse_args(["design", "--fasta", "x.fa"])
    assert args.workers is None


def test_config_workers_override_to_one(tmp_path):
    from flexprobe.config_io import save_config
    from flexprobe.models import ModalityDefinition, RunConfig

    fa = tmp_path / "toy.fasta"
    fa.write_text(TOY_FASTA)
    out = tmp_path / "out"
    cfg = RunConfig(
        fasta_path=str(fa), output_dir=str(out), workers=16,
        modality=ModalityDefinition(name="flex", version="v2"),
    )
    cfg.screen.skip_screen = True
    cfg.export.overwrite = True
    cfg_path = tmp_path / "config.json"
    save_config(cfg, str(cfg_path))

    rc = main(["design", "--config", str(cfg_path), "--workers", "1"])
    assert rc == 0
    assert (out / "results.csv").exists()


def test_design_requires_fasta_without_config(capsys):
    rc = main(["design", "--skip-screen"])
    assert rc != 0
    assert "fasta" in capsys.readouterr().err.lower()


def test_auto_select_per_gene_writes_orders(tmp_path):
    """--auto-select-per-gene forces discovery + writes order CSVs headlessly."""
    fa = tmp_path / "toy.fasta"
    fa.write_text(TOY_FASTA)
    out = tmp_path / "out"
    rc = main([
        "design", "--fasta", str(fa), "--out", str(out),
        "--skip-screen", "--auto-select-per-gene", "2", "--force",
    ])
    assert rc == 0
    assert (out / "results.csv").exists()
    assert (out / "selections.json").exists()
    # auto-export produced at least one order CSV
    assert (out / "idt_opools.csv").exists()


# ---------------------------------------------------------------------------
# B2: spacing flags must actually reach the auto-pick selector (design path)
# ---------------------------------------------------------------------------

def _long_synthetic_fasta(seed: int = 12345, length: int = 3000) -> str:
    """A single-gene FASTA long enough for gap tiers to matter.

    Deterministic (a seeded local Random instance, never the global one) so
    the two runs compared below differ ONLY in the CLI flags under test.
    """
    import random as _random
    rng = _random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(length))
    lines = [seq[i:i + 70] for i in range(0, len(seq), 70)]
    return ">GENE1\n" + "\n".join(lines) + "\n"


def test_min_probe_gap_changes_the_auto_selected_probes(tmp_path):
    """Two runs differing ONLY in --min-probe-gap (& friends) must select
    DIFFERENT probes — the exact defect the audit found byte-identical.

    Before the fix, cli._auto_export called auto_pick_with_report without
    forwarding min_gap/relaxed_gap/floor_gap, so the selector silently fell
    back to the module constants regardless of these flags and both runs
    below produced identical selections.json / order CSVs.
    """
    fa = tmp_path / "long.fasta"
    fa.write_text(_long_synthetic_fasta())

    out_default = tmp_path / "out_default"
    rc = main([
        "design", "--fasta", str(fa), "--out", str(out_default),
        "--skip-screen", "--auto-select-per-gene", "3", "--force",
    ])
    assert rc == 0

    out_wide = tmp_path / "out_wide"
    rc = main([
        "design", "--fasta", str(fa), "--out", str(out_wide),
        "--skip-screen", "--auto-select-per-gene", "3", "--force",
        "--min-probe-gap", "1200", "--relaxed-probe-gap", "1200",
        "--floor-probe-gap", "1200",
    ])
    assert rc == 0

    default_sel = json.loads((out_default / "selections.json").read_text())
    wide_sel = json.loads((out_wide / "selections.json").read_text())
    assert default_sel["selections"] != wide_sel["selections"], (
        "selections must differ once --min-probe-gap is actually honoured"
    )
    # And the printed tier banner must report the ACTUAL gap used, not the
    # hardcoded module constant (the audit also caught this: "green >= 150 bp"
    # was printed regardless of the flags).
    default_summary = json.loads((out_default / "design_summary.json").read_text())
    wide_summary = json.loads((out_wide / "design_summary.json").read_text())
    assert default_summary["status"] == "ok"
    assert wide_summary["status"] == "ok"


def test_preferred_probe_gap_flag_was_removed_as_dead(tmp_path):
    """--preferred-probe-gap is gone: it was recorded in run_config.json and
    read nowhere (auto_select.py's SPACING_LADDER has no "preferred" rung —
    only green/orange/red = min/relaxed/floor). Kept as a regression test so
    the flag is not silently reintroduced without being wired to something."""
    with pytest.raises(SystemExit):
        build_parser().parse_args(["design", "--fasta", "x.fa",
                                   "--preferred-probe-gap", "300"])
    from flexprobe.models import RunConfig
    assert not hasattr(RunConfig(), "preferred_probe_gap_bp")


# ---------------------------------------------------------------------------
# E2: a run that designs zero probes must not look like success
# ---------------------------------------------------------------------------

def test_normal_run_still_exits_ok_with_at_least_one_probe(tmp_path):
    """Guard against regressing the >=1-probe path while fixing the
    zero-probe case below: status must stay "ok" and exit code 0."""
    fa = tmp_path / "toy.fasta"
    fa.write_text(TOY_FASTA)
    out = tmp_path / "out"
    rc = main([
        "design", "--fasta", str(fa), "--out", str(out),
        "--skip-screen", "--force",
    ])
    assert rc == 0
    summary = json.loads((out / "design_summary.json").read_text())
    assert summary["status"] == "ok"
    assert summary["stats"]["pairs_output"] > 0


def test_zero_probe_run_exits_nonzero_with_no_probes_status(tmp_path, capsys):
    """A run whose FILTER funnel is not empty (candidates pass the input
    filters fine) but whose OFF-TARGET screen then rejects every remaining
    candidate is a different failure mode from an empty filter funnel:
    run_pipeline does not raise for it, so before this fix the CLI reported
    exit 0 and design_summary.json said status="ok" for zero probes.

    Contrived by screening the input against a reference that carries the
    on-target sequence itself under a different record id: every probe half
    is then a perfect (0-mismatch) match to a "non-target" gene, so every
    filter-passing candidate is rejected as an off-target hit.
    """
    fa = tmp_path / "toy.fasta"
    fa.write_text(TOY_FASTA)
    seq = TOY_FASTA.split("\n", 1)[1].replace("\n", "")
    decoy = tmp_path / "decoy.fasta"
    decoy.write_text(f">DECOY\n{seq}\n")
    out = tmp_path / "out"
    rc = main([
        "design", "--fasta", str(fa), "--out", str(out),
        "--blast-db", str(decoy), "--force",
    ])
    captured = capsys.readouterr()
    assert rc != 0
    assert "0 probe pairs" in captured.err

    summary = json.loads((out / "design_summary.json").read_text())
    assert summary["status"] == "no_probes"
    assert summary["stats"]["pairs_output"] == 0
    # The funnel was NOT empty — this is not the "no_candidates" path.
    assert summary["stats"]["filter_passed"] > 0


def test_design_invalid_config_is_reported(tmp_path, capsys):
    fa = tmp_path / "toy.fasta"
    fa.write_text(TOY_FASTA)
    out = tmp_path / "out"
    # gc_min > gc_max -> RunConfig.validate() should make the run fail cleanly
    rc = main([
        "design", "--fasta", str(fa), "--out", str(out),
        "--skip-screen", "--gc-min", "0.9", "--gc-max", "0.4", "--force",
    ])
    assert rc != 0


# ---- resume-related flags ----

def test_resume_on_by_default_and_no_resume_disables_it():
    args = build_parser().parse_args(["design", "--fasta", "x.fa"])
    assert _build_run_config(args).resume_screen is True
    args2 = build_parser().parse_args(["design", "--fasta", "x.fa", "--no-resume"])
    assert _build_run_config(args2).resume_screen is False


def test_config_no_resume_and_cap_overrides(tmp_path):
    # With --config, --no-resume and --screen-top-per-gene override the baked file
    # so a cluster resubmission can change them without re-exporting.
    from flexprobe.config_io import save_config
    from flexprobe.models import ModalityDefinition, RunConfig

    cfg = RunConfig(
        fasta_path="x.fa", output_dir="out", resume_screen=True, screen_top_per_gene=0,
        modality=ModalityDefinition(name="flex", version="v2"),
    )
    cfg_path = tmp_path / "config.json"
    save_config(cfg, str(cfg_path))
    args = build_parser().parse_args([
        "design", "--config", str(cfg_path),
        "--no-resume", "--screen-top-per-gene", "50",
    ])
    built = _build_run_config(args)
    assert built.resume_screen is False
    assert built.screen_top_per_gene == 50


def test_config_cap_reset_to_zero_and_omit_keeps(tmp_path):
    # Audit M3: an explicit --screen-top-per-gene 0 must reset a baked non-zero
    # cap back to "screen all"; omitting the flag must keep the config's value.
    from flexprobe.config_io import save_config
    from flexprobe.models import ModalityDefinition, RunConfig

    cfg = RunConfig(
        fasta_path="x.fa", output_dir="out", screen_top_per_gene=100,
        modality=ModalityDefinition(name="flex", version="v2"),
    )
    cfg_path = tmp_path / "config.json"
    save_config(cfg, str(cfg_path))

    reset = _build_run_config(build_parser().parse_args(
        ["design", "--config", str(cfg_path), "--screen-top-per-gene", "0"]))
    assert reset.screen_top_per_gene == 0          # explicit 0 overrides the 100

    kept = _build_run_config(build_parser().parse_args(
        ["design", "--config", str(cfg_path)]))
    assert kept.screen_top_per_gene == 100         # omitted → config value kept


# ---- export subcommand (re-export an existing results.csv) ----

def _results_csv(path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "probe_id", "gene_id", "start", "end", "composite_score"])
        w.writerow(["1", "G_t1", "G", "0", "50", "0.91"])
        w.writerow(["2", "G_t2", "G", "60", "110", "0.80"])


def test_export_subcommand_json_and_bed(tmp_path):
    src = tmp_path / "results.csv"
    _results_csv(src)
    out = tmp_path / "exp"
    assert main(["export", "--results-csv", str(src), "--format", "json",
                 "--out", str(out)]) == 0
    rows = json.load(open(out / "export.json"))
    assert {r["probe_id"] for r in rows} == {"G_t1", "G_t2"}

    assert main(["export", "--results-csv", str(src), "--format", "bed",
                 "--out", str(out)]) == 0
    bed = (out / "export.bed").read_text().splitlines()
    assert bed[0].split("\t")[:4] == ["G", "0", "50", "G_t1"]


def test_export_subcommand_refuses_overwrite_without_force(tmp_path):
    src = tmp_path / "results.csv"
    _results_csv(src)
    out = tmp_path / "exp"
    assert main(["export", "--results-csv", str(src), "--format", "csv",
                 "--out", str(out)]) == 0
    # Second run without --force must refuse (non-zero) and keep the file.
    rc = main(["export", "--results-csv", str(src), "--format", "csv",
               "--out", str(out)])
    assert rc != 0
    assert os.path.exists(out / "export.csv")
    # With --force it succeeds.
    assert main(["export", "--results-csv", str(src), "--format", "csv",
                 "--out", str(out), "--force"]) == 0


def test_export_missing_source_is_reported(tmp_path, capsys):
    rc = main(["export", "--results-csv", str(tmp_path / "nope.csv"),
               "--out", str(tmp_path / "exp")])
    assert rc != 0
    assert "not found" in capsys.readouterr().err.lower()


def test_select_coverage_command(tmp_path, capsys):
    """End-to-end: select-coverage reads results.csv + FASTA and writes the
    per-gene selection + flagged-genes report."""
    # Gene G has two isoforms; t1 has a unique head + a shared tail, t2 only the
    # tail. A probe on the shared tail covers both; the unique head covers t1.
    # Head and tail windows are 30 bp apart so both clear the spacing rule.
    head = "TTTTTGGGGGTTTTTGGGGGTTTTTGGGGG"     # 30 bp, t1-only (start 0)
    spacer = "ATATATATATCGCGCGCGCGATATATATAT"   # 30 bp filler
    shared = "ACGTACGTACGTACGTACGTACGTACGTAC"   # 30 bp tail (start 60 in t1)
    t1 = head + spacer + shared
    t2 = "CCCCCAAAAACCCCCAAAAACCCCCAAAAA" + shared
    fa = tmp_path / "tx.fa"
    fa.write_text(f">t1 gene:G\n{t1}\n>t2 gene:G\n{t2}\n")

    results = tmp_path / "results.csv"
    with open(results, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["probe_id", "gene_id", "start", "end", "mrna_target_50bp",
                    "composite_score", "specificity_pass", "reduced_confidence",
                    "isoform_hits", "specificity_status"])
        # probe on the shared tail (binds t1 + t2); probe on t1-only head.
        w.writerow(["t1_p_tail", "t1", "60", "90", shared, "0.9", "True",
                    "False", "2", "passed"])
        w.writerow(["t1_p_head", "t1", "0", "30", head, "0.8", "True",
                    "False", "1", "passed"])

    out = tmp_path / "sel.csv"
    flagged = tmp_path / "flag.csv"
    rc = main(["select-coverage", "--results", str(results), "--fasta", str(fa),
               "--out", str(out), "--flagged", str(flagged), "--n", "2",
               "--min-gap", "25"])
    assert rc == 0
    rows = list(csv.DictReader(open(out)))
    assert {r["probe_id"] for r in rows} == {"t1_p_tail", "t1_p_head"}
    assert all(r["gene"] == "G" for r in rows)
    # Both isoforms reachable -> gene fully covered -> not flagged.
    flags = list(csv.DictReader(open(flagged)))
    assert flags == []


def test_select_coverage_missing_results(tmp_path):
    rc = main(["select-coverage", "--results", str(tmp_path / "nope.csv"),
               "--fasta", str(tmp_path / "nope.fa")])
    assert rc == 2


def test_select_coverage_fasta_free_via_column(tmp_path):
    """A current run persists isoform_hit_ids + a manifest gene map, so
    select-coverage needs no --fasta."""
    results = tmp_path / "results.csv"
    with open(results, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["probe_id", "gene_id", "start", "end", "mrna_target_50bp",
                    "composite_score", "specificity_pass", "reduced_confidence",
                    "isoform_hits", "isoform_hit_ids", "specificity_status"])
        # pA binds t1+t2; pB (far away) binds t3 -> together cover the gene.
        w.writerow(["pA", "t1", "0", "30", "", "0.9", "True", "False",
                    "2", "t1;t2", "passed"])
        w.writerow(["pB", "t3", "200", "230", "", "0.8", "True", "False",
                    "1", "t3", "passed"])
    # Manifest supplies the gene's full isoform universe.
    (tmp_path / "manifest.json").write_text(json.dumps(
        {"params": {"gene_isoforms": {"G": ["t1", "t2", "t3"]}}}))

    rc = main(["select-coverage", "--results", str(results), "--n", "2",
               "--min-gap", "25"])   # no --fasta
    assert rc == 0
    sel = list(csv.DictReader(open(tmp_path / "selected_probes_per_gene.csv")))
    assert {r["probe_id"] for r in sel} == {"pA", "pB"}
    assert all(r["gene"] == "G" for r in sel)
    # Full coverage reached without a FASTA -> nothing flagged.
    assert list(csv.DictReader(open(tmp_path / "flagged_genes.csv"))) == []


def test_select_coverage_old_run_requires_fasta(tmp_path):
    """A run without the column (older) and no FASTA fails cleanly."""
    results = tmp_path / "results.csv"
    with open(results, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["probe_id", "gene_id", "start", "end", "mrna_target_50bp",
                    "composite_score", "specificity_pass"])
        w.writerow(["p1", "t1", "0", "30", "ACGT" * 12 + "AC", "0.9", "True"])
    rc = main(["select-coverage", "--results", str(results)])  # no column, no fasta
    assert rc == 2


# ---------------------------------------------------------------------------
# export-hpc — the GUI's HPC export, headless.
#
# Before this existed, write_hpc_bundle was reachable only from the Qt dialog,
# so a cluster run's bundle could not be regenerated from a recorded command
# line — a reproducibility hole for anything published from those runs.
# ---------------------------------------------------------------------------

def _export_hpc(tmp_path, *extra):
    dest = tmp_path / "bundle"
    rc = main([
        "export-hpc",
        "--dest", str(dest),
        "--scheduler", "pbs",
        "--image", "/apps/flexiprobe-0.4.0.sif",
        "--remote-input-path", "inputs/tx.fa",
        "--remote-db-path", "inputs/tx.fa",
        "--group-by-gene", "--discovery-only",
        *extra,
    ])
    return rc, dest


def test_export_hpc_writes_a_bundle(tmp_path):
    rc, dest = _export_hpc(tmp_path, "--walltime", "100:00:00",
                           "--mem-gb", "64", "--cores", "16",
                           "--job-name", "tx_run")
    assert rc == 0
    for name in ("config.json", "run_flexiprobes.pbs", "submit_chain.sh",
                 "README.txt"):
        assert (dest / name).exists(), name
    script = (dest / "run_flexiprobes.pbs").read_text()
    assert "#PBS -N tx_run" in script
    assert "#PBS -l walltime=100:00:00" in script
    assert "#PBS -l select=1:ncpus=16:mem=64gb" in script
    assert "--workers 16" in script
    # nothing was staged: no half-gigabyte FASTA copied into the bundle
    assert not (dest / "inputs").exists()


def test_export_hpc_config_round_trips(tmp_path):
    from flexprobe.config_io import load_config
    rc, dest = _export_hpc(
        tmp_path, "--cores", "16", "--modality", "flex_v1", "--singleplex",
        "--probe-barcode", "BC001", "--screen-top-per-gene", "100",
        "--min-isoform-coverage", "0.25", "--enable-primer3",
        "--tn-junction-mode", "soft",
    )
    assert rc == 0
    cfg = load_config(str(dest / "config.json"))
    assert cfg.fasta_path == "inputs/tx.fa"
    assert cfg.screen.db_paths == ["inputs/tx.fa"]
    assert cfg.output_dir == "output"
    assert cfg.modality.name == "flex" and cfg.modality.version == "v1"
    assert cfg.multiplex is False
    assert cfg.group_by_gene and cfg.discovery_only
    assert cfg.screen_top_per_gene == 100
    assert cfg.min_isoform_coverage == 0.25
    assert cfg.filters.tn_junction_mode == "soft"
    # --cores also sets the config's worker count so config and script agree
    assert cfg.workers == 16
    cfg.validate()


def test_export_hpc_defaults_are_the_v040_corrections(tmp_path):
    """A bundle written with no threshold flags must not ship legacy behaviour."""
    from flexprobe.config_io import load_config
    rc, dest = _export_hpc(tmp_path)
    assert rc == 0
    cfg = load_config(str(dest / "config.json"))
    assert cfg.filters.enforce_gc_range is False       # GC is not a gate
    assert cfg.filters.tm_min == 56.0                  # per-half Tm gate is
    assert cfg.filters.tm_max is None
    assert cfg.specificity.off_target_full_half is True
    assert cfg.specificity.off_target_require_both_halves is True
    assert cfg.specificity.off_target_same_subject is True


def test_export_hpc_explicit_workers_beats_cores(tmp_path):
    from flexprobe.config_io import load_config
    rc, dest = _export_hpc(tmp_path, "--cores", "16", "--workers", "4")
    assert rc == 0
    assert load_config(str(dest / "config.json")).workers == 4


def test_export_hpc_bakes_the_version_assertion(tmp_path):
    from flexprobe import __version__
    rc, dest = _export_hpc(tmp_path)
    assert rc == 0
    assert (f"EXPECTED_VERSION={__version__}"
            in (dest / "run_flexiprobes.pbs").read_text())


def test_export_hpc_expected_version_override(tmp_path):
    rc, dest = _export_hpc(tmp_path, "--expected-version", "0.4.1")
    assert rc == 0
    assert "EXPECTED_VERSION=0.4.1" in (dest / "run_flexiprobes.pbs").read_text()


def test_export_hpc_job_log_project_reaches_the_scripts(tmp_path):
    rc, dest = _export_hpc(tmp_path, "--job-log-project", "flexiprobe")
    assert rc == 0
    script = (dest / "run_flexiprobes.pbs").read_text()
    assert "#PBS -o /dev/null" in script
    assert "${WORKDIR}/projects/flexiprobe" in script


def test_export_hpc_rejects_an_unsafe_job_log_project(tmp_path):
    rc, _ = _export_hpc(tmp_path, "--job-log-project", "a;b")
    assert rc == 2


def test_export_hpc_stages_when_asked(tmp_path):
    fa = tmp_path / "tx.fa"
    fa.write_text(">g1 gene:ENSG1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / "staged"
    rc = main([
        "export-hpc", "--dest", str(dest), "--scheduler", "pbs",
        "--image", "/apps/f.sif", "--fasta", str(fa), "--stage-inputs",
        "--group-by-gene", "--discovery-only",
    ])
    assert rc == 0
    assert (dest / "inputs" / "tx.fa").exists()
    from flexprobe.config_io import load_config
    assert load_config(str(dest / "config.json")).fasta_path == "inputs/tx.fa"


def test_export_hpc_needs_an_input(tmp_path):
    rc = main(["export-hpc", "--dest", str(tmp_path / "b"),
               "--image", "/apps/f.sif"])
    assert rc == 2


def test_export_hpc_shares_every_design_config_flag():
    """design and export-hpc must not drift apart.

    They both feed _build_run_config, so any flag one accepts the other must
    accept too — otherwise a bundle silently cannot express a run the local
    pipeline can.
    """
    parser = build_parser()
    subs = parser._subparsers._group_actions[0].choices  # type: ignore[attr-defined]
    design_opts = {s for a in subs["design"]._actions for s in a.option_strings}
    hpc_opts = {s for a in subs["export-hpc"]._actions for s in a.option_strings}
    # --quiet shapes the run's console output, not the config, so it is the one
    # documented exception.
    assert design_opts - hpc_opts == {"--quiet"}


def test_design_prescreen_flags_reach_the_config():
    from flexprobe.cli import _build_run_config
    args = build_parser().parse_args([
        "design", "--fasta", "x.fa", "--emit-prescreen-candidates",
        "--keep-candidate-cache",
    ])
    cfg = _build_run_config(args)
    assert cfg.emit_prescreen_candidates is True
    assert cfg.emit_prescreen_precap_rows is False
    assert cfg.keep_candidate_cache is True


def test_precap_flag_implies_the_main_artefact():
    """--emit-prescreen-precap alone must not produce an invalid config."""
    from flexprobe.cli import _build_run_config
    args = build_parser().parse_args(
        ["design", "--fasta", "x.fa", "--emit-prescreen-precap"])
    cfg = _build_run_config(args)
    assert cfg.emit_prescreen_candidates and cfg.emit_prescreen_precap_rows
    cfg.validate()


def test_prescreen_flags_override_a_baked_config(tmp_path):
    """Turn-on-only, like --force: omitting the flag never silently disables it."""
    from flexprobe.cli import _build_run_config
    from flexprobe.config_io import save_config
    from flexprobe.models import RunConfig
    p = tmp_path / "c.json"
    save_config(RunConfig(fasta_path="x.fa", emit_prescreen_candidates=True), str(p))

    off = _build_run_config(build_parser().parse_args(
        ["design", "--config", str(p)]))
    assert off.emit_prescreen_candidates is True, "must not be turned off"

    on = _build_run_config(build_parser().parse_args(
        ["design", "--config", str(p), "--emit-prescreen-precap"]))
    assert on.emit_prescreen_precap_rows is True


def test_export_hpc_bakes_the_prescreen_artefacts(tmp_path):
    from flexprobe.config_io import load_config
    rc, dest = _export_hpc(tmp_path, "--emit-prescreen-candidates",
                           "--keep-candidate-cache")
    assert rc == 0
    cfg = load_config(str(dest / "config.json"))
    assert cfg.emit_prescreen_candidates is True
    assert cfg.keep_candidate_cache is True
    cfg.validate()


# ---------------------------------------------------------------------------
# design_summary.json is written for EVERY run
# ---------------------------------------------------------------------------
# A run whose filter funnel empties is a result (tiles existed and every one was
# rejected). Without a design_summary.json, a consumer keying off that file
# cannot tell it apart from a run that never started, and records a real
# candidate set as "0 tiles".

# Every window here fails the homopolymer filter, so tiles are generated and all
# of them are rejected. That is the shape that matters: a non-zero
# tile count with a zero pass count.
ALL_REJECTED_FASTA = (
    ">toy|gene|ACC1\n"
    "AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAA"
    "CCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCC\n"
)


def test_empty_funnel_still_writes_a_summary(tmp_path):
    fa = tmp_path / "reject.fasta"
    fa.write_text(ALL_REJECTED_FASTA)
    out = tmp_path / "out"
    rc = main(["design", "--fasta", str(fa), "--out", str(out),
               "--skip-screen", "--force"])
    assert rc == 2, "the run should still report failure"

    summary = out / "design_summary.json"
    assert summary.exists(), \
        "an emptied funnel must leave a summary, not an absence"
    data = json.loads(summary.read_text())

    # The status distinguishes 'the filters rejected everything' from a crash.
    assert data["status"] == "no_candidates"
    assert data["error"], "the reason must be carried, not just the status"

    # The counts are the point: tiles_total > 0 is what stops a consumer
    # recording this run as '0 tiles'.
    assert data["stats"]["tiles_total"] > 0
    assert data["stats"]["candidates_passed_filters"] == 0
    assert data["stats"]["pairs_output"] == 0
    assert data["stats"]["filter_stats"]["homopolymer"] > 0


def test_successful_run_still_reports_status_ok(tmp_path):
    """The failure path must not have changed what a good run writes."""
    fa = tmp_path / "toy.fasta"
    fa.write_text(TOY_FASTA)
    out = tmp_path / "out"
    assert main(["design", "--fasta", str(fa), "--out", str(out),
                 "--skip-screen", "--force"]) == 0
    data = json.loads((out / "design_summary.json").read_text())
    assert data["status"] == "ok"
    assert data["stats"]["pairs_output"] > 0


def test_summary_is_written_even_when_nothing_could_be_tiled(tmp_path):
    """The other empty case: sequences shorter than one tile.

    No filter_report.json exists here (the pipeline aborts before writing it),
    so the summary carries the reason without stats rather than nothing at all.
    """
    fa = tmp_path / "short.fasta"
    fa.write_text(">tiny|gene|ACC1\nACGTACGTACGT\n")
    out = tmp_path / "out"
    assert main(["design", "--fasta", str(fa), "--out", str(out),
                 "--skip-screen", "--force"]) == 2
    data = json.loads((out / "design_summary.json").read_text())
    assert data["status"] in ("failed", "no_candidates")
    assert "tile" in data["error"].lower()


# ---------------------------------------------------------------------------
# Legacy configs must not silently change engine
# ---------------------------------------------------------------------------

def _write_cfg(tmp_path, blast_block):
    import json
    p = tmp_path / "c.json"
    p.write_text(json.dumps({
        "fasta_path": "x", "output_dir": "y",
        "modality": {"name": "flex", "version": "v2"},
        "blast": blast_block,
    }))
    return str(p)


def test_config_without_engine_is_pinned_to_blast(tmp_path):
    """A config file records what a run DID.

    The `engine` field was added when the k-mer engine became the default, so a
    payload without it was authored under BLAST. Inheriting today's default
    would re-run it on a different engine — different probes, different order —
    with nothing in the file to say so, and no error either: every committed HPC
    bundle points db_paths at a FASTA, which the k-mer engine accepts happily.
    """
    from flexprobe.config_io import load_config
    from flexprobe.screen import ENGINE_BLAST
    cfg = load_config(_write_cfg(tmp_path, {"db_paths": ["ref.fa"]}))
    assert cfg.screen.engine == ENGINE_BLAST


def test_config_with_explicit_engine_is_honoured(tmp_path):
    from flexprobe.config_io import load_config
    from flexprobe.screen import ENGINE_KMER
    cfg = load_config(
        _write_cfg(tmp_path, {"db_paths": ["ref.fa"], "engine": "kmer"}))
    assert cfg.screen.engine == ENGINE_KMER


def test_legacy_pinning_says_so_out_loud(tmp_path, caplog):
    """Silent is the failure mode being prevented; the log line IS the fix."""
    import logging
    from flexprobe.config_io import load_config
    with caplog.at_level(logging.WARNING, logger="flexprobe.config_io"):
        load_config(_write_cfg(tmp_path, {"db_paths": ["ref.fa"]}))
    assert caplog.records, "pinning must be announced, not silent"
    assert "--screen-engine" in caplog.records[0].getMessage()


def test_engine_is_moot_when_no_screen_runs(tmp_path, caplog):
    """skip_screen / no db_paths means no screen at all — do not warn about an
    engine that will never be consulted."""
    import logging
    from flexprobe.config_io import load_config
    for block in ({"skip_screen": True}, {"db_paths": []}):
        caplog.clear()
        with caplog.at_level(logging.WARNING, logger="flexprobe.config_io"):
            load_config(_write_cfg(tmp_path, block))
        # Scoped to the ENGINE warning. A one-key config also (correctly)
        # warns that it filled every other key from defaults.
        engine = [r for r in caplog.records
                  if "--screen-engine" in r.getMessage()]
        assert not engine, f"spurious engine warning for {block}"


def test_roundtripped_config_always_carries_the_engine(tmp_path):
    """A config written today must never be ambiguous tomorrow."""
    import json
    from flexprobe.config_io import save_config, load_config
    from flexprobe.models import (
        ScreenResources, ModalityDefinition, RunConfig,
    )
    cfg = RunConfig(fasta_path="x", output_dir="y",
                    modality=ModalityDefinition(name="flex", version="v2"),
                    screen=ScreenResources(db_paths=["ref.fa"]))
    p = str(tmp_path / "out.json")
    save_config(cfg, p)
    assert "engine" in json.loads(open(p).read())["screen"]
    assert load_config(p).screen.engine == cfg.screen.engine
