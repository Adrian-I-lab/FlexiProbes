"""Tests for the HPC bundle writer and result importer."""
import json
import os
import tarfile
import zipfile

import pytest

from flexprobe.config_io import load_config
from flexprobe.hpc_export import (
    ABORT_SENTINEL, HpcResources, import_hpc_output, write_hpc_bundle,
)
from flexprobe.models import ScreenResources, RunConfig


def _config(fasta: str) -> RunConfig:
    return RunConfig(
        fasta_path=fasta,
        output_dir="/local/laptop/out",
        group_by_gene=True,
        screen=ScreenResources(db_paths=["/laptop/human_tx"]),
    )


def _resources() -> HpcResources:
    return HpcResources(
        job_name="wt", walltime="12:00:00", mem_gb=128, cores=32,
        partition="bio", account="acct1",
    )


def test_write_bundle_slurm(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">g1 gene:ENSG1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / "bundle"
    write_hpc_bundle(
        _config(str(fa)), scheduler="slurm", resources=_resources(),
        image="/apps/flexiprobe.sif", dest_dir=str(dest), stage_inputs=True,
    )
    script = (dest / "run_flexiprobes.slurm").read_text()
    assert script.startswith("#!/bin/bash")
    assert "#SBATCH --time=12:00:00" in script
    assert "#SBATCH --mem=128G" in script
    assert "#SBATCH --cpus-per-task=32" in script
    assert "#SBATCH --partition=bio" in script
    assert "#SBATCH --account=acct1" in script
    assert ("flexiprobes design --config config.json --workers 32 "
            "--screen-workers 32 --force") in script
    assert "/apps/flexiprobe.sif" in script
    assert "$SLURM_SUBMIT_DIR" in script
    # input was staged; config path rewritten to be bundle-relative
    assert (dest / "inputs" / "in.fa").exists()
    cfg = load_config(str(dest / "config.json"))
    assert cfg.fasta_path == os.path.join("inputs", "in.fa")
    assert cfg.output_dir == "output"


def test_write_bundle_pbs(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">g1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / "bundle"
    write_hpc_bundle(
        _config(str(fa)), scheduler="pbs", resources=_resources(),
        image="/apps/flexiprobe.sif", dest_dir=str(dest), stage_inputs=False,
    )
    script = (dest / "run_flexiprobes.pbs").read_text()
    assert "#PBS -N wt" in script
    assert "#PBS -l walltime=12:00:00" in script
    assert "#PBS -l select=1:ncpus=32:mem=128gb" in script
    assert "#PBS -q bio" in script
    assert "#PBS -A acct1" in script
    assert "$PBS_O_WORKDIR" in script
    assert (dest / "README.txt").exists()
    # not staged -> no inputs dir, fasta path left as-is
    assert not (dest / "inputs").exists()


def test_gene_aware_bundle_includes_selection_step(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">g1 gene:ENSG1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / "bundle"
    write_hpc_bundle(  # _config sets group_by_gene=True
        _config(str(fa)), scheduler="pbs", resources=_resources(),
        image="/apps/flexiprobe.sif", dest_dir=str(dest), stage_inputs=True,
    )
    script = (dest / "run_flexiprobes.pbs").read_text()
    # The coverage-selection step runs after design, before the DONE sentinel.
    assert "flexiprobes select-coverage" in script
    assert "selected_probes_per_gene.csv" in script
    assert "flagged_genes.csv" in script
    assert script.index("select-coverage") < script.index('touch "$DONE"')
    # README documents the extra outputs.
    readme = (dest / "README.txt").read_text()
    assert "selected_probes_per_gene.csv" in readme
    assert "flagged_genes.csv" in readme


def test_non_gene_aware_bundle_omits_selection_step(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">g1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / "bundle"
    cfg = RunConfig(fasta_path=str(fa), output_dir="/local/out",
                    group_by_gene=False,
                    screen=ScreenResources(db_paths=["/laptop/db"]))
    write_hpc_bundle(
        cfg, scheduler="slurm", resources=_resources(),
        image="/apps/flexiprobe.sif", dest_dir=str(dest), stage_inputs=False,
    )
    script = (dest / "run_flexiprobes.slurm").read_text()
    assert "select-coverage" not in script


def test_db_reuses_staged_input_when_same_file(tmp_path):
    # Whole-transcriptome: input FASTA also serves as the off-target DB. A DB
    # entry that is the same file (a separate copy, same name + size) must be
    # repointed at the staged inputs/ copy so the bundle is self-contained.
    content = ">g1 gene:ENSG1\n" + "ACGT" * 30 + "\n"
    fa = tmp_path / "tx.fa"
    fa.write_text(content)
    db = tmp_path / "dbdir" / "tx.fa"          # different location, same name+size
    db.parent.mkdir()
    db.write_text(content)
    cfg = RunConfig(
        fasta_path=str(fa), output_dir="x", group_by_gene=True,
        screen=ScreenResources(db_paths=[str(db)]),
    )
    dest = tmp_path / "bundle"
    write_hpc_bundle(
        cfg, scheduler="slurm", resources=_resources(),
        image="img.sif", dest_dir=str(dest), stage_inputs=True,
    )
    out = load_config(str(dest / "config.json"))
    assert out.screen.db_paths == [os.path.join("inputs", "tx.fa")]
    # README should say the DB is bundled, not referenced on the cluster.
    readme = (dest / "README.txt").read_text()
    assert "BLAST DB is bundled" in readme
    assert "self-contained" in readme


def test_db_left_alone_when_different_file(tmp_path):
    # A genuinely separate DB (different basename) is NOT rewritten.
    fa = tmp_path / "tx.fa"
    fa.write_text(">g1\n" + "ACGT" * 30 + "\n")
    cfg = RunConfig(
        fasta_path=str(fa), output_dir="x",
        screen=ScreenResources(db_paths=["/scratch/host_genome"]),
    )
    dest = tmp_path / "bundle"
    write_hpc_bundle(
        cfg, scheduler="slurm", resources=_resources(),
        image="img.sif", dest_dir=str(dest), stage_inputs=True,
    )
    out = load_config(str(dest / "config.json"))
    assert out.screen.db_paths == ["/scratch/host_genome"]
    assert "referenced on the cluster" in (dest / "README.txt").read_text()


def test_remote_db_override(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">g1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / "bundle"
    write_hpc_bundle(
        _config(str(fa)), scheduler="slurm", resources=_resources(),
        image="img.sif", dest_dir=str(dest), stage_inputs=True,
        remote_db_paths=["/cluster/scratch/human_tx"],
    )
    cfg = load_config(str(dest / "config.json"))
    assert cfg.screen.db_paths == ["/cluster/scratch/human_tx"]


def test_image_path_is_shell_quoted(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">g1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / "bundle"
    write_hpc_bundle(
        _config(str(fa)), scheduler="slurm", resources=_resources(),
        image="/apps/my images/flexiprobe.sif", dest_dir=str(dest),
        stage_inputs=True,
    )
    script = (dest / "run_flexiprobes.slurm").read_text()
    # space in the path must be quoted so the shell sees one token
    assert "CONTAINER='/apps/my images/flexiprobe.sif'" in script


def test_stage_missing_fasta_raises(tmp_path):
    cfg = _config(str(tmp_path / "does_not_exist.fa"))
    with pytest.raises(FileNotFoundError):
        write_hpc_bundle(
            cfg, scheduler="slurm", resources=_resources(),
            image="img.sif", dest_dir=str(tmp_path / "b"), stage_inputs=True,
        )


def test_write_bundle_rejects_bad_scheduler(tmp_path):
    with pytest.raises(ValueError):
        write_hpc_bundle(
            _config("x.fa"), scheduler="lsf", resources=_resources(),
            image="img.sif", dest_dir=str(tmp_path / "b"),
        )


def test_does_not_mutate_caller_config(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">g1\n" + "ACGT" * 30 + "\n")
    cfg = _config(str(fa))
    write_hpc_bundle(
        cfg, scheduler="slurm", resources=_resources(),
        image="img.sif", dest_dir=str(tmp_path / "b"), stage_inputs=True,
        remote_db_paths=["/cluster/db"],
    )
    # caller's config untouched
    assert cfg.output_dir == "/local/laptop/out"
    assert cfg.fasta_path == str(fa)
    assert cfg.screen.db_paths == ["/laptop/human_tx"]


# ---- resumable chaining ----

def _write_bundle(tmp_path, scheduler):
    fa = tmp_path / "in.fa"
    fa.write_text(">g1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / "bundle"
    write_hpc_bundle(
        _config(str(fa)), scheduler=scheduler, resources=_resources(),
        image="/apps/flexiprobe.sif", dest_dir=str(dest), stage_inputs=True,
    )
    return dest


@pytest.mark.parametrize("scheduler,script_name,submit,dep", [
    ("pbs", "run_flexiprobes.pbs", "qsub", "-W depend=afterany:"),
    ("slurm", "run_flexiprobes.slurm", "sbatch --parsable", "--dependency=afterany:"),
])
def test_submit_chain_script_written(tmp_path, scheduler, script_name, submit, dep):
    dest = _write_bundle(tmp_path, scheduler)
    chain = dest / "submit_chain.sh"
    assert chain.exists()
    # executable
    assert os.access(str(chain), os.X_OK)
    body = chain.read_text()
    assert body.startswith("#!/bin/bash")
    assert script_name in body
    assert submit in body          # submits the first link
    assert dep in body             # links subsequent jobs with afterany
    assert "K=" in body            # configurable number of links


def test_run_script_has_done_guard_and_resumes(tmp_path):
    dest = _write_bundle(tmp_path, "pbs")
    script = (dest / "run_flexiprobes.pbs").read_text()
    # idempotent: a finished run no-ops via the DONE sentinel
    assert ".flexiprobes_DONE" in script
    assert 'if [ -f "$DONE" ]' in script
    # resume is communicated and checkpoints are named
    assert ".candidate_cache" in script
    assert ".screen_cache" in script
    # the sentinel is only written on success, BEFORE packaging the results
    done_idx = script.index('touch "$DONE"')
    tar_idx = script.index("tar -czf flexiprobes_output.tar.gz")
    assert done_idx < tar_idx
    # success is gated on the run's exit code, not unconditional
    assert 'rc=$?' in script
    assert 'if [ "$rc" -ne 0 ]' in script


@pytest.mark.parametrize("scheduler,script_name", [
    ("pbs", "run_flexiprobes.pbs"),
    ("slurm", "run_flexiprobes.slurm"),
])
def test_job_log_is_per_job_not_a_shared_filename(tmp_path, scheduler, script_name):
    """Every chain link must log to its own file.

    Regression guard for a confirmed data loss: with a fixed ``#PBS -o
    flexiprobes.log``, the leftover chain links that wake after the run
    completes, see the DONE sentinel and exit each overwrote that file with
    their two-line no-op output — so both the mouse and rat whole-transcriptome
    runs (July 2026) ended up with a 68-byte log written hours after the run
    finished, and their per-stage timings were unrecoverable. PBS spools stdout
    and copies it into place after the job exits, so only the directive can fix
    this.
    """
    dest = _write_bundle(tmp_path, scheduler)
    script = (dest / script_name).read_text()
    directives = [ln for ln in script.splitlines()
                  if ln.startswith("#PBS -o") or ln.startswith("#SBATCH --output")]

    if scheduler == "pbs":
        # PBS has no %j token, so the only safe option is to emit no -o at all
        # and fall through to the default <job_name>.o<jobid>.
        assert directives == [], f"PBS must not pin stdout to a filename: {directives}"
        assert "#PBS -j oe" in script          # stdout+stderr still merged
    else:
        # Slurm can name the file, but only with the per-job %j token in it.
        assert directives, "Slurm bundle emitted no --output directive at all"
        for line in directives:
            assert "%j" in line, (
                f"log target is a fixed filename, links will clobber: {line}")


def test_readme_documents_where_logs_go(tmp_path):
    readme = (_write_bundle(tmp_path, "pbs") / "README.txt").read_text()
    assert "Where the logs go" in readme
    assert "<job_name>.o<jobid>" in readme
    # the manifest is named as the durable fallback for the timings
    assert "manifest.json" in readme


# ---- chain budget: walltime PER LINK x link count ----------------------------
#
# The link count was hardcoded at 12 inside the generated submit_chain.sh and
# settable from neither the CLI nor the GUI, while "walltime" read as a
# per-RUN figure everywhere it was shown. Both halves of the budget are now
# explicit, so both are pinned here.

def test_chain_defaults_are_short_links_not_one_long_job(tmp_path):
    """A long walltime is held until after a maintenance window rather than
    started, so the default link is short and the chain is what covers the run.
    """
    from flexprobe.hpc_export import (
        DEFAULT_CHAIN_LINKS, DEFAULT_WALLTIME, HpcResources,
    )
    assert DEFAULT_WALLTIME == "24:00:00"
    res = HpcResources()
    assert res.walltime == DEFAULT_WALLTIME
    assert res.links == DEFAULT_CHAIN_LINKS == 12


@pytest.mark.parametrize("scheduler", ["slurm", "pbs"])
def test_link_count_reaches_the_generated_chain_script(tmp_path, scheduler):
    fa = tmp_path / "in.fa"
    fa.write_text(">g1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / f"b_{scheduler}"
    res = HpcResources(job_name="wt", walltime="24:00:00", mem_gb=28,
                       cores=8, links=30)
    write_hpc_bundle(
        _config(str(fa)), scheduler=scheduler, resources=res,
        image="/apps/f.sif", dest_dir=str(dest), stage_inputs=True,
    )
    chain = (dest / "submit_chain.sh").read_text()
    # The default K the submitter uses when given no argument.
    assert 'K="${1:-30}"' in chain
    # ...and the budget is spelled out, since neither number means anything alone.
    assert "30 links x 24:00:00" in chain


def test_per_link_walltime_is_the_scheduler_directive(tmp_path):
    """The walltime in the job script is ONE link's, not the whole run's."""
    fa = tmp_path / "in.fa"
    fa.write_text(">g1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / "b"
    res = HpcResources(job_name="wt", walltime="24:00:00", mem_gb=28,
                       cores=8, links=12)
    write_hpc_bundle(
        _config(str(fa)), scheduler="pbs", resources=res,
        image="/apps/f.sif", dest_dir=str(dest), stage_inputs=True,
    )
    assert "#PBS -l walltime=24:00:00" in (dest / "run_flexiprobes.pbs").read_text()
    readme = (dest / "README.txt").read_text()
    assert "PER LINK" in readme
    # The reason, not just the fact — it is why nobody should "fix" this by
    # lengthening the walltime.
    assert "maintenance window" in readme
    assert "submit 12 linked jobs" in readme


def test_readme_documents_chain(tmp_path):
    dest = _write_bundle(tmp_path, "pbs")
    readme = (dest / "README.txt").read_text()
    assert "submit_chain.sh" in readme
    assert ".candidate_cache" in readme
    assert "resume" in readme.lower()


# ---- importer ----

def _make_output(dirpath, with_orders=False):
    os.makedirs(dirpath, exist_ok=True)
    with open(os.path.join(dirpath, "results.csv"), "w") as fh:
        fh.write("rank,probe_id,gene_id\n1,p1,g1\n")
    with open(os.path.join(dirpath, "manifest.json"), "w") as fh:
        json.dump({"params": {"modality": "flex_v2", "discovery_only": True}}, fh)
    if with_orders:
        with open(os.path.join(dirpath, "idt_opools.csv"), "w") as fh:
            fh.write("pool,seq\n")


def test_import_from_folder_ready_for_review(tmp_path):
    out = tmp_path / "job_output"
    _make_output(str(out))
    # `dest_root` is unused on the directory-import path below, but
    # `import_hpc_output` computes its default (`_imported_root()`) BEFORE
    # branching on directory-vs-archive, and that default's side effect is
    # `os.makedirs(~/.flexprobe/imported)` on the REAL developer machine.
    # Passing an explicit tmp_path keeps this test from touching real $HOME.
    info = import_hpc_output(str(out), dest_root=str(tmp_path / "imported"))
    assert info["output_dir"] == os.path.abspath(str(out))
    assert info["status"] == "ready_for_review"
    assert info["modality_label"] == "flex_v2"


def test_import_without_manifest_still_flagged(tmp_path):
    # Only results.csv copied back (no manifest) — must still be badged HPC.
    out = tmp_path / "job_output"
    out.mkdir()
    (out / "results.csv").write_text("rank,probe_id,gene_id\n1,p1,g1\n")
    info = import_hpc_output(str(out), dest_root=str(tmp_path / "imported"))
    assert info["stats"].get("imported") is True
    assert info["status"] == "ready_for_review"


def test_import_folder_with_orders_is_complete(tmp_path):
    out = tmp_path / "job_output"
    _make_output(str(out), with_orders=True)
    assert import_hpc_output(
        str(out), dest_root=str(tmp_path / "imported"))["status"] == "complete"


def test_import_from_targz(tmp_path):
    out = tmp_path / "job_output"
    _make_output(str(out))
    arch = tmp_path / "flexiprobes_output.tar.gz"
    with tarfile.open(arch, "w:gz") as tf:
        tf.add(str(out), arcname=".")
    info = import_hpc_output(str(arch), dest_root=str(tmp_path / "imported"))
    assert os.path.exists(os.path.join(info["output_dir"], "results.csv"))


def test_import_from_zip(tmp_path):
    out = tmp_path / "job_output"
    _make_output(str(out))
    arch = tmp_path / "out.zip"
    with zipfile.ZipFile(arch, "w") as zf:
        zf.write(os.path.join(str(out), "results.csv"), "results.csv")
        zf.write(os.path.join(str(out), "manifest.json"), "manifest.json")
    info = import_hpc_output(str(arch), dest_root=str(tmp_path / "imported"))
    assert os.path.exists(os.path.join(info["output_dir"], "results.csv"))


def test_import_missing_results_raises(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError):
        import_hpc_output(str(empty), dest_root=str(tmp_path / "imported"))


def test_import_unsupported_archive_raises(tmp_path):
    bad = tmp_path / "thing.rar"
    bad.write_text("nope")
    with pytest.raises(ValueError):
        import_hpc_output(str(bad), dest_root=str(tmp_path / "imported"))


# ---- container version assertion ----
#
# The .sif bakes the engine in at build time, so a stale image runs the bundle's
# config with the WRONG engine and produces plausible, wrong results under the
# right filenames. Nothing downstream can catch that, so the job script has to.

def _apptainer_shim(bin_dir, reported_version):
    """Write a fake ``apptainer`` on PATH that reports ``reported_version``.

    Any invocation other than ``... flexiprobes version`` succeeds silently, so
    the rest of the script (design, select-coverage, packaging) runs as a no-op
    and we observe only what the assertion did.
    """
    bin_dir.mkdir(parents=True, exist_ok=True)
    shim = bin_dir / "apptainer"
    shim.write_text(
        "#!/bin/bash\n"
        'for a in "$@"; do\n'
        '  if [ "$a" = "version" ]; then\n'
        f'    echo "flexiprobes {reported_version}"\n'
        "    exit 0\n"
        "  fi\n"
        "done\n"
        "exit 0\n"
    )
    shim.chmod(0o755)
    return bin_dir


def _run_job_script(dest, bin_dir):
    import subprocess
    env = dict(os.environ)
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    env["PBS_O_WORKDIR"] = str(dest)
    return subprocess.run(
        ["bash", str(dest / "run_flexiprobes.pbs")],
        env=env, capture_output=True, text=True, cwd=str(dest),
    )


def test_script_asserts_container_version(tmp_path):
    dest = _write_bundle(tmp_path, "pbs")
    script = (dest / "run_flexiprobes.pbs").read_text()
    from flexprobe import __version__
    assert f"EXPECTED_VERSION={__version__}" in script
    assert "flexiprobes version" in script
    # the check must run BEFORE any pipeline work
    assert script.index("EXPECTED_VERSION=") < script.index("flexiprobes design")


def test_version_mismatch_aborts_the_job(tmp_path):
    """The assertion must FIRE, not merely be present in the text."""
    dest = _write_bundle(tmp_path, "pbs")
    proc = _run_job_script(dest, _apptainer_shim(tmp_path / "bin", "0.3.10"))
    assert proc.returncode == 3, proc.stdout + proc.stderr
    assert "VERSION MISMATCH" in proc.stderr
    assert "0.3.10" in proc.stderr
    # and it aborted before doing anything: no success sentinel, no archive
    assert not (dest / ".flexiprobes_DONE").exists()
    assert not (dest / "flexiprobes_output.tar.gz").exists()


def test_matching_version_proceeds(tmp_path):
    """Control for the test above: the right version is not blocked."""
    from flexprobe import __version__
    dest = _write_bundle(tmp_path, "pbs")
    proc = _run_job_script(dest, _apptainer_shim(tmp_path / "bin", __version__))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "Container engine version OK" in proc.stdout
    assert (dest / ".flexiprobes_DONE").exists()


def test_unreadable_container_version_aborts(tmp_path):
    """A container that can't report a version is also a refusal, not a warning."""
    dest = _write_bundle(tmp_path, "pbs")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    shim = bin_dir / "apptainer"
    shim.write_text("#!/bin/bash\nexit 1\n")   # e.g. the .sif isn't there
    shim.chmod(0o755)
    proc = _run_job_script(dest, bin_dir)
    assert proc.returncode == 3, proc.stdout + proc.stderr
    assert "could not read the engine version" in proc.stderr.lower()
    # It must also stop the CHAIN, not just this link. A preflight failure is
    # deterministic, so without the sentinel every remaining link reproduces it
    # — which is what happened on 18 Aug (12 jobs in 70 seconds).
    assert (dest / ".flexiprobes_ABORT").exists(), proc.stdout + proc.stderr


def test_expected_version_is_overridable(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">g1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / "bundle"
    write_hpc_bundle(
        _config(str(fa)), scheduler="pbs", resources=_resources(),
        image="/apps/flexiprobe.sif", dest_dir=str(dest), stage_inputs=True,
        expected_version="1.2.3",
    )
    assert "EXPECTED_VERSION=1.2.3" in (dest / "run_flexiprobes.pbs").read_text()
    assert "flexiprobes 1.2.3" in (dest / "README.txt").read_text()


def test_readme_leads_with_the_container_rebuild(tmp_path):
    readme = (_write_bundle(tmp_path, "pbs") / "README.txt").read_text()
    assert "Rebuild the container FIRST" in readme
    # the Lustre-quota gotcha must be in the build recipe, not folklore
    assert "APPTAINER_TMPDIR" in readme
    assert "APPTAINER_CACHEDIR" in readme
    # and it must precede the submit instructions
    assert (readme.index("Rebuild the container FIRST")
            < readme.index("Submit (recommended"))


# ---- pre-screen candidate artefacts are declared in the bundle README ----

def _bundle_with(tmp_path, **cfg_kw):
    fa = tmp_path / "in.fa"
    fa.write_text(">g1 gene:ENSG1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / "bundle"
    cfg = _config(str(fa))
    for k, v in cfg_kw.items():
        setattr(cfg, k, v)
    write_hpc_bundle(
        cfg, scheduler="pbs", resources=_resources(),
        image="/apps/flexiprobe.sif", dest_dir=str(dest), stage_inputs=True,
    )
    return dest


def test_readme_says_no_prescreen_capture_by_default(tmp_path):
    readme = (_bundle_with(tmp_path) / "README.txt").read_text()
    assert "does NOT retain the pre-screen candidate set" in readme


def test_readme_declares_postcap_capture(tmp_path):
    readme = (_bundle_with(tmp_path, emit_prescreen_candidates=True)
              / "README.txt").read_text()
    assert "CAPTURE MODE: post-cap" in readme
    assert "candidates_prescreen.csv.gz" in readme
    assert "candidates_precap.csv.gz" not in readme


def test_readme_declares_precap_capture(tmp_path):
    """Six months later, nobody should have to guess which set a run kept."""
    readme = (_bundle_with(tmp_path, emit_prescreen_candidates=True,
                           emit_prescreen_precap_rows=True)
              / "README.txt").read_text()
    assert "CAPTURE MODE: pre-cap" in readme
    assert "candidates_precap.csv.gz" in readme
    assert "prescreen_transcripts.csv.gz" in readme
    assert "prescreen_index.json" in readme


# --- chain abort on deterministic failure (13 Aug regression) --------------
# `afterany` links start regardless of the previous exit status. A reproducible
# crash therefore got retried once per link: 12 jobs, ~1h36m of 16-core
# walltime, zero progress. The generated script must classify the failure.

def _pbs_body():
    from flexprobe.hpc_export import _script_body
    return _script_body("pbs", "../flexiprobe-0.5.1.sif", 16)


def test_script_aborts_chain_on_non_signal_exit():
    body = _pbs_body()
    assert ABORT_SENTINEL in body
    # rc <= 128 means the pipeline chose to exit -> deterministic -> abort
    assert '[ "$rc" -le 128 ]' in body


def test_script_no_ops_when_abort_sentinel_present():
    body = _pbs_body()
    assert 'if [ -f "$ABORT" ]; then' in body
    assert "Chain ABORTED by an earlier link" in body


def test_script_still_resumes_signal_kills():
    """Walltime and OOM kills must keep resuming — that is the whole point of
    the chain, and over-tightening this would break long runs."""
    body = _pbs_body()
    assert "retryable — walltime/OOM/node fault" in body


def test_script_retries_a_worker_oom_despite_a_low_exit_code():
    """A worker OOM is invisible in the exit status, and must still resume.

    With ``--screen-workers`` the OOM killer takes a worker, not the parent.
    The parent raises ``BrokenProcessPool``, unwinds cleanly and exits rc=1 —
    a resource fault wearing the exit status of a logic bug. Unhandled, one
    such kill aborts every remaining link of a chain.
    """
    body = _pbs_body()
    # The link's output has to be captured, or there is nothing to classify on.
    assert 'LINK_LOG=' in body
    assert 'tee "$LINK_LOG"' in body
    assert "rc=${PIPESTATUS[0]}" in body
    # ... and the low-exit-code branch must consult it before aborting.
    assert 'BrokenProcessPool' in body
    assert 'MemoryError' in body
    assert '[ "$rc" -le 128 ] && grep -qE' in body


def test_script_still_aborts_a_genuine_deterministic_failure():
    """The counterweight to the test above: widening the retry test must not
    turn a reproducible failure into an infinite retry loop."""
    body = _pbs_body()
    # The abort branch survives, as the *else* of the resource-fault check.
    assert 'elif [ "$rc" -le 128 ]; then' in body
    assert "Deterministic failure" in body
    # Matching on "Traceback" would have swallowed real bugs, so it must not
    # appear in the retry condition.
    cond = [l for l in body.splitlines() if '-le 128 ] && grep -qE' in l]
    assert len(cond) == 1, cond
    assert "Traceback" not in cond[0]


def test_script_has_progress_based_stall_guard():
    """Independent of the exit code and of any error text: a link that
    advanced no checkpoint did no work."""
    body = _pbs_body()
    assert "progress_fingerprint" in body
    assert 'before_fp="$(progress_fingerprint)"' in body
    assert '[ "$stalled" -ge 2 ]' in body


def test_generated_pbs_script_is_valid_shell(tmp_path):
    import subprocess
    p = tmp_path / "run.sh"
    p.write_text("#!/bin/bash\n" + _pbs_body())
    r = subprocess.run(["bash", "-n", str(p)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_missing_container_runtime_aborts_the_chain(tmp_path):
    """No apptainer on PATH is the single most likely batch-job failure: a job
    does not inherit the modules you loaded interactively. It must stop the
    chain, not be retried once per link."""
    dest = _write_bundle(tmp_path, "pbs")
    empty_bin = tmp_path / "emptybin"
    empty_bin.mkdir()
    proc = _run_job_script(dest, empty_bin)
    assert proc.returncode != 0
    assert "No container runtime" in proc.stderr
    assert (dest / ".flexiprobes_ABORT").exists(), proc.stdout + proc.stderr
    # and the message must name the way out
    assert "APPTAINER_MODULE" in (dest / ".flexiprobes_ABORT").read_text()


def test_abort_sentinel_makes_the_next_link_a_noop(tmp_path):
    dest = _write_bundle(tmp_path, "pbs")
    (dest / ".flexiprobes_ABORT").write_text("earlier link said no\n")
    empty_bin = tmp_path / "emptybin2"
    empty_bin.mkdir()
    proc = _run_job_script(dest, empty_bin)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "Chain ABORTED by an earlier link" in proc.stdout


def test_submit_chain_forwards_the_module_to_every_link():
    from flexprobe.hpc_export import _submit_chain_script
    body = _submit_chain_script("pbs")
    assert "APPTAINER_MODULE" in body
    assert "-v APPTAINER_MODULE=$APPTAINER_MODULE" in body


# ---- container module: the commonest way a chain dies -----------------------
#
# A batch job does not inherit interactively-loaded modules, so a chain without
# its container module aborts seconds into link 1. The escape hatch
# (APPTAINER_MODULE=... ./submit_chain.sh) is not something anyone remembers to
# type, so the module name is recorded in the bundle.

def _bundle_with_module(tmp_path, module, name="b"):
    fa = tmp_path / "in.fa"
    fa.write_text(">g1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / name
    write_hpc_bundle(
        _config(str(fa)), scheduler="pbs", resources=HpcResources(),
        image="/apps/f.sif", dest_dir=str(dest), stage_inputs=True,
        apptainer_module=module,
    )
    return dest


def test_module_is_baked_into_the_job_script(tmp_path):
    dest = _bundle_with_module(tmp_path, "apptainer/1.4.1")
    script = (dest / "run_flexiprobes.pbs").read_text()
    assert "BUNDLE_APPTAINER_MODULE=apptainer/1.4.1" in script
    # ...and the environment still wins, so one bundle works on two clusters.
    assert 'APPTAINER_MODULE="${APPTAINER_MODULE:-$BUNDLE_APPTAINER_MODULE}"' \
        in script
    assert "module load apptainer/1.4.1" in (dest / "README.txt").read_text()


def test_no_module_says_so_rather_than_looking_configured(tmp_path):
    dest = _bundle_with_module(tmp_path, "", name="c")
    script = (dest / "run_flexiprobes.pbs").read_text()
    assert "BUNDLE_APPTAINER_MODULE=''" in script
    chain = (dest / "submit_chain.sh").read_text()
    assert "no container module named" in chain
    assert "module avail apptainer" in (dest / "README.txt").read_text()


def test_unversioned_module_is_tried_before_giving_up(tmp_path):
    """A site retiring apptainer/1.4.1 must not abort a chain whose only real
    problem is a stale string in a config file."""
    script = (_bundle_with_module(tmp_path, "apptainer/9.9.9")
              / "run_flexiprobes.pbs").read_text()
    assert "for _m in apptainer singularity; do" in script
    # The abort message has to explain the batch/interactive distinction, which
    # is the actual misunderstanding behind the failure.
    assert "DOES NOT INHERIT THE MODULES YOU LOADED" in script


def test_submit_chain_refuses_to_queue_into_an_aborted_bundle(tmp_path):
    """The sentinel makes every link no-op, so submitting K jobs on top of one
    is K wasted submissions and a log dir that looks like a fresh failure."""
    import subprocess
    dest = _bundle_with_module(tmp_path, "apptainer/1.4.1", name="d")
    (dest / ".flexiprobes_ABORT").write_text(
        "No container runtime: ...\nJob: 34817504.pbs01\n")
    proc = subprocess.run(
        ["./submit_chain.sh"], cwd=str(dest), capture_output=True, text=True,
        env={"PATH": "/usr/bin:/bin"},
    )
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert "REFUSING TO SUBMIT" in proc.stderr
    # It quotes the previous failure, so the reason is on screen, not in a file
    # the user has to know to open.
    assert "34817504.pbs01" in proc.stderr
    assert "rm .flexiprobes_ABORT" in proc.stderr
    # And nothing was queued: qsub is not even on the PATH we gave it.
    assert "link 1/" not in proc.stdout
