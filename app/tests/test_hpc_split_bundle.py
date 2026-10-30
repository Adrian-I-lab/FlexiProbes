"""The split prepare/screen HPC bundle: two job types, one chain, one config.

The unsplit bundle asks the scheduler for the MAXIMUM of what any stage needs,
for the whole walltime, and on a whole transcriptome those maxima come from two
different stages. The same screening work peaks far higher in a link that
tiled first than in one that resumed from a checkpoint, because CPython does
not return freed pages to the OS, so the screen forks its workers off a parent
still holding the tiling high-water mark. An unsplit job sized for that peak
is easily killed by a memory or walltime limit.

These tests pin the four things that make the split safe rather than merely
smaller, each of which fails silently if it regresses:

* the screen job REFUSES to start without a prepared checkpoint, instead of
  quietly re-tiling inside an allocation sized for a job that does not tile,
* the prepare job never writes ``.flexiprobes_DONE`` (which would strand the
  bundle looking finished),
* every failure-handling property of the unsplit script survives in both, and
* the unsplit script itself is untouched.
"""
import os
import subprocess

import pytest

from flexprobe.config_io import load_config
from flexprobe.constants import PREPARED_MARKER
from flexprobe import __version__ as engine_version
from flexprobe.hpc_export import (
    ABORT_SENTINEL, DEFAULT_PREPARE_CORES, DEFAULT_PREPARE_LINKS,
    DEFAULT_PREPARE_MEM_GB, DEFAULT_PREPARE_WALLTIME, DONE_SENTINEL,
    HpcResources, PREPARE_SCRIPT_STEM, SCREEN_SCRIPT_STEM, write_hpc_bundle,
)
from flexprobe.models import RunConfig, ScreenResources


def _config(fasta: str, **kw) -> RunConfig:
    return RunConfig(
        fasta_path=fasta,
        output_dir="/local/laptop/out",
        group_by_gene=True,
        screen=ScreenResources(db_paths=["/laptop/human_tx"], **kw),
    )


def _bundle(tmp_path, scheduler="pbs", *, split=True, cfg_kw=None,
            prepare_resources=None, resources=None):
    fa = tmp_path / "in.fa"
    fa.write_text(">g1 gene:ENSG1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / "bundle"
    write_hpc_bundle(
        _config(str(fa), **(cfg_kw or {})),
        scheduler=scheduler,
        resources=resources or HpcResources(
            job_name="wt", walltime="24:00:00", mem_gb=48, cores=8, links=12),
        image="../flexiprobe-0.6.2.sif",
        dest_dir=str(dest),
        split_prepare_screen=split,
        prepare_resources=prepare_resources,
    )
    return dest


def _ext(scheduler):
    return "slurm" if scheduler == "slurm" else "pbs"


# ---------------------------------------------------------------------------
# shape of the bundle
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("scheduler", ["pbs", "slurm"])
def test_split_bundle_writes_two_scripts_and_no_single_one(tmp_path, scheduler):
    """The unsplit script must be ABSENT, not merely unused.

    A leftover run_flexiprobes.<ext> beside the two is a trap: submitting it
    does the whole run in one process at whichever resource request happens to
    be in it, which is exactly what this bundle exists to prevent.
    """
    dest = _bundle(tmp_path, scheduler)
    ext = _ext(scheduler)
    assert (dest / f"{PREPARE_SCRIPT_STEM}.{ext}").exists()
    assert (dest / f"{SCREEN_SCRIPT_STEM}.{ext}").exists()
    assert not (dest / f"run_flexiprobes.{ext}").exists()
    assert (dest / "submit_chain.sh").exists()
    assert (dest / "config.json").exists()


def test_unsplit_bundle_is_unchanged(tmp_path):
    """The default is still one job with the name it always had."""
    dest = _bundle(tmp_path, "pbs", split=False)
    assert (dest / "run_flexiprobes.pbs").exists()
    assert not (dest / f"{PREPARE_SCRIPT_STEM}.pbs").exists()
    assert not (dest / f"{SCREEN_SCRIPT_STEM}.pbs").exists()


def test_the_two_jobs_ask_for_different_resources(tmp_path):
    dest = _bundle(tmp_path, "pbs")
    prep = (dest / f"{PREPARE_SCRIPT_STEM}.pbs").read_text()
    scr = (dest / f"{SCREEN_SCRIPT_STEM}.pbs").read_text()

    assert (f"#PBS -l select=1:ncpus={DEFAULT_PREPARE_CORES}:"
            f"mem={DEFAULT_PREPARE_MEM_GB}gb") in prep
    assert f"#PBS -l walltime={DEFAULT_PREPARE_WALLTIME}" in prep
    assert "#PBS -l select=1:ncpus=8:mem=48gb" in scr
    assert "#PBS -l walltime=24:00:00" in scr
    # ...and they are distinguishable in the queue, which matters because the
    # two logs land as <job_name>.o<jobid> in the same directory.
    assert "#PBS -N wt-prepare" in prep
    assert "#PBS -N wt" in scr


def test_resource_figures_are_parameters_not_literals(tmp_path):
    """Every number in the split is overridable from the generator."""
    dest = _bundle(
        tmp_path, "pbs",
        prepare_resources=HpcResources(
            job_name="custom-prep", walltime="01:02:03", mem_gb=7, cores=3,
            links=5),
        resources=HpcResources(
            job_name="custom", walltime="04:05:06", mem_gb=9, cores=2, links=4),
    )
    prep = (dest / f"{PREPARE_SCRIPT_STEM}.pbs").read_text()
    scr = (dest / f"{SCREEN_SCRIPT_STEM}.pbs").read_text()
    chain = (dest / "submit_chain.sh").read_text()
    assert "#PBS -l select=1:ncpus=3:mem=7gb" in prep
    assert "#PBS -l walltime=01:02:03" in prep
    assert "--workers 3" in prep
    assert "#PBS -l select=1:ncpus=2:mem=9gb" in scr
    assert "#PBS -l walltime=04:05:06" in scr
    assert "--workers 2 --screen-workers 2" in scr
    assert "P=5" in chain
    assert 'K="${1:-4}"' in chain


# ---------------------------------------------------------------------------
# what each job runs
# ---------------------------------------------------------------------------

def test_prepare_job_runs_prepare_only_and_does_not_finish_the_run(tmp_path):
    dest = _bundle(tmp_path, "pbs")
    prep = (dest / f"{PREPARE_SCRIPT_STEM}.pbs").read_text()

    assert "--prepare-only" in prep
    # No screen workers: nothing forks, because the run stops before the screen.
    # Checked on the invocation line itself — the surrounding comments are
    # phase-specific prose and must not make this assertion vacuous.
    design_line = [l for l in prep.splitlines()
                   if l.lstrip().startswith("flexiprobes design")
                   or "$CONTAINER" in l and "flexiprobes design" in l]
    assert len(design_line) == 1, design_line
    assert "--screen-workers" not in design_line[0]
    assert "--prepare-only" in design_line[0]
    # It must NOT finish the run: no DONE sentinel, no tarball, no per-gene
    # selection (there is no results.csv to select from).
    assert 'touch "$DONE"' not in prep
    assert "tar -czf" not in prep
    assert "select-coverage" not in prep


def test_screen_job_finishes_the_run(tmp_path):
    dest = _bundle(tmp_path, "pbs")
    scr = (dest / f"{SCREEN_SCRIPT_STEM}.pbs").read_text()

    assert "--prepare-only" not in scr
    assert "--workers 8 --screen-workers 8" in scr
    assert 'touch "$DONE"' in scr
    assert "tar -czf flexiprobes_output.tar.gz" in scr
    assert "select-coverage" in scr, "gene-aware runs still get the selection step"


def test_screen_job_refuses_to_start_without_the_prepare_marker(tmp_path):
    """The single most important guard in the split.

    Started without a checkpoint this job does NOT fail fast on its own: it
    tiles the transcriptome inside a 48 GB allocation chosen for a job that
    does not tile, and dies on memory or walltime hours later.
    """
    dest = _bundle(tmp_path, "pbs")
    scr = (dest / f"{SCREEN_SCRIPT_STEM}.pbs").read_text()

    assert f'PREPARED="output/{PREPARED_MARKER}"' in scr
    assert 'if [ ! -f "$PREPARED" ]; then' in scr
    assert "No prepared candidate checkpoint" in scr
    # It aborts the CHAIN, not just the link: every screen link would hit the
    # same wall, and the fix is to resubmit the prepare half.
    assert "abort_chain" in scr.split('if [ ! -f "$PREPARED" ]')[1][:600]
    # ...and the check runs before any pipeline work.
    assert scr.index('if [ ! -f "$PREPARED" ]') < scr.index("flexiprobes design")


def test_prepare_job_no_ops_when_already_prepared(tmp_path):
    """submit_chain.sh queues a spare prepare link as a free retry, and
    `afterany` starts it whatever happened before. Without this guard link 2
    would re-tile the transcriptome for nothing."""
    dest = _bundle(tmp_path, "pbs")
    prep = (dest / f"{PREPARE_SCRIPT_STEM}.pbs").read_text()
    assert 'if [ -f "$PREPARED" ]; then' in prep
    assert "already prepared" in prep


def test_prepare_job_fails_loudly_if_it_exits_0_without_the_marker(tmp_path):
    """Exit 0 and no marker is a bug, and it must stop the chain HERE, where
    the cause is still in this link's log — not two jobs later when the screen
    job refuses for a reason it cannot explain."""
    dest = _bundle(tmp_path, "pbs")
    prep = (dest / f"{PREPARE_SCRIPT_STEM}.pbs").read_text()
    # rsplit: `rm -f "$STALL"` also appears inside the rc!=0 branch above.
    tail = prep.rsplit('rm -f "$STALL"', 1)[1]
    assert 'if [ ! -f "$PREPARED" ]; then' in tail
    assert "abort_chain" in tail


# ---------------------------------------------------------------------------
# every safety property of the unsplit script survives in both
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("stem", [PREPARE_SCRIPT_STEM, SCREEN_SCRIPT_STEM])
def test_both_scripts_keep_the_chain_safety_machinery(tmp_path, stem):
    """The DONE no-op, the ABORT stop, the stall counter and the worker-OOM
    classifier were each paid for in dead cluster jobs (13 Aug, 18 Aug,
    27 Aug, 1 Sep). Splitting the job must not drop any of them."""
    dest = _bundle(tmp_path, "pbs")
    body = (dest / f"{stem}.pbs").read_text()

    assert f'DONE="{DONE_SENTINEL}"' in body
    assert 'if [ -f "$DONE" ]; then' in body
    assert f'ABORT="{ABORT_SENTINEL}"' in body
    assert "Chain ABORTED by an earlier link" in body
    assert "progress_fingerprint" in body
    assert '[ "$stalled" -ge 2 ]' in body
    assert '[ "$rc" -le 128 ] && grep -qE' in body
    assert "BrokenProcessPool" in body
    assert 'elif [ "$rc" -le 128 ]; then' in body
    assert "Deterministic failure" in body
    assert "retryable — walltime/OOM/node fault" in body
    assert "EXPECTED_VERSION=" in body
    assert body.index("EXPECTED_VERSION=") < body.index("flexiprobes design")
    assert "No container runtime" in body


@pytest.mark.parametrize("name", [
    f"{PREPARE_SCRIPT_STEM}.pbs", f"{SCREEN_SCRIPT_STEM}.pbs",
    "submit_chain.sh",
])
def test_generated_split_scripts_are_valid_shell(tmp_path, name):
    dest = _bundle(tmp_path, "pbs")
    r = subprocess.run(["bash", "-n", str(dest / name)],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


# ---------------------------------------------------------------------------
# the guards, EXECUTED
# ---------------------------------------------------------------------------
# String matching proves the text is in the file. These run the generated
# scripts under a fake `apptainer` and prove the shell logic does what the text
# says, in both directions — a guard that fires on every input is measuring
# nothing.

def _apptainer_shim(bin_dir, version=None):
    """A fake `apptainer` that reports a matching engine version and otherwise
    succeeds silently, so the pipeline invocation is a no-op and only the
    script's own logic is under test.

    The reported version DEFAULTS TO THE LIVE ENGINE VERSION. It was a
    hard-coded literal, which meant every release bump silently turned these
    tests from "the gate passes on a matching container" into "the gate fires
    on a mismatched one" — the opposite assertion, still green-looking until
    someone read it. Pass an explicit version to test the mismatch direction.
    """
    if version is None:
        version = engine_version
    bin_dir.mkdir(parents=True, exist_ok=True)
    shim = bin_dir / "apptainer"
    shim.write_text(
        "#!/bin/bash\n"
        'for a in "$@"; do\n'
        '  if [ "$a" = "version" ]; then\n'
        f'    echo "flexiprobes {version}"\n'
        "    exit 0\n"
        "  fi\n"
        "done\n"
        "exit 0\n"
    )
    shim.chmod(0o755)
    return bin_dir


def _run_script(dest, bin_dir, name):
    env = dict(os.environ)
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    env["PBS_O_WORKDIR"] = str(dest)
    return subprocess.run(["bash", str(dest / name)], env=env,
                          capture_output=True, text=True, cwd=str(dest))


def test_screen_script_actually_aborts_without_the_marker(tmp_path):
    dest = _bundle(tmp_path, "pbs")
    shim = _apptainer_shim(tmp_path / "bin")
    proc = _run_script(dest, shim, f"{SCREEN_SCRIPT_STEM}.pbs")

    assert proc.returncode == 4, proc.stdout + proc.stderr
    assert (dest / ABORT_SENTINEL).exists()
    assert "No prepared candidate checkpoint" in (dest / ABORT_SENTINEL).read_text()


def test_screen_script_proceeds_once_the_marker_exists(tmp_path):
    """The other direction. Without this the test above would pass just as well
    against a script that refused unconditionally."""
    dest = _bundle(tmp_path, "pbs")
    shim = _apptainer_shim(tmp_path / "bin")
    (dest / "output").mkdir(exist_ok=True)
    (dest / "output" / PREPARED_MARKER).write_text('{"candidates": 3}\n')

    proc = _run_script(dest, shim, f"{SCREEN_SCRIPT_STEM}.pbs")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert not (dest / ABORT_SENTINEL).exists()
    assert "Prepared checkpoint found" in proc.stdout
    # The shim makes `flexiprobes design` a silent success, so the script runs
    # to the end and finishes the run exactly as it would in production.
    assert (dest / DONE_SENTINEL).exists()


def test_prepare_script_no_ops_when_the_marker_exists(tmp_path):
    dest = _bundle(tmp_path, "pbs")
    shim = _apptainer_shim(tmp_path / "bin")
    (dest / "output").mkdir(exist_ok=True)
    (dest / "output" / PREPARED_MARKER).write_text('{"candidates": 3}\n')

    proc = _run_script(dest, shim, f"{PREPARE_SCRIPT_STEM}.pbs")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "already prepared" in proc.stdout
    assert not (dest / DONE_SENTINEL).exists(), \
        "the prepare job must never mark the run finished"


def test_prepare_script_aborts_when_the_pipeline_writes_no_marker(tmp_path):
    """The shim makes `flexiprobes design` exit 0 without doing anything, which
    is exactly the failure this guard is for: a success that produced no
    checkpoint. It must stop here, not two jobs later."""
    dest = _bundle(tmp_path, "pbs")
    shim = _apptainer_shim(tmp_path / "bin")

    proc = _run_script(dest, shim, f"{PREPARE_SCRIPT_STEM}.pbs")
    assert proc.returncode == 5, proc.stdout + proc.stderr
    assert (dest / ABORT_SENTINEL).exists()
    assert "did not write" in (dest / ABORT_SENTINEL).read_text()
    assert not (dest / DONE_SENTINEL).exists()


# ---------------------------------------------------------------------------
# the chain
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("scheduler,dep", [
    ("pbs", "-W depend=afterany:"),
    ("slurm", "--dependency=afterany:"),
])
def test_chain_submits_prepare_then_screen_with_afterany(
        tmp_path, scheduler, dep):
    dest = _bundle(tmp_path, scheduler)
    chain = (dest / "submit_chain.sh").read_text()
    ext = _ext(scheduler)
    assert f'PREPARE_SCRIPT="{PREPARE_SCRIPT_STEM}.{ext}"' in chain
    assert f'SCREEN_SCRIPT="{SCREEN_SCRIPT_STEM}.{ext}"' in chain
    assert dep in chain
    # Order matters: the prepare loop must come first.
    assert chain.index('submit "$PREPARE_SCRIPT"') < \
        chain.index('submit "$SCREEN_SCRIPT"')
    assert f"P={DEFAULT_PREPARE_LINKS}" in chain


def test_chain_still_refuses_to_queue_into_an_aborted_bundle(tmp_path):
    dest = _bundle(tmp_path, "pbs")
    chain = (dest / "submit_chain.sh").read_text()
    assert f'if [ -f "{ABORT_SENTINEL}" ]; then' in chain
    assert "REFUSING TO SUBMIT" in chain


# ---------------------------------------------------------------------------
# one config, two jobs
# ---------------------------------------------------------------------------

def test_prepare_only_is_never_baked_into_the_bundle_config(tmp_path):
    """Both jobs read the SAME config.json; the flag is per-job.

    Baked in, the screen job would prepare again and never screen. The exporter
    therefore clears it rather than trusting the caller's config.
    """
    cfg = RunConfig(fasta_path=str(tmp_path / "in.fa"), output_dir="o",
                    prepare_only=True,
                    screen=ScreenResources(db_paths=["/db"]))
    (tmp_path / "in.fa").write_text(">g1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / "b"
    write_hpc_bundle(cfg, scheduler="pbs", resources=HpcResources(),
                     image="../x.sif", dest_dir=str(dest),
                     split_prepare_screen=True)
    assert load_config(str(dest / "config.json")).prepare_only is False
    # The prepare job supplies it on the command line instead.
    assert "--prepare-only" in (dest / f"{PREPARE_SCRIPT_STEM}.pbs").read_text()


@pytest.mark.parametrize("chunk", [2500, 5000])
def test_query_chunk_reaches_the_bundle_config(tmp_path, chunk):
    """Part of the same change: a bundle can request a different screen chunk
    without anyone editing constants.py."""
    dest = _bundle(tmp_path, "pbs", cfg_kw={"query_chunk": chunk})
    assert load_config(str(dest / "config.json")).screen.query_chunk == chunk


def test_query_chunk_default_is_unchanged(tmp_path):
    """The default must stay 5000 — every measurement quoted in this change
    was taken at 5000, and moving it silently would invalidate them."""
    from flexprobe.constants import BLAST_QUERY_CHUNK
    assert BLAST_QUERY_CHUNK == 5000
    dest = _bundle(tmp_path, "pbs")
    assert load_config(str(dest / "config.json")).screen.query_chunk == 5000


# ---------------------------------------------------------------------------
# the README a human reads months later
# ---------------------------------------------------------------------------

def test_readme_explains_the_two_jobs(tmp_path):
    dest = _bundle(tmp_path, "pbs")
    readme = (dest / "README.txt").read_text()
    assert "TWO JOBS, NOT ONE" in readme
    assert f"{PREPARE_SCRIPT_STEM}.pbs" in readme
    assert f"{SCREEN_SCRIPT_STEM}.pbs" in readme
    assert PREPARED_MARKER in readme
    assert "status \"prepared\"" in readme
    assert "REFUSES to start" in readme
    # It must not tell the reader to submit a script the bundle does not have.
    assert "qsub run_flexiprobes.pbs" not in readme


def test_unsplit_readme_says_nothing_about_two_jobs(tmp_path):
    dest = _bundle(tmp_path, "pbs", split=False)
    readme = (dest / "README.txt").read_text()
    assert "TWO JOBS, NOT ONE" not in readme
