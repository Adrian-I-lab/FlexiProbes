"""Real-time PBS job logs under <project root>/hpc_job_logs/ (--job-log-project).

PBS spools a job's stdout on the compute node and copies it into place only
when the job ENDS, so a running link has nothing readable. The only live output
was .flexiprobes_link.log: a dotfile inside the bundle, with no timestamp or job
id in its name, overwritten by every link of the chain.

With ``job_log_project`` set, each PBS script discards the spooled copy
(``#PBS -o /dev/null``) and its first command redirects everything into
``$WORKDIR/projects/<project>/hpc_job_logs/<YYYYMMDD_HHMMSS>_<job>.<jobid>.log``.

These tests EXECUTE the generated scripts under a fake ``apptainer``, because
the property that matters (a log that exists, in the right place, while the
job runs, and the retry classifier still reading its own copy) is shell
behaviour, not text.
"""
import os
import re
import subprocess

import pytest

from flexprobe import __version__ as engine_version
from flexprobe.constants import PREPARED_MARKER
from flexprobe.hpc_export import (
    ABORT_SENTINEL, DONE_SENTINEL, HpcResources, PREPARE_SCRIPT_STEM,
    SCREEN_SCRIPT_STEM, write_hpc_bundle,
)
from flexprobe.models import RunConfig, ScreenResources

PROJECT = "flexiprobe"
SCREEN = f"{SCREEN_SCRIPT_STEM}.pbs"
PREPARE = f"{PREPARE_SCRIPT_STEM}.pbs"


def _bundle(tmp_path, *, split=True, project=PROJECT, scheduler="pbs",
            name="bundle"):
    fa = tmp_path / "in.fa"
    fa.write_text(">g1 gene:ENSG1\n" + "ACGT" * 30 + "\n")
    dest = tmp_path / name
    write_hpc_bundle(
        RunConfig(fasta_path=str(fa), output_dir="/x", group_by_gene=True,
                  screen=ScreenResources(db_paths=["/laptop/tx"])),
        scheduler=scheduler,
        resources=HpcResources(job_name="wt", walltime="24:00:00", mem_gb=48,
                               cores=8, links=12),
        image="../f.sif", dest_dir=str(dest),
        split_prepare_screen=split, job_log_project=project,
    )
    return dest


def _shim(bin_dir):
    """Fake apptainer: reports the live engine version; `design` prints
    $SHIM_TEXT and exits $SHIM_RC; anything else succeeds silently."""
    bin_dir.mkdir(parents=True, exist_ok=True)
    shim = bin_dir / "apptainer"
    shim.write_text(
        "#!/bin/bash\n"
        'for a in "$@"; do\n'
        f'  if [ "$a" = "version" ]; then echo "flexiprobes {engine_version}"; exit 0; fi\n'
        '  if [ "$a" = "design" ]; then echo "${SHIM_TEXT:-}"; exit "${SHIM_RC:-0}"; fi\n'
        "done\n"
        "exit 0\n"
    )
    shim.chmod(0o755)
    return bin_dir


def _run(dest, bin_dir, script, **env_kw):
    env = {k: v for k, v in os.environ.items() if k not in ("WORKDIR", "PROJ")}
    env["PATH"] = f"{bin_dir}:{env['PATH']}"
    env["PBS_O_WORKDIR"] = str(dest)
    env.update({k: str(v) for k, v in env_kw.items()})
    return subprocess.run(["bash", str(dest / script)], env=env,
                          capture_output=True, text=True, cwd=str(dest))


def _prepared(dest):
    (dest / "output").mkdir(exist_ok=True)
    (dest / "output" / PREPARED_MARKER).write_text('{"candidates": 3}\n')


def _logs(root):
    d = root / "hpc_job_logs"
    return sorted(d.iterdir()) if d.is_dir() else []


# ---------------------------------------------------------------------------
# shape
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("split,scripts", [
    (True, [PREPARE, SCREEN]), (False, ["run_flexiprobes.pbs"]),
])
def test_header_discards_the_spooled_log(tmp_path, split, scripts):
    dest = _bundle(tmp_path, split=split)
    for name in scripts:
        text = (dest / name).read_text()
        assert "#PBS -j oe" in text
        assert "#PBS -o /dev/null" in text


@pytest.mark.parametrize("split", [True, False])
def test_off_by_default(tmp_path, split):
    """Negative control: without a project the scripts carry none of it."""
    dest = _bundle(tmp_path, split=split, project="")
    for p in dest.glob("*.pbs"):
        text = p.read_text()
        assert "#PBS -o" not in text
        assert "hpc_job_logs" not in text
        assert "exec >" not in text


@pytest.mark.parametrize("name", [PREPARE, SCREEN])
def test_redirect_is_the_first_thing_after_set(tmp_path, name):
    """Anything that runs above the exec prints to /dev/null and is lost."""
    lines = (_bundle(tmp_path) / name).read_text().splitlines()
    set_idx = lines.index("set -uo pipefail")
    # Code only: the block's comments name BASH_SOURCE to explain avoiding it.
    code = [ln for ln in lines[set_idx + 1:]
            if ln.strip() and not ln.lstrip().startswith("#")]
    exec_at = code.index('exec >"${LOG_FILE}" 2>&1')
    before = code[:exec_at]
    # Only the log-location logic may precede the exec.
    assert before[0].startswith("PROJECT_ROOT=")
    assert len(before) <= 10, before
    assert exec_at < code.index('cd "$PBS_O_WORKDIR"')
    # Located from the environment, never from the script's own path: PBS runs
    # a spooled copy.
    assert not any("BASH_SOURCE" in ln or "$0" in ln for ln in before), before


def test_readme_and_screen_gate_point_at_the_new_log(tmp_path):
    dest = _bundle(tmp_path)
    readme = (dest / "README.txt").read_text()
    assert f"$WORKDIR/projects/{PROJECT}/hpc_job_logs/" in readme
    assert "<job_name>.o<jobid>" not in readme
    screen = (dest / SCREEN).read_text()
    assert "hpc_job_logs/<YYYYMMDD_HHMMSS>_<job_name>-prepare.<jobid>.log" in screen


@pytest.mark.parametrize("bad", ["a b", "x;rm -rf", "$(id)", "p/q", ""])
def test_project_name_is_validated(tmp_path, bad):
    if bad == "":
        _bundle(tmp_path, project=bad)       # empty means off, not invalid
        return
    with pytest.raises(ValueError, match="job_log_project"):
        _bundle(tmp_path, project=bad)


def test_slurm_refuses_rather_than_silently_ignoring(tmp_path):
    with pytest.raises(ValueError, match="PBS-only"):
        _bundle(tmp_path, scheduler="slurm")


# ---------------------------------------------------------------------------
# executed
# ---------------------------------------------------------------------------

def test_log_lands_under_workdir_project_with_stripped_array_id(tmp_path):
    dest = _bundle(tmp_path)
    _prepared(dest)
    wd = tmp_path / "wd"
    proc = _run(dest, _shim(tmp_path / "bin"), SCREEN,
                WORKDIR=wd, PBS_JOBID="12345[7].pbs01")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    # Everything went to the file, nothing to the (discarded) job stdout.
    assert proc.stdout == "" and proc.stderr == ""
    logs = _logs(wd / "projects" / PROJECT)
    assert len(logs) == 1, logs
    assert re.fullmatch(r"\d{8}_\d{6}_wt\.12345_7\.pbs01\.log", logs[0].name)
    text = logs[0].read_text()
    assert f"[INFO] Bundle: {dest}" in text
    assert "Prepared checkpoint found" in text
    assert "[WARNING]" not in text
    assert (dest / DONE_SENTINEL).exists()
    # The fallback location was not used.
    assert not (dest / "hpc_job_logs").exists()


def test_proj_overrides_workdir(tmp_path):
    dest = _bundle(tmp_path)
    _prepared(dest)
    proj = tmp_path / "elsewhere"
    proc = _run(dest, _shim(tmp_path / "bin"), SCREEN,
                WORKDIR=tmp_path / "wd", PROJ=proj, PBS_JOBID="1.srv")
    assert proc.returncode == 0
    assert len(_logs(proj)) == 1
    assert _logs(tmp_path / "wd" / "projects" / PROJECT) == []


def test_unset_workdir_falls_back_into_the_bundle_and_says_so(tmp_path):
    """A ${WORKDIR:?} here would kill every link before it could log anywhere."""
    dest = _bundle(tmp_path)
    _prepared(dest)
    proc = _run(dest, _shim(tmp_path / "bin"), SCREEN, PBS_JOBID="2.srv")
    assert proc.returncode == 0, proc.stdout + proc.stderr
    logs = _logs(dest)
    assert len(logs) == 1
    assert "[WARNING] Job log could not go under" in logs[0].read_text()


def test_each_link_gets_its_own_log(tmp_path):
    """The old live log was overwritten by every link; these must not be."""
    dest = _bundle(tmp_path)
    _prepared(dest)
    bin_dir = _shim(tmp_path / "bin")
    wd = tmp_path / "wd"
    assert _run(dest, bin_dir, SCREEN, WORKDIR=wd, PBS_JOBID="10.srv").returncode == 0
    assert _run(dest, bin_dir, SCREEN, WORKDIR=wd, PBS_JOBID="11.srv").returncode == 0
    logs = _logs(wd / "projects" / PROJECT)
    assert [p.name.split("_wt.")[1] for p in logs] == ["10.srv.log", "11.srv.log"]
    assert "Prepared checkpoint found" in logs[0].read_text()
    assert "Run already complete" in logs[1].read_text()


def test_worker_oom_is_still_retried_under_the_redirect(tmp_path):
    """The classifier greps .flexiprobes_link.log; the redirect must not break it."""
    dest = _bundle(tmp_path)
    _prepared(dest)
    wd = tmp_path / "wd"
    proc = _run(dest, _shim(tmp_path / "bin"), SCREEN, WORKDIR=wd,
                PBS_JOBID="20.srv", SHIM_RC=1,
                SHIM_TEXT="concurrent.futures.process.BrokenProcessPool: died")
    assert proc.returncode == 1
    assert not (dest / ABORT_SENTINEL).exists()
    log = _logs(wd / "projects" / PROJECT)[0].read_text()
    assert "BrokenProcessPool" in log            # the pipeline output is in the job log
    assert "worker death" in log                 # ...and the classifier's verdict
    assert "BrokenProcessPool" in (dest / ".flexiprobes_link.log").read_text()


def test_genuine_failure_still_aborts_under_the_redirect(tmp_path):
    """Control for the test above: a real bug must not become retryable."""
    dest = _bundle(tmp_path)
    _prepared(dest)
    wd = tmp_path / "wd"
    proc = _run(dest, _shim(tmp_path / "bin"), SCREEN, WORKDIR=wd,
                PBS_JOBID="21.srv", SHIM_RC=1,
                SHIM_TEXT="RuntimeError: refusing to apply a partial chunk")
    assert proc.returncode == 1
    assert (dest / ABORT_SENTINEL).exists()
    assert "Deterministic failure" in _logs(wd / "projects" / PROJECT)[0].read_text()


def test_prepare_job_logs_under_its_own_name(tmp_path):
    dest = _bundle(tmp_path)
    wd = tmp_path / "wd"
    _run(dest, _shim(tmp_path / "bin"), PREPARE, WORKDIR=wd, PBS_JOBID="30.srv")
    logs = _logs(wd / "projects" / PROJECT)
    assert len(logs) == 1
    assert re.fullmatch(r"\d{8}_\d{6}_wt-prepare\.30\.srv\.log", logs[0].name)
