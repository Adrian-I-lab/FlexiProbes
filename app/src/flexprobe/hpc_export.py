# flexiprobes — probe design for 10x Genomics Flex v1 / Flex v2 / Visium.
# Copyright (C) 2026 Adrian Ilich and contributors
#
# This file is part of flexiprobes.
#
# flexiprobes is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# flexiprobes is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Export a run as a self-contained HPC bundle and re-import its results.

The whole-transcriptome method is too memory-heavy for a laptop but runs fine on
a cluster. Rather than tie execution to the GUI's in-process thread, this module
lets the app write a bundle the user copies to an HPC, edits, and submits:

    <bundle>/
      config.json              serialized RunConfig (paths relative to bundle)
      inputs/<fasta>           input FASTA (only when staged into the bundle)
      output/                  the job writes results here
      run_flexiprobes.slurm    (or .pbs) — scheduler directives + container run
      README.txt               how to edit, submit, and bring results back

The job runs the headless CLI inside an Apptainer/Singularity container:
``flexiprobes design --config config.json``. When it finishes the user copies the
output back (folder or archive) and :func:`import_hpc_output` registers it so the
GUI can review/select probes exactly as for a local run.
"""
from __future__ import annotations

import os
import re
import shlex
import shutil
import tarfile
import uuid
import zipfile
from dataclasses import dataclass
from typing import List, Optional

from . import __version__
from .config_io import config_to_dict, save_config
from .constants import (
    MIN_PROBE_GAP_BP, MIN_PROBE_GAP_FLOOR_BP, PREPARED_MARKER,
    RELAXED_PROBE_GAP_BP,
)
from .models import RunConfig

SCHEDULERS = ("slurm", "pbs")

# Files that, if present in an output dir, mean the run is already finalized
# (probes picked + order CSVs written) rather than awaiting manual review.
_ORDER_FILES = ("idt_opools.csv", "cellranger_probe_set.csv")
RESULTS_FILE = "results.csv"


# Chain defaults. WALLTIME IS PER LINK — the total resumable budget is
# walltime x links, and both halves have to be read together.
#
# 24 h rather than a long single link. A long walltime request queues behind
# everything: schedulers backfill short jobs into gaps, and a request that
# would still be running when a maintenance window opens is simply held until
# after it, which can be a week. A killed link is not lost work — every stage
# checkpoints, and the next link resumes from where it stopped — so the only
# thing a longer link buys is fewer restarts.
DEFAULT_WALLTIME = "24:00:00"
# 12 x 24 h = 288 h of resumable budget, which leaves wide headroom over a
# whole-transcriptome run with the k-mer screen.
DEFAULT_CHAIN_LINKS = 12
DEFAULT_MEM_GB = 28
DEFAULT_CORES = 8

# ---------------------------------------------------------------------------
# Split prepare/screen chain (--split-prepare-screen)
# ---------------------------------------------------------------------------
# One job that does everything has to request the MAXIMUM of what any stage
# needs, for the whole walltime. On a whole transcriptome those maxima come
# from two different stages, and one of them is the reason a single job tends
# to run out of memory.
#
# The same screening work costs several times more memory in a process that
# did the tiling first, because CPython does not return freed pages to the OS
# and the screen then forks its workers off that inflated parent. Splitting
# the run lets each half ask for what it needs:
#
#   prepare  ncpus=8  mem=96gb  walltime=6:00:00   x 2 links
#   screen   ncpus=8  mem=48gb  walltime=24:00:00  x 12 links
#
# These defaults are FlexiProbes sizing choices from whole-transcriptome runs
# (human and mouse). None is a proven ceiling; a larger genome or a different
# cluster may need more.
#
# PREPARE_MEM_GB = 96. The tile+filter peak is the largest memory demand in a
#   run and varies a lot across species (mouse is far above human), so this
#   is the largest observed peak plus headroom. Do not trim it.
#
# PREPARE_CORES = 8. The tile stage is single-threaded and short. The FILTER
#   stage is the expensive half and it is PARALLEL, using several cores on
#   average at --workers 8. With too few cores the slowest filter can exceed
#   the prepare walltime, so this matches the worker count. (A stage that
#   resumed from a checkpoint and filtered nothing reports a near-zero core
#   count. That is noise, not a measurement.)
#
# PREPARE_WALLTIME = 6:00:00. About twice the slowest tile+filter wall clock
#   seen on a whole transcriptome, the margin a stage with a large spread
#   across species deserves.
#
# PREPARE_LINKS = 2. One free retry for a node fault. Prepare is not resumable
#   from a partial state — a link killed mid-tiling starts over — so more than
#   two is buying nothing.
#
# SCREEN_MEM_GB = 48. Set from links that resumed from a checkpoint and
#   therefore measured the screen WITHOUT the tiling high-water mark, at 8
#   screen workers (the parallelism this asks for), with over 2x headroom. A
#   job that cannot be placed is worse than one that runs slowly: very large
#   memory requests can sit unschedulable, and 48 GB places on far more nodes
#   than 128 GB.
#
# SCREEN_CORES = 8 and SCREEN_WALLTIME/LINKS at 24 h x 12, the same as the
#   unsplit chain.
DEFAULT_PREPARE_MEM_GB = 96
DEFAULT_PREPARE_CORES = 8
DEFAULT_PREPARE_WALLTIME = "06:00:00"
DEFAULT_PREPARE_LINKS = 2
DEFAULT_SCREEN_MEM_GB = 48
DEFAULT_SCREEN_CORES = 8
DEFAULT_SCREEN_WALLTIME = DEFAULT_WALLTIME
DEFAULT_SCREEN_LINKS = DEFAULT_CHAIN_LINKS

# Which half of a split chain a generated script is. PHASE_SINGLE is the
# unsplit script, and its output must stay byte-identical so that a bundle
# exported earlier can be regenerated exactly.
PHASE_SINGLE = "single"
PHASE_PREPARE = "prepare"
PHASE_SCREEN = "screen"

PREPARE_SCRIPT_STEM = "run_flexiprobes_prepare"
SCREEN_SCRIPT_STEM = "run_flexiprobes_screen"

# Environment-module name that puts `apptainer` on PATH, baked into the job
# script at export time.
#
# A BATCH JOB DOES NOT INHERIT THE MODULES YOU LOADED INTERACTIVELY, and this
# is the single most common way a resubmission dies: the chain exits in about
# a second, writes the ABORT sentinel, and every remaining link no-ops, so a
# whole chain is lost.
#
# The escape hatch (APPTAINER_MODULE=... ./submit_chain.sh) is exactly the
# kind of thing nobody remembers under time pressure. Recording
# it in the bundle makes remembering unnecessary; the environment variable
# still overrides it, so a bundle moved to another cluster is not stuck with
# the wrong name.
DEFAULT_APPTAINER_MODULE = ""

# Real-time PBS job logs (``job_log_project``). OFF by default.
#
# PBS writes a job's stdout to a spool area on the compute node and copies it
# into place only when the job ENDS, so a running link has nothing readable —
# precisely when a log is wanted. With a project name set, the PBS script
# discards that spooled copy (``#PBS -o /dev/null``) and its first command
# redirects everything into
#
#     $WORKDIR/projects/<project>/hpc_job_logs/<YYYYMMDD_HHMMSS>_<job>.<jobid>.log
#
# which grows while the job runs. Opt-in because it assumes a site layout
# ($WORKDIR/projects/<name>) that other clusters do not have, and because the
# unsplit script's default output must stay byte-identical.
JOB_LOG_DIR = "hpc_job_logs"
_JOB_LOG_PROJECT_RE = re.compile(r"^[A-Za-z0-9._-]+$")


@dataclass
class HpcResources:
    """Scheduler resource request, collected from the GUI form."""
    job_name: str = "flexiprobes"
    # PER LINK, not per run. submit_chain.sh chains `links` copies of the same
    # script with `afterany`, and the run resumes from its checkpoints, so the
    # budget a user actually gets is walltime x links. A short link is usually
    # the better trade: it schedules sooner, it is what fits between cluster
    # maintenance windows, and a walltime kill costs at most one link's worth
    # of un-checkpointed work.
    walltime: str = DEFAULT_WALLTIME
    mem_gb: int = DEFAULT_MEM_GB
    cores: int = DEFAULT_CORES
    # How many links submit_chain.sh queues by default. Over-provisioning is
    # free: once the run finishes, every remaining link sees the DONE sentinel
    # and exits in about a second.
    links: int = DEFAULT_CHAIN_LINKS
    partition: str = ""      # SLURM partition / PBS queue
    account: str = ""        # SLURM --account / PBS -A (allocation)


# ---------------------------------------------------------------------------
# Script templates
# ---------------------------------------------------------------------------

def _slurm_header(res: HpcResources) -> str:
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name={res.job_name}",
        f"#SBATCH --time={res.walltime}",
        f"#SBATCH --mem={res.mem_gb}G",
        f"#SBATCH --cpus-per-task={res.cores}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        "#SBATCH --output=flexiprobes-%j.log",
    ]
    if res.partition:
        lines.insert(2, f"#SBATCH --partition={res.partition}")
    if res.account:
        lines.insert(2, f"#SBATCH --account={res.account}")
    return "\n".join(lines)


def _pbs_header(res: HpcResources, realtime_log: bool = False) -> str:
    lines = [
        "#!/bin/bash",
        f"#PBS -N {res.job_name}",
        f"#PBS -l walltime={res.walltime}",
        f"#PBS -l select=1:ncpus={res.cores}:mem={res.mem_gb}gb",
        # NO fixed "-o flexiprobes.log": every link of the submit_chain.sh chain
        # runs this same script, so a fixed name means the leftover links — which
        # wake after the run finishes, see the DONE sentinel and exit in a second
        # — each copy their two-line no-op output over the real run's log,
        # destroying the stage-by-stage progress record. PBS spools stdout and
        # copies it into place AFTER the job exits, so the script itself cannot
        # prevent this. It has to be handled here, in the directive.
        # Omitting -o gives PBS's default per-job "<job_name>.o<jobid>" in the
        # submit directory, which is what the Slurm header already does via %j.
        "#PBS -j oe",
    ]
    if realtime_log:
        # /dev/null is NOT a shared filename in the sense above: nothing is
        # written to it, so a leftover link has nothing to overwrite. The log
        # the reader wants is the per-link file _realtime_log_block opens.
        lines.append("#PBS -o /dev/null")
    if res.partition:
        lines.insert(2, f"#PBS -q {res.partition}")
    if res.account:
        lines.insert(2, f"#PBS -A {res.account}")
    return "\n".join(lines)


DONE_SENTINEL = ".flexiprobes_DONE"
# A deterministic crash must stop the chain. `afterany` links start regardless
# of the previous exit status, so without this a reproducible bug is retried
# once per link: on 13 Aug one TypeError burned 12 jobs and ~1h36m of 16-core
# walltime for zero progress. The chain resumes resource kills and only those.
ABORT_SENTINEL = ".flexiprobes_ABORT"
STALL_COUNTER = ".flexiprobes_STALL"
# Links that may end having advanced no checkpoint before the chain gives up.
# 1 is a free retry for a transient stall; the 2nd identical no-op is a pattern.
MAX_STALLED_LINKS = 2


def _log_filename_token(job_name: str) -> str:
    """The job name as it appears inside a log filename: shell- and path-safe."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", job_name) or "flexiprobes"


def _realtime_log_block(project: str, job_name: str) -> str:
    """Redirect the whole PBS script into a real-time log under the project root.

    Must sit FIRST in the body, directly under ``set``: the header sends PBS's
    own copy to /dev/null, so anything that fails above the ``exec`` is lost.
    That cost is also why a missing ``$WORKDIR`` falls back to a log inside the
    bundle instead of stopping — a ``${WORKDIR:?}`` here would kill every link
    of the chain before it could write a line anywhere, and before
    ``abort_chain`` exists to stop the rest.

    ``$0``/``BASH_SOURCE`` are useless for locating anything: PBS runs a
    spooled copy of the script, not the file that was submitted.
    """
    template = r"""
# ---------------------------------------------------------------------------
# Real-time job log. KEEP THIS BLOCK FIRST, directly under `set`.
#
# PBS spools a job's output on the compute node and copies it back only when
# the job ENDS, so a running link would have nothing to read. The header sends
# that copy to /dev/null (#PBS -o /dev/null) and the exec below sends
# everything this script prints, from here on, into a file that grows live:
#
#     $PROJECT_ROOT/@LOGDIR@/<YYYYMMDD_HHMMSS>_@JOBNAME@.<jobid>.log
#
# One file per link, so no link overwrites another, and the timestamp prefix
# sorts the folder: `ls "$PROJECT_ROOT/@LOGDIR@" | tail -1` is the newest.
#
# PROJECT_ROOT comes from the environment, never from this file's location:
# PBS runs a spooled COPY of the script, so $0 and BASH_SOURCE point at the
# copy. PROJ wins if set; otherwise $WORKDIR/projects/@PROJECT@.
#
# Anything that fails ABOVE the exec is lost, because stdout is /dev/null until
# then. So an unset WORKDIR or an unwritable project root falls back to
# @LOGDIR@/ inside this bundle, and the log's first lines say so.
PROJECT_ROOT="${PROJ:-${WORKDIR:+${WORKDIR}/projects/@PROJECT@}}"
LOG_DIR="${PROJECT_ROOT:+${PROJECT_ROOT}/@LOGDIR@}"
LOG_FALLBACK_FROM=""
if [ -z "$LOG_DIR" ] || ! mkdir -p "$LOG_DIR"; then
    LOG_FALLBACK_FROM="${LOG_DIR:-(neither PROJ nor WORKDIR is set)}"
    LOG_DIR="${PBS_O_WORKDIR:-$PWD}/@LOGDIR@"
    mkdir -p "$LOG_DIR"
fi
# Array jobs have ids like 12345[0].server; brackets make awkward filenames.
_JOBID="${PBS_JOBID:-nojobid}"; _JOBID="${_JOBID//\[/_}"; _JOBID="${_JOBID//\]/}"
LOG_FILE="${LOG_DIR}/$(date +%Y%m%d_%H%M%S)_@JOBNAME@.${_JOBID}.log"
exec >"${LOG_FILE}" 2>&1
echo "[INFO] Job ${PBS_JOBID:-none} starting at $(date)"
# Every bundle in the project shares this folder, so name the bundle.
echo "[INFO] Bundle: ${PBS_O_WORKDIR:-$PWD}"
echo "[INFO] Log: ${LOG_FILE}"
if [ -n "$LOG_FALLBACK_FROM" ]; then
    echo "[WARNING] Job log could not go under ${LOG_FALLBACK_FROM};"
    echo "[WARNING] writing it inside the bundle instead. Export WORKDIR or PROJ."
fi
"""
    return (template
            .replace("@PROJECT@", project)
            .replace("@JOBNAME@", _log_filename_token(job_name))
            .replace("@LOGDIR@", JOB_LOG_DIR))


def _select_block(fasta_rel: str, n_per_gene: int, min_gap: int,
                  relaxed_gap: int = RELAXED_PROBE_GAP_BP,
                  floor_gap: int = MIN_PROBE_GAP_FLOOR_BP) -> str:
    """Post-pipeline step: pick the top-N spaced probes per gene maximising
    isoform coverage. Emitted only for gene-aware (``group_by_gene``) runs.
    Non-fatal — the full candidate list (results.csv) is the primary output."""
    fasta_q = shlex.quote(fasta_rel)
    return f"""
# Top-{n_per_gene}-per-gene selection (isoform-coverage-maximising). Reads the
# full candidate list (output/results.csv) plus the input FASTA and writes
# output/selected_probes_per_gene.csv + output/flagged_genes.csv. Non-fatal: the
# full candidate list is the primary deliverable, so a failure here only warns.
echo "Selecting top-{n_per_gene} probes per gene (maximising isoform coverage)..."
set +e
"$RUNNER" exec --bind "$PWD" "$CONTAINER" \\
    flexiprobes select-coverage \\
        --results output/results.csv \\
        --fasta {fasta_q} \\
        --out output/selected_probes_per_gene.csv \\
        --flagged output/flagged_genes.csv \\
        --n {n_per_gene} --min-gap {min_gap} \\
        --relaxed-gap {relaxed_gap} --floor-gap {floor_gap}
sel_rc=$?
set -e
if [ "$sel_rc" -ne 0 ]; then
    echo "WARNING: per-gene selection failed (rc=$sel_rc). results.csv is intact;"
    echo "re-run 'flexiprobes select-coverage' or pick probes in the app instead."
fi
"""


def _version_assert_block(expected_version: str) -> str:
    """Abort the job unless the container's engine is the version this bundle
    was written for.

    The ``.sif`` bakes the engine in at build time (``pip install
    /opt/flexiprobe-src`` in flexiprobe.def), so a stale image on the cluster is
    invisible: the job runs happily and produces results from whatever engine
    was frozen into it. That is the single most dangerous failure mode of this
    workflow — an image built before 0.4.0 carries the pre-0.4 GC gate and the
    pre-0.4 off-target rule, so it silently reproduces the OLD, wrong numbers
    under a config that says otherwise. Nothing downstream would flag it: the
    config, the filenames and the import path are all identical.

    So the check is a hard gate, not a warning, and it runs before any work.
    """
    ver_q = shlex.quote(expected_version)
    return f"""
# ---------------------------------------------------------------------------
# Container version assertion — DO NOT REMOVE.
#
# The .sif has the engine baked into it (flexiprobe.def: pip install
# /opt/flexiprobe-src). An image built from an older checkout will run this
# bundle's config without complaint and return results from the OLD engine —
# wrong numbers, right filenames. Rebuild the image whenever the engine version
# changes; see "Rebuild the container" in README.txt.
#
# If you deliberately want a different engine version, edit EXPECTED_VERSION
# below and say so in your methods. There is no environment-variable bypass:
# skipping this check is exactly the mistake it exists to prevent.
EXPECTED_VERSION={ver_q}
IMAGE_VERSION="$("$RUNNER" exec "$CONTAINER" flexiprobes version 2>/dev/null \\
    | tail -n 1 | awk '{{print $NF}}')"
if [ -z "$IMAGE_VERSION" ]; then
    ABORT_RC=3 abort_chain "Could not read the engine version from the \
container (tried: $RUNNER exec $CONTAINER flexiprobes version). Check that \
CONTAINER points at a real flexiprobes image, and that the .sif sits in the \
PARENT of this bundle directory — it is referenced relatively."
fi
if [ "$IMAGE_VERSION" != "$EXPECTED_VERSION" ]; then
    echo "==========================================================" >&2
    echo "ERROR: CONTAINER VERSION MISMATCH — REFUSING TO RUN." >&2
    echo "  bundle expects flexiprobes: $EXPECTED_VERSION" >&2
    echo "  container provides:         $IMAGE_VERSION" >&2
    echo "" >&2
    echo "  The .sif bakes the engine in at build time, so this image" >&2
    echo "  would silently produce results from the WRONG engine under" >&2
    echo "  this bundle's config. Rebuild the image from the $EXPECTED_VERSION" >&2
    echo "  source tree (see 'Rebuild the container' in README.txt) and" >&2
    echo "  resubmit. Do not delete this check." >&2
    echo "==========================================================" >&2
    ABORT_RC=3 abort_chain "Container version mismatch: bundle expects \
$EXPECTED_VERSION, image provides $IMAGE_VERSION."
fi
echo "Container engine version OK: flexiprobes $IMAGE_VERSION"
"""


def _phase_gate_block(phase: str, realtime_log: bool = False) -> str:
    """The extra precondition each half of a split chain checks before working.

    Empty for ``PHASE_SINGLE``, so the unsplit script this module has always
    written is unchanged byte for byte.

    PREPARE gets an idempotence guard. ``submit_chain.sh`` queues more than one
    prepare link (one free retry for a node fault) and they are chained with
    ``afterany``, so link 2 wakes up whether or not link 1 succeeded. Without
    this it would re-tile the transcriptome for nothing.

    SCREEN gets a REFUSAL. This is the important one. The screen job's whole
    premise is that it resumes from a checkpoint somebody else built, and its
    allocation is sized on that basis — 48 GB, against the 96 GB the tiling
    needs. Started without a prepared checkpoint it would not fail cleanly: it
    would happily tile the transcriptome itself, inside an allocation chosen
    for a job that does not tile, and die on memory or walltime some hours
    later. So it stops at second zero and says why. This is an ABORT rather
    than a plain exit because every screen link would hit the same wall, and
    the fix is to resubmit the prepare half, not to queue more screen links.
    """
    if phase == PHASE_PREPARE:
        return f"""
PREPARED="output/{PREPARED_MARKER}"

# Already prepared by an earlier link. submit_chain.sh queues a spare prepare
# link as a free retry against a node fault, and `afterany` starts it whatever
# happened before, so this is the normal path for link 2 of 2 — not an error.
if [ -f "$PREPARED" ]; then
    echo "Candidates already prepared ($PREPARED present) — nothing to do."
    echo "The screen job resumes from output/.candidate_cache/."
    exit 0
fi
"""
    if phase == PHASE_SCREEN:
        prepare_log = (
            f"({JOB_LOG_DIR}/<YYYYMMDD_HHMMSS>_<job_name>-prepare.<jobid>.log \
under the project root)" if realtime_log else
            "(<job_name>-prepare.o<jobid>)")
        return f"""
PREPARED="output/{PREPARED_MARKER}"

# REFUSE to start without a prepared checkpoint. This job asks the scheduler
# for a screen-sized allocation on the explicit assumption that ingest, tiling
# and filtering have already happened in the prepare job. Run it anyway and it
# does NOT fail fast — it re-tiles the whole transcriptome inside an allocation
# that was never sized for it and dies on memory or walltime hours later, which
# is exactly the failure this split exists to remove.
if [ ! -f "$PREPARED" ]; then
    ABORT_RC=4 abort_chain "No prepared candidate checkpoint: $PREPARED is \
missing, so the prepare job has not finished successfully. This screen job is \
sized for a run that RESUMES from that checkpoint and must not tile the \
transcriptome itself. Check the prepare job's log \
{prepare_log}, fix what stopped it, then: rm -f \
{ABORT_SENTINEL} && ./submit_chain.sh"
fi
echo "Prepared checkpoint found:"
sed 's/^/    /' "$PREPARED"
"""
    return ""


def _design_invocation(phase: str, cores: int) -> str:
    """The `flexiprobes design` command line for this phase.

    PREPARE passes ``--prepare-only`` and no ``--screen-workers``: nothing
    forks, because the run stops before the screen. ``--workers`` still matters
    and is not 1. The FILTER stage is parallel and uses several cores on
    average, which is why the prepare job asks for more than one
    core in the first place.
    """
    if phase == PHASE_PREPARE:
        return (f"    flexiprobes design --config config.json "
                f"--prepare-only --workers {cores} --force")
    return (f"    flexiprobes design --config config.json "
            f"--workers {cores} --screen-workers {cores} --force")


def _success_tail(phase: str, select_block: str) -> str:
    """What a link does after the pipeline exits 0.

    PREPARE writes NO ``.flexiprobes_DONE`` and packages nothing. The marker it
    reports on was written by the pipeline itself, inside output/, and says
    only that the front half finished — see ``constants.PREPARED_MARKER`` for
    why using DONE here would strand the bundle looking complete and stop the
    screen job from ever running.
    """
    if phase == PHASE_PREPARE:
        return f"""
# Prepared, NOT finished. No {DONE_SENTINEL} and no tarball: this run has
# designed no probes and there is no results.csv. The pipeline wrote
# output/{PREPARED_MARKER} itself once the candidate checkpoint and the
# pre-screen artefacts were safely on disk; the screen job gates on it.
if [ ! -f "$PREPARED" ]; then
    ABORT_RC=5 abort_chain "flexiprobes design --prepare-only exited 0 but \
did not write $PREPARED. Treating that as a bug rather than a success: the \
screen job would refuse to start and the chain would stop anyway, so stop \
here where the cause is still in this log."
fi
echo "Prepare complete. The screen job will resume from output/.candidate_cache/."
"""
    return f"""{select_block}
# Success: mark done (so later chained links no-op) and package the results.
# Remove any stale archive from a previous submission first.
touch "$DONE"
rm -f flexiprobes_output.tar.gz
tar -czf flexiprobes_output.tar.gz -C output .
echo "Done. Copy flexiprobes_output.tar.gz back and import it in flexiprobes."
"""


def _run_notes(phase: str, cores: int) -> tuple:
    """The two comment blocks around the pipeline invocation.

    Phase-specific because they describe the parallelism, and the prepare job
    has none of the kind the original text is about. Leaving the unsplit
    wording in the prepare script would have it explain ``--screen-workers`` in
    a job that deliberately does not pass it — the sort of stale comment that
    later gets "fixed" by adding the flag back.

    ``PHASE_SINGLE`` returns the original text verbatim, blank lines and the
    two sentences that got interleaved in it included, so the unsplit script
    stays byte-identical.
    """
    if phase == PHASE_PREPARE:
        run_note = (
            "# Run the FRONT HALF only: ingest, tiling and filtering.\n"
            "# --prepare-only writes the pre-screen candidate checkpoint\n"
            "# (output/.candidate_cache/) plus the pre-screen artefacts and then\n"
            "# exits, BEFORE the off-target screen. No probes are designed here and\n"
            "# there is no results.csv; the screen job produces those.\n"
            "# --workers parallelises the FILTER stage, which is where this job's\n"
            "# time goes (it uses several cores on average on a whole transcriptome).\n"
            "# Tiling itself is single-threaded but costs under 90 s.\n"
            "# --bind \"$PWD\" makes the bundle (config.json, inputs/, output/) visible\n"
            "# inside the container. Add more --bind paths if your input FASTA lives\n"
            "# elsewhere. Do NOT delete output/ between attempts.\n"
        )
        thread_note = (
            f"# One source of parallelism only. --workers uses {cores} processes, and\n"
            f"# each would otherwise start its own BLAS thread pool — {cores}x{cores}\n"
            f"# threads on {cores} cores, which is slower than running serially.\n"
        )
        return run_note, thread_note
    if phase == PHASE_SCREEN:
        run_note = (
            "# Run the BACK HALF: resume the candidate checkpoint the prepare job\n"
            "# left in output/.candidate_cache/, then screen, rank and export. The\n"
            "# resume is what keeps this job's memory down — a process that tiled\n"
            "# first carries that high-water mark for the rest of its life, and the\n"
            "# screen forks its workers off it.\n"
            f"# --screen-workers fans the OFF-TARGET SCREEN across {cores} processes,\n"
            "# which is where essentially all of this job's walltime goes. It is NOT\n"
            "# in any cache fingerprint, so changing it resumes an in-flight run\n"
            "# rather than invalidating it: a chain already part-way through can be\n"
            "# edited to a different value and resubmitted. Every screen chunk is\n"
            "# checkpointed (output/.screen_cache/), so a job killed by walltime\n"
            "# resumes where it stopped on the next chained link — specificity is\n"
            "# scored chunk-by-chunk so memory stays bounded.\n"
            "# --bind \"$PWD\" makes the bundle (config.json, inputs/, output/) visible\n"
            "# inside the container. Add more --bind paths if your input FASTA or\n"
            "# BLAST DB live elsewhere. Do NOT delete output/ between attempts.\n"
        )
        thread_note = (
            f"# One source of parallelism only. --screen-workers forks {cores} workers,\n"
            f"# and each would otherwise start its own BLAS thread pool — {cores}x{cores}\n"
            f"# threads on {cores} cores, which is slower than running serially.\n"
        )
        return run_note, thread_note
    # PHASE_SINGLE: kept verbatim, including the order of its sentences, so
    # the unsplit script stays byte-identical.
    run_note = (
        "# Run the pipeline. It checkpoints the pre-screen candidate set\n"
        "# --screen-workers fans the OFF-TARGET SCREEN across that many processes, which\n"
        "# on a whole-transcriptome run is where essentially all the walltime goes\n"
        "# (--workers only covers tiling+filtering). It is NOT in any cache fingerprint,\n"
        "# so changing it resumes an in-flight run rather than invalidating it: a chain\n"
        "# already part-way through can be edited to a different value and resubmitted.\n"
        "# (output/.candidate_cache/) and every BLAST chunk (output/.screen_cache/), so a\n"
        "# job killed by walltime resumes where it stopped on the next chained link —\n"
        "# specificity is scored chunk-by-chunk so memory stays bounded.\n"
        "# --bind \"$PWD\" makes the bundle (config.json, inputs/, output/) visible inside\n"
        "# the container. Add more --bind paths if your input FASTA or BLAST DB live\n"
        "# elsewhere on the cluster filesystem. Do NOT delete output/ between attempts.\n"
    )
    thread_note = (
        f"# One source of parallelism only. --screen-workers forks {cores} workers, and\n"
        f"# each would otherwise start its own BLAS thread pool — {cores}x{cores} threads\n"
        f"# on {cores} cores, which is slower than running serially.\n"
    )
    return run_note, thread_note


def _script_body(scheduler: str, image: str, cores: int, *,
                 select_fasta: Optional[str] = None,
                 select_n: int = 3,
                 select_min_gap: int = MIN_PROBE_GAP_BP,
                 select_relaxed_gap: int = RELAXED_PROBE_GAP_BP,
                 select_floor_gap: int = MIN_PROBE_GAP_FLOOR_BP,
                 expected_version: str = __version__,
                 phase: str = PHASE_SINGLE,
                 apptainer_module: str = DEFAULT_APPTAINER_MODULE,
                 job_log_project: str = "",
                 job_name: str = "flexiprobes") -> str:
    """The body of a generated job script.

    ``phase`` selects which of the three scripts this is. PHASE_SINGLE is the
    original, unsplit job and its output is unchanged byte for byte, so an
    exported bundle can be regenerated exactly. The failure handling it encodes
    (the DONE no-op, the ABORT stop, the stall counter, the worker-OOM
    classifier) covers faults that otherwise kill cluster jobs. The two split
    phases reuse all of it and differ only in three places: an extra
    precondition (``_phase_gate_block``), the design command line
    (``_design_invocation``) and what happens on success (``_success_tail``).

    ``job_log_project`` (PBS only) prepends ``_realtime_log_block``. Empty, the
    body is unchanged byte for byte.
    """
    realtime_log = bool(job_log_project) and scheduler == "pbs"
    log_block = (
        _realtime_log_block(job_log_project, job_name) if realtime_log else "")
    workdir = "$SLURM_SUBMIT_DIR" if scheduler == "slurm" else "$PBS_O_WORKDIR"
    # Quoted: it lands inside a shell assignment in the generated script.
    module_q = shlex.quote(apptainer_module or "")
    # shlex.quote the user-supplied image path so a path with spaces or shell
    # metacharacters can't break (or inject into) the generated script.
    image_q = shlex.quote(image)
    # Gene-aware runs append a coverage-selection step after a successful run so
    # the bundle ships both the full candidate list and the top-N picks. Never
    # on the prepare phase: there is no results.csv to select from yet.
    select_block = (
        _select_block(select_fasta, select_n, select_min_gap,
                      select_relaxed_gap, select_floor_gap)
        if (select_fasta and phase != PHASE_PREPARE) else ""
    )
    version_block = _version_assert_block(expected_version)
    phase_gate = _phase_gate_block(phase, realtime_log)
    design_cmd = _design_invocation(phase, cores)
    success_tail = _success_tail(phase, select_block)
    run_note, thread_note = _run_notes(phase, cores)
    return f"""
set -uo pipefail
{log_block}
# Move to the bundle directory the job was submitted from.
cd "{workdir}"

DONE="{DONE_SENTINEL}"

# Resumable chain: each chained link runs this same script. If a previous link
# already finished the run, this one is a no-op — so the tail of the chain costs
# nothing once the work is done. (Delete this file + output/ to force a rerun.)
if [ -f "$DONE" ]; then
    echo "Run already complete ($DONE present) — nothing to do."
    echo "(This is a leftover chain link. The finished run's log is a separate"
    echo " per-job file — see 'Where the logs go' in README.txt.)"
    exit 0
fi

ABORT="{ABORT_SENTINEL}"
STALL="{STALL_COUNTER}"

# A previous link hit a fault that resuming cannot fix. Every remaining link
# no-ops rather than reproducing it, so the chain stops costing allocation the
# moment the diagnosis is in. Fix the cause, delete $ABORT, resubmit.
if [ -f "$ABORT" ]; then
    echo "=========================================================="
    echo "Chain ABORTED by an earlier link — not retrying."
    echo ""
    cat "$ABORT"
    echo ""
    echo "Resuming cannot clear this: the same input reproduces it. Fix the"
    echo "cause, 'rm $ABORT', and resubmit the chain. Checkpoints in output/"
    echo "are intact, so the rerun picks up where the run stopped."
    echo "=========================================================="
    exit 0
fi

# Cheap, monotone summary of checkpoint state. Both components only grow as the
# run advances, so an unchanged fingerprint across a whole link means that link
# produced nothing — the signature of a crash loop, and independent of whatever
# the error text happens to say.
progress_fingerprint() {{
    local n b
    n=$(find output -type f 2>/dev/null | wc -l | tr -d " ")
    b=$(du -sk output 2>/dev/null | cut -f1)
    printf "%s %s" "${{n:-0}}" "${{b:-0}}"
}}

# Stop the whole chain, not just this link. Every caller below is a fault that
# an identical resubmission cannot fix — a missing runtime, the wrong image, a
# deterministic crash — so the remaining links must no-op instead of
# reproducing it. This is defined BEFORE the preflight checks on purpose: the
# first version of this script classified only post-pipeline failures, so a
# preflight exit still burned every link in the chain (observed 18 Aug, 12 jobs
# in 70 seconds on a missing `module load`).
abort_chain() {{
    {{
      echo "$@"
      echo "Job: ${{PBS_JOBID:-${{SLURM_JOB_ID:-unknown}}}}"
    }} > "$ABORT"
    echo "ERROR: $* — aborting the chain; remaining links will no-op." >&2
    exit "${{ABORT_RC:-1}}"
}}

# Locate the container runtime.
#
# A BATCH JOB DOES NOT INHERIT THE MODULES YOU LOADED INTERACTIVELY. If
# apptainer is only on PATH after a `module load` at your site — which is
# usual — name the module here, or pass it at submit time:
#
#     APPTAINER_MODULE=apptainer/1.4.1 ./submit_chain.sh
#
# The bundle's own baked-in name is below; the environment variable overrides
# it, so the same bundle works on a cluster that spells the module differently.
BUNDLE_APPTAINER_MODULE={module_q}
APPTAINER_MODULE="${{APPTAINER_MODULE:-$BUNDLE_APPTAINER_MODULE}}"
RUNNER="$(command -v apptainer || command -v singularity || true)"
if [ -z "$RUNNER" ] && [ -n "$APPTAINER_MODULE" ]; then
    echo "No container runtime on PATH; trying: module load $APPTAINER_MODULE"
    module load "$APPTAINER_MODULE" 2>/dev/null || true
    RUNNER="$(command -v apptainer || command -v singularity || true)"
fi
# Last resort before giving up: the UNVERSIONED module name. Sites almost
# always provide `apptainer` as an alias for their current default, and a
# version bump that retires `apptainer/1.4.1` would otherwise abort a chain
# whose only real problem is a stale string in a config file.
if [ -z "$RUNNER" ]; then
    for _m in apptainer singularity; do
        echo "Still no runtime; trying: module load $_m"
        module load "$_m" 2>/dev/null || true
        RUNNER="$(command -v apptainer || command -v singularity || true)"
        [ -n "$RUNNER" ] && {{ APPTAINER_MODULE="$_m"; break; }}
    done
fi
if [ -z "$RUNNER" ]; then
    abort_chain "No container runtime: neither 'apptainer' nor 'singularity' is \
on PATH$([ -n "$APPTAINER_MODULE" ] && echo ", and 'module load $APPTAINER_MODULE' did not provide one"), \
and a bare 'module load apptainer' / 'module load singularity' did not either. \
A BATCH JOB DOES NOT INHERIT THE MODULES YOU LOADED INTERACTIVELY, so loading \
it in your shell before submitting does nothing. Find the name with \
'module avail apptainer' on a login node, then either re-export the bundle \
with --apptainer-module NAME or resubmit as: \
APPTAINER_MODULE=NAME ./submit_chain.sh"
fi
echo "Container runtime: $RUNNER${{APPTAINER_MODULE:+ (module $APPTAINER_MODULE)}}"

CONTAINER={image_q}
{version_block}
mkdir -p output
{phase_gate}
{run_note}before_fp="$(progress_fingerprint)"
set +e
{thread_note}export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# The link's output is TEED to a file as well as to the job log, because the
# retry classifier below has to READ it. PIPESTATUS keeps rc as the pipeline's
# own exit status rather than tee's, which is always 0.
LINK_LOG=".flexiprobes_link.log"
"$RUNNER" exec --bind "$PWD" "$CONTAINER" \
{design_cmd} 2>&1 | tee "$LINK_LOG"
rc=${{PIPESTATUS[0]}}
set -e
after_fp="$(progress_fingerprint)"

if [ "$rc" -ne 0 ]; then
    # Did this link advance any checkpoint at all?
    if [ "$after_fp" = "$before_fp" ]; then
        stalled=$(( $(cat "$STALL" 2>/dev/null || echo 0) + 1 ))
        echo "$stalled" > "$STALL"
    else
        stalled=0
        rm -f "$STALL"
    fi

    # A signal (rc > 128) is the scheduler or the kernel killing us: walltime,
    # OOM, node fault. Those are exactly the faults a resume fixes. Anything
    # else is the pipeline deciding to exit — a bug, bad config, missing input
    # — and is reproducible by construction, so the chain must not retry it.
    #
    # THE EXIT STATUS ALONE IS NOT ENOUGH. With --screen-workers the OOM
    # killer takes a WORKER, not the parent. The parent sees its pool die,
    # raises concurrent.futures.process.BrokenProcessPool, unwinds cleanly and
    # exits rc=1. An out-of-memory kill therefore arrives wearing the exit
    # status of a logic bug. PBS appends its own "cgroup/OOM" line only AFTER
    # this script has decided, so that is no help either. Trusting rc alone
    # would abort the rest of the chain on a fault a resume would clear.
    #
    # So a low exit status is checked against the log before it is trusted. The
    # distinction is a pool that DIED (external kill — retryable) against an
    # exception a worker RAISED (a real screen failure — not retryable). Those
    # produce different strings, which is why this matches the pool and
    # allocation failures specifically rather than "Traceback".
    #
    # The match must stay this narrow, so a real bug is never retried forever.
    if [ "$rc" -le 128 ] && grep -qE 'BrokenProcessPool|MemoryError|Cannot allocate memory|std::bad_alloc' "$LINK_LOG" 2>/dev/null; then
        echo "Pipeline exited rc=$rc, but the log shows a worker death or an"
        echo "allocation failure, not a deterministic error. Treating it as"
        echo "retryable so the next link resumes from the last checkpoint:"
        grep -nE 'BrokenProcessPool|MemoryError|Cannot allocate memory|std::bad_alloc' "$LINK_LOG" | tail -3
        echo "If two consecutive links stop at the SAME checkpoint, the stall"
        echo "guard below aborts the chain. That is the signal that this"
        echo "memory/worker combination genuinely cannot run this chunk, and"
        echo "that the fix is more memory or fewer --screen-workers."
    elif [ "$rc" -le 128 ]; then
        ABORT_RC=$rc abort_chain "Deterministic failure: flexiprobes exited \
rc=$rc (not a signal). The pipeline chose to exit, so every link would fail \
identically. See the log above."
    fi

    # Retryable, but produced nothing — twice. Resuming a link that cannot
    # advance just re-buys the same walltime, so stop. This is also the guard
    # that catches a memory request which is simply too small: an OOM that
    # always lands on the same chunk makes no progress, twice, and stops here.
    if [ "$stalled" -ge {MAX_STALLED_LINKS} ]; then
        ABORT_RC=$rc abort_chain "No progress in $stalled consecutive links \
(last rc=$rc). Checkpoint state was identical before and after each attempt, \
so resuming is not advancing the run — more walltime will not help."
    fi

    echo "Pipeline exited rc=$rc (retryable — walltime/OOM/node fault). Checkpoints"
    echo "kept and progress was made, so the next chained link will resume."
    echo "If this was the last link, resubmit the chain."
    exit "$rc"
fi

rm -f "$STALL"
{success_tail}""".lstrip("\n")


def _submit_chain_script(scheduler: str,
                         default_links: int = DEFAULT_CHAIN_LINKS,
                         walltime: str = DEFAULT_WALLTIME,
                         apptainer_module: str = DEFAULT_APPTAINER_MODULE) -> str:
    """A small wrapper that submits a chain of K linked jobs.

    Each link depends on the previous one with ``afterany`` so it starts whether
    or not the previous link succeeded; the per-job ``DONE`` guard makes finished
    links no-op. This is the robust, cluster-friendly form of "resubmit until it
    finishes" — no qsub-from-compute-node required.
    """
    env_flag = "--export=ALL," if scheduler == "slurm" else "-v"
    module_q = shlex.quote(apptainer_module or "")
    if scheduler == "slurm":
        script, submit_one, submit_dep = (
            "run_flexiprobes.slurm",
            'sbatch --parsable $ENV_ARGS "$SCRIPT"',
            'sbatch --parsable $ENV_ARGS --dependency=afterany:"$prev" "$SCRIPT"',
        )
    else:
        script, submit_one, submit_dep = (
            "run_flexiprobes.pbs",
            'qsub $ENV_ARGS "$SCRIPT"',
            'qsub $ENV_ARGS -W depend=afterany:"$prev" "$SCRIPT"',
        )
    return f"""#!/bin/bash
# Submit a resumable chain of K jobs (default {default_links}). Each link resumes
# from the last checkpoint; once the run finishes, the remaining links detect the
# {DONE_SENTINEL} sentinel and exit immediately, so over-provisioning K is free.
#
# THE BUDGET IS K x THE PER-LINK WALLTIME. This bundle asks for {walltime} per
# link, so the default K gives {default_links} x {walltime}. Raise K, not the
# walltime in {script}: a short link is backfilled into gaps the scheduler
# cannot give a long one, and it is what fits before a maintenance window —
# a job whose walltime would overrun the window is held until after it.
#
#   ./submit_chain.sh          # {default_links} links x {walltime}
#   ./submit_chain.sh 24       # 24 links (longer total resumable wall-time)
set -euo pipefail

K="${{1:-{default_links}}}"
SCRIPT="{script}"

# A batch job does not inherit your interactive environment, so if the container
# runtime only appears after a `module load` at your site, pass the module name
# through to every link:
#
#     APPTAINER_MODULE=apptainer/1.4.1 ./submit_chain.sh
#
ENV_ARGS=""
if [ -n "${{APPTAINER_MODULE:-}}" ]; then
    ENV_ARGS="{env_flag} APPTAINER_MODULE=$APPTAINER_MODULE"
    echo "Forwarding APPTAINER_MODULE=$APPTAINER_MODULE to every link."
elif [ -n {module_q} ]; then
    echo "Container module baked into the job script: {apptainer_module}"
else
    # Neither source names one. That is not necessarily wrong — some sites put
    # apptainer on the default PATH — but it IS the documented commonest way a
    # chain dies, so say it before spending K submissions rather than after.
    echo "NOTE: no container module named (bundle or APPTAINER_MODULE)."
    echo "      If apptainer is not on a BATCH job's PATH at your site, every"
    echo "      link will abort in about a second. A batch job does not inherit"
    echo "      the modules you loaded in this shell. Check with:"
    echo "          module avail apptainer"
    echo "      then resubmit as: APPTAINER_MODULE=<name> ./submit_chain.sh"
fi

# Refuse to queue into an aborted bundle. The sentinel makes EVERY link no-op,
# so submitting K jobs on top of one is K wasted submissions and a log
# directory that looks like a fresh failure. The previous run wrote the reason
# into the file — show it, and stop.
if [ -f "{ABORT_SENTINEL}" ]; then
    echo "REFUSING TO SUBMIT: {ABORT_SENTINEL} is present, so every link would" >&2
    echo "no-op. The previous attempt stopped for this reason:" >&2
    echo >&2
    sed 's/^/    /' "{ABORT_SENTINEL}" >&2
    echo >&2
    echo "Fix the cause, then: rm {ABORT_SENTINEL} && ./submit_chain.sh" >&2
    echo "Checkpoints in output/ are intact; the rerun resumes from them." >&2
    exit 1
fi

prev=""
for i in $(seq 1 "$K"); do
    if [ -z "$prev" ]; then
        prev=$({submit_one})
    else
        prev=$({submit_dep})
    fi
    echo "link $i/$K submitted: $prev"
done
echo "Chain submitted. It self-continues across walltime kills until the run"
echo "completes; spare links no-op once {DONE_SENTINEL} appears."
""".lstrip("\n")


def _split_submit_chain_script(
    scheduler: str,
    *,
    prepare_links: int = DEFAULT_PREPARE_LINKS,
    screen_links: int = DEFAULT_SCREEN_LINKS,
    prepare_walltime: str = DEFAULT_PREPARE_WALLTIME,
    screen_walltime: str = DEFAULT_SCREEN_WALLTIME,
    prepare_mem_gb: int = DEFAULT_PREPARE_MEM_GB,
    prepare_cores: int = DEFAULT_PREPARE_CORES,
    screen_mem_gb: int = DEFAULT_SCREEN_MEM_GB,
    screen_cores: int = DEFAULT_SCREEN_CORES,
    apptainer_module: str = DEFAULT_APPTAINER_MODULE,
) -> str:
    """Submit the split chain: P prepare links, then K screen links.

    ONE chain, not two submissions. The screen links depend on the last prepare
    link with the same ``afterany`` the unsplit chain uses, so they queue
    immediately and start whatever the prepare half did — and refuse to work
    when it did not finish (see ``_phase_gate_block``). That ordering matters
    on a busy queue: submitting the screen half only after the prepare half
    succeeds would put it at the back of the queue hours later, and the point
    of the split is to finish sooner, not to add a manual step.

    Everything the unsplit chain guarantees still holds. A finished run's
    leftover links no-op on the DONE sentinel, an ABORT stops every remaining
    link, and the stall counter catches a link that advanced nothing.
    """
    env_flag = "--export=ALL," if scheduler == "slurm" else "-v"
    module_q = shlex.quote(apptainer_module or "")
    ext = "slurm" if scheduler == "slurm" else "pbs"
    prepare_script = f"{PREPARE_SCRIPT_STEM}.{ext}"
    screen_script = f"{SCREEN_SCRIPT_STEM}.{ext}"
    if scheduler == "slurm":
        submit_one = 'sbatch --parsable $ENV_ARGS "$1"'
        submit_dep = 'sbatch --parsable $ENV_ARGS --dependency=afterany:"$2" "$1"'
    else:
        submit_one = 'qsub $ENV_ARGS "$1"'
        submit_dep = 'qsub $ENV_ARGS -W depend=afterany:"$2" "$1"'
    return f"""#!/bin/bash
# Submit the SPLIT chain: {prepare_links} prepare link(s), then K screen links
# (default {screen_links}), every link chained with `afterany`.
#
# WHY TWO JOB TYPES. One job has to request the maximum of what any stage needs,
# for the whole walltime. Those maxima come from two different stages:
#
#   {prepare_script}
#       ingest + tile + filter, then STOP. Writes the candidate checkpoint,
#       the pre-screen artefacts and output/{PREPARED_MARKER}.
#       ncpus={prepare_cores} mem={prepare_mem_gb}gb walltime={prepare_walltime}
#       Memory-hungry: this is the run's largest memory peak.
#
#   {screen_script}
#       resumes from that checkpoint and runs the off-target screen, ranking,
#       chemistry and export. Produces results.csv and the tarball.
#       ncpus={screen_cores} mem={screen_mem_gb}gb walltime={screen_walltime}
#       The same screening work needs several times less memory here than
#       in a process that tiled first.
#
# THE BUDGET IS K x {screen_walltime} FOR THE SCREEN HALF. Raise K, not the
# walltime: a short link is backfilled into gaps a long one cannot use, and a
# job whose walltime would overrun a maintenance window is held until after it.
#
#   ./submit_chain.sh          # {prepare_links} prepare + {screen_links} screen links
#   ./submit_chain.sh 24       # {prepare_links} prepare + 24 screen links
set -euo pipefail

K="${{1:-{screen_links}}}"
P={prepare_links}
PREPARE_SCRIPT="{prepare_script}"
SCREEN_SCRIPT="{screen_script}"

# A batch job does not inherit your interactive environment, so if the container
# runtime only appears after a `module load` at your site, pass the module name
# through to every link:
#
#     APPTAINER_MODULE=apptainer/1.4.1 ./submit_chain.sh
#
ENV_ARGS=""
if [ -n "${{APPTAINER_MODULE:-}}" ]; then
    ENV_ARGS="{env_flag} APPTAINER_MODULE=$APPTAINER_MODULE"
    echo "Forwarding APPTAINER_MODULE=$APPTAINER_MODULE to every link."
elif [ -n {module_q} ]; then
    echo "Container module baked into the job scripts: {apptainer_module}"
else
    echo "NOTE: no container module named (bundle or APPTAINER_MODULE)."
    echo "      If apptainer is not on a BATCH job's PATH at your site, every"
    echo "      link will abort in about a second. A batch job does not inherit"
    echo "      the modules you loaded in this shell. Check with:"
    echo "          module avail apptainer"
    echo "      then resubmit as: APPTAINER_MODULE=<name> ./submit_chain.sh"
fi

# Refuse to queue into an aborted bundle. The sentinel makes EVERY link no-op,
# so submitting on top of one is a wasted submission and a log directory that
# looks like a fresh failure. The previous run wrote the reason into the file.
if [ -f "{ABORT_SENTINEL}" ]; then
    echo "REFUSING TO SUBMIT: {ABORT_SENTINEL} is present, so every link would" >&2
    echo "no-op. The previous attempt stopped for this reason:" >&2
    echo >&2
    sed 's/^/    /' "{ABORT_SENTINEL}" >&2
    echo >&2
    echo "Fix the cause, then: rm {ABORT_SENTINEL} && ./submit_chain.sh" >&2
    echo "Checkpoints in output/ are intact; the rerun resumes from them." >&2
    exit 1
fi

submit() {{      # $1 = script, $2 = dependency job id ("" for none)
    if [ -z "$2" ]; then
        {submit_one}
    else
        {submit_dep}
    fi
}}

prev=""
for i in $(seq 1 "$P"); do
    prev=$(submit "$PREPARE_SCRIPT" "$prev")
    echo "prepare link $i/$P submitted: $prev"
done
for i in $(seq 1 "$K"); do
    prev=$(submit "$SCREEN_SCRIPT" "$prev")
    echo "screen link $i/$K submitted: $prev"
done
echo
echo "Split chain submitted: $P prepare + $K screen links."
echo "The screen links are queued NOW and start after the prepare half"
echo "whatever it did — if it did not finish, the first screen link stops the"
echo "chain immediately rather than re-tiling in a screen-sized allocation."
echo "Spare links no-op once {DONE_SENTINEL} appears."
""".lstrip("\n")


def _rebuild_block(image: str, expected_version: str) -> str:
    """The container-rebuild warning + recipe.

    Deliberately the FIRST thing in the README and phrased as an instruction,
    not a note: a stale ``.sif`` is the one failure mode of this workflow that
    produces plausible-looking wrong answers instead of an error.
    """
    return f"""Rebuild the container FIRST — this is not optional
--------------------------------------------------
This bundle expects flexiprobes {expected_version}.

The .sif does NOT read your source tree at run time. flexiprobe.def bakes the
engine into the image (`pip install /opt/flexiprobe-src`), so an image built
from an older checkout keeps running the OLD engine — with this bundle's config,
this bundle's filenames, and this bundle's output layout. The results would look
completely normal and be wrong. If a .sif already exists at

    {image}

assume it is stale and rebuild it, unless you built it from a {expected_version}
checkout yourself.

The job script asserts this for you: it runs `flexiprobes version` inside the
container and refuses to start (exit 3) on any mismatch. Do not delete that
check — read the error and rebuild.

Build it (once per engine version), from the repository's app/ directory:

    cd /path/to/FlexiProbes/app

    # Apptainer's build + cache dirs MUST be on node-local storage. On a Lustre
    # working filesystem the build blows through the directory quota part-way
    # and dies with confusing "no space left"/quota errors — the temporary
    # layers are many GB. /tmp on the build node is the right place.
    export APPTAINER_TMPDIR=/tmp/$USER/apptainer-tmp
    export APPTAINER_CACHEDIR=/tmp/$USER/apptainer-cache
    mkdir -p "$APPTAINER_TMPDIR" "$APPTAINER_CACHEDIR"

    apptainer build packaging/container/flexiprobe.sif \\
        packaging/container/flexiprobe.def

    # Verify before you copy it anywhere:
    apptainer exec packaging/container/flexiprobe.sif flexiprobes version
    #   -> must print exactly: flexiprobes {expected_version}

Then put the image where this bundle expects it ({image}), or edit the
CONTAINER=... line in the job script.

"""


def _prescreen_outputs_block(capture: str) -> str:
    """Describe the pre-screen candidate artefacts this bundle will produce.

    Stated in the bundle's own README, not just in the code, because "did this
    run keep the candidates it considered, or only the ones that survived?" is
    the question that decides whether the rank statistics can be corrected — and
    the previous runs answered it wrongly and silently.
    """
    if capture == "none":
        return """
  (This run does NOT retain the pre-screen candidate set. results.csv holds only
   the probes that PASSED the off-target screen, so any statistic computed from
   it is conditioned on that outcome. Re-export the bundle with
   --emit-prescreen-candidates if you need the denominator.)"""
    common = """
  output/candidates_prescreen.csv.gz     Every candidate the off-target screen
                                        actually saw (i.e. after the
                                        screen_top_per_gene cap). The
                                        denominator for any claim about the
                                        screen — results.csv is the survivors.
  output/prescreen_cap_summary.csv.gz    Per gene: how many candidates passed
                                        the chemistry filters BEFORE the cap vs
                                        how many the cap passed to BLAST.
  output/prescreen_index.json            States which candidate set this run
                                        captured, the schema of each file, and
                                        the reconstruction rules."""
    if capture == "post_cap":
        return common + """

  CAPTURE MODE: post-cap. The per-gene cap ran before BLAST and only its
  survivors were recorded row by row; model the cap from the summary above."""
    return common + """
  output/candidates_precap.csv.gz       EVERY window that passed the chemistry
                                        filters, before the cap. Compact
                                        schema (see prescreen_index.json):
                                        dictionary-encoded transcript index,
                                        start, pre-screen score, junction base.
  output/prescreen_transcripts.csv.gz    The transcript dictionary for the file
                                        above.

  CAPTURE MODE: pre-cap. Both selection steps applied before BLAST — the
  chemistry filters and the per-gene cap — are fully reconstructable from these
  files, so a null model does not have to assume either of them away."""


def _readme(scheduler: str, image: str, staged_input: bool,
            db_paths: List[str], gene_aware: bool = False,
            expected_version: str = __version__,
            prescreen_capture: str = "none",
            links: int = DEFAULT_CHAIN_LINKS,
            walltime: str = DEFAULT_WALLTIME,
            split: bool = False,
            prepare_res: Optional["HpcResources"] = None,
            screen_res: Optional["HpcResources"] = None,
            apptainer_module: str = DEFAULT_APPTAINER_MODULE,
            job_log_project: str = "") -> str:
    module_note = (
        f"This bundle loads the container runtime with "
        f"`module load {apptainer_module}`.\n"
        if apptainer_module else
        "This bundle names NO container module, so it relies on apptainer\n"
        "already being on a batch job's PATH. If it is not, every link aborts\n"
        "in about a second. Find the name with `module avail apptainer` on a\n"
        "login node and resubmit as:\n"
        "    APPTAINER_MODULE=<name> ./submit_chain.sh\n"
    )
    submit = "sbatch run_flexiprobes.slurm" if scheduler == "slurm" \
        else "qsub run_flexiprobes.pbs"
    cancel_cmd = "scancel" if scheduler == "slurm" else "qdel"
    script = "run_flexiprobes.slurm" if scheduler == "slurm" \
        else "run_flexiprobes.pbs"
    ext = "slurm" if scheduler == "slurm" else "pbs"
    qsub = "sbatch" if scheduler == "slurm" else "qsub"
    prepare_script = f"{PREPARE_SCRIPT_STEM}.{ext}"
    screen_script = f"{SCREEN_SCRIPT_STEM}.{ext}"
    split_note = ""
    if split:
        # In a split bundle there is no run_flexiprobes.<ext> at all. Naming a
        # script that does not exist would be worse than useless here: the
        # reader would go looking for it, and the two that DO exist are not
        # interchangeable with it or with each other.
        script = f"{prepare_script} / {screen_script}"
        submit = (f"{qsub} {prepare_script}   "
                  f"# then, after it finishes: {qsub} {screen_script}")
        pr = prepare_res or HpcResources()
        sr = screen_res or HpcResources()
        split_note = f"""
TWO JOBS, NOT ONE — read this before submitting
-----------------------------------------------
This bundle is SPLIT. The run is done by two different job scripts with two
different resource requests, chained into one submission:

  1. {prepare_script}
     ncpus={pr.cores} mem={pr.mem_gb}gb walltime={pr.walltime}, x{pr.links} link(s)
     Runs ingest + tiling + filtering with `--prepare-only`, writes
     output/.candidate_cache/, the pre-screen artefacts and
     output/{PREPARED_MARKER}, then STOPS. It designs no probes: there is no
     results.csv when it finishes, and design_summary.json says
     status "prepared". That is success, not failure.

  2. {screen_script}
     ncpus={sr.cores} mem={sr.mem_gb}gb walltime={sr.walltime}, x{sr.links} links
     Resumes from that checkpoint and runs the off-target screen, ranking,
     chemistry and export. This is the job that produces results.csv and
     flexiprobes_output.tar.gz, and the one that writes .flexiprobes_DONE.
     It REFUSES to start if output/{PREPARED_MARKER} is absent, because it is
     sized for a run that does not tile.

WHY. On a whole transcriptome, the same screening work peaks several times
higher in a process that tiled first than in one that resumed from a
checkpoint. CPython does
not return freed pages to the OS, so the screen forks its workers off a parent
still holding the tiling high-water mark. Splitting lets the memory-hungry half
ask for a lot of memory for a few hours, and the long half ask for much less
for a day at a time.

`./submit_chain.sh` queues both halves in the right order in one command. Use
it rather than submitting the two scripts by hand.
"""
    # Relative DB paths (e.g. inputs/<fasta>) live inside the bundle; absolute
    # ones must already exist on the cluster and be bound into the container.
    bundled_dbs = [p for p in db_paths if not os.path.isabs(p)]
    cluster_dbs = [p for p in db_paths if os.path.isabs(p)]
    if not db_paths:
        db_note = "  - No BLAST database configured (off-target screening skipped)."
    else:
        parts = []
        if bundled_dbs:
            parts.append(
                "  - BLAST DB is bundled (the input FASTA doubles as the off-target\n"
                "    DB): " + ", ".join(bundled_dbs) + "\n"
                "    The container builds the index on first run — nothing to stage."
            )
        if cluster_dbs:
            parts.append(
                "  - BLAST DB(s) referenced on the cluster:\n      "
                + "\n      ".join(cluster_dbs) + "\n"
                "    Make sure these exist on the cluster and add a matching\n"
                "    --bind for their directory to the exec line in the script."
            )
        db_note = "\n".join(parts)
    input_note = (
        "  - The input FASTA is bundled under inputs/ — nothing to stage."
        if staged_input else
        "  - The input FASTA is referenced by the path in config.json; ensure it\n"
        "    exists on the cluster (edit config.json's \"fasta_path\" if needed)."
    )
    self_contained = staged_input and not cluster_dbs
    if job_log_project and scheduler == "pbs":
        log_root = f"$WORKDIR/projects/{job_log_project}/{JOB_LOG_DIR}"
        log_where = f"""Each job writes its OWN log file, IN REAL TIME, to the project's log folder:

    {log_root}/<YYYYMMDD_HHMMSS>_<job_name>.<jobid>.log

(If PROJ is set in the job's environment, it replaces
$WORKDIR/projects/{job_log_project}.) PBS's own spooled copy is discarded with
#PBS -o /dev/null, because PBS only copies it back once a job has ENDED. Watch
the newest log with:

    tail -f "{log_root}/$(ls "{log_root}" | tail -1)"

Every bundle in the project shares that folder, so each log names its bundle
directory on its second line. If WORKDIR is unset inside the job, the log goes
to {JOB_LOG_DIR}/ inside this bundle instead and its first lines say so.

.flexiprobes_link.log in this directory is NOT the job log. It holds only the
current link's pipeline output, for the retry classifier, and each link
overwrites it."""
    else:
        log_where = """Each job writes its OWN log file, named per job id, in this directory:

    """ + ("flexiprobes-<jobid>.log       (Slurm, --output=flexiprobes-%j.log)"
           if scheduler == "slurm" else
           "<job_name>.o<jobid>           (PBS default, stdout+stderr merged)")
    log_note = ("""Where the logs go
-----------------
""" + log_where + """

The chain submits several links and they all run this same script, so you get
one log per link. Read them in job-id order: the link that did the real work is
the big one; the leftovers are two-liners saying the run was already complete.

They are deliberately NOT written to one shared filename — a leftover link
finishing after the run would overwrite it and destroy the stage-by-stage
progress record (which is the only place the per-stage narration lives).

Per-stage elapsed times are also recorded in output/manifest.json, under
"elapsed_seconds", so the timings survive even if the logs are lost. That
covers the design pipeline AND the per-gene selection step that follows it
(recorded as "select_coverage"), so "total" is comparable to the walltime
your scheduler reports for the job.
""")
    selection_outputs = ("""
  output/selected_probes_per_gene.csv   Top-N spaced probes per gene, chosen to
                                        maximise isoform coverage (gene-aware
                                        runs). Written automatically after a
                                        successful run.
  output/flagged_genes.csv              Genes that couldn't reach full isoform
                                        coverage or N spaced probes, with the
                                        missing isoforms listed.""" if gene_aware else "")
    return f"""flexiprobes — HPC job bundle
============================

This bundle runs the flexiprobes probe-design pipeline on an HPC cluster, then
lets you import the results back into the desktop app for review and selection.

Contents
--------
  config.json            The exact run configuration from the app.
  {script:<22} The {scheduler.upper()} batch script (resumable; safe to resubmit).
  submit_chain.sh        Submits a self-continuing chain of these jobs.
  inputs/                Staged input files (if any).
  output/               The job writes results here (created on run).
  README.txt             This file.
{split_note}
{_rebuild_block(image, expected_version)}Copy this bundle to your cluster
--------------------------------
{"This whole folder is self-contained — copy it across and submit from inside it:"
 if self_contained else
 "Copy this whole folder across, then sort out the cluster paths noted below:"}

    rsync -avP <this-folder>/ you@cluster:/path/to/run_dir/

Then `cd` into the copied folder on the cluster before submitting.

Before you submit
-----------------
1. Rebuild the container at {expected_version} (see the section above) and put it at:
       {image}
   Edit the CONTAINER=... line in {script} if it lives elsewhere.
{input_note}
{db_note}
2. Adjust the scheduler directives (walltime, memory, cores, partition/queue,
   account) at the top of {script} to match your cluster's limits.
3. Container module. A BATCH JOB DOES NOT INHERIT THE MODULES YOU LOADED
   INTERACTIVELY, so `module load apptainer` in your shell before submitting
   does nothing for the job. {module_note}

Submit (recommended: the resumable chain)
-----------------------------------------
    ./submit_chain.sh          # submit {links} linked jobs (override: ./submit_chain.sh 24)

This submits several jobs linked with `afterany` dependencies. Each link resumes
from the last checkpoint; once the run finishes, the leftover links detect the
.flexiprobes_DONE sentinel and exit instantly — so picking a generous link count
is free. This is the hands-off way to outlast any single walltime limit.

The walltime in {script} is PER LINK, not for the whole run. This bundle asks
for {walltime} x {links} links = the total resumable budget. If the run needs
longer, raise the LINK COUNT rather than the walltime: a short job is backfilled
into scheduling gaps a long one cannot use, and — the reason that matters most —
a job whose walltime would run past a scheduled maintenance window is not
started at all, it simply sits queued until the window closes.

Or submit a single job (resubmit it yourself if it hits walltime):
    {submit}

How resume works
----------------
The run checkpoints two things as it goes, so a job killed by walltime (or a dead
node) continues where it stopped on the next attempt:
  - output/.candidate_cache/  the tile+filter+cap result — a resumed run skips
                              straight to BLAST instead of re-tiling everything.
  - output/.screen_cache/      every finished BLAST chunk; specificity is scored
                              chunk-by-chunk so memory stays bounded (no end-of-
                              BLAST blow-up).
Both are written atomically and removed automatically once the run completes, so
they never end up in the results archive. Do NOT delete output/ between attempts,
or you throw the progress away. To force a clean rerun, delete BOTH output/ and
the .flexiprobes_DONE sentinel, then resubmit.

If a link keeps dying immediately, cancel the rest of the chain with `{cancel_cmd}`.

{log_note}
Outputs
-------
  output/results.csv                    Every probe candidate that passed all
                                        filters + specificity (the full list).
  output/manifest.json                  Run provenance: the complete config,
                                        input + BLAST DB checksums, engine
                                        commit, tool versions, stage timings.
  output/run_config.json                The exact config this run executed, in
                                        the form `flexiprobes design --config`
                                        takes — so the run rides home inside the
                                        tarball and can be re-executed
                                        verbatim.{selection_outputs}
{_prescreen_outputs_block(prescreen_capture)}

When it finishes
----------------
The job writes flexiprobes_output.tar.gz in this directory. Copy it back to your
machine, then in the app choose  Tools -> Import HPC run...  and select either
that archive or the unpacked output/ folder. The run then appears in the app and
opens in the review dialog like any local run.
"""


# ---------------------------------------------------------------------------
# Bundle writer
# ---------------------------------------------------------------------------

def _capture_mode(cfg: RunConfig) -> str:
    """Which pre-screen candidate set this config will preserve.

    Mirrors the same derivation in ``pipeline.run_pipeline`` — kept as one
    expression in each place rather than shared through an import, because the
    bundle README is written long before (and often on a different machine
    from) the run that produces the artefacts.
    """
    from .prescreen_export import CAPTURE_NONE, CAPTURE_POSTCAP, CAPTURE_PRECAP
    if not cfg.emit_prescreen_candidates:
        return CAPTURE_NONE
    return CAPTURE_PRECAP if cfg.emit_prescreen_precap_rows else CAPTURE_POSTCAP


def write_hpc_bundle(
    config: RunConfig,
    *,
    scheduler: str,
    resources: HpcResources,
    image: str,
    dest_dir: str,
    stage_inputs: bool = False,
    remote_input_path: Optional[str] = None,
    remote_db_paths: Optional[List[str]] = None,
    expected_version: str = __version__,
    apptainer_module: str = DEFAULT_APPTAINER_MODULE,
    split_prepare_screen: bool = False,
    prepare_resources: Optional[HpcResources] = None,
    job_log_project: str = "",
) -> str:
    """Write a self-contained HPC bundle into ``dest_dir``. Returns ``dest_dir``.

    The serialized config has its paths rewritten so the job is relocatable:
    ``output_dir`` becomes ``output`` (relative to the bundle), and the input
    FASTA becomes ``inputs/<name>`` when staged or ``remote_input_path`` when the
    file already lives on the cluster.

    BLAST databases are referenced on the cluster by default (``remote_db_paths``)
    — they are too large to stage. The one exception: when the input FASTA is
    staged AND a configured DB IS that same FASTA (the whole-transcriptome case,
    where the transcriptome is both input and off-target DB), that DB entry is
    repointed at the staged ``inputs/<name>`` copy so the bundle is
    self-contained. ``remote_db_paths``, if given, overrides everything.

    ``split_prepare_screen`` writes the two-job form instead of one:
    ``run_flexiprobes_prepare.<ext>`` (ingest + tile + filter, ``--prepare-only``)
    and ``run_flexiprobes_screen.<ext>`` (resume from the checkpoint, screen,
    rank, export), chained by ``submit_chain.sh``. ``resources`` then describes
    the SCREEN job and ``prepare_resources`` the prepare job; the module
    defaults for each are derived from measured per-stage peaks, see
    ``DEFAULT_PREPARE_MEM_GB`` and friends.

    OFF by default, deliberately. The GUI's HPC export produces the unsplit
    script, and the output is byte-identical either way, so turning the split on
    is an explicit choice a caller makes rather than something a version bump
    does to an existing workflow. A split bundle carries a prepare and a
    screen script.

    NOTE that ``config.prepare_only`` is NOT baked into config.json for the
    split form, and must not be: both jobs read the SAME config.json, and the
    prepare job adds ``--prepare-only`` on its command line. Baking it in would
    leave the screen job preparing again and never screening.

    ``job_log_project`` (PBS only) makes every job script log in real time to
    ``$WORKDIR/projects/<job_log_project>/hpc_job_logs/`` instead of PBS's
    end-of-job spooled copy. See ``JOB_LOG_DIR``. Empty (the default) leaves the
    scripts exactly as before.
    """
    if scheduler not in SCHEDULERS:
        raise ValueError(f"scheduler must be one of {SCHEDULERS}, got {scheduler!r}")
    if job_log_project:
        # Raised rather than ignored: a caller who asked for real-time logs and
        # silently got none would find out only when a running job had no log.
        if scheduler != "pbs":
            raise ValueError(
                "job_log_project is PBS-only; Slurm already writes --output "
                "in real time")
        # It lands unquoted inside a shell parameter expansion and in paths.
        if not _JOB_LOG_PROJECT_RE.match(job_log_project):
            raise ValueError(
                f"job_log_project must match {_JOB_LOG_PROJECT_RE.pattern}, "
                f"got {job_log_project!r}")

    os.makedirs(dest_dir, exist_ok=True)

    # Work on a shallow copy so we don't mutate the caller's config.
    cfg = config_from_copy(config)
    cfg.output_dir = "output"

    staged_rel: Optional[str] = None
    staged_size: Optional[int] = None
    if stage_inputs and cfg.fasta_path:
        if not os.path.exists(cfg.fasta_path):
            # Don't ship a bundle whose config points at a staged file that was
            # never created — the cluster job would fail with a confusing error.
            raise FileNotFoundError(
                f"cannot stage input FASTA into bundle: {cfg.fasta_path} not found"
            )
        inputs_dir = os.path.join(dest_dir, "inputs")
        os.makedirs(inputs_dir, exist_ok=True)
        base = os.path.basename(cfg.fasta_path)
        staged_size = os.path.getsize(cfg.fasta_path)
        shutil.copy2(cfg.fasta_path, os.path.join(inputs_dir, base))
        staged_rel = os.path.join("inputs", base)
        cfg.fasta_path = staged_rel
    elif remote_input_path:
        cfg.fasta_path = remote_input_path

    if remote_db_paths is not None:
        cfg.screen.db_paths = list(remote_db_paths)
    elif staged_rel is not None:
        # Whole-transcriptome runs use the same transcriptome FASTA as both the
        # input and the off-target BLAST DB. When we've just staged that FASTA,
        # point any DB entry that IS that same file at the staged copy so the
        # bundle is self-contained — the user copies one folder, no cluster DB
        # path to fix. Match on basename + size (cheap, no hashing) so we only
        # reuse a genuine same-file, not an unrelated DB that shares a name.
        base = os.path.basename(staged_rel)
        cfg.screen.db_paths = [
            staged_rel
            if (os.path.basename(p) == base
                and os.path.isfile(p)
                and os.path.getsize(p) == staged_size)
            else p
            for p in cfg.screen.db_paths
        ]

    # `prepare_only` is a per-JOB command-line flag, never a bundle setting:
    # both scripts read this one file. Cleared rather than trusted, so a caller
    # who set it on the config they exported cannot produce a bundle whose
    # screen job silently prepares again and never screens.
    cfg.prepare_only = False
    save_config(cfg, os.path.join(dest_dir, "config.json"))

    ext = "slurm" if scheduler == "slurm" else "pbs"
    # Gene-aware (isoform) runs get the extra top-N-per-gene selection step,
    # which needs the input FASTA path inside the bundle/cluster (the same path
    # the design step reads).
    select_fasta = cfg.fasta_path if cfg.group_by_gene and cfg.fasta_path else None

    def _write_script(name: str, res: HpcResources, phase: str) -> None:
        body = _script_body(
            scheduler, image, res.cores, select_fasta=select_fasta,
            select_min_gap=cfg.min_probe_gap_bp,
            select_relaxed_gap=cfg.relaxed_probe_gap_bp,
            select_floor_gap=cfg.floor_probe_gap_bp,
            expected_version=expected_version,
            phase=phase,
            apptainer_module=apptainer_module,
            job_log_project=job_log_project,
            job_name=res.job_name)
        header = (_pbs_header(res, realtime_log=bool(job_log_project))
                  if scheduler == "pbs" else _slurm_header(res))
        path = os.path.join(dest_dir, name)
        with open(path, "w") as fh:
            fh.write(header + "\n\n" + body)
        os.chmod(path, 0o755)

    prep_res: Optional[HpcResources] = None
    if split_prepare_screen:
        prep_res = prepare_resources or HpcResources(
            job_name=f"{resources.job_name}-prepare",
            walltime=DEFAULT_PREPARE_WALLTIME,
            mem_gb=DEFAULT_PREPARE_MEM_GB,
            cores=DEFAULT_PREPARE_CORES,
            links=DEFAULT_PREPARE_LINKS,
            partition=resources.partition,
            account=resources.account,
        )
        _write_script(f"{PREPARE_SCRIPT_STEM}.{ext}", prep_res, PHASE_PREPARE)
        _write_script(f"{SCREEN_SCRIPT_STEM}.{ext}", resources, PHASE_SCREEN)
        # No run_flexiprobes.<ext> is written at all. A leftover unsplit script
        # beside these two is a trap: submitting it would do the whole run in
        # one process at whichever resource request happened to be in it, which
        # is precisely what this bundle exists to stop.
        chain = _split_submit_chain_script(
            scheduler,
            prepare_links=prep_res.links,
            screen_links=resources.links,
            prepare_walltime=prep_res.walltime,
            screen_walltime=resources.walltime,
            prepare_mem_gb=prep_res.mem_gb,
            prepare_cores=prep_res.cores,
            screen_mem_gb=resources.mem_gb,
            screen_cores=resources.cores,
            apptainer_module=apptainer_module,
        )
    else:
        _write_script(f"run_flexiprobes.{ext}", resources, PHASE_SINGLE)
        chain = _submit_chain_script(
            scheduler, resources.links, resources.walltime, apptainer_module)

    # Chain submitter — the one-command, resumable way to run the job.
    chain_path = os.path.join(dest_dir, "submit_chain.sh")
    with open(chain_path, "w") as fh:
        fh.write(chain)
    os.chmod(chain_path, 0o755)

    with open(os.path.join(dest_dir, "README.txt"), "w") as fh:
        fh.write(_readme(scheduler, image, stage_inputs, cfg.screen.db_paths,
                         gene_aware=bool(cfg.group_by_gene),
                         expected_version=expected_version,
                         prescreen_capture=_capture_mode(cfg),
                         links=resources.links,
                         walltime=resources.walltime,
                         split=split_prepare_screen,
                         prepare_res=prep_res,
                         screen_res=resources,
                         apptainer_module=apptainer_module,
                         job_log_project=job_log_project))

    return dest_dir


def config_from_copy(config: RunConfig) -> RunConfig:
    """Deep-ish copy of a RunConfig via its JSON round-trip (avoids aliasing the
    nested dataclasses we're about to mutate)."""
    from .config_io import config_from_dict
    return config_from_dict(config_to_dict(config))


# ---------------------------------------------------------------------------
# Result importer
# ---------------------------------------------------------------------------

def _looks_like_hex_id(name: str) -> bool:
    """True for the 12-char hex names we generate for archive extraction dirs."""
    return len(name) == 12 and all(c in "0123456789abcdef" for c in name)


def _imported_root() -> str:
    base = os.environ.get("FLEXPROBE_HOME") or os.path.join(
        os.path.expanduser("~"), ".flexprobe"
    )
    root = os.path.join(base, "imported")
    os.makedirs(root, exist_ok=True)
    return root


def _find_results_dir(root: str) -> Optional[str]:
    """Return the directory under ``root`` that holds results.csv (root first,
    then a shallow search), or None."""
    if os.path.exists(os.path.join(root, RESULTS_FILE)):
        return root
    for dirpath, _dirs, files in os.walk(root):
        if RESULTS_FILE in files:
            return dirpath
    return None


_MANIFEST_FILE = "manifest.json"
_RUN_CONFIG_FILE = "run_config.json"

#: How far below a picked folder :func:`find_run_folders` looks. The usual
#: layout is ``<base>/<project>/<run>`` (two levels); one more absorbs an
#: unpacked bundle's ``<bundle>/output`` or an extra grouping folder.
RUN_SCAN_DEPTH = 3

#: Written into ``output/`` only by the bundle's select step
#: (``flexiprobes select-coverage``, see ``_select_block``). The desktop app
#: never writes it.
_CLUSTER_ONLY_FILES = ("selected_probes_per_gene.csv",)


def find_run_folders(root: str,
                     max_depth: int = RUN_SCAN_DEPTH,
                     ) -> tuple[list[str], list[str]]:
    """Finished run folders at or below ``root``: ``(runs, bare)``.

    A run folder holds BOTH ``results.csv`` and ``manifest.json``; those are
    ``runs``. A folder with a ``results.csv`` and no manifest is returned in
    ``bare`` so the caller can say why it was passed over, rather than either
    registering an unidentifiable file or dropping it without a word.

    Looks ``max_depth`` levels below ``root`` and no further. Never descends
    into a run folder (its only subfolders are checkpoint caches), a hidden
    folder, or a symlinked folder, so a link loop cannot hang the scan. Nothing
    here reads results.csv, which reaches 1 GB on a transcriptome.
    """
    runs: list[str] = []
    bare: list[str] = []

    def _visit(path: str, depth: int) -> None:
        try:
            with os.scandir(path) as it:
                entries = sorted(it, key=lambda e: e.name)
        except OSError:
            return
        files = {e.name for e in entries if e.is_file()}
        if RESULTS_FILE in files:
            if _MANIFEST_FILE in files:
                runs.append(os.path.abspath(path))
                return
            bare.append(os.path.abspath(path))
        if depth >= max_depth:
            return
        for e in entries:
            if e.name.startswith(".") or not e.is_dir(follow_symlinks=False):
                continue
            _visit(e.path, depth + 1)

    _visit(root, 0)
    return runs, bare


def cluster_evidence(output_dir: str, manifest: Optional[dict] = None) -> str:
    """Why ``output_dir`` looks like a cluster job's output, or ``""``.

    This decides the "HPC" badge for a folder brought in by "Import runs…"
    (the "Import HPC run…" dialog badges HPC unconditionally, because the user
    said where it came from). In order:

    1. ``selected_probes_per_gene.csv`` is present. Only the bundle's select
       step writes it.
    2. ``manifest.json`` records ``select_coverage`` under ``elapsed_seconds``,
       which the same step folds in, so this survives the CSV being deleted.
    3. The folder is the ``output/`` of an unpacked bundle in place: its parent
       holds ``.flexiprobes_DONE``, which only a finished job script writes.

    Anything else reads as a desktop run. ``select-coverage`` run by hand on a
    workstation also satisfies 1 and 2; the badge is a hint, and the manifest
    stays the record of where a run was made.
    """
    for name in _CLUSTER_ONLY_FILES:
        if os.path.exists(os.path.join(output_dir, name)):
            return f"{name} is present"
    elapsed = (manifest or {}).get("elapsed_seconds")
    if isinstance(elapsed, dict) and "select_coverage" in elapsed:
        return "manifest.json records a select_coverage stage"
    parent = os.path.dirname(os.path.abspath(output_dir))
    if os.path.exists(os.path.join(parent, DONE_SENTINEL)):
        return f"the enclosing bundle holds {DONE_SENTINEL}"
    return ""


def _read_json_object(path: str) -> Optional[dict]:
    """A JSON object from ``path``, or None when absent, unreadable or not a dict."""
    import json
    try:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return None
    return data if isinstance(data, dict) else None


def import_hpc_output(src: str, dest_root: Optional[str] = None) -> dict:
    """Register completed HPC output and return a dict describing the run.

    ``src`` may be a directory or an archive (.zip / .tar.gz / .tgz). Archives
    (and directories that aren't usable in place) are extracted under
    ``~/.flexprobe/imported/<id>/``. Raises ``FileNotFoundError`` if no
    ``results.csv`` is found.

    Returns: ``{output_dir, name, modality_label, status, stats}`` ready to
    populate a ``RunRecord``.
    """
    if not os.path.exists(src):
        raise FileNotFoundError(f"{src} does not exist")

    dest_root = dest_root or _imported_root()

    if os.path.isdir(src):
        results_dir = _find_results_dir(src)
        if results_dir is None:
            raise FileNotFoundError(f"no {RESULTS_FILE} found under {src}")
    else:
        target = os.path.join(dest_root, uuid.uuid4().hex[:12])
        os.makedirs(target, exist_ok=True)
        try:
            _extract_archive(src, target)
            results_dir = _find_results_dir(target)
        except Exception:
            shutil.rmtree(target, ignore_errors=True)
            raise
        if results_dir is None:
            # Don't leave an orphaned extraction dir behind on a bad archive.
            shutil.rmtree(target, ignore_errors=True)
            raise FileNotFoundError(
                f"no {RESULTS_FILE} found in archive {os.path.basename(src)}"
            )

    return _describe_run(results_dir)


def _extract_archive(src: str, target: str) -> None:
    lower = src.lower()
    if lower.endswith(".zip"):
        # ZipFile.extractall strips leading slashes and ".." components since
        # Python 3.6.2, so member paths can't escape `target`.
        with zipfile.ZipFile(src) as zf:
            zf.extractall(target)
    elif lower.endswith((".tar.gz", ".tgz", ".tar")):
        # filter="data" (py>=3.12) rejects unsafe members; fall back when absent.
        with tarfile.open(src) as tf:
            try:
                tf.extractall(target, filter="data")
            except TypeError:
                tf.extractall(target)
    else:
        raise ValueError(
            f"unsupported archive type: {os.path.basename(src)} "
            f"(expected a folder, .zip, .tar.gz, or .tgz)"
        )


def describe_run_folder(output_dir: str) -> dict:
    """Everything a run folder says about itself, for registering it as a run.

    Shared by "Import HPC run…" and "Import runs…", so the two cannot drift.
    Reads only small sidecars (manifest.json, design_summary.json, the order
    files' existence); never results.csv.

    Returns ``{output_dir, name, modality_label, status, stats}`` (the original
    contract) plus the identity and provenance keys ``manifest_readable``,
    ``manifest_timestamp``, ``fasta_sha256``, ``fasta_path``, ``app_version``,
    ``elapsed_seconds`` (the manifest's total, or None), ``has_run_config`` and
    ``cluster_evidence`` (see :func:`cluster_evidence`).
    """
    from .summary import SUMMARY_NAME

    name = os.path.basename(os.path.normpath(output_dir)) or "Imported HPC run"
    modality_label = ""
    stats: dict = {}
    identity: dict = {
        "manifest_readable": False, "manifest_timestamp": "",
        "fasta_sha256": "", "fasta_path": "", "app_version": "",
        "elapsed_seconds": None,
    }
    manifest = None
    if os.path.exists(os.path.join(output_dir, _MANIFEST_FILE)):
        manifest = _read_json_object(os.path.join(output_dir, _MANIFEST_FILE))
    if manifest is not None:
        params = manifest.get("params", {}) or {}
        from .select_export import modality_label as _modality_label
        modality_label = _modality_label(params) or modality_label
        stats["discovery_only"] = bool(params.get("discovery_only"))
        inputs = manifest.get("inputs", {}) or {}
        fasta = inputs.get("fasta", "") if isinstance(inputs, dict) else ""
        # Prefer a human-friendly name from the input FASTA over a random
        # extraction-dir hex when importing an archive.
        if fasta and _looks_like_hex_id(name):
            name = os.path.splitext(os.path.basename(fasta))[0] or name
        elapsed = manifest.get("elapsed_seconds")
        total = elapsed.get("total") if isinstance(elapsed, dict) else None
        identity.update({
            "manifest_readable": True,
            "manifest_timestamp": str(manifest.get("timestamp") or ""),
            "fasta_sha256": (str(inputs.get("fasta_sha256") or "")
                             if isinstance(inputs, dict) else ""),
            "fasta_path": str(fasta or ""),
            "app_version": str(manifest.get("app_version") or ""),
            "elapsed_seconds": (float(total)
                                if isinstance(total, (int, float)) else None),
        })

    # design_summary.json carries the run's counts (pairs_output and friends).
    # The home table reads those off ``stats``, and an import that ignored this
    # file left the Pairs column blank for every imported run. The manifest's
    # own reading of discovery_only wins, since it is what the run executed.
    summary = _read_json_object(os.path.join(output_dir, SUMMARY_NAME))
    summary_stats = (summary or {}).get("stats")
    if isinstance(summary_stats, dict):
        stats = {**summary_stats, **stats}
    # ``imported`` must be set unconditionally, because it drives the import
    # badge in the home-page table. A run with results.csv but no manifest.json (e.g.
    # only the results file was copied back) is still an imported run.
    stats["imported"] = True

    has_orders = any(
        os.path.exists(os.path.join(output_dir, f)) for f in _ORDER_FILES
    )
    status = "complete" if has_orders else "ready_for_review"

    return {
        "output_dir": os.path.abspath(output_dir),
        "name": name,
        "modality_label": modality_label,
        "status": status,
        "stats": stats,
        "has_run_config": os.path.exists(
            os.path.join(output_dir, _RUN_CONFIG_FILE)),
        "cluster_evidence": cluster_evidence(output_dir, manifest),
        **identity,
    }


#: The name this function had before it was made public; kept for callers.
_describe_run = describe_run_folder
