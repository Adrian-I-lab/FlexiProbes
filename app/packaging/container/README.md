# flexiprobes HPC container

The whole-transcriptome method (gene-aware tiling plus an off-target screen
against a full reference transcriptome) is too memory-heavy for a laptop. The desktop app can
instead **export an HPC job bundle** (New query → *Export for HPC…*) that runs
the headless `flexiprobes` CLI on a cluster inside this container, then lets you
**import the results back** (Tools → *Import HPC run…*).

This directory holds the container definition. Build it for your cluster's
architecture, and **rebuild it for every FlexiProbes version**: each job
bundle checks the engine version inside the image and refuses to start
(exit 3) on a mismatch. The resulting image path is what you enter in the
export form.

Run every build command below from the repository's `app/` directory (the
`%files` and `COPY` paths are relative to it). Start by reading the version
from `pyproject.toml`, so the image is named after the release you are
building:

```bash
cd app
VER=$(grep -m1 '^version' pyproject.toml | cut -d'"' -f2)
```

## What's inside

- `flexiprobe.def` — Apptainer/Singularity definition (recommended for HPC).
- `Dockerfile` — OCI build, for sites that prefer Docker → Apptainer conversion.

Both install the `flexiprobes` CLI plus the same vendored **NCBI BLAST+ 2.17.0**
(`blastn`, `makeblastdb`, `blastdbcmd`) and **MAFFT 7.526** the desktop app
ships, copied from `packaging/vendor/`. A cluster run therefore uses the same
tool versions as a desktop run. The build fails if the BLAST binaries' SHA-256 hashes or either tool's
version do not match.

## Build with Apptainer (recommended)

```bash
apptainer build packaging/container/flexiprobe.sif packaging/container/flexiprobe.def
```

Copy `flexiprobe.sif` to the cluster (or build it there) and point the export
form's *Container image* field at its path.

## Build with Docker, then convert

```bash
docker build -f packaging/container/Dockerfile -t flexiprobe:$VER .
apptainer build flexiprobe.sif docker-daemon://flexiprobe:$VER
```

## Build locally for HPC (from macOS / a machine without Apptainer)

Apptainer is **Linux-only** — it can't run natively on macOS. So on a Mac you
build the *Docker* image locally and convert it to a `.sif` on the cluster (or in
a Linux VM). Two gotchas:

- **Target the cluster's CPU architecture.** Apple Silicon builds `arm64` by
  default, but most HPCs are `x86_64`. Pass `--platform linux/amd64` or the image
  won't run on the cluster. (The build then runs under emulation — slower, but
  fine.)
- **The build needs internet** (it pulls the `python:3.11-slim` base image, Python
  packages, and `libgomp1`), but the resulting tarball is self-contained, so the HPC does
  **not** need internet.

**On your Mac** (needs Docker Desktop), from the `app/` directory, export to a tarball so no registry is required:

```bash
docker build --platform linux/amd64 \
    -f packaging/container/Dockerfile -t flexiprobe:$VER .
docker save flexiprobe:$VER -o flexiprobe-$VER.tar
```

**Copy to the HPC and convert** to a `.sif`:

```bash
scp flexiprobe-$VER.tar you@hpc:/scratch/you/

# on the HPC:
module load apptainer/1.4.1            # any 1.x is fine
export APPTAINER_CACHEDIR=$SCRATCH/.apptainer/cache   # keep layers off your home quota
export APPTAINER_TMPDIR=$SCRATCH/.apptainer/tmp
apptainer build flexiprobe.sif docker-archive://flexiprobe-<version>.tar   # VER is not set on the HPC
# if it complains about privileges, add --fakeroot
```

Then point the export form's *Container image* field at
`/scratch/you/flexiprobe.sif`.

*Alternative (HPC has registry access):* `docker push` the image to a registry
and run `apptainer build flexiprobe.sif docker://<registry>/flexiprobe:<version>` on
the cluster instead of the tarball steps.

*Alternative (produce the `.sif` entirely on the Mac):* run `apptainer build`
inside a Linux VM such as [Lima](https://lima-vm.org/) (`brew install lima`),
then `scp` the finished `flexiprobe.sif` over.

## Verify

```bash
apptainer exec flexiprobe.sif flexiprobes version
apptainer exec flexiprobe.sif blastn -version
```

## How the job uses it

The exported `run_flexiprobes.slurm` / `.pbs` script does, in effect:

```bash
apptainer exec --bind "$PWD" flexiprobe.sif \
    flexiprobes design --config config.json --workers <cores> --force
```

`--bind "$PWD"` makes the bundle (`config.json`, `inputs/`, `output/`) visible
inside the container. If your input FASTA or BLAST database live elsewhere on the
cluster filesystem, add more `--bind` paths in the script.
