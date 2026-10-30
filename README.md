<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="app/packaging/assets/flexiprobes-wordmark-dark.svg">
    <img src="app/packaging/assets/flexiprobes-wordmark-light.svg" alt="FlexiProbes" width="420">
  </picture>
</p>

FlexiProbes is a custom probe design tool for
**10x Genomics Flex v1 / Flex v2 / Visium spatial**. FlexiProbes takes a target
FASTA, tiles it, filters and off-target-screens the candidates, ranks them, and
exports orderable IDT oPool and Cell Ranger probe-set CSVs.

> [!NOTE]
> The FlexiProbes manuscript will be posted on bioRxiv soon. Until then,
> please cite the software itself, using [`CITATION.cff`](CITATION.cff) or
> **Cite this repository** on GitHub.

Requires Python **3.11 or newer** to run from source. The desktop binaries
bundle their own Python and need nothing installed. 

## Install the desktop app

Download the archive for your platform from the
[**Releases** page](https://github.com/Adrian-I-lab/FlexiProbes/releases),
check it against the `SHA256SUMS` file published with each release, and open it.

The app opens the GUI. Drag a target FASTA onto the window to design probes. See
`app/README.md` for the full feature tour (manual review, auto-pick,
homology-mode consensus).

| Platform | Asset | What is inside |
|---|---|---|
| macOS (Apple Silicon) | `flexiprobes-macos-arm64.zip` | `flexiprobes.app` |
| Windows (x64) | `flexiprobes-windows-x64.zip` | `flexiprobes/flexiprobes.exe` (onedir) |
| Linux (x64) | `flexiprobes-linux-x64.tar.gz` | `flexiprobes/flexiprobes` (onedir) |

Download `SHA256SUMS` into the same folder as the archive, then:

```bash
shasum -a 256 -c SHA256SUMS --ignore-missing          # macOS / Linux
certutil -hashfile flexiprobes-windows-x64.zip SHA256  # Windows: compare with SHA256SUMS
```

On Linux the CLI (`flexiprobes design …`) is the usual way to use FlexiProbes,
and both the CLI and the GUI also run from source (see **Run from source**
below).

**Nothing to install.** Each build bundles its own Python, Qt, **BLAST+**
(`blastn`, `makeblastdb`, `blastdbcmd`) and **MAFFT**, so the
`--screen-engine blast` path and multi-sequence alignment work out of the box
alongside the default k-mer engine.

The macOS bundle is **ad-hoc signed, not notarized**. On first launch Gatekeeper
may refuse it: right-click → Open, or
`xattr -dr com.apple.quarantine flexiprobes.app`.

## Run from source (CLI / GUI)

Python **3.11+**. From the repository root:

```bash
cd app
python3.11 -m venv .venv && source .venv/bin/activate
pip install -e ".[gui,dev]" -c constraints.txt   # -c pins the tested versions
flexiprobes-gui                                   # desktop GUI
flexiprobes design --fasta my_genes.fa --out out  # CLI
pytest                                            # headless test suite
```

`constraints.txt` is the reproducible-install lockfile: it holds the exact
dependency versions a release was built and tested against. `pyproject.toml`
declares looser ranges so CI catches dependency drift. Use `-c constraints.txt`
whenever you want to reproduce a result rather than test against the latest.

BLAST+ is bundled only in the packaged app. From source, `--screen-engine blast`
needs it on `PATH` (`brew install blast`, `apt-get install ncbi-blast+`) or
`$FLEXPROBE_BLAST_DIR` pointed at a folder containing `blastn`. The default
k-mer engine needs nothing.

## Run on an HPC cluster

Whole-transcriptome designs are too memory-heavy for a laptop. In the desktop
app, configure the run then choose **New query → Export for HPC…** to write a
self-contained PBS/SLURM job bundle that runs the pipeline on a cluster inside an
Apptainer/Singularity container. When the job finishes, bring the results back
with **Tools → Import HPC run…**. How to build the container here:
[`app/packaging/container/README.md`](app/packaging/container/README.md).

## Layout

```
FlexiProbes/
├── app/                     # application source
│   ├── src/flexprobe/       #   core engine + CLI  (entry: flexprobe.cli:main)
│   ├── src/flexprobe_gui/   #   PySide6 desktop GUI (entry: flexprobe_gui.app:main)
│   ├── tests/               #   headless pytest suite
│   ├── packaging/           #   PyInstaller spec, build.py, container recipes
│   └── README.md            #   full feature and option reference
└── .github/workflows/ci.yml # tests, type-check, lint, desktop builds
```

## Citing FlexiProbes

See [`CITATION.cff`](CITATION.cff), or use **Cite this repository** on GitHub.

The FlexiProbes app was developed by Adrian Ilich and collaborators at the Systems 
Virology Lab, QIMR Berghofer.

## Licence

**GPL-3.0-or-later**. See [`LICENSE`](LICENSE), [`app/NOTICE`](app/NOTICE) and
[`app/THIRD-PARTY-LICENSES/`](app/THIRD-PARTY-LICENSES). Research software, not a
medical device.
