# flexiprobes

Probe design for **10x Genomics Flex v1 / Flex v2 / Visium spatial**. Takes a
target FASTA, tiles it, filters and off-target-screens candidates, ranks them,
and exports orderable IDT oPool and Cell Ranger probe-set CSVs. The screen is
an exhaustive k-mer search by default; `--screen-engine blast` selects BLAST+.

The probe chemistry (handles, Probe Barcodes, pCS1, poly-A) comes from the 10x
Genomics tech notes CG000621 Rev D and CG000839 Rev A, and the test suite
reproduces the EGFP example in CG000839 Tables 4 & 5 byte-for-byte. Some
screening thresholds are FlexiProbes' own heuristics rather than 10x
specifications (for example the 9 bp junction-duplex gate and the 56 °C Tm
floor). `src/flexprobe/constants.py` marks which is which.

---

## Install

Python **3.11 or newer**. Run these from this directory (`app/`):

```bash
python3.11 -m venv .venv && source .venv/bin/activate
pip install -e ".[gui,dev]" -c constraints.txt   # core + PySide6 GUI + pytest
```

`-c constraints.txt` pins the exact dependency versions this release was built
and tested against. `pyproject.toml` declares looser ranges so CI catches
dependency drift across the Python matrix; use the constraints file whenever you
want to reproduce a result rather than test against the latest.

To install the desktop binary instead of running from source, see the
[repository README](../README.md#install-the-desktop-app) — it comes from the
GitHub Release, not from this tree.

BLAST+ is only needed for `--screen-engine blast`; the default k-mer engine
needs nothing. The **packaged desktop app bundles BLAST+**, so end users install
nothing. When running **from source**, install it (or point
`$FLEXPROBE_BLAST_DIR` at a folder containing `blastn`):

```bash
brew install blast              # or: apt-get install ncbi-blast+
```

## Run

```bash
flexiprobes-gui                  # desktop GUI
# or
python -m flexprobe_gui
```

CLI:

```bash
# Basic design (filter-only smoke test)
flexiprobes design --fasta my_genes.fa --out out --modality flex_v2 --skip-screen --force

# Whole-transcriptome: k-mer off-target screen + parallelise across processes
flexiprobes design --fasta tx.fa --out out --off-target-ref human_tx.fa --workers 8

# Isoform-collapse: probes shared by a quorum of isoforms (0 = all)
flexiprobes design --fasta gene_isoforms.fa --out out \
    --multi-strain --strategy constitutive --isoform-min 0

# Don't reject a target's own isoforms as off-target
flexiprobes design --fasta gene.fa --out out --off-target-ref human_tx.fa \
    --extra-on-target NM_000546.6 --extra-on-target NM_001126112.3

# Whole-transcriptome, auto-pick top N per gene with no manual review
flexiprobes design --fasta tx.fa --out out \
    --off-target-ref human_tx.fa --auto-select-per-gene 5

# Whole-transcriptome speedup: only screen the best 100 candidates per gene
flexiprobes design --fasta tx.fa --out out \
    --off-target-ref human_tx.fa --screen-top-per-gene 100

# The BLAST+ engine instead (needed for genome-scale references)
flexiprobes build-db --input human_tx.fa --db-prefix ~/blast_dbs/human_tx
flexiprobes design --fasta tx.fa --out out --screen-engine blast \
    --off-target-ref ~/blast_dbs/human_tx

# Run from a config exported by the GUI's "Export for HPC…" (paths relocatable)
flexiprobes design --config config.json --workers 16 --force

flexiprobes validate-blast --db ~/blast_dbs/human_tx --out /tmp/check
flexiprobes version
```

Run `flexiprobes design -h` for the full flag list (GC bounds, N/homopolymer,
TN-junction, multiplex/barcodes, primer3, multi-strain strategies, etc.).

To prepare input, copy a template from [`data/templates/`](data/templates/)
and see [`data/README.md`](data/README.md) (also in the GUI under
**Tools → Input data guide**). Then drag your FASTA onto the window.

In the New-query dialog, the off-target screen (k-mer engine) and primer3
Tm/hairpin scoring are **on by default**; the chemistry-derived filter/tile values and the
composite-score ranking weights live in a foldable **Advanced options** section
(collapsed by default) since they rarely need changing.

---

## What it does

```
ingest → tile → filter → screen → specificity → rank → chemistry → export
```

- **Tile**: sliding 50-bp windows (25 bp LHS half + 25 bp RHS half). Aligned /
  consensus FASTAs with `-` gaps are accepted (gap-containing tiles are skipped,
  coordinates preserved). Multi-gene runs can fan out across worker processes
  (`--workers N`).
- **Filter funnel**: per-half Tm floor (56 °C default — the duplex-stability
  gate), ≤1 N per half, homopolymer runs, optional TN ligation-junction (A at
  mRNA position 25). GC per half is computed and ranked on but does not gate by default
  (`--enforce-gc-range` makes the 44–72% band a hard gate). A
  `filter_report.json` records why tiles were rejected (viewable via the
  **Filters** button).
- **Off-target screen** (optional, multi-DB). The default engine is an
  **exhaustive k-mer search**: every probe half is compared to every candidate
  subject position across all 25 positions, so the mismatch count is the number
  of mismatches rather than a local-alignment proxy. `--screen-engine blast`
  selects the previous path — real `blastn-short`, chunked with progress — and
  is required for a reference larger than the k-mer index can address
  (~2.1 Gbases, i.e. a genome). **BLAST+ (`blastn` + `makeblastdb`) ships inside
  the desktop app** — no install needed; `$FLEXPROBE_BLAST_DIR` overrides for a
  custom build.
  Pass DB prefixes *or* raw FASTAs (auto-indexed, SHA-256 cached). On/off-target
  identity & length floors are configurable. If no DB is configured, probes are
  passed forward but tagged `specificity_status = not_screened` so a "pass" is
  never confused with "never screened". Accessions listed as *extra on-target*
  (e.g. the target gene's other isoforms) are excused from off-target rejection.
- **Rank**: a transparent, **non-calibrated preference score** — a weighted sum
  of GC proximity (to the ≈58% centre of the 44–72% band, *not* 50%), GC
  balance, specificity margin (distance from a passing probe's closest near-miss
  off-target to the disqualification threshold), homopolymer penalty and central
  placement, plus an optional isoform-coverage term in gene-aware mode. Weights
  are overridable per run and the six sub-scores are written to `results.csv` so
  any rank is auditable. It is a *preference* ordering, not a probability of
  probe success. primer3 Tm, hairpin, self-/hetero-dimer
  and per-half Tm-balance adjustments are on by default (`--no-primer3` turns
  them off) (per-candidate, failure-isolated); reported
  Tm uses primer3's default salt conditions and is a **relative** indicator, not
  the in-assay Flex buffer.
- **Chemistry**: assembles LHS/RHS oligos per modality (handles, Probe
  Barcodes, pCS1, poly-A).
- **Export**: `results.csv` (every candidate), then order CSVs.

### Modes

| Mode | What it gives you |
|------|-------------------|
| **Manual review** (default, tile-based) | Discovery writes `results.csv`; you pick final probes in the review dialog (visual map, spacing check, per-sequence columns), then order CSVs are written from your selection. The map colours each probe by its closest neighbour — green ≥150 bp, orange 40–150, red <40 — and a **Spacing** box can forbid the tighter tiers so they are never selectable. |
| **Auto-pick top N** | Whole-transcriptome convenience: greedily picks the top-N probes per gene, spread as widely as the transcript allows, walking a 150 → 40 → 15 bp ladder and stopping at the first rung that finds 3 probes. Reports probes per tier and exports without manual review. |
| **Maximise isoform coverage** | Gene-aware runs: picks up to N spaced probes per gene that together bind as many of the gene's isoforms as possible, flagging genes that can't be fully covered. Uses the persisted `isoform_hit_ids` column, so no input FASTA is needed at selection time. |
| **Homology mode (consensus)** | One FASTA = N related sequences of one fragment. Builds a consensus (IUPAC at variants; Biopython pairwise for indels) and designs either **conserved** probes (bind every sequence) or **sequence-unique** probes (distinguish one sequence by ≥K mismatches). Any ambiguous base in a designed probe defaults to the first input sequence, so order files contain no N's. `results.csv` gains per-sequence mismatch + ΔTm columns. |
| **Isoform-collapse (constitutive)** | One FASTA = N isoforms of a gene. Designs probes from regions shared by a **quorum** of isoforms (`--isoform-min`, 0 = all), so one probe binds every covered isoform. The strict "all isoforms" case equals conserved mode. For many isoforms / large exon indels, pre-align (MAFFT/MUSCLE) and pass equal-length records. |
| **Gene-aware / whole-transcriptome** (`--group-by-gene`) | One FASTA = a whole transcriptome (all isoforms of many genes). Groups transcripts by their Ensembl `gene:` header field so sibling isoforms aren't off-targets, and **scores each probe on isoform coverage** (fraction of its gene's isoforms it binds). Soft scoring rewards pan-isoform probes; `--min-isoform-coverage` optionally hard-drops low-coverage ones. Unlike isoform-collapse, no alignment/per-gene prep — it scales to all genes at once. Use the same transcriptome as input + BLAST DB. **Speed:** off-target screening every candidate dominates runtime — set **Screen top N / gene** (`--screen-top-per-gene`, e.g. 100) to screen only the best N windows per gene by the screen-free pre-score and cut screening volume ~10–50× with no quality tradeoff. |
| **Flex v1 multiplex** | Choose how many samples you are multiplexing; the run uses Probe Barcodes BC001..BC00N. Every selected probe is repeated for each barcode and each pool is named after its barcode (`BC001`, `BC008`, …). |
| **HPC export / import** | For whole-transcriptome runs too memory-heavy for a laptop. **New query → Export for HPC…** writes a PBS/SLURM job bundle (config + script + staged input) that runs the pipeline on a cluster inside an Apptainer/Singularity container. When the input transcriptome doubles as the off-target DB, the bundle is **self-contained** — `rsync` the whole folder over and submit, no paths to edit. Bring the results back with **Tools → Import HPC run…** (folder or `.tar.gz`/`.zip`) — imported runs are badged **HPC** and open in the review dialog like local ones. Build the container once from [`packaging/container/`](packaging/container/README.md). |

### Projects & bulk

- Group runs into **projects**; drag a folder or many FASTAs to queue one run
  each (sequential), all sharing the dialog's settings.
- **Bulk export** combines every confirmed selection in a project into one set
  of order CSVs with run-prefixed pool names.
- **Bulk select** on the home page: tick the checkbox on any completed runs (or
  **Select all** to tick the whole page) and then **Move to project…** or
  **Delete selected** in one action. Deleting can optionally remove the runs'
  output folders from disk too.
- **Manage projects…** supports the same multi-select: tick projects (or
  **Select all**) and **Delete selected** to remove several at once. Rename and
  edit-description act on a single ticked project. Deleting a project never
  touches its runs — they stay on disk and just lose their project filter.

Persistent state lives in `~/.flexprobe/` (`runs.json`, `projects.json`) and in
each run's output dir (`results.csv`, `manifest.json`, `selections.json`,
`filter_report.json`, `consensus.fa`).

---

## Output columns

`results.csv` has one row per ranked tile: the original mRNA window, its RC,
the two 25-bp halves, both full oligos, GC for the full 50-mer **and** each half
(`gc_content_full`, `gc_lhs_half`, `gc_rhs_half`; the review table shows all
three as percentages), Tm (LHS/RHS/full; `n/a` when a sequence has ambiguous
bases), homopolymer, TN-junction, specificity pass + **status**, isoform
coverage (`isoform_coverage` / `isoform_hits` / `isoform_total`; `1/1/1.0`
when gene grouping is off or the input has no `gene:` header fields), composite score, and — for homology-mode /
isoform runs — per-sequence mismatch / ΔTm columns. The per-strain `est_dTm`
(and the `likely_binds` flag) are **comparative** rank aids from a simple
per-mismatch heuristic (≈−2 °C/mismatch, extra near the junction), **not**
absolute nearest-neighbour thermodynamics — use them to compare strains within a
probe, not as predicted melting temperatures.

Additional columns:

- **Off-target evidence** (per passing probe): `offtarget_nearmiss_count`, and
  the closest near-miss off-target's `best_offtarget_pident` /
  `best_offtarget_bitscore` / `best_offtarget_evalue`. The pass/fail decision
  itself stays on the 10x-sourced identity+length rule; these report the
  evidence behind a "passed" so specificity is interpretable (the manifest also
  records each BLAST DB's size for the multiple-testing context).
- **Score sub-components**: `score_gc_proximity`, `score_gc_balance`,
  `score_specificity_margin`, `score_homopolymer_penalty`,
  `score_central_placement`, `score_isoform_coverage` — the six 0–1 sub-scores
  whose weighted, normalised sum is `composite_score` (a **non-calibrated**
  preference score).
- **Thermodynamics** (0 with `--no-primer3`): `tm_balance` (|Tm(LHS) −
  Tm(RHS)|), `self_dimer_dg`, `hetero_dimer_dg` — relative indicators.

### Specificity gate reference

Numeric defaults live in `src/flexprobe/constants.py`, which cites each to its
source; they are surfaced on `FilterThresholds` (`src/flexprobe/models.py`), the
CLI, and the New query dialog. Each gate has an opt-out or a control.

| Gate | Default | Opt-out / control |
|---|---|---|
| TN ligation junction (mRNA index 25 == `A`) | **on**, `strict_a` | `--no-require-tn-junction`; `--tn-junction-mode {auto,strict_a,soft,off}` |
| Per-half Tm floor | **56.0 °C** | `--tm-min` (negative disables) |
| Per-half Tm ceiling | none | `--tm-max` |
| GC band per half | **ranking only** (44–72%) | `--enforce-gc-range` makes it a hard gate; `--gc-min` / `--gc-max` |
| Off-target: halves that must fail | **both** | `--off-target-either-half` |
| Off-target: mismatch counting | **full 25 nt**, ≤ 4 disqualifies | `--off-target-local-alignment` (local alignment: ≥ 80% id over ≥ 20 bp) |
| Off-target search engine | **k-mer (exhaustive, end-to-end)** | `--screen-engine blast` (required for genome-scale references) |
| Off-target: junction-spanning duplex | **on**, 9 bp | `--off-target-junction-duplex N`; `--no-off-target-junction-duplex` |
| Off-target: subject scoping | **per off-target gene** | `--no-off-target-same-subject` (score over the whole probe) |
| Reference composition QC | **ALT/PATCH excluded** | `--no-exclude-alt-scaffolds` |
| N bases per half | 1 | `--max-n-per-half` |
| Homopolymer run | reject > 4 | `--max-homopolymer-run` |

The Tm floor is a **relative** figure — primer3 `calc_tm` under the PCR-primer
salt conditions in `constants.py`, because 10x does not publish the Flex
hybridisation buffer. It is calibrated in the same units it is compared
against, so it is a valid screen, but it must not be read as "these probes melt
at 56 °C in the Flex assay". If you know your buffer, adjust `TM_MV_CONC` …
`TM_DNA_CONC` and re-derive the floor.

### Defaults worth knowing

- **TN junction** is **on** by default. It is the **hard** adenine gate: a tile is
  rejected unless the mRNA base opposite the LHS probe's 3′ end is `A`.
  Disable it with `--no-require-tn-junction` — defensible, since only ~51% of
  the real 10x v1 panel carries the junction adenine, but state it in your
  methods. `--tn-junction-mode soft` relaxes it to "any base but C" for
  Flex v1/v2 only — 10x's panels avoid C at that position ~99% of the time —
  and is refused for Visium, whose tech note states the adenine flatly.
  Passing `--require-tn-junction` together with `--tn-junction-mode soft` is
  rejected as a contradiction. `junction_base` and `tn_junction_ok` are written
  to `results.csv` for every probe whatever gate ran.
- **GC per half** does not reject by default (`--enforce-gc-range` applies the
  hard 44–72% band); the **per-half Tm floor** (`--tm-min`, default 56.0 °C) is the
  duplex-stability gate instead. 8.36% of 10x's own human-panel halves fall
  outside their published GC band while none falls below 56.32 °C. The Tm is
  **relative** — primer3 under PCR-primer salt, because 10x does not publish
  the Flex hybridisation buffer — so it ranks and gates consistently but is not
  the in-assay melting temperature.
- **Off-target screening** rejects a pair only when *both* halves are within 4
  mismatches of an off-target, counted across the full 25-mer, per
  CG000621/CG000839 ("at least five mismatches in at least one of the LHS or
  RHS probes"). `--off-target-either-half` and `--off-target-local-alignment` apply a
  stricter interpretation.
- The **homopolymer** default rejects runs longer than 4 — *stricter* than 10x, which gives no explicit length.
- **BLAST DB header convention**: on-/off-target classification matches the
  subject id against the input record's id using accession/`|`-delimited
  boundaries (not `_`), so a DB whose subjects use `_`-delimited compound ids
  (e.g. `NC_045512`) can misclassify. For transcriptome runs, label the DB by
  gene id (`ENSG|ENST`) so sibling isoforms are on-target and other genes
  off-target.

---

## Tests

```bash
pytest            # the default suite; headless, no BLAST required.
                  # GUI tests auto-skip if PySide6 is absent; primer3
                  # Tm/dimer tests skip if primer3 is absent.
pytest -m slow    # the slow exhaustive ground-truth oracle tests
```

`addopts = -m 'not slow'` in `pyproject.toml` deselects the `slow` tests
from a plain `pytest` run. They are the only check that compares the specificity
screen against ground truth rather than against itself. They take seconds, and CI
runs them as a separate step, so run them too if you touch the specificity
path.

CI (the repository root's `.github/workflows/ci.yml`) runs pytest + `pytest -m
slow` + mypy + pyflakes across the Python matrix, and builds the macOS
`.app`, the Windows `.exe` and the Linux binary (`build-macos`,
`build-windows`, `build-linux`). On a `v*` tag push the binaries are uploaded, and the `release` job
publishes them, with a `SHA256SUMS` file, as a GitHub Release.

## Package a desktop binary

PyInstaller is not part of the `dev` extra, so install it into the same
environment first:

```bash
pip install pyinstaller
python packaging/build.py --clean              # PyInstaller → dist/flexiprobes(.app)
python packaging/build.py --clean --sign -     # + ad-hoc codesign (macOS)
```

PyInstaller is host-specific and does not cross-compile: build the macOS `.app`
on macOS, the Windows `flexiprobes.exe` (onedir) on Windows, and the Linux
binary on Linux. The artifact always lands in the **repository-root** `dist/`,
never `app/dist/`, whichever directory you launch the build from.

Releases include macOS, Windows and Linux binaries. On Linux the CLI is the
usual way to run FlexiProbes, from the binary or from source.

The build bundles the BLAST+ and MAFFT binaries under `packaging/vendor/`
for the build host's platform (see `packaging/vendor/blast/README.md` and
`packaging/vendor/mafft/README.md`). If a platform's folder is empty, the app
still builds, but without that tool.

### Code signing and notarisation (macOS)

To distribute a macOS build beyond your own machine, sign it with a Developer
ID and notarise it:

```bash
python packaging/build.py --clean --sign "Developer ID Application: NAME (TEAMID)"
ditto -c -k --keepParent dist/flexiprobes.app flexiprobes.zip
xcrun notarytool submit flexiprobes.zip --apple-id you@example.com \
      --team-id TEAMID --password APP_SPECIFIC_PW --wait
xcrun stapler staple dist/flexiprobes.app
```

`packaging/entitlements.plist` grants the hardened-runtime exceptions a
PyInstaller app needs.

## Layout

```
app/
  pyproject.toml
  LICENSE  NOTICE  THIRD-PARTY-LICENSES/   # GPL-3.0 + third-party breakdown
  src/
    flexprobe/          # core library + CLI
      blast.py          # BLAST discovery (bundled-first) + run + DB build
      chemistry.py      # per-modality oligo assembly (10x-sourced constants)
      consensus.py      # homology-mode consensus + binding scoring
      pipeline.py       # orchestrator
      select_export.py  # selection → order CSVs (+ bulk + multi-barcode)
      writers/
    flexprobe_gui/      # PySide6 GUI (home, new-query, review, filter-report)
  tests/
  data/                 # input-format templates + guide (no sample sequences)
    packaging/            # PyInstaller spec + build.py
    container/          # Docker / Apptainer recipes for HPC runs
    vendor/blast/       # bundled NCBI BLAST+ binaries (per-platform; public domain)
    vendor/mafft/       # bundled MAFFT portable packages (per-platform; BSD)
```

## License

**GPL-3.0-or-later** (copyleft inherited from primer3-py; PySide6/Qt used under
its GPL-3.0 option). Bundled NCBI BLAST+ binaries are public domain, and
bundled MAFFT is BSD-licensed. See
[`LICENSE`](LICENSE), [`NOTICE`](NOTICE), and
[`THIRD-PARTY-LICENSES/`](THIRD-PARTY-LICENSES). Research software — **not a
medical device**; no patent rights granted.
