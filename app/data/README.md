# Preparing your input data

flexiprobes takes **FASTA** files. What goes in each file depends on the mode.
Copy a template from [`templates/`](templates/) and replace the placeholder
sequences with yours. Lines starting with `#` are ignored, so the template
notes can stay in or be deleted.

Accepted extensions: `.fa`, `.fasta`, `.fna`, `.ffn`, `.frn`.
Sequences may be upper- or lower-case; IUPAC ambiguity codes (e.g. `N`) are
allowed but reduce probe yield (an `N` in a probe half can't have a Tm and
may fail filters).

---

## 1. Standard mode — Flex v1 / Flex v2 / Visium

**One file, one record per gene.** Each sequence is the sense-strand
mRNA/cDNA (usually the CDS) you want probes against. Put as many genes in the
file as you like (a whole transcriptome is fine — use *Auto-pick top N* mode).

```
>GENE_A   optional description
ATG...sequence...
>GENE_B
ATG...sequence...
```

Template: [`templates/single_gene.template.fa`](templates/single_gene.template.fa)

---

## 2. Homology mode (consensus)

**One file per fragment.** Every record in the file is the *same* fragment
from a *different* related sequence (e.g. strain or isoform). Put the **actual
individual sequences** in — flexiprobes aligns them, builds a consensus, then
designs **conserved** probes (bind every sequence) or **sequence-unique** probes
(distinguish one sequence).

Wherever the sequences disagree the consensus carries an ambiguous base; in a
designed probe that position **defaults to the first record** in the FASTA (the
reference/default sequence), so the exported Cell Ranger / IDT order files
contain only concrete A/C/G/T — never an `N`.

For a 3-segment virus with 2 sequences, make **3 files** (one per segment),
each with 2 records, and drag all three in together (they become 3 runs):

```
OROV_L.fa     >seqA ...   >seqB ...
OROV_M.fa     >seqA ...   >seqB ...
OROV_S.fa     >seqA ...   >seqB ...
```

Records should be roughly the same length. Records of different lengths are aligned first to handle indels: with MAFFT
(bundled in the desktop app), or MUSCLE if no MAFFT is found, falling back to
the built-in progressive aligner. You can also pre-align and supply
equal-length records.

Template: [`templates/multi_strain_segment.template.fa`](templates/multi_strain_segment.template.fa)

### Isoform-collapse (constitutive) mode

Same input shape as multi-strain — **one file, one record per isoform** of a
gene — but choose the **Constitutive (isoform-collapse)** strategy. flexiprobes
designs probes from regions shared by a *quorum* of isoforms ("Min isoforms
shared", `--isoform-min`; 0 = all), so a single probe detects every covered
isoform. Use it to make pan-isoform probes against constitutive (always-included)
exons. The strict "all isoforms" case is identical to *conserved* mode.

Aligned/gapped FASTAs are accepted (gap columns from skipped exons are handled).
The built-in aligner is progressive pairwise — best for a few closely related
isoforms; for many isoforms or large exon insertions/deletions, pre-align with
MAFFT/MUSCLE and supply equal-length, already-aligned records.

---

## Optional: off-target screening

> The default engine is an **exhaustive k-mer search** that reads reference
> FASTAs directly and needs no extra software. The `--screen-engine blast`
> alternative uses BLAST+, which the desktop app bundles (`blastn`,
> `makeblastdb`, `blastdbcmd`). `FLEXPROBE_BLAST_DIR` overrides which BLAST+ is
> used; running from source uses a system install.
>
> **Privacy:** screening is fully offline. NCBI BLAST+ reports usage stats to
> NCBI by default, and the app runs it with that **disabled**. Your references
> and FASTAs are read locally and nothing is uploaded.

The screen is an **off-target screen**: the references you supply are
sequences a probe must **NOT** hit. Each 25-bp half is compared against every
reference, with mismatches counted across the full 25 nt. A probe is rejected
when **both** halves sit within **4 mismatches** of the same off-target
(CG000621/CG000839: at least five mismatches in at least one half). An on-target hit is *not*
required — the probe is tiled from your target, so on-target binding is given.

> ⚠ **Do not put your design target in the off-target references** — it would screen out
> your own probes. Screen against the host transcriptome (e.g. GENCODE) and any
> other organisms you want to exclude.

To screen, point the run at one or more reference FASTAs (flexiprobes indexes
each FASTA the first time and caches the index):

| Role | File | Purpose |
|------|------|---------|
| Host off-target | host transcriptome FASTA (e.g. GENCODE) | reject probes that bind host RNA |
| Cross-reactivity | all *other* (non-target) organisms in one FASTA | reject probes shared with non-targets |

To confirm coverage/cross-reactivity *within* a viral panel (a probe should hit
its own sequences but not others), use **homology mode** (records = the related
sequences) — those sequences are then recognised — or list a target's other
isoforms under **On-target isoforms** (`--extra-on-target`) to excuse them from
the screen.

Add the FASTAs under **Off-target references** in the New-query dialog, or
with `--off-target-ref` on the CLI (repeatable).

For the BLAST engine you can instead build a reusable database once and pass
its prefix:

```
flexiprobes build-db --input human_tx.fa --db-prefix ~/blast_dbs/human_tx
flexiprobes design --fasta tx.fa --out out --screen-engine blast \
    --off-target-ref ~/blast_dbs/human_tx
```

> If no off-target reference is given, probes are still produced but every row in
> `results.csv` is tagged `specificity_status = not_screened` so you know they
> were **not** off-target checked.

**Isoforms vs off-targets.** If your off-target reference contains other transcript isoforms of the target gene, a valid probe will hit them and look non-specific.
List those accessions as **On-target isoforms** in the New-query dialog (or
`--extra-on-target ACCESSION`, repeatable) so they're excused from off-target
rejection without being required.
