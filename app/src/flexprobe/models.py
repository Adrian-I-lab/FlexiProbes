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

"""Core data models for flexprobe."""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
import datetime

from .constants import (
    GC_MIN_DEFAULT, GC_MAX_DEFAULT, MAX_HOMOPOLYMER_RUN_DEFAULT,
    TILE_LENGTH_DEFAULT, GROUP_BY_GENE_DEFAULT,
    TM_MIN_DEFAULT, TM_MAX_DEFAULT,
    OFF_TARGET_PIDENT_MIN, OFF_TARGET_ALN_MIN,
    ON_TARGET_PIDENT_MIN, ON_TARGET_ALN_MIN,
    OFF_TARGET_MIN_MISMATCH, OFF_TARGET_FULL_HALF_DEFAULT,
    OFF_TARGET_REQUIRE_BOTH_HALVES_DEFAULT, OFF_TARGET_SAME_SUBJECT_DEFAULT,
    OFF_TARGET_JUNCTION_DUPLEX_BP, OFF_TARGET_JUNCTION_DUPLEX_MIN_BP,
    SCREEN_ENGINE_DEFAULT,
    STRAIN_UNIQUE_PER_HALF_DEFAULT,
    TN_MODE_AUTO, TN_MODES,
    DEGENERATE_JUNCTION_CLEARANCE, PROBE_ROLE_PLAIN,
    BLAST_QUERY_CHUNK, BLAST_MAX_TARGET_SEQS,
    MIN_PROBE_GAP_BP, MIN_PROBE_GAP_FLOOR_BP,
    RELAXED_PROBE_GAP_BP,
    JUNCTION_WINDOW,
    REFERENCE_STRAND_DEFAULT,
)
from .errors import ConfigError
# Kept next to the classifier it governs rather than in constants: the default
# is only defensible alongside the measurement that justifies it, and
# reference_qc imports nothing from this package, so there is no cycle.
from .reference_qc import EXCLUDE_ALT_SCAFFOLDS_DEFAULT


@dataclass
class ModalityDefinition:
    """Describes an assay modality (Flex v1/v2, Visium spatial)."""
    name: str = "default"
    version: str = "v1"
    spatial: bool = False


@dataclass
class FilterThresholds:
    """Thresholds applied during the filter funnel.

    Defaults are sourced from ``flexprobe.constants`` (the single source of
    truth, cited to the 10x tech notes) so a bare ``FilterThresholds()`` matches
    the documented Flex spec — GC 44-72% per half, <=1 N per half — rather than
    drifting from the CLI/GUI defaults.
    """
    gc_min: float = GC_MIN_DEFAULT          # 0.44 (CG000839)
    gc_max: float = GC_MAX_DEFAULT          # 0.72 (CG000839)
    # GC is a RANKING term by default, not a gate. 10x's own panels put 8.4% of
    # their halves outside 44-72% while never once dropping below 56.3 C per
    # half, so the duplex-stability constraint they really honour is thermal.
    # Set True to restore the hard GC band (reproduces pre-0.4 runs).
    enforce_gc_range: bool = False
    # Per-half Tm gate (deg C), the replacement duplex-stability filter.
    # RELATIVE, not an in-assay melting temperature — primer3 calc_tm under the
    # PCR-primer salt conditions in constants.py, because 10x does not publish
    # the Flex hybridisation buffer. Floor derived from 10x's shipped panels
    # (observed minimum 56.32 C human / 56.35 C mouse). Set tm_min=None to
    # disable the gate; tm_max=None means no upper bound (the default — 10x's
    # panels show no upper boundary).
    tm_min: Optional[float] = TM_MIN_DEFAULT
    tm_max: Optional[float] = TM_MAX_DEFAULT
    # Degenerate-base budget per 25 nt half. Counts EVERY non-ACGT base — 'N'
    # and every IUPAC ambiguity code (constants.DEGENERATE_BASES) — so a
    # consensus window carrying a two-way code like 'R' is measured, not waved
    # through. Named max_n_per_half for continuity with the CLI flag and with
    # existing run_config.json files, although it counts more than 'N'.
    max_n_per_half: int = 1
    max_homopolymer_run: int = MAX_HOMOPOLYMER_RUN_DEFAULT  # reject runs > 4
    reject_junction_base_n: bool = True
    # No degenerate base within this many bp either side of the ligation
    # junction (constants.DEGENERATE_JUNCTION_CLEARANCE, CG000621 Rev D). A
    # check of only the two bases straddling the nick would admit ambiguity at
    # the neighbouring positions, where the ligase is least tolerant of one.
    degenerate_junction_clearance: int = DEGENERATE_JUNCTION_CLEARANCE
    # Same band, concrete bases. ``reject_junction_base_n`` above catches an
    # IUPAC ambiguity code sitting in the clearance band; this catches a probe
    # that is pure ACGT but carries a REAL mismatch there against a record it
    # is supposed to bind. CG000621 Rev D asks for both in one sentence:
    #
    #   "Avoid common single nucleotide polymorphisms (SNPs) and potential
    #    mismatches at the ligation junction. Refer to the UCSC Genome Browser
    #    and the Single Nucleotide Polymorphism Database (dbSNP). If avoiding
    #    SNPs is not possible, SNPs and mismatches should be at least four bp
    #    away from the ligation junction."
    #
    # The two are separate switches because they are separately evidenced: a
    # degenerate code is visible in the probe string itself, while a concrete
    # mismatch is only visible against the per-record alignment
    # (``ProbeCandidate.strain_binding``), which exists only on a multi-strain
    # run. Turning one off must not silently turn the other off with it.
    #
    # WHICH records count as "supposed to bind" is NOT decided here — it is
    # mode-dependent and is passed in per run; see
    # ``filters.on_target_record_ids``. Without that scope the gate is inert,
    # which is the safe direction: in strain-unique mode a junction mismatch
    # against the OTHER strains is the entire point of the design.
    reject_junction_mismatch: bool = True
    # Half-width (bp) of the ligation-junction zone, used for SCORING rather
    # than placement: the per-strain dTm penalty doubles inside it
    # (consensus._binding_score_for_window) and off-target hits report how many
    # of their mismatches fall in it (specificity.hit_junction_mismatches).
    #
    # Separate from degenerate_junction_clearance above, which is a hard
    # placement rule, because they answer different questions — "may an
    # ambiguous base sit here?" versus "how much does a mismatch here matter?".
    # Both default to CG000621's four bp, and most users should move them
    # together; they are separable because only one of them rejects windows.
    junction_window: int = JUNCTION_WINDOW
    # TN junction (mRNA index 25 = the base opposite the LHS probe's 3' end).
    # Master switch: True = enforce the junction adenine, False = no
    # junction-identity gate at all.
    #
    # Defaults to True, so the library and CLI agree with the GUI modality
    # presets, which also set it True. A bare FilterThresholds() therefore
    # applies the same junction gate as every other entry point. True is also
    # the conservative reading of both tech notes: CG000621 states the adenine
    # as a requirement, and CG000839,
    # having relaxed it, still concludes "for a simplified and robust design
    # process, 10x Genomics recommends using the T at the specified position".
    require_tn_junction: bool = True
    # Which gate to apply. "auto" (default) follows require_tn_junction above.
    # The gate is HARD ("strict_a", mRNA[25] == 'A') wherever it is engaged;
    # "soft" (reject only 'C') is an explicit opt-in available for Flex v1/v2
    # only — asking for it on Visium is a ConfigError. An explicit mode wins
    # over require_tn_junction, which is a default rather than a signal of
    # intent. See constants.TN_MODE_* for the verbatim citations and the panel
    # base composition behind the soft mode.
    tn_junction_mode: str = TN_MODE_AUTO


@dataclass
class SpecificityThresholds:
    """BLAST hit identity/length floors for on- and off-target classification.

    Defaults implement the 10x off-target rule as written (CG000621 Rev D /
    CG000839 Rev A): *"Matches to off-target genes should have at least five
    mismatches in at least one of the LHS or RHS probes to prevent efficient
    hybridization."* — i.e. mismatches counted over the FULL 25 nt half, and the
    pair rejected only when BOTH halves are too close to an off-target. On-target
    classification is unchanged (>=90% over 21+ bp).

    The two booleans below exist so runs made under the stricter, non-spec
    legacy interpretation can be reproduced exactly. Both default to the
    10x-faithful behaviour.
    """
    off_target_pident_min: float = OFF_TARGET_PIDENT_MIN
    off_target_aln_min: int = OFF_TARGET_ALN_MIN
    on_target_pident_min: float = ON_TARGET_PIDENT_MIN
    on_target_aln_min: int = ON_TARGET_ALN_MIN
    # Count mismatches over the whole probe half (True, the rule as written)
    # rather than over a >=off_target_aln_min bp local alignment (False, legacy).
    off_target_full_half: bool = OFF_TARGET_FULL_HALF_DEFAULT
    off_target_min_mismatch: int = OFF_TARGET_MIN_MISMATCH
    # Reject a pair only when BOTH halves fail the off-target gate (True, the
    # rule as written — RTL needs both halves to hybridise before the ligase can
    # make a library molecule). False restores "a hit on EITHER half rejects".
    off_target_require_both_halves: bool = OFF_TARGET_REQUIRE_BOTH_HALVES_DEFAULT
    # Scope the disjunction per off-target GENE, as 10x's sentence says: reject
    # only when both halves match the SAME subject. ON by default, on panel
    # evidence (see constants). False restores the whole-probe reading. Has no effect when off_target_require_both_halves is
    # False.
    off_target_same_subject: bool = OFF_TARGET_SAME_SUBJECT_DEFAULT
    # Junction-spanning clean duplex, in bp REQUIRED ON EACH SIDE of the nick
    # before a non-target gene counts as ligatable. ON by default at 9 bp, which
    # is the argmax of a measured separation curve, not a chosen round number —
    # see constants.OFF_TARGET_JUNCTION_DUPLEX_BP for the table.
    #
    # THIS ONE IS OURS, not 10x's. Every other field in this dataclass applies
    # a published 10x sentence (or restores a superseded reading of one); this
    # adds a criterion 10x does not state, and it only ever rejects MORE than
    # their count rule does. 0 disables it and restores exactly the published
    # rule; any other value below OFF_TARGET_JUNCTION_DUPLEX_MIN_BP is refused
    # rather than clamped, because the scan cannot see below that floor.
    off_target_junction_duplex_bp: int = OFF_TARGET_JUNCTION_DUPLEX_BP


@dataclass
class ScreenResources:
    """Paths and settings for BLAST databases."""
    db_paths: List[str] = field(default_factory=list)
    min_blast_version: str = "2.12.0"
    # Probe queries per blastn invocation. Bounds peak memory and gives
    # per-chunk progress on a run that would otherwise be one silent giant
    # search; it also sets the resume granularity, since a killed run restarts
    # at a chunk boundary. Larger = fewer process launches, bigger TSVs.
    query_chunk: int = BLAST_QUERY_CHUNK
    # Max subjects blastn may return per query half. BLAST+ defaults to 500 and
    # applies that cap HEURISTICALLY DURING the search — it is not guaranteed to
    # keep the best hits — so at the default a 25-mer against a transcriptome
    # can silently drop qualifying off-targets and let a non-specific probe
    # pass. Lower it only if you know why; the pipeline warns when any query
    # actually reaches the cap, which means possible truncation.
    max_target_seqs: int = BLAST_MAX_TARGET_SEQS
    skip_screen: bool = False
    use_cached_hits: bool = False
    cached_hits_path: Optional[str] = None
    # Which search engine backs the off-target screen (see screen.make_screen).
    # "kmer" by default; see constants.SCREEN_ENGINE_DEFAULT for why, and for
    # why the engine also changes ranking.
    engine: str = SCREEN_ENGINE_DEFAULT
    # Drop ALT haplotype and assembly PATCH scaffold records from the off-target
    # reference. ON by default: they are 7.6% of GRCh38 cdna.all's bases and
    # produced 100% of its spurious perfect-50bp off-target calls, because an
    # ALT/PATCH record is a near-verbatim copy of a primary region — so a good
    # probe is rejected for matching a second copy of its own gene. A no-op for
    # any reference without Ensembl location fields, which is every viral panel
    # and most custom FASTAs. k-mer engine only (BLAST searches a prebuilt
    # database, not the FASTA). See flexprobe.reference_qc.
    exclude_alt_scaffolds: bool = EXCLUDE_ALT_SCAFFOLDS_DEFAULT
    # Which strand(s) of the off-target reference to search: "both" (default)
    # or "plus" (skip minus-strand/reverse-complement
    # matches — a probe cannot hybridise to them, so counting one is a
    # phantom off-target). "plus" is only correct for a sense-oriented
    # cDNA/transcriptome reference; it is WRONG for a genome or a
    # negative-sense viral reference, neither of which this codebase can
    # currently tell apart from a sense cDNA file by inspection, so it stays
    # an explicit opt-in rather than the default. See
    # constants.REFERENCE_STRAND_DEFAULT for the full rationale.
    reference_strand: str = REFERENCE_STRAND_DEFAULT
    # Rebuild a k-mer index whose builder_sha256 does not match this
    # build's kmer_index.py, instead of stopping on
    # KmerIndexStampMismatchError. Off by default: a stale index is refused
    # loudly, never reused or rebuilt without being asked. k-mer engine only.
    rebuild_kmer_index: bool = False


@dataclass
class ExportSettings:
    """Controls which output writers are enabled."""
    write_idt_opools: bool = True
    write_cellranger_probe_set: bool = True
    write_results_table: bool = True
    write_manifest: bool = True
    overwrite: bool = False


@dataclass
class RunConfig:
    """Top-level configuration for a flexprobe pipeline run."""
    fasta_path: str = ""
    output_dir: str = "output"
    # Legacy header handling: truncate each FASTA record id at the first '|'.
    # OFF by default — an NCBI-style header set (``>gi|123|ref|NM_000546.6|``)
    # then collapses onto the id ``gi`` and the records silently overwrite one
    # another in the pipeline's {id: seq} dict. Opt in only for inputs where the
    # first pipe-field really is the intended, unique id.
    split_header_on_pipe: bool = False
    modality: ModalityDefinition = field(default_factory=ModalityDefinition)
    filters: FilterThresholds = field(default_factory=FilterThresholds)
    screen: ScreenResources = field(default_factory=ScreenResources)
    export: ExportSettings = field(default_factory=ExportSettings)
    tile_length: int = TILE_LENGTH_DEFAULT
    tile_step: int = 1
    # Multiplex vs singleplex — controls which RHS handle is used. Relevant
    # for Flex v1 and Flex v2; ignored for Visium.
    multiplex: bool = True
    # Probe Barcode labels. Single-element list = singleplex / one-barcode
    # multiplex; multi-element list = Flex v1 multi-sample multiplex, where
    # the RHS oligo carries a {BC} placeholder during discovery and is
    # expanded once per barcode at export time.
    probe_barcodes: List[str] = field(default_factory=lambda: ["BC001"])
    # Discovery mode: write only results.csv + manifest.json. The user reviews
    # the ranked candidates in the GUI and picks the final probes; the order
    # CSVs (idt_opools, cellranger_probe_set) are written later by
    # ``flexprobe.select_export.export_selected_probes`` based on that selection.
    discovery_only: bool = False
    # Whole-transcriptome convenience: if > 0 and discovery_only is True, the
    # run manager auto-picks the top-N probes per gene (by composite score,
    # kept >= MIN_PROBE_GAP_BP apart) and writes the order CSVs without review.
    auto_select_per_gene: int = 0
    # --- Multi-strain consensus mode ----------------------------------
    # When True, every record in the input FASTA is treated as a different
    # strain of the same target fragment. The pipeline builds a per-position
    # consensus, restricts probe candidates to either conserved or strain-
    # unique windows, and annotates each candidate with per-strain mismatch
    # info in results.csv.
    multi_strain: bool = False
    # Strategy:
    #   "conserved"     — probes identical across EVERY record (pan-strain)
    #   "strain_unique" — probes that distinguish ``strain_unique_target``
    #   "constitutive"  — isoform-collapse: probes shared by a QUORUM of records
    #                     (>= isoform_min_count); designed from the shared
    #                     sequence so one probe binds every covered isoform.
    consensus_strategy: str = "conserved"
    # Constitutive (isoform-collapse) quorum. 0 = require all records (strict);
    # else the minimum number of records a window must be shared by.
    isoform_min_count: int = 0
    # Which input record (0-based) is the consensus TEMPLATE for "conserved" and
    # "constitutive": it supplies the alignment coordinate system and, more
    # visibly, the base that every ambiguous consensus position resolves to when
    # the orderable probe is written out (see consensus.resolve_to_reference).
    #
    # This is not cosmetic. At a variant position the emitted probe carries the
    # template's base, so a run only ever yields the template's "flavour" of each
    # window; the other records' equally valid pan-strain probes are simply never
    # generated. On a two-strain input the two templates can give quite
    # different probe sets, and only their union recovers every pan-strain
    # probe that binds both strains. One run can be half the answer, which is
    # why the GUI can fan a run out across every record.
    #
    # Implemented by rotating the record list so the chosen record leads, i.e.
    # exactly equivalent to reordering the input FASTA. Ignored by
    # "strain_unique", where ``strain_unique_target`` already picks the sequence.
    consensus_reference: int = 0
    # When strategy == "strain_unique", which strain (by index into the
    # input FASTA records, 0-based) to design distinguishing probes for.
    strain_unique_target: int = 0
    # Minimum mismatches between the target strain and EVERY other strain
    # for a window to qualify as strain-unique. Defaults to 5 to match the
    # 10X off-target ≥5 mismatches rule (CG000621 / CG000839).
    strain_unique_min_mismatches: int = 5
    # Count those mismatches within a SINGLE 25 nt half (True, the reading the
    # ≥5 is borrowed from and the one specificity.py uses) rather than over
    # the whole 50 nt window (False, legacy). Mismatches split across the
    # junction leave both halves able to bind the non-target strain, so the
    # probe still hybridises and ligates there and is not discriminating.
    # See STRAIN_UNIQUE_PER_HALF_DEFAULT in constants.py.
    strain_unique_per_half: bool = STRAIN_UNIQUE_PER_HALF_DEFAULT
    # WITHIN a half, mismatches count the same whether they sit at the nick or
    # out at the probe's terminal end — but ligation only cares about the nick.
    # This reuses the off-target junction-duplex scan (same function, same
    # default 9 bp, same "0 disables" convention as
    # ``SpecificityThresholds.off_target_junction_duplex_bp`` above) to veto a
    # window whose per-half count clears the threshold but whose mismatches are
    # all far enough from the nick that a non-target strain would still present
    # a sealing duplex there — the strain-unique mirror of the off-target gate.
    # Ignored when strain_unique_per_half is False (see consensus.py). Only
    # ever removes windows the count rule already granted; never contradicts
    # 10x's published rule for the same reason the off-target twin does not.
    strain_unique_junction_duplex_bp: int = OFF_TARGET_JUNCTION_DUPLEX_BP
    # "conserved" strategy tolerance. DEFAULT 1: a window qualifies if each
    # 25 nt half carries at most one variant position, provided none falls
    # within ``filters.degenerate_junction_clearance`` bp of the ligation
    # junction. Set 0 for the strict rule, where a window qualifies only if
    # every record has identical bases across it.
    #
    # 1 is the default because 0 is not a neutral choice. It is a severe one
    # whose output reads as biology: on a real strain pair, exact identity can
    # leave a whole segment with no pan-strain probe once the chemistry filters
    # run, which says more about the rule than about the virus. 1 admits
    # windows carrying at most one variant position PER 25 nt half, provided no
    # variant falls within ``filters.degenerate_junction_clearance`` bp of the
    # ligation junction. The window is then emitted as a degenerate consensus
    # parent plus one concrete variant per record, and typically yields many
    # times more designable windows than exact identity does.
    conserved_max_degenerate_per_half: int = 1
    # Probe spacing used by the per-gene auto-pick (three-rung ladder: green
    # `min`, orange `relaxed`, red `floor` — see auto_select.SPACING_LADDER).
    # `min` is the design target: a hard floor below which a pair is never
    # chosen except by falling through to a looser rung. There is no
    # "preferred" rung in this ladder — that concept belongs to the SEPARATE
    # isoform-coverage selector (`select-coverage` / `isoform_select.coverage_pick`,
    # its own `--preferred-gap` flag), which is why this dataclass carries no
    # `preferred_probe_gap_bp` field: on the `design` path such a field would
    # have no tier to drive.
    min_probe_gap_bp: int = MIN_PROBE_GAP_BP
    # Fallback spacing for genes that cannot yield PROBE_GAP_RELAX_MIN_PROBES
    # probes at `min`. A value >= `min` skips this rung only; the `floor` rung
    # still runs unless it is also >= `min` (auto_select.ladder_from keeps a
    # rung only when it is tighter than the one above). Probes chosen under
    # it are always flagged, never silently emitted at the tighter spacing.
    relaxed_probe_gap_bp: int = RELAXED_PROBE_GAP_BP
    # Last rung. Used only when tighter than the rung above it, so a value
    # >= `relaxed` drops it and refuses red-tier spacing, accepting fewer
    # probes on the genes that would otherwise have needed it.
    floor_probe_gap_bp: int = MIN_PROBE_GAP_FLOOR_BP
    # Provenance only. The pipeline is fully deterministic (no RNG / sampling),
    # so this seed drives nothing at runtime; it is recorded in the manifest
    # purely to stamp the run. Do not read it as evidence of stochasticity.
    seed: int = 42
    # ON by default. It adds primer3 Tm, hairpin, homodimer and heterodimer
    # terms to the composite score, and it is the only way the degenerate
    # consensus aggregates (hairpin_dg_worst, tm_balance_max, the Tm/GC ranges)
    # get populated at all — a parent cannot be scored directly, so it carries
    # its variants' worst case or nothing.
    #
    # The cost is small: well under a millisecond per candidate, so seconds for
    # a viral panel and minutes for a whole transcriptome, against an
    # off-target screen measured in hours. Turning it off saves a negligible
    # fraction of a run's wall clock and leaves the aggregate columns blank.
    primer3_scoring: bool = True
    # Resumable BLAST: checkpoint each off-target BLAST chunk to disk
    # (``<output>/.screen_cache/``) so a run killed by an HPC walltime can be
    # resubmitted and skip the chunks it already finished. The checkpoints are
    # removed automatically on a fully successful run. On by default; turn off
    # for throwaway local runs where the extra per-chunk writes aren't wanted.
    resume_screen: bool = True
    # Stop after ingest + tile + filter. The candidate checkpoint and the
    # pre-screen artefacts are written, `PREPARED_MARKER` is dropped in the
    # output directory, and the run exits BEFORE the off-target screen. No
    # screen, no ranking, no chemistry, no export, no results.csv.
    #
    # NOT `screen.skip_screen`, which is a different request with a different
    # answer. `skip_screen` says "design probes without an off-target screen":
    # it runs the whole funnel and produces final, orderable output whose
    # specificity_status is `not_screened`. This says "do the front half of
    # this run now and the screen in a separate process later", and its output
    # directory is INCOMPLETE by construction — `design_summary.json` records
    # status `prepared`, and the manifest carries `prepare_only: true`.
    #
    # WHY IT EXISTS. Python does not return freed pages to the OS, so a
    # process that tiled a transcriptome carries the tiling high-water mark for
    # the rest of its life, and the screen then forks 8 workers off that
    # inflated parent. On whole-transcriptome runs, a process that did the
    # tiling peaks at roughly twice the memory of one that resumed from the
    # checkpoint and never tiled. Splitting the run at this line lets the two
    # halves ask the scheduler for what they each actually need instead of one
    # job asking for the maximum of both for its whole walltime.
    #
    # Requires `resume_screen` (the checkpoint is the whole deliverable) and
    # refuses to combine with skip_screen, cached hits or multi_strain — see
    # `validate`, which explains each.
    prepare_only: bool = False
    # WHAT the resume checkpoint stores. Only affects the streaming
    # (whole-transcriptome) path; the non-streaming path always needs raw hits.
    #
    #   "verdicts" (default) — store the ~11 scored fields per candidate.
    #       Size is INDEPENDENT of hit volume, which is the only thing that
    #       scales at transcriptome scale: a human transcriptome averages
    #       thousands of hits per probe, so raw hits run to gigabytes per chunk
    #       and terabytes per run. As verdicts the same chunk is a few MB.
    #       COST: a verdict is only valid for the thresholds that produced it, so
    #       changing `specificity` invalidates the cache and forces a re-BLAST.
    #
    #   "hits" — store the raw BLAST hits.
    #       Scoring-agnostic, so you can re-score at a new `specificity`
    #       threshold for free without re-running BLAST — worth it while tuning
    #       thresholds on a SMALL run. Also marginally smaller below ~2-3
    #       hits/probe (verdicts carry a fixed ~322 B/probe overhead) and keeps
    #       the crash-safer save-then-score order.
    #       COST: size scales linearly with hits/probe. Do NOT use at
    #       whole-transcriptome scale, where the cache outgrows ordinary disks.
    #
    # Switching modes changes the on-disk format, so it is part of the cache
    # fingerprint: flipping this invalidates an existing cache rather than
    # misreading it.
    screen_cache_mode: str = "verdicts"
    # Parallelism for the tile+filter stage on multi-gene (whole-transcriptome)
    # runs. 1 = the serial streaming path (default, unchanged). >1 fans genes
    # out across that many worker processes and merges the results
    # deterministically. Ignored for single-gene and multi-strain runs.
    workers: int = 1
    # Parallelism for the OFF-TARGET SCREEN, which is a different stage from
    # `workers` above and the one that dominates a whole-transcriptome run,
    # where tile+filter is a small fraction of the total. 1 = the serial path
    # (default).
    #
    # >1 splits each BLAST chunk into that many contiguous candidate slices and
    # searches AND SCORES each slice in a worker process, so only the ~11
    # VERDICT_FIELDS per probe cross the process boundary — never the hits. That
    # is what makes this memory-neutral: the hits for a chunk can run to
    # gigabytes, and W workers each hold 1/W of them rather than W copies.
    # Returning raw hits instead would multiply both the pickle cost and the
    # peak.
    #
    # Only the k-mer engine uses it. The BLAST engine's scoring can need
    # SubjectWindows, whose blastdbcmd prefetch is per-chunk shared state that a
    # per-slice split would fragment into one spawn per slice — the exact cost
    # SubjectWindows.prefetch exists to avoid.
    #
    # NOT part of any cache fingerprint (see candidate_cache.FINGERPRINT_IGNORED
    # and screen_cache.fingerprint): it changes who does the arithmetic, never
    # the answer, so raising it must resume an in-flight run rather than
    # invalidate its checkpoints.
    screen_workers: int = 1
    # Off/on-target BLAST stringency (see SpecificityThresholds).
    specificity: SpecificityThresholds = field(default_factory=SpecificityThresholds)
    # Isoform awareness: accessions (e.g. other transcript isoforms of the
    # target gene) that a valid probe legitimately hits and so must NOT be
    # counted as off-target. They are excused from off-target rejection without
    # being required by the strain rule. Matched by exact id or header boundary.
    extra_on_target_ids: List[str] = field(default_factory=list)
    # Whole-transcriptome speedup: if > 0, only the top-N candidates per gene
    # (by the BLAST-free pre-score — GC, balance, homopolymer, central
    # placement) are off-target screened; the rest are dropped before BLAST.
    # Cuts BLAST volume ~10-50x on a whole transcriptome while keeping tile_step
    # resolution where it matters, since you only select a few probes per gene
    # anyway. 0 = screen every candidate (default, unchanged behaviour).
    screen_top_per_gene: int = 0
    # Gene-aware (isoform) mode for whole-transcriptome runs. When True, the
    # pipeline groups the input transcripts by gene — parsed from Ensembl-style
    # ``gene:ENSG...`` header fields — so a probe's sibling isoforms count as
    # on-target (not off-target), and each probe is scored on isoform coverage
    # (the fraction of its gene's transcripts it binds). ON by default;
    # switched off for the run, with a log line, when no
    # input record carries a gene: field, so unlabelled input is unchanged.
    # Ignored in multi-strain mode (which has its own consensus logic).
    group_by_gene: bool = GROUP_BY_GENE_DEFAULT
    # Optional hard gate (0 = off): with group_by_gene, drop probes whose
    # isoform_coverage is below this fraction. Only applies when BLAST ran.
    # e.g. 0.8 keeps only probes binding >= 80% of their gene's isoforms.
    min_isoform_coverage: float = 0.0
    # --- pre-screen candidate set as an analysis artefact -------------------
    # results.csv holds only the candidates that SURVIVED the off-target screen,
    # so any statistic computed from it is conditioned on the outcome. A correct
    # null needs the set that was CONSIDERED. The pipeline materialises that set
    # (it is what goes into BLAST) but has only ever kept it as an internal,
    # schema-versioned resume cache that is deleted on success — unusable as
    # evidence. When True, write it out properly: output/candidates_prescreen.csv.gz
    # (one row per candidate BLAST screened, i.e. post-cap) plus
    # output/prescreen_cap_summary.csv.gz (per gene, how many passed filters
    # before the screen_top_per_gene cap vs how many the cap kept).
    # Off by default: it is a large extra file that most runs don't need.
    emit_prescreen_candidates: bool = False
    # Additionally stream EVERY filter-passing window, before the cap, to
    # output/candidates_precap.csv.gz. Strictly more information for the null —
    # the cap is itself an outcome-correlated selection step — and 20-50x the
    # rows: hundreds of millions at whole-transcriptome scale under a soft
    # junction gate. Requires emit_prescreen_candidates.
    emit_prescreen_precap_rows: bool = False
    # Keep <output>/.candidate_cache/ after a successful run instead of deleting
    # it. Belt-and-braces backup for the artefacts above (different code path,
    # different format), at the cost of leaving the resume checkpoint in the
    # output tree. It is dot-prefixed, so it still stays out of the tarball.
    keep_candidate_cache: bool = False
    # ---- null-model calibration of the off-target gate ----
    # Screen dinucleotide-shuffled copies of a bounded sample of this run's
    # windows through the same gate, and write the resulting false-positive rate
    # to output/null_model.json. A specificity gate whose rejection rate on
    # composition-matched noise is unknown cannot be defended in a methods
    # section. On test panels that rate varies widely between runs, which is
    # not a detail a reader can infer.
    #
    # Cost is bounded by null_model_max_windows regardless of run size: at the
    # defaults it is at most 200 x 20 = 4,000 extra BLAST queries whether the
    # run screened 40 windows or 35 million.
    #
    # 0 shuffles disables it. Off by default, because turning it on changes
    # every run's runtime and output surface. `--null-model` / `--null-model-shuffles` opt in.
    null_model_shuffles: int = 0
    null_model_max_windows: int = 200
    null_model_seed: int = 42
    # Optional overrides for the composite-ranking weights. Keys are a subset of
    # RankingWeights' fields (gc_proximity, gc_balance, specificity_margin,
    # homopolymer_penalty, central_placement); unspecified keys keep the default.
    ranking_weights: Optional[Dict[str, float]] = None

    def validate(self) -> None:
        """Validate the configuration, raising ``ConfigError`` on the first
        problem. Cheap, side-effect-free, and called at the top of
        ``run_pipeline`` so a bad config fails fast with a clear message
        instead of deep in a stage (or, worse, silently producing nonsense).
        """
        # Imported here, not at module scope: `screen` imports this module for
        # ProbeCandidate, so a top-level import would be a cycle. The names are
        # the engine ids `screen.make_screen` resolves.
        from .screen import ENGINE_BLAST as _ENGINE_BLAST, ENGINE_KMER as _ENGINE_KMER

        errs: List[str] = []
        if not (0.0 <= self.filters.gc_min <= self.filters.gc_max <= 1.0):
            errs.append(
                f"GC range must satisfy 0 <= gc_min ({self.filters.gc_min}) "
                f"<= gc_max ({self.filters.gc_max}) <= 1.0 (fractions, not %)."
            )
        if self.filters.tn_junction_mode not in TN_MODES:
            errs.append(
                f"tn_junction_mode must be one of {list(TN_MODES)}, got "
                f"{self.filters.tn_junction_mode!r}."
            )
        tmin, tmax = self.filters.tm_min, self.filters.tm_max
        if tmin is not None and tmax is not None and tmin > tmax:
            errs.append(
                f"Tm range must satisfy tm_min ({tmin}) <= tm_max ({tmax}) "
                f"(degrees C)."
            )
        if self.specificity.off_target_min_mismatch < 1:
            errs.append("off_target_min_mismatch must be >= 1.")
        _jd = self.specificity.off_target_junction_duplex_bp
        if _jd < 0:
            errs.append(
                f"off_target_junction_duplex_bp must be >= 0 (got {_jd}); "
                "0 disables the gate."
            )
        elif 0 < _jd < OFF_TARGET_JUNCTION_DUPLEX_MIN_BP:
            # Refused, not clamped: the junction scan seeds on an exact 12-mer
            # centred on the nick, so it cannot see a site with fewer than 6
            # clean bp per side. Accepting a smaller number would run the gate
            # at 6 while the config, the manifest and the Methods section all
            # said something else.
            errs.append(
                f"off_target_junction_duplex_bp must be 0 (off) or "
                f">= {OFF_TARGET_JUNCTION_DUPLEX_MIN_BP} (got {_jd}): the "
                "junction scan seeds on an exact 12-mer spanning the nick and "
                "cannot measure a clean run shorter than half of it. Below "
                "that the number would also be chance homology — a 12-mer "
                "occurs ~24 times by chance in a transcriptome."
            )
        # The gate is implemented ONLY in the k-mer engine. Refused, not warned
        # about and not silently dropped: the config, the manifest and the
        # Methods section would all still say the gate was at 9 while nothing
        # measured it, which is exactly the kind of unearned claim this
        # validation exists to stop. See `blast.py`'s module docstring for why
        # the BLAST path cannot produce the junction records the gate needs.
        # Only when a screen runs: under --skip-screen no gate is measured or
        # claimed by either engine, so refusing there would only block a
        # harmless config.
        if (_jd > 0 and self.screen.engine == _ENGINE_BLAST
                and not self.screen.skip_screen):
            errs.append(
                f"off_target_junction_duplex_bp is {_jd} but --screen-engine is "
                f"{_ENGINE_BLAST!r}, which cannot measure it: the gate needs "
                "junction-seeded records, and a BLAST hit list provably cannot "
                "contain the sites it exists to catch (a non-target gene "
                "matching perfectly across the nick has 6+ mismatches in each "
                "half, so blastn reports no alignment there at all). Either run "
                f"--screen-engine {_ENGINE_KMER} (the default), or pass "
                "--no-off-target-junction-duplex (or "
                "--off-target-junction-duplex 0) to turn the gate off "
                "explicitly and accept the weaker screen."
            )
        if self.filters.max_n_per_half < 0:
            errs.append("max_n_per_half must be >= 0.")
        if self.filters.max_homopolymer_run < 1:
            errs.append("max_homopolymer_run must be >= 1.")
        if not (1 <= self.filters.junction_window <= max(1, self.tile_length // 2)):
            errs.append(
                f"junction_window must be between 1 and the {self.tile_length // 2} "
                f"bp half (got {self.filters.junction_window})."
            )
        if self.tile_length < 2:
            errs.append(f"tile_length must be >= 2 (got {self.tile_length}).")
        elif self.tile_length % 2 != 0:
            errs.append(
                f"tile_length must be even so it splits into two equal halves "
                f"(got {self.tile_length})."
            )
        else:
            # Chemistry assembly drops every tile that is not exactly the
            # chemistry's window (pipeline Stage 7), so any other length would
            # validate and then yield ZERO probes with no error. Refuse it
            # here, where the user can still fix it.
            try:
                from .chemistry import get_chemistry
                _chem_len = get_chemistry(self.modality).tile_length
            except Exception:       # an unknown modality is reported elsewhere
                _chem_len = None
            if _chem_len is not None and self.tile_length != _chem_len:
                errs.append(
                    f"tile_length {self.tile_length} cannot be assembled: the "
                    f"{self.modality.name} {self.modality.version} chemistry "
                    f"needs {_chem_len} bp windows, and every other length "
                    f"yields no probes."
                )
        if self.tile_step < 1:
            errs.append(f"tile_step must be >= 1 (got {self.tile_step}).")
        if self.conserved_max_degenerate_per_half < 0:
            # consensus would otherwise clamp a negative value to 0 silently.
            errs.append(
                f"conserved_max_degenerate_per_half must be >= 0 (got "
                f"{self.conserved_max_degenerate_per_half})."
            )
        if self.workers < 1:
            errs.append(f"workers must be >= 1 (got {self.workers}).")
        if self.screen_workers < 1:
            errs.append(
                f"screen_workers must be >= 1 (got {self.screen_workers})."
            )
        if self.consensus_strategy not in ("conserved", "strain_unique", "constitutive"):
            errs.append(
                f"consensus_strategy must be 'conserved', 'strain_unique' or "
                f"'constitutive' (got {self.consensus_strategy!r})."
            )
        if self.isoform_min_count < 0:
            errs.append("isoform_min_count must be >= 0 (0 = all records).")
        if self.consensus_reference < 0:
            errs.append("consensus_reference must be >= 0 (index into the "
                        "input FASTA records).")
        if not (0.0 <= self.min_isoform_coverage <= 1.0):
            errs.append(
                f"min_isoform_coverage must be between 0.0 and 1.0 "
                f"(got {self.min_isoform_coverage})."
            )
        if self.multi_strain and self.consensus_strategy == "strain_unique":
            if self.strain_unique_target < 0:
                errs.append("strain_unique_target must be >= 0.")
            if self.strain_unique_min_mismatches < 1:
                errs.append("strain_unique_min_mismatches must be >= 1.")
            # In per-half mode the count is capped by the half, not the tile, so
            # a threshold above the half length can never be met and every run
            # would silently return zero windows.
            half = self.tile_length // 2
            if (self.strain_unique_per_half
                    and self.strain_unique_min_mismatches > half):
                errs.append(
                    f"strain_unique_min_mismatches "
                    f"({self.strain_unique_min_mismatches}) exceeds the "
                    f"{half} bp probe half it is counted over in per-half mode; "
                    f"no window can qualify. Lower it to <= {half}, or pass "
                    f"--strain-unique-whole-window to count over the whole "
                    f"{self.tile_length} bp window."
                )
        _sujd = self.strain_unique_junction_duplex_bp
        if _sujd < 0:
            errs.append(
                f"strain_unique_junction_duplex_bp must be >= 0 "
                f"(got {_sujd}); 0 disables the gate."
            )
        elif 0 < _sujd < OFF_TARGET_JUNCTION_DUPLEX_MIN_BP:
            # Same floor as off_target_junction_duplex_bp and for the same
            # reason: junction_clean_duplex cannot see a run shorter than this.
            errs.append(
                f"strain_unique_junction_duplex_bp must be 0 (off) or "
                f">= {OFF_TARGET_JUNCTION_DUPLEX_MIN_BP} (got {_sujd}): see "
                "off_target_junction_duplex_bp above for why."
            )
        # NB: probe_barcodes may legitimately be empty — e.g. Flex v2, where the
        # GUI doesn't select a per-sample barcode and the pipeline applies the
        # chemistry default. So it is intentionally NOT validated as non-empty.
        if self.ranking_weights is not None:
            valid_keys = {
                "gc_proximity", "gc_balance", "specificity_margin",
                "homopolymer_penalty", "central_placement", "isoform_coverage",
                "mismatch_placement",
            }
            bad = set(self.ranking_weights) - valid_keys
            if bad:
                errs.append(
                    f"ranking_weights has unknown key(s) {sorted(bad)}; "
                    f"valid keys are {sorted(valid_keys)}."
                )
        if self.auto_select_per_gene < 0:
            errs.append("auto_select_per_gene must be >= 0.")
        if self.screen_top_per_gene < 0:
            errs.append("screen_top_per_gene must be >= 0 (0 = screen all).")
        if self.emit_prescreen_precap_rows and not self.emit_prescreen_candidates:
            errs.append(
                "emit_prescreen_precap_rows needs emit_prescreen_candidates: the "
                "pre-cap dump is an addition to the pre-screen artefacts, not a "
                "replacement for them."
            )
        # --prepare-only's entire deliverable is the candidate checkpoint that a
        # later screen job resumes from. Each of these four combinations would
        # produce a run that exits having written no checkpoint at all, i.e. a
        # job that burned its walltime for nothing. Refused rather than warned
        # about, because the failure is silent: the marker would still be
        # written, and the screen job would then tile the transcriptome again
        # inside its own smaller allocation.
        if self.prepare_only:
            if not self.resume_screen:
                errs.append(
                    "prepare_only needs resume_screen: the candidate checkpoint "
                    "IS the deliverable, and resume_screen=False is what turns "
                    "checkpointing off (see pipeline._load_candidate_checkpoint). "
                    "Drop --no-resume."
                )
            if self.screen.skip_screen:
                errs.append(
                    "prepare_only and skip_screen are contradictory requests. "
                    "skip_screen designs probes WITHOUT an off-target screen and "
                    "produces final output; prepare_only stops before the screen "
                    "and produces no probes, so a later job can run it. Pick one."
                )
            if self.screen.use_cached_hits:
                errs.append(
                    "prepare_only with use_cached_hits has nothing to prepare: "
                    "pre-computed hits bypass the screen entirely, and the "
                    "candidate checkpoint is only written on the path that will "
                    "screen (pipeline._load_candidate_checkpoint gates on "
                    "screen_hits_by_probe being absent)."
                )
            if self.multi_strain:
                errs.append(
                    "prepare_only cannot split a multi_strain run: the consensus "
                    "alignment is rebuilt on every invocation and the candidate "
                    "checkpoint is deliberately not used there, so the screen job "
                    "would redo all of the tiling inside its own allocation. Run "
                    "a multi-strain design as a single job."
                )
        if errs:
            raise ConfigError("Invalid RunConfig:\n  - " + "\n  - ".join(errs))


@dataclass
class ProbeCandidate:
    """A single candidate probe target site on the mRNA (50-bp window).

    ``sequence`` is the 50-bp sense-strand window from the input FASTA (the
    "original gene/genome sequence"). Use the helpers in ``flexprobe.chemistry``
    to derive the actual LHS/RHS oligo strings.
    """
    gene_id: str = ""
    probe_id: str = ""
    sequence: str = ""             # 50-bp mRNA target window, sense
    start: int = 0
    end: int = 0
    gc_content: float = 0.0
    gc_left: float = 0.0
    gc_right: float = 0.0
    homopolymer_max: int = 0
    passed_filters: bool = True
    filter_fail_reasons: List[str] = field(default_factory=list)
    specificity_pass: bool = False
    specificity_margin: float = 0.0
    # Quantitative off-target evidence for a passing probe, so a
    # "passed" verdict carries the numbers behind it: the closest near-miss
    # off-target's identity / BLAST bitscore / e-value, and how many near-miss
    # off-targets cleared the alignment-length floor. All 0 when nothing came
    # close. The pass/fail decision itself still uses the 10x-sourced
    # pident+length rule (bitscore/e-value are reported, not the gate).
    offtarget_nearmiss_count: int = 0
    best_offtarget_pident: float = 0.0
    best_offtarget_bitscore: float = 0.0
    best_offtarget_evalue: float = 0.0
    # Explicit provenance for the specificity verdict so a "pass" can't be
    # silently confused with "never screened". One of:
    #   passed               — BLAST ran, on-target hit, no off-target
    #   failed_off_target    — BLAST ran, an off-target hit cleared thresholds
    #   failed_junction_duplex — the 10x COUNT rule passed this pair, and our
    #                   stricter junction gate rejected it: a non-target gene
    #                   presents a clean matched duplex >= off_target_junction_
    #                   duplex_bp on BOTH sides of the ligation nick, so it can
    #                   hybridise AND ligate. Deliberately NOT folded into
    #                   failed_off_target: one status is the vendor's rule and
    #                   the other is ours, and a reader has to be able to tell
    #                   which sentence rejected a probe.
    #   no_on_target_hit     — BLAST ran but no on-target hit found
    #   not_screened         — no BLAST DB configured / BLAST skipped
    #   blast_failed         — BLAST errored; treated as unscreened
    specificity_status: str = "not_screened"
    # HOW the off-target mismatch counts behind that verdict were made, so a run
    # can never be read as exactly-recounted when it was only bounded. One of:
    #   exact         — every off-target hit was recounted against the real
    #                   subject bases retrieved from the BLAST DB
    #   conservative  — at least one subject window could not be retrieved, so
    #                   its unaligned query positions were counted as MATCHES.
    #                   The counts are LOWER BOUNDS: the screen may over-reject,
    #                   never over-pass.
    #   legacy_local  — off_target_full_half=False; the legacy local-alignment
    #                   rule ran and no full-half count was made
    #   not_screened  — no BLAST hits were scored (BLAST skipped or failed)
    # See specificity.full_half_mismatches for why the distinction matters: a
    # count inferred from HSP geometry alone could only ever OVERSTATE
    # specificity.
    offtarget_recount_basis: str = "not_screened"
    # How much of a probe half BLAST actually ALIGNED, at its worst across the
    # non-target hits behind the verdict (0-1; NaN when there were no non-target
    # hits at all, which is the strongest pass and must not read as coverage 0).
    #
    # The companion to the basis above. `offtarget_recount_basis` says WHETHER
    # unexamined positions had to be reconstructed from the database;
    # this says HOW MUCH of the half was unexamined to begin with. A probe that
    # passed with every hit spanning all 25 nt is a different kind of pass from
    # one that passed on a 9 nt core plus a blastdbcmd lookup, and this column
    # is what tells them apart in results.csv.
    offtarget_min_qcov: float = float("nan")
    # Profile of the CLOSEST disqualifying off-target: how many mismatches it
    # sits at over the full half, and how many of those fall next to the
    # ligation junction (-1 = unknown, e.g. the HSP never reached that zone).
    #
    # The junction count reads BACKWARDS from the on-target case and that is the
    # point. Against a strain you mean to bind, a junction mismatch is fatal —
    # the ligase will not seal across it. Against an off-target it is
    # PROTECTIVE for the same reason. Four mismatches bunched at the junction
    # make a harmless off-target; four at the outer ends barely destabilise it,
    # so it hybridises, ligates, and puts reads on the wrong gene. The flat
    # ">= 5 in one half" rule scores those two identically.
    #
    # Reported, not gated: the gate stays the vendor's flat count so it remains
    # traceable to CG000621/CG000839. This is the evidence for a reader (or a
    # future weighted rank) to see past it.
    offtarget_worst_mismatches: int = -1
    offtarget_worst_junction_mm: int = -1
    # The junction-spanning CLEAN DUPLEX against the worst non-target gene: the
    # unbroken run of matched bases immediately 5' of the ligation nick and
    # immediately 3' of it, measured in CONTIGUOUS reference sequence so both
    # halves are necessarily adjacent on one molecule in one gene. That
    # adjacency is the whole point — it is the condition under which the ligase
    # can actually seal, which the flat per-half mismatch count never checks.
    #
    # This is the evidence behind the ``failed_junction_duplex`` verdict, and it
    # is reported whether or not the gate is enabled — like
    # offtarget_worst_junction_mm above, a reader should be able to see the
    # number even in a run that did not act on it.
    #
    # THREE READINGS, and they are different things:
    #   -1  the scan did not run (BLAST engine, screen skipped, hits supplied
    #       from a cache written before the scan existed). NOT zero: unknown
    #       must never render as "measured, and nothing found".
    #    0  the scan ran and found no non-target site at all. Its detection
    #       floor is constants.OFF_TARGET_JUNCTION_DUPLEX_MIN_BP per side, so
    #       this reads "below 6 on at least one side", not "exactly 0".
    #   >0  measured. Always >= 6 on both sides, for the same reason: values of
    #       1-5 cannot appear, because a site that short is not enumerated.
    offtarget_clean_left_bp: int = -1
    offtarget_clean_right_bp: int = -1
    # ---- isoform coverage (group-by-gene / whole-transcriptome mode) ----
    # When group_by_gene is on, the input transcripts are grouped by gene and a
    # probe is scored on how many of its gene's isoforms it binds.
    #   isoform_total    — number of transcripts in this probe's gene
    #   isoform_hits     — how many of them the probe binds (BLAST on-target)
    #   isoform_coverage — isoform_hits / isoform_total
    # Defaults (1 / 1 / 1.0) mean "not grouped" or a single-isoform gene, so
    # they have no effect on ranking outside gene-aware runs.
    isoform_total: int = 1
    isoform_hits: int = 1
    isoform_coverage: float = 1.0
    # The actual transcript ids this probe binds (the gene's isoforms whose
    # sequence the probe's target matches on-target). Persisted to results.csv
    # as a ';'-joined ``isoform_hit_ids`` column so downstream selection knows
    # *which* isoforms each probe covers without re-reading the input FASTA.
    isoform_hit_ids: List[str] = field(default_factory=list)
    composite_score: float = 0.0
    # Optional primer3 thermodynamics (only populated when primer3_scoring is
    # on). ΔG in cal/mol (more negative = stronger structure); NaN
    # when not computed. Relative indicators, not in-assay values.
    # All three are computed on the target-binding 25-mers actually synthesised
    # (target_LHS / target_RHS), NOT on the mRNA window — see ranking.py.
    tm_balance: float = 0.0          # |Tm(target_LHS) - Tm(target_RHS)| in °C
    self_dimer_dg: float = 0.0       # worst homodimer ΔG of the two 25-mers
    hetero_dimer_dg: float = 0.0     # target_LHS-vs-target_RHS heterodimer ΔG
    # Per-component 0-1 sub-scores that make up composite_score, so a rank is
    # auditable (which axis drove it). Populated by compute_composite_score.
    # NOTE: composite_score is a non-calibrated *preference* score (weighted sum
    # of these), not a probability of probe success.
    score_components: Dict[str, float] = field(default_factory=dict)
    rank: int = 0
    reduced_confidence: bool = False
    # ---- chemistry-derived (filled in by pipeline after assembly) ----
    probe_rc: str = ""            # 50-bp RC of mRNA target — what the probe contains
    target_lhs: str = ""          # 25 bp, goes into LHS oligo
    target_rhs: str = ""          # 25 bp, goes into RHS oligo
    lhs_oligo: str = ""           # full LHS oligo (with 5' handle)
    rhs_oligo: str = ""           # full RHS oligo (/5Phos/, target_RHS, 3' handle)
    # ---- ligation-junction identity ----
    # ``junction_base`` is the mRNA base at index 25 of the 50-bp window — the
    # base opposite the LHS probe's 3'-most nucleotide. ``tn_junction_ok`` is
    # the derived convenience boolean (junction_base == 'A', i.e. a true TN
    # motif). BOTH are populated by the filter funnel for every candidate that
    # reaches it, whatever TN gate is in force, so ranking can always prefer a
    # real TN junction — that is what keeps A-junction probes on top when the
    # soft ('not C') gate admits G and T. ``junction_base`` is '' only for a
    # candidate shorter than one half, or one never passed through filters.
    junction_base: str = ""
    tn_junction_ok: bool = False
    # Melting temperatures (deg C). Computed via primer3.calc_tm where
    # available; NaN if primer3 is missing or the sequence contains any
    # non-ACGT base (primer3 rejects N / IUPAC ambiguity codes). Returning
    # NaN rather than 0.0 keeps "couldn't compute" distinguishable from
    # a genuinely low Tm in results.csv.
    tm_target_lhs: float = float("nan")
    tm_target_rhs: float = float("nan")
    tm_full_50: float = float("nan")
    # ---- Multi-strain mode only ----
    # Filled in by the consensus pipeline. Each entry is a small dict
    # ``{strain_id, total_mismatches, lhs_mismatches, rhs_mismatches,
    # junction_mismatches, mismatch_positions, estimated_dTm, likely_binds}``
    # — one entry per strain in the input alignment. Empty in single-
    # sequence mode.
    strain_binding: List[Dict[str, Any]] = field(default_factory=list)
    # Convenience scalars for sorting / ranking in multi-strain mode.
    min_strain_mismatches: int = 0
    max_strain_mismatches: int = 0
    worst_strain_dTm: float = 0.0
    # The input records THIS window was designed to cover, when that is a subset
    # of the alignment and the run's mode alone cannot say which subset.
    #
    # Populated only by constitutive (isoform-collapse) mode, where each window
    # is the shared sequence of its largest isoform group and the isoforms
    # outside that group are not targets:
    # ``consensus.constitutive_start_positions`` returns the covered indices per
    # window and the pipeline maps them to record ids here. Empty everywhere
    # else, which is what keeps every other path on its pre-existing behaviour.
    #
    # *** POPULATING THIS NARROWS A SAFETY GATE. *** It is the per-candidate
    # scope for the junction-mismatch rule (filters.candidate_on_target_records)
    # and it OVERRIDES the run-level scope, so a path that sets it on a probe
    # required to bind more records than it lists would silently stop gating
    # those records. Set it only where the covered set is the whole set of
    # records the probe must bind.
    covered_record_ids: List[str] = field(default_factory=list)
    # ---- Degenerate (consensus) probes: the parent / variant model ----
    # See constants.PROBE_ROLE_*. A single-sequence run, and any window that is
    # identical across every record, leaves all of these at their defaults and
    # is unaffected by the parent / variant model.
    #   probe_role           — "", "consensus" (degenerate parent) or "variant"
    #   parent_probe_id      — the consensus probe a variant resolves from
    #   resolves_to          — record ids whose sequence this variant matches
    #                          (more than one when records resolve identically)
    #   degenerate_positions — 0-based offsets into the 50-mer carrying IUPAC
    #   iupac_codes          — the code at each of those offsets, same order
    #   n_variants           — distinct variants under a consensus parent
    probe_role: str = PROBE_ROLE_PLAIN
    parent_probe_id: str = ""
    resolves_to: List[str] = field(default_factory=list)
    degenerate_positions: List[int] = field(default_factory=list)
    iupac_codes: List[str] = field(default_factory=list)
    n_variants: int = 0
    # Aggregates over a consensus parent's variants, NaN on every other row.
    #
    # A degenerate parent cannot be scored directly — primer3 rejects non-ACGT
    # input — but its variants can, and they all go into the same pool. What
    # decides whether that pool works is therefore its WORST member, not an
    # absent number: leaving these blank would hide exactly the probes carrying
    # the most structural risk from anyone sorting the review table by hairpin
    # or Tm. So the parent carries the range (Tm, GC) or the worst case
    # (hairpin, heterodimer, Tm balance) across its own variants, and
    # ``reduced_confidence`` on a parent means "aggregated from variants",
    # not "unscored".
    #
    # This is also where the IUPAC widening pays off: an 'R' has two
    # resolutions and a narrow range, an 'N' has four and a wide one, so the
    # spread itself reports how much the ambiguity costs.
    tm_lhs_min: float = float("nan")
    tm_lhs_max: float = float("nan")
    tm_rhs_min: float = float("nan")
    tm_rhs_max: float = float("nan")
    gc_lhs_min: float = float("nan")
    gc_lhs_max: float = float("nan")
    gc_rhs_min: float = float("nan")
    gc_rhs_max: float = float("nan")
    hairpin_dg: float = float("nan")         # per-variant, worst of the two halves
    hairpin_dg_worst: float = float("nan")   # parent aggregate over variants
    hetero_dimer_dg_worst: float = float("nan")
    tm_balance_max: float = float("nan")


@dataclass
class ProbePair:
    """A paired probe (LHS + RHS) for Flex-style assays.

    Both halves point to the same underlying tile-level ``ProbeCandidate``
    (held on ``lhs``); ``rhs`` carries the per-strand oligo string for the
    benefit of writers that need it.
    """
    lhs: ProbeCandidate = field(default_factory=ProbeCandidate)
    rhs: ProbeCandidate = field(default_factory=ProbeCandidate)
    probe_barcode_label: str = ""    # e.g. "BC001"
    probe_barcode_seq: str = ""      # e.g. "ACTTTAGG"
    pair_score: float = 0.0


@dataclass
class RunManifest:
    """Records provenance for reproducibility."""
    # Empty by default; the manifest writer fills in flexprobe.__version__.
    app_version: str = ""
    timestamp: str = field(
        default_factory=lambda: datetime.datetime.now(
            datetime.timezone.utc
        ).isoformat()
    )
    # Identity of every input file: the path AND its content hash, so "which
    # FASTA produced Table 1" is answerable after the file has been moved,
    # renamed, or edited. A bare path is a claim about a filesystem that may
    # no longer exist by the time anyone asks.
    inputs: Dict[str, str] = field(default_factory=dict)
    # The COMPLETE RunConfig this run executed (see pipeline.run_pipeline). Not
    # a curated subset: "what GC range and off-target identity produced this?"
    # must be answerable from the outputs alone.
    params: Dict[str, Any] = field(default_factory=dict)
    # Which build of the engine ran. app_version alone is too coarse — several
    # commits ship under one version, and validation runs are often made from a
    # working tree. Keys: 'commit', 'dirty', 'describe'; absent when the source
    # isn't a git checkout (e.g. an installed wheel or the packaged app).
    code_version: Dict[str, str] = field(default_factory=dict)
    # External tool + library versions actually used by this run (blastn,
    # primer3, python, ...). Manuscript Methods sections need these; a value is
    # simply absent when the tool couldn't be found or probed.
    tool_versions: Dict[str, str] = field(default_factory=dict)
    # Wall-clock seconds per pipeline stage (ingest, tile, filter, blast,
    # specificity, rank, chemistry, export) and 'total'. A whole-transcriptome
    # run takes days, and runtime is worth comparing across runs, so it is
    # recorded rather than left to the job log.
    elapsed_seconds: Dict[str, float] = field(default_factory=dict)
    # SHA-256 of each BLAST database actually screened against, keyed by DB
    # prefix. "Specificity passed" is only interpretable against the database
    # that produced it, and a DB prefix is a path on a cluster that gets
    # rebuilt — the hash is what makes the claim checkable.
    db_sha256s: Dict[str, str] = field(default_factory=dict)
    # What each of those references CONTAINED, keyed the same way: record and
    # base counts split primary vs ALT/PATCH, the gene-biotype breakdown, and
    # whether non-primary records were excluded from the search. The hash above
    # IDENTIFIES a reference; this CHARACTERISES it, and the two answer
    # different questions — 7.6% of GRCh38 cdna.all's bases are near-verbatim
    # copies of its own primary assembly, which is why two runs can cite the
    # same hash and still have incomparable rejection rates. See
    # flexprobe.reference_qc.
    reference_composition: Dict[str, Any] = field(default_factory=dict)
    # WHICH ALIGNER BUILT THE CONSENSUS, for multi-strain runs only.
    # ConsensusResult.aligner is also emitted as a progress event, but `_emit`
    # passes only (stage, msg, pct) to the timeline writer, whose record is a
    # fixed ten keys, so without this field the value would reach only a log
    # line and the GUI callback. A Methods section cannot cite a log line, and
    # the bundled MAFFT and the Biopython pairwise fallback do not produce the
    # same alignment, so the choice is load-bearing.
    #
    # ABSENT, not empty, on a single-species run. _stage_multi_strain_tiling
    # never executes there, so there is no consensus and no aligner to name;
    # `params['multi_strain']` is False in that case and says so unambiguously.
    # Do not read an absent key here as a failure to record.
    consensus: Dict[str, Any] = field(default_factory=dict)
    seed: int = 42   # provenance only; the pipeline is deterministic (no RNG)
    # Keys the run's config did NOT set, filled from this build's defaults
    # when it was loaded from a file (config_io.filled_from_defaults). None
    # (written as null) when the config was built in memory.
    config_filled_from_defaults: Optional[List[str]] = None
