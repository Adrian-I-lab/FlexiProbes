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

"""Central home for tunable constants, each cited to its source.

Keeping these in one place (rather than scattered across modules) makes it
obvious which numbers come from a 10x Genomics tech note vs. an internal
heuristic, and gives a single edit point when a threshold needs tuning.
"""
from __future__ import annotations

# ---------------------------------------------------------------------------
# BLAST
# ---------------------------------------------------------------------------
# Minimum supported BLAST+ version. 2.12 is the first release whose
# blastn-short output format matches what parse_hits_tsv expects.
MIN_BLAST_VERSION = (2, 12, 0)

# How many probe queries to put in a single blastn invocation. Larger DBs
# (human transcriptome) make a single giant query slow with no feedback;
# chunking gives per-chunk progress and bounds peak memory.
BLAST_QUERY_CHUNK = 5000

# Max off-target subjects BLAST may return per query half. BLAST+ defaults to
# 500 and applies the cap *heuristically during* the search — it is NOT
# guaranteed to keep the best hits — so with the default a 25-mer against a
# transcriptome-scale DB can silently drop qualifying off-targets and let a
# non-specific probe pass. We set it high and warn if a query actually reaches
# the cap (possible truncation → re-run with a larger value). Internal safety
# choice, not a 10x parameter.
BLAST_MAX_TARGET_SEQS = 100000


# Which engine backs the off-target screen. "kmer" is the default; "blast" is
# still selectable (--screen-engine blast) and is required for a reference the
# k-mer index cannot address. See kmer_index.MAX_INDEXABLE_BASES.
#
# Defined as a literal rather than imported from ``screen`` because this module
# imports nothing from the package; ``test_screen_abstraction`` asserts it is a
# member of ``screen.ENGINES`` so the two cannot drift apart.
#
# WHY k-mer. blastn reports the maximum-scoring LOCAL segment while the 10x rule
# is stated over the whole 25 nt half, so BLAST's verdicts rest on hits it may
# never have reported. Where the engines disagree, checking against exhaustive
# brute force (not against either engine) shows BLAST erring BOTH ways: it
# certifies non-specific probes clean (unsafe) and discards specific ones
# (yield). Under BLAST, `offtarget_min_qcov` can sit well below 1.0, meaning the
# verdict rests on reconstructed positions. Under the k-mer engine that column
# is 1.0 and `offtarget_recount_basis` is `exact`, both by construction rather
# than by successful retrieval.
#
# THE ENGINE ALSO CHANGES RANKING, substantially. Top-N selections from the two
# engines can be near-disjoint. The cause is the 0.30-weight
# `specificity_margin` term: under BLAST it is derived from trimmed-HSP
# identity, which the 9 nt pushdown floor pins at ~100% for nearly every probe,
# so it is frequently a single constant carrying no information. Under the
# k-mer engine the term discriminates between probes. Runs made with different
# engines are not rank-comparable; `screen_engine` in results.csv and the
# manifest is what tells them apart.
SCREEN_ENGINE_DEFAULT = "kmer"


# Which strand(s) of the off-target reference the screen searches.
#
# A probe half is single-stranded DNA that hybridises to a sense-strand
# mRNA/cDNA molecule. A hit that only matches the reference's MINUS strand is
# sequence the probe is chemically incapable of binding — see
# ``chemistry.py`` (the probe sequence is RC(target)) — so counting it as an
# off-target is a "phantom" rejection of a probe that is actually fine. Such
# hits account for a large share of the junction-duplex gate's chance-driven
# firings at the default 9 bp threshold.
#
# THE ASSUMPTION THAT MAKES "plus" CORRECT, AND WHY IT IS NOT THE DEFAULT.
# "plus" is only right when the reference is a SENSE-ORIENTED
# transcriptome/cDNA — every Ensembl ``cdna.all``, GENCODE ``transcripts`` or
# Cell Ranger reference this tool has been validated against. It is WRONG for
# two reference kinds this tool also explicitly supports and neither
# ``reference_qc`` nor anything else in the codebase can currently tell apart
# from a sense cDNA file by inspection:
#   1. A GENOME FASTA. README.md documents ``--screen-engine blast`` as
#      "required for a reference larger than the k-mer index can address
#      (~2.1 Gbases -- i.e. a genome)" -- i.e. genome-scale references are a
#      real, supported use case, not a hypothetical. On a genome, a
#      minus-strand hit is a gene on the opposite strand that really can
#      transcribe the off-target -- a GENUINE off-target, not a phantom one.
#   2. A negative-sense viral reference deposited in vRNA (genomic) sense,
#      for example an orthobunyavirus such as OROV. There the
#      mRNA is the reverse complement of the deposited record, so restricting
#      to "plus" would make the PHANTOM hit invisible and the REAL one is the
#      minus-strand hit that "plus" would discard.
# ``reference_qc.py`` classifies assembly (primary/ALT-PATCH) and gene
# biotype, never orientation, so there is no existing signal to gate on
# automatically. Defaulting to "plus" would therefore silently convert a
# fail-SAFE bug (discards a good probe; costs yield) into a fail-OPEN one
# (certifies a probe against a real off-target; costs the experiment) for
# anyone running a genome or negative-sense-viral reference. "both" is
# therefore the default; "plus" is an explicit opt-in for a
# reference the caller asserts is sense-oriented cDNA (``--reference-strand
# plus`` / ``ScreenResources.reference_strand``).
REFERENCE_STRAND_BOTH = "both"
REFERENCE_STRAND_PLUS = "plus"
REFERENCE_STRANDS = (REFERENCE_STRAND_BOTH, REFERENCE_STRAND_PLUS)
REFERENCE_STRAND_DEFAULT = REFERENCE_STRAND_BOTH


# ---------------------------------------------------------------------------
# Specificity (off-/on-target classification of BLAST hits)
# ---------------------------------------------------------------------------
# Source, VERBATIM and COMPLETE — CG000621 Rev D (p.3, Visium section; p.12,
# Chromium Flex section) and CG000839 Rev A (p.2), identical wording in all
# three places:
#
#   "Matches to off-target genes should have at least five mismatches in at
#    least one of the LHS or RHS probes to prevent efficient hybridization."
#
# Two things follow from the full sentence, and both are lost if it is quoted
# truncated at "...at least five mismatches":
#
#  1. The clause "in at least one of the LHS or RHS probes" makes the criterion
#     a *disjunction*. RTL chemistry only generates a library molecule when BOTH
#     halves hybridise and the nick between them is sealed by the ligase. One
#     half carrying >=5 mismatches is enough to abolish the product, so such a
#     pair emits no false signal and must NOT be rejected. See
#     ``off_target_require_both_halves`` in specificity.evaluate_specificity.
#  2. "five mismatches" is counted over the whole 25 nt probe half, not over
#     whatever sub-window a local aligner happens to report. Counting over a
#     20 nt local alignment instead (the legacy OFF_TARGET_ALN_MIN form) rejects a
#     25-mer whose 6 mismatches are clustered at its ends but which still
#     contains a clean 20 nt core. See ``off_target_full_half`` below.
#
# Number of mismatches over the FULL probe half at which an off-target match
# stops being able to hybridise efficiently. A half is disqualifying at
# <= OFF_TARGET_MIN_MISMATCH - 1 (i.e. <= 4) mismatches.
OFF_TARGET_MIN_MISMATCH = 5

# Default ON: count mismatches over the full 25 nt half (rule as written).
# Set False to restore the legacy local-alignment form below.
OFF_TARGET_FULL_HALF_DEFAULT = True

# Default ON: a pair is rejected only when BOTH halves fail the off-target
# gate (rule as written). Set False to restore the legacy "either half
# disqualifies" behaviour.
OFF_TARGET_REQUIRE_BOTH_HALVES_DEFAULT = True

# Default ON: a pair is rejected only when both halves have a disqualifying
# match to the SAME off-target gene.
#
# 10x writes "matches to off-target GENES", which scopes the disjunction per
# off-target gene rather than over the probe as a whole: an LHS bound to gene X
# and an RHS bound to gene Y are never adjacent on one molecule, so the ligase
# has nothing to seal and neither gene yields a product. Rejecting such a pair
# discards a probe that cannot generate false signal by the RTL mechanism.
#
# The default rests on evidence, not on the wording alone. On 10x's public
# Flex v1 human panel (all 53,444 `included=TRUE` and all 1,136
# `included=FALSE` probes, 54,580 in total), screened exhaustively against
# GENCODE v32:
#
#                        shipped rejected   excluded rejected (n=1,136)   sep
#   OFF (per-probe)           17.6%                82.9%                +65.3
#   ON  (per-gene)             6.8%                76.6%                +69.7
#
# i.e. per-gene scoping retains 5,774 additional probes that 10x itself ships
# (10.8% of the panel) at a cost of 72 of the 1,136 excluded probes (6.3%).
# Separation, which is how well the screen reproduces 10x's judgement, rises
# from +65.3 to +69.7 points.
#
# THE EXCLUDED SET IS NOT ONE PILE. 623 of the 1,136 carry a `DEPRECATED_`
# gene-id prefix, meaning the gene was withdrawn from the annotation rather
# than that the probe failed anything. The prefix is stripped before the
# lookup, so all 623 resolve, but they remain a different population. Over the
# other 513 the per-probe rule rejects 90.8% and the per-gene rule 88.9%, so
# per-gene costs 10 of those 513 (1.9%) and separation rises +73.2 -> +82.0.
# Both denominators are defensible.
#
# Do not narrow this to per-LOCUS scoping: requiring both halves to be close to
# the same off-target SITE (not merely the same gene) rejects 6.4% of shipped
# and 75.0% of excluded probes, a separation of +68.6. That is 1.1 points BELOW
# per-gene on the full panel, and 0.9 points below it on the 513-probe subset.
# Narrowing from gene to site costs agreement and gains none.
OFF_TARGET_SAME_SUBJECT_DEFAULT = True

# ---------------------------------------------------------------------------
# Junction-spanning clean duplex — OURS, NOT 10x's
# ---------------------------------------------------------------------------
# Everything above this block is 10x's rule, applied as written. This block is
# not: 10x publishes no criterion of this kind, and nothing here can be cited to
# CG000621 or CG000839. flexiprobes is "10x's rules, made slightly stricter in
# places"; this is one of the places, and it is labelled as ours so a reader
# can tell which rules came from the vendor and which did not.
#
# WHAT IT ADDS. 10x's sentence is a pure COUNT — "at least five mismatches in at
# least one of the LHS or RHS probes" — and says nothing about WHERE the
# mismatches sit. RTL only writes a library molecule when both halves hybridise
# AND the ligase seals the nick between them, and the ligase reads the bases
# either side of the nick, not the ends of the oligos. A pair whose mismatches
# against some other gene are all terminal therefore still presents a long,
# perfectly matched duplex ACROSS the junction: it hybridises, it ligates, and
# it puts reads on the wrong gene, while satisfying the count rule comfortably.
# This gate measures that duplex directly — the unbroken run of matched bases
# immediately 5' and immediately 3' of the nick, in CONTIGUOUS reference
# sequence — and rejects when both sides are long enough to seal.
#
# It can only ever reject MORE than the count rule, never less, so it cannot
# contradict the vendor spec. A probe it rejects is one 10x's own sentence would
# have shipped; ``specificity_status`` says ``failed_junction_duplex`` rather
# than ``failed_off_target`` so the two populations stay separable in results.csv.
#
# WHY 9, MEASURED, NOT CHOSEN. 10x's public Human Chromium Transcriptome Probe
# Set v1.1.0 against Ensembl GRCh38 cdna.all, primary assembly only,
# readthrough/ALT/PATCH artefacts excluded. Requiring >= S clean bp on BOTH
# sides of the nick against a non-target gene:
#
#   S      catches 10x's OWN rejects      fires on what 10x SHIPS   separation
#   7            97.0%                          53.4%                 poor
#   9            87.0%                           0.6%                0.864
#   11           80.5%                           0.2%                0.803
#
# 9 maximises the separation between the two populations. 7 is useless — it
# fires on half of everything 10x ships. 11 and above trade real sensitivity for
# almost nothing. The value is the argmax of a measured curve on 10x's public
# panel, not a judgement call. It is still a FlexiProbes heuristic, not a 10x
# specification.
#
# PER SIDE, NOT TOTAL. 9+9 is rejected and 15+3 is not, even though the second
# has the longer total run. The ligase needs a stable duplex on EACH side of the
# nick; three matched bases on one side do not hold it there. A total-length
# form would have scored those two identically and would have been measuring
# homology rather than ligatability.
#
# 0 disables the gate. Any other value below
# OFF_TARGET_JUNCTION_DUPLEX_MIN_BP is REFUSED rather than silently clamped —
# see that constant.
OFF_TARGET_JUNCTION_DUPLEX_BP = 9

# Floor below which the gate cannot be honestly applied, and therefore refused.
#
# The scan seeds on an EXACT k-mer centred on the nick (screen.KmerScreen.K = 12,
# so 6 bases each side). That seeding is complete — every site with >= 6 clean bp
# per side is enumerated, and a site with fewer is not a site this metric can
# see. Asking for a threshold below 6 would therefore run a gate at 6 while
# reporting it as the requested number, which is the kind of silent
# disagreement between a knob and its effect that this codebase refuses
# elsewhere (see the recount basis). It is also below the noise floor: a
# specific L-mer occurs ~403e6 / 4**L times by chance in a transcriptome, so at
# L=12 (6+6) that is ~24 hits per probe from chance alone. Anything shorter
# measures the size of the reference, not the probe.
OFF_TARGET_JUNCTION_DUPLEX_MIN_BP = 6

# Ceiling on the number of reference sites the junction scan will measure for
# one probe, per reference. A specific 12-mer occurs ~403e6 / 4**12 = 24 times
# by chance in a transcriptome, so 2,000 is ~80x what chance can produce: a
# probe that reaches it does not have 2,000 homologues, it has a junction that
# sits inside a repeat element.
#
# On reaching it the scan reports the probe as UNMEASURED (-1 in the two
# clean-bp columns) and warns, rather than measuring the first 2,000 sites and
# reporting that as the worst. Reference position order is not danger order, so
# a truncated maximum is not a maximum, and a number that silently means
# "the worst of an arbitrary subset" is worse than an admitted gap. Such a
# probe keeps whatever verdict 10x's count rule gives it.
JUNCTION_DUPLEX_MAX_SITES = 2000

# Legacy local-alignment form, retained so prior results can be reproduced
# (``off_target_full_half=False``): >=5 mismatches over a 25 bp half was
# translated as <80% identity over a >=20 bp local alignment.
OFF_TARGET_PIDENT_MIN = 80.0
OFF_TARGET_ALN_MIN = 20
ON_TARGET_PIDENT_MIN = 90.0
ON_TARGET_ALN_MIN = 21

# Slack (bp) allowed when pushing the full-half off-target floor down into
# BLAST itself. A gapped HSP spends query positions on gap columns, so its
# reported identity/coverage understate how much of the 25-mer it really
# covers; without this margin BLAST would drop HSPs the Python classifier
# would have counted, making the screen silently *less* sensitive. Mirrors the
# same 3 bp margin used by _derive_qcov_hsp_perc.
OFF_TARGET_PUSHDOWN_GAP_MARGIN = 3

# The scoring this screen asks blastn for, and the mismatch cost that follows
# from it. The cost is DERIVED, not typed, because a scoring parameter living
# in two places is how the two copies drift apart.
# blast.run_blastn passes these, pipeline._off_target_pushdown consumes the
# cost, and the tests import all three rather than restating them.
#
# blastn-short's own defaults are reward +1 / penalty -3, giving cost 4. This
# screen overrides them: at -3 a 25-mer four mismatches from a subject scores
# 21 - 12 = 9 at full length, so a short exact core outscores the full
# alignment and BLAST reports the core, which cannot answer a criterion
# defined across the whole half.
BLASTN_REWARD = 1
BLASTN_PENALTY = -1
#
# This is the constant that makes the off-target length pushdown correct.
# BLAST reports the maximum-scoring LOCAL segment, so a half sitting m
# mismatches from a subject is reported with length >= half_len - cost*m: the
# full-length alignment scores half_len - cost*m, and any segment's length is
# at least its score. Using the identity-derived `keep` (21/25) as a length
# floor instead would assume the HSP spans the half, and would discard exactly
# the trimmed cores the full-half recount exists to catch. See
# pipeline._off_target_pushdown and specificity.full_half_mismatches.
#
# Direction of danger, so a future edit knows which way is safe: a cost that is
# too SMALL lowers the length floor, keeps hits that need not be kept, and is
# merely slower. A cost that is too LARGE raises the floor and can discard a
# disqualifying hit, which is a silent false negative. If in doubt, go low.
#
# Checked by brute force over all 15,276 arrangements of up to 4 mismatches in
# a 25-mer, scoring every local segment and taking the maximum: at penalty -1
# the shortest such segment is exactly 17, which equals the floor this yields,
# and 0 arrangements are dropped. The same brute force with penalty -3 against
# a floor of 17 drops 3,899, so the check does distinguish the two and is not
# merely agreeing with itself.
BLASTN_SHORT_MISMATCH_COST = BLASTN_REWARD - BLASTN_PENALTY


# ---------------------------------------------------------------------------
# Off-target recount: the only screen
# ---------------------------------------------------------------------------
# There is deliberately no cheaper "bounded" mode that keeps a BLAST length
# floor of 21 instead of 9. Such a floor is not merely blind in the 3-5
# mismatch band. It FALSELY REJECTS, and by a lot: checked against exhaustive
# <=4-mismatch ground truth (a pigeonhole q-gram index, with no heuristic and
# nothing reconstructed), it discards a large share of usable probes.
#
# The mechanism is the floor itself: at 21 the pushdown demands an HSP covering
# >=72% of the half, so where a subject's best local alignment is a short exact
# core, BLAST reports a stretched and shifted alignment instead. The recount then
# rebuilds the flanks at the wrong offset and the hit scores far closer than it
# is. There is no configuration that makes that sound, so the option does not
# exist: a knob whose every setting but one is known to be wrong is a trap, not
# a choice.
#
# The correct screen costs far more hits, memory and time than a floor of 21
# (orders of magnitude more hits). It is not tunable: a 9 nt exact match is
# common by chance, and demanding 100% identity still leaves most of those hits.
# The price is paid for the bases BLAST never looked at, and nothing else buys
# them.

# Warn before a run commits to a large screen. 20,000 halves is ~10,000 probes —
# far above any viral panel, far below a whole-transcriptome chunk.
RECOUNT_EXACT_WARN_HALVES = 20000


# ---------------------------------------------------------------------------
# Filter funnel defaults
# ---------------------------------------------------------------------------
# Source: CG000839 Rev A p.2 (and CG000621 Rev D p.2 / p.12, identical) —
# "GC content should be between 44 - 72% for each 25 bp probe half."
#
# This band is NOT a hard gate by default; it is available as an opt-in filter
# (``FilterThresholds.enforce_gc_range``) and otherwise feeds ranking only.
# Reason: the band is a proxy for duplex stability written for mammalian
# transcriptomes (10x probe halves have median GC 0.52), and 10x's own panels
# do not obey it. Measured over the 53,459 included probes of 10x's public Chromium Human
# Transcriptome Probe Set v1.1.0 (106,918 halves), with the primer3 conditions
# below:
#   • 8,113 halves (7.59%) sit BELOW GC 0.44 and 827 above 0.72 — 8,940 (8.36%)
#     violate 10x's own published band, down to GC 0.28;
#   • yet NOT ONE half of the panel falls below 56.32 C, and the histogram is
#     empty in every bin under it. Mouse v1.1.1 gives 56.35 C.
# A hard GC band therefore discards windows 10x itself would ship, and the
# constraint 10x actually respects is thermodynamic. On AT-rich targets (for
# example bunyavirus segments at GC 0.31-0.42) the band discards most windows,
# far more than the Tm floor below does.
GC_MIN_DEFAULT = 0.44
GC_MAX_DEFAULT = 0.72

# Per-half Tm floor (deg C) — the hard duplex-stability gate, replacing the GC
# band. Derived empirically (not published by 10x): the observed minimum
# per-half Tm across 10x's shipped panels is 56.32 C (human v1.1.0 / v2.0.0 /
# v1.0) and 56.35 C (mouse v1.1.1), so 56.0 admits 100% of 10x's own probes
# while still excluding thermodynamically weak windows.
#
# *** THIS Tm IS RELATIVE, NOT AN IN-ASSAY MELTING TEMPERATURE. *** It is
# primer3 calc_tm under the PCR-primer salt conditions in TM_MV_CONC..
# TM_DNA_CONC below, because 10x does not publish the Flex hybridisation
# buffer. The floor is calibrated in exactly the same units as the number it is
# compared against, so it is a valid *relative* screen; it must not be read as
# "these probes melt at 56 C in the Flex assay". Anyone who knows their buffer
# should adjust the salt constants and re-derive this floor with them.
TM_MIN_DEFAULT = 56.0

# No upper Tm gate by default: 10x's panels run to 78.7 C per half with no
# boundary visible in the distribution, so there is nothing to calibrate one
# against. Set FilterThresholds.tm_max to impose one.
TM_MAX_DEFAULT = None

# "Avoid homopolymer repeats" — no explicit run length given; 4 is a
# conservative reject-if-longer threshold.
MAX_HOMOPOLYMER_RUN_DEFAULT = 4


# ---------------------------------------------------------------------------
# TN ligation junction
# ---------------------------------------------------------------------------
# The mRNA base opposite the LHS probe's 3'-most nucleotide (mRNA position 26
# of the 50 bp window, index 25). The spec is chemistry-dependent:
#
#   • CG000621 Rev D, Visium section (p.3) AND Chromium Flex section (p.12),
#     identical wording — "The 25th nucleotide of the probe (3' most nucleotide
#     of the LHS probe) must be a T. The opposing nucleotide in the target RNA
#     must be an A."
#   • CG000839 Rev A (GEM-X Flex v2, p.2) downgrades it to a recommendation —
#     "It is recommended to have a TN ligation junction, where T is the 25th
#     nucleotide of the probe (3' most nucleotide of the LHS probe). In this
#     case, the opposing nucleotide in the target RNA would be an A." and then,
#     explicitly: "While preferred, the TN motif is not always strictly
#     required, and other nucleotide motifs can also function effectively.
#     However, for a simplified and robust design process, 10x Genomics
#     recommends using the T at the specified position."
#
# Empirically 10x does not obey its own "must" even for Flex v1. mRNA-side base
# at index 25, complement of probe_seq[24], over the shipped panels:
#
#   panel                    n        A       C       G       T     not-C
#   human v1.1.0 (Flex v1)   54,580   51.30%   1.00%  36.04%  11.67%  99.00%
#   human v2.0.0 (Flex v2)   54,018   50.96%   1.01%  36.25%  11.78%  98.99%
#   mouse v1.1.1 (Flex v1)   55,538   52.16%   0.71%  34.91%  12.22%  99.29%
#   rat   v2.0.0 (Flex v2)   62,506   50.34%   1.25%  35.71%  12.69%  98.75%
#
# "Must be A" holds ~51% of the time in 10x's own Flex v1 panel; the invariant
# they actually respect is "not C" (~99%).
#
# That finding does NOT become the default. The junction gate, when engaged, is
# the strict adenine for every chemistry — CG000839's own conclusion is "for a
# simplified and robust design process, 10x Genomics recommends using the T at
# the specified position", and it is the conservative choice. The relaxation is
# available as an explicit opt-in (``soft``) for the Flex chemistries only,
# because CG000839's "not always strictly required" language is scoped to
# GEM-X Flex v2 and CG000621's Visium section states the requirement flatly
# with no such qualifier.
TN_MODE_STRICT_A = "strict_a"   # HARD (default): reject unless mRNA[25] == 'A'
TN_MODE_SOFT = "soft"           # opt-in, Flex only: reject only mRNA[25] == 'C'
TN_MODE_OFF = "off"             # no junction-identity gate at all
TN_MODE_AUTO = "auto"           # resolve from the run's chemistry + legacy flag
TN_MODES = (TN_MODE_AUTO, TN_MODE_STRICT_A, TN_MODE_SOFT, TN_MODE_OFF)


# ---------------------------------------------------------------------------
# Tm calculation conditions (primer3 calc_tm)
# ---------------------------------------------------------------------------
# primer3's defaults are PCR-primer conditions, NOT the 10x Flex
# hybridisation/ligation buffer (which 10x doesn't publish). Reported Tm is
# therefore a *relative* indicator for ranking/uniformity, not the in-assay
# melting temperature. Exposed here so a user who knows their buffer can adjust
# them in one place. Values below == primer3 defaults.
TM_MV_CONC = 50.0      # monovalent cations (mM)
TM_DV_CONC = 1.5       # divalent cations (mM)
TM_DNTP_CONC = 0.6     # dNTP (mM)
TM_DNA_CONC = 50.0     # oligo/DNA (nM)

# ---------------------------------------------------------------------------
# Multi-strain ΔTm heuristics (flexprobe.consensus)
# ---------------------------------------------------------------------------
# Comparative, not absolute. Source rationale in consensus.py docstring.
DTM_PER_MISMATCH = -2.0          # baseline penalty per internal mismatch
# Retained as a named constant (re-exported by consensus.py, and referred to by
# name in specificity.py's explanation of the OFF-target reading). The
# ON-target ΔTm does not use it; it uses the graded ramp below instead, which
# reproduces this step's endpoint exactly (see MISMATCH_WEIGHT_AT_NICK).
DTM_JUNCTION_BONUS = -4.0        # extra penalty within JUNCTION_WINDOW
JUNCTION_WINDOW = 4              # bp each side of the ligation junction (CG000621)

# ---------------------------------------------------------------------------
# Graded positional weighting of an ON-TARGET mismatch.
# ---------------------------------------------------------------------------
#
# INTERNAL heuristic, NOT 10x-sourced — exactly like the primer3 nudge in
# ranking.py. CG000621 Rev D says only that mismatches "should be at least four
# bp away from the ligation junction"; the *shape* of the penalty either side of
# that line is ours, and is a rule of thumb, not a thermodynamic model.
#
# WHY A RAMP AND NOT A STEP.
# A step that charges ``DTM_JUNCTION_BONUS`` to every mismatch inside
# JUNCTION_WINDOW and to none outside it makes, on a 50-mer, a mismatch at
# position 28 cost -6.0 C and one at position 29 cost -2.0 C (a 3x cliff for one
# base of movement at an arbitrary boundary), and then FLAT -2.0 C across the
# remaining 21 bases of the half. Nothing then prefers a mismatch at the outer
# end of a half, where the duplex frays anyway and the ligase never sees it,
# over one sitting in the middle of that half. Duplex destabilisation is not a step
# function of position; it falls off smoothly with distance from the nick, so a
# ramp is both the more honest shape and the one that can actually rank the
# in-between cases the step could not distinguish.
#
# THE ENDPOINTS MATCH THE STEP.
# The weights below multiply DTM_PER_MISMATCH (-2.0), so a mismatch at the nick
# costs 3.0 * -2.0 = -6.0 C, precisely
# DTM_PER_MISMATCH + DTM_JUNCTION_BONUS, and one at the probe end costs
# 1.0 * -2.0 = -2.0 C, precisely the flat cost. Only the interior differs.
# Nothing here changes any mismatch COUNT; counts are reported unweighted.
MISMATCH_WEIGHT_AT_NICK = 3.0    # d = 0   -> -6.0 C (== per-mismatch + junction)
MISMATCH_WEIGHT_AT_END = 1.0     # d = max -> -2.0 C (== per-mismatch alone)


def mismatch_distance_from_nick(position: int, half_length: int) -> int:
    """Distance in bases from the ligation nick to ``position``.

    ``position`` is a 0-based offset into the target window; ``half_length`` is
    where the LHS half ends and the RHS half begins. The nick sits *between*
    ``half_length - 1`` and ``half_length``, so both of those bases are at
    distance 0 and the two outermost bases of the window are at distance
    ``half_length - 1`` (24 on the standard 50-mer).

    Out-of-window positions are clamped rather than rejected: the callers get
    positions from alignment arithmetic and a defensive clamp is cheaper than a
    raise they would only have to catch.
    """
    d_max = half_length - 1
    if d_max <= 0:
        return 0
    d = (half_length - 1 - position) if position < half_length \
        else (position - half_length)
    return min(max(d, 0), d_max)


def mismatch_position_weight(position: int, half_length: int) -> float:
    """Multiplier on ``DTM_PER_MISMATCH`` for a mismatch at ``position``.

    Linear from ``MISMATCH_WEIGHT_AT_NICK`` at distance 0 down to
    ``MISMATCH_WEIGHT_AT_END`` at the probe ends::

        w(d) = 3.0 + (1.0 - 3.0) * d / (half_length - 1)

    On the standard 50-mer (half_length 25, d_max 24) that is
    ``w(d) = 3.0 - d / 12``: -6.0 C at the nick, -4.0 C halfway out, -2.0 C at
    either end. See the block comment above for why the two endpoints are the
    values they are.

    A degenerate ``half_length <= 1`` (no interior to grade) returns the
    at-nick weight, the conservative end of the range.
    """
    d_max = half_length - 1
    if d_max <= 0:
        return MISMATCH_WEIGHT_AT_NICK
    d = mismatch_distance_from_nick(position, half_length)
    span = MISMATCH_WEIGHT_AT_END - MISMATCH_WEIGHT_AT_NICK
    return MISMATCH_WEIGHT_AT_NICK + span * (d / d_max)


# ---------------------------------------------------------------------------
# Degenerate (consensus) probes — parent / variant model
# ---------------------------------------------------------------------------
# A pan-strain window that is not identical across every input record carries
# IUPAC ambiguity at the positions where the records disagree. Such a window is
# emitted as a "consensus" candidate that is NEVER orderable, plus one
# "variant" candidate per input record — each a concrete ACGT resolution of the
# same window against that record's own sequence.
#
# Resolutions are taken PER RECORD, not as the cartesian product of the IUPAC
# codes: the product would generate sequences no input organism carries, while
# the per-record resolution is exactly the oligo that binds that strain or
# isoform perfectly. Two records resolving identically collapse to one variant
# that lists both.
PROBE_ROLE_PLAIN = ""            # single-sequence run, or an exact-match window
PROBE_ROLE_CONSENSUS = "consensus"   # degenerate parent — never exported
PROBE_ROLE_VARIANT = "variant"       # concrete ACGT child — the orderable one

# Ambiguity must sit clear of the ligation junction. CG000621 Rev D asks for
# SNPs and mismatches at least four bp away from it, and the ligase simply will
# not seal across a junction-proximal mispair — so a degenerate base inside this
# window would make the probe's behaviour on the non-reference record
# unpredictable in the one place the chemistry cannot tolerate it. Reuses
# JUNCTION_WINDOW so the tolerance rule and the per-strain ΔTm penalty above
# describe the same region of the probe.
DEGENERATE_JUNCTION_CLEARANCE = JUNCTION_WINDOW

# Every non-ACGT base counts against the per-half degenerate budget, not just
# literal 'N'. The consensus builder emits M/R/W/S/Y/K/V/H/D/B wherever the
# disagreeing bases are all in ACGT and only falls back to 'N' for a gap, so a
# budget that counted 'N' alone would let a two-way ambiguity through unmeasured
# — and it is the two-way codes that dominate a real strain alignment.
DEGENERATE_BASES = frozenset("NMRWSYKVHDB")


# ---------------------------------------------------------------------------
# Strain-unique discrimination (flexprobe.consensus)
# ---------------------------------------------------------------------------
# The strain-unique gate asks the mirror image of the off-target question:
# instead of "can this probe still hybridise + ligate on a gene it should not
# hit?", it asks "will this probe FAIL to hybridise + ligate on the strains it
# should not hit?". Both are answered by the same physics, so both must be
# counted the same way — over a probe HALF, not over the whole 50-mer.
#
# ``RunConfig.strain_unique_min_mismatches`` defaults to 5 precisely because it
# is taken from the 10x off-target rule (CG000621 Rev D / CG000839 Rev A), whose
# full sentence scopes that number PER HALF:
#
#   "Matches to off-target genes should have at least five mismatches in at
#    least one of the LHS or RHS probes to prevent efficient hybridization."
#
# Counting the 5 over the whole 50 nt window instead accepts windows that the
# engine's OWN off-target classifier says will still hybridise and ligate. A
# window with 3 LHS + 2 RHS mismatches against the other strain totals 5, but
# neither half reaches 5, so ``evaluate_specificity`` would (correctly) call
# that pair a live off-target hit — while the strain-unique gate simultaneously
# called it strain-discriminating. The two verdicts cannot both be right, and
# the per-half one is the one backed by the source and by RTL chemistry.
#
# Default ON (the consistent, spec-faithful reading). Set False for the legacy
# whole-window count, which reproduces panels designed under that reading.
# Same escape-hatch pattern as OFF_TARGET_FULL_HALF_DEFAULT above.
STRAIN_UNIQUE_PER_HALF_DEFAULT = True

# per_half fixes WHICH half a strain's mismatches must fall in; it says
# nothing about WHERE in that half. A raw >= 5 count treats a mismatch at the
# nick identically to one at the probe's outer end, but ligation only cares
# about the nick — 5 mismatches bunched at the far end of a half leave the
# nick-spanning region a perfectly matched duplex, and the off-target engine's
# OWN ``junction_clean_duplex`` says a duplex that long across the nick still
# seals and ligates (see OFF_TARGET_JUNCTION_DUPLEX_BP above). So a window can
# satisfy the per-half count while the strain it is supposed to exclude would
# still make a library molecule from it — over-certification, not a rounding
# error.
#
# ``strain_unique_start_positions`` (consensus.py) closes this by reusing that
# SAME off-target scan, unmodified, pointed at the other strain instead of a
# non-target gene: a window only certifies discrimination against a given
# strain when it both reaches ``strain_unique_min_mismatches`` in a half AND
# that strain's window does not present a sealing duplex around the nick. It
# can only remove windows the count rule already granted, so it never
# contradicts 10x's published sentence — exactly the relationship the
# off-target gate already has to the same sentence. Reusing the function and
# the threshold rather than inventing a second one keeps "does this half
# actually ligate here" answered by one piece of arithmetic on both sides of
# the strain/off-target mirror, instead of two that could drift apart.
#
# The over-certification this prevents is reproducible on synthetic input.
# Threshold lives in ``RunConfig.strain_unique_junction_duplex_bp``
# (default OFF_TARGET_JUNCTION_DUPLEX_BP, "0 disables", the same convention as
# its off-target twin). Only meaningful when STRAIN_UNIQUE_PER_HALF_DEFAULT (or
# the per-run override) is True; the legacy whole-window reading is untouched
# so panels designed under it keep reproducing bit-for-bit.


# ---------------------------------------------------------------------------
# Probe geometry
# ---------------------------------------------------------------------------
TILE_LENGTH_DEFAULT = 50         # 25 bp LHS half + 25 bp RHS half
# Group input transcripts by their Ensembl 'gene:' header field (isoform-aware
# specificity). ON by default. A no-op, with a log line, on input that carries
# no 'gene:' fields: see pipeline._ingest.
GROUP_BY_GENE_DEFAULT = True
HALF_LENGTH_DEFAULT = 25

# Spacing (bp) between two chosen probes on the same target. Auto-pick and the
# review-dialog overlap check enforce the minimum so selected probes stay spread
# out along the target rather than clustering on a single hot spot; selection
# prefers the larger gap when it can still hit its coverage / probe-count goal.
#
# TARGET, not a hard floor. Selection walks a three-rung ladder (see
# ``auto_select.SPACING_LADDER``) and stops at the first rung that yields
# ``PROBE_GAP_RELAX_MIN_PROBES`` probes, so a gene only ever pays for the
# tightness it actually needs. Each rung pairs a spacing floor with the terminal
# fraction of the transcript it will avoid, because both have to give at once on
# a gene that is simply too short to satisfy either:
#
#   tier    gap        terminal avoided   meaning
#   GREEN   >= 150 bp   15% each end      as designed
#   ORANGE   40-150 bp  15% each end      tight; flagged
#   RED      15-40 bp   10% each end      last resort; flagged prominently
#
# Nothing below MIN is ever silently accepted — every probe carries its tier
# into the reports, the exports and the review map.
#
# 150 comes from 10x's own predesigned panels rather than from a tech note —
# they publish no spacing rule. Mapping their panels into TRANSCRIPT coordinates
# (genomic distance is meaningless: two probes either side of an intron are
# adjacent on the mRNA) gives the median gap between consecutive probes on a
# gene:
#
#   human v1.0.1 vs GENCODE v32   52,170 probes placed   median 151 bp
#   mouse v1.1.1 vs GENCODE vM33  54,652 probes placed   median 175 bp
#
# Two independent public panels against two independently matched annotations,
# so 150 is a replicated figure and not a quirk of one species. It is still a
# FlexiProbes heuristic, not a 10x specification.
#
# 10x's own CLOSEST pair per gene is much tighter than that median (human median
# 71 bp, mouse 84; 70.9% / 66.2% of genes carry a pair under 150 bp), so 150 as
# a hard floor would be *stricter* than the commercial panels. That is why it is
# the top of a ladder rather than a floor: the lower rungs keep the probe count,
# the tier keeps the tightness visible.
MIN_PROBE_GAP_BP = 150

# ORANGE tier floor — the second rung. Below this two 50-mers are close enough
# to compete for the same locally-melted stretch of transcript, and the pair
# starts to stop behaving like two independent measurements.
#
# 40 bp is 10x's own p20: 20.2% of consecutive gaps in the human panel and 17.9%
# in the mouse one fall below it. So it is the point past which their panels stop
# spacing routinely and start spacing exceptionally — the right place for the
# boundary between "tight but ordinary" and "worth looking at".
RELAXED_PROBE_GAP_BP = 40

# RED tier floor — the third and last rung, and the absolute floor: selection
# never emits a pair closer than this, it emits fewer probes instead. 15 bp is
# roughly one and a half helical turns of B-form duplex between two footprints:
# enough that the two 25 nt halves flanking each nick are not sterically
# competing for the same bases, and not much more.
#
# Stated honestly, this is OUR floor and not 10x's: 8.5% of consecutive gaps in
# their human panel and 7.4% in the mouse one are tighter than 15 bp, so they do
# go further. We stop here because below it the two probes stop being separable
# evidence, and because a gene needing that rung is a gene whose transcript
# cannot hold three separated probes at all — at which point the honest options
# are "three crowded probes, loudly flagged" or "fewer probes", which is exactly
# what the GUI's spacing filter lets a user choose between.
MIN_PROBE_GAP_FLOOR_BP = 15

# Relaxation trigger: a gene that cannot yield this many probes on one rung
# drops to the next. Three is the mode of every 10x predesigned panel
# (15,309 / 19,074 human v1.0.1 genes carry exactly 3 probes) and the point
# below which per-gene signal stops averaging usefully.
PROBE_GAP_RELAX_MIN_PROBES = 3

# --- Spacing tiers -------------------------------------------------------- #
# One vocabulary for the whole app: the selector records a tier per gene, the
# exports carry it per probe, the review map colours bars by it and the GUI
# filter is expressed in it. Ordered so a HIGHER number is better spacing, which
# lets callers write ``tier >= GAP_TIER_ORANGE`` for "at least this good".
GAP_TIER_RED = 0
GAP_TIER_ORANGE = 1
GAP_TIER_GREEN = 2
GAP_TIER_NAMES = {GAP_TIER_RED: "red", GAP_TIER_ORANGE: "orange",
                  GAP_TIER_GREEN: "green"}
# Floor gap for each tier, best first — the ladder's gap column.
GAP_TIER_FLOORS = {GAP_TIER_GREEN: MIN_PROBE_GAP_BP,
                   GAP_TIER_ORANGE: RELAXED_PROBE_GAP_BP,
                   GAP_TIER_RED: MIN_PROBE_GAP_FLOOR_BP}


def gap_tier(gap: "int | float | None") -> int:
    """Which spacing tier a gap in bp falls in.

    ``None`` (a gene with a single probe, so no pair and no spacing) counts as
    GREEN: there is nothing crowded about one probe, and reporting it as a
    spacing problem would bury the genes that have one. A gap below
    ``MIN_PROBE_GAP_FLOOR_BP`` is still RED rather than a fourth tier — the
    selector cannot produce one, but a hand-picked selection can, and the point
    of the colour is "too close", which it is.
    """
    if gap is None:
        return GAP_TIER_GREEN
    if gap >= MIN_PROBE_GAP_BP:
        return GAP_TIER_GREEN
    if gap >= RELAXED_PROBE_GAP_BP:
        return GAP_TIER_ORANGE
    return GAP_TIER_RED

# Soft TIE-BREAK for the coverage-maximising selector only. NOT a rung of the
# ladder, NOT a maximum, and NOT larger-than-150-so-it-wins: 150
# (``MIN_PROBE_GAP_BP``) is the design target, and this never overrides it.
#
# Used in exactly one place — ``isoform_select._greedy``, as the middle element
# of the lexicographic key ``(new isoforms covered, meets_pref, score)``. It
# therefore only separates candidates that already add EQUAL isoform coverage,
# and can never cost a gene a probe. The default top-N selector
# (``auto_select.auto_pick_per_gene``) does not take a preferred gap at all, so
# on a normal design run this value never enters.
#
# In particular this is NOT the floor of the adaptive spread target.
# ``auto_select.spread_gap_for`` floors at ``min_gap`` (i.e. MIN_PROBE_GAP_BP)
# and deliberately does NOT floor at PREFERRED_PROBE_GAP_BP: on a 1 kb
# transcript three 50-mers at 300 bp gaps need 750 nt, 75% of it, so a 300 bp
# floor would FORCE probes into the terminal 15% that the placement score is
# simultaneously pushing them out of.
PREFERRED_PROBE_GAP_BP = 300


# ---------------------------------------------------------------------------
# Placement along the transcript (ranking.central_placement)
# ---------------------------------------------------------------------------
# Fraction of the transcript at EACH end that the placement term discourages.
# Inside the remaining middle 70% the term is flat, so ranking expresses no
# preference between one central position and another and the *spacing* target
# is what spreads probes out. A linear tent peaking at the exact midpoint would
# instead pull every top-ranked window on a long transcript into a narrow band
# around the centre.
#
# 0.15 is where 10x themselves sit. On the length-matched lncRNA subset of 10x's
# public human v1.0.1 panel (no CDS, so the whole transcript is fair game and the
# design's own positional preference is visible undistorted) 14.6% of probes
# fall in the first 15% and 7.3% in the last 15% — under the 30% a uniform
# design would put there, but far from zero. Protein-coding transcripts of the
# same length instead show 20.5% / 1.6%, which is 3'UTR avoidance rather than a
# positional rule; mouse protein-coding transcripts give 19.2% / 1.5% and an
# identical median position (0.297), so this replicates across species.
END_EXCLUSION_FRAC = 0.15

# Terminal fraction avoided on the RED rung. A gene only reaches that rung
# because its transcript cannot hold three separated probes, and on such a
# transcript the middle 70% is often simply too small to contain them — so the
# terminal band has to give at the same time the spacing does, or the rung buys
# nothing. Loosening to 10% while the gap drops to 15 bp keeps the two knobs
# turning together instead of one silently blocking the other.
#
# 0.10 is where 10x's own panels sit at their loosest: 14.6% of human and 11.1%
# of mouse probes fall in the terminal 10%, against 20% for a uniform design. So
# even the last rung stays inside what the commercial panels do rather than
# inventing new territory.
END_EXCLUSION_FRAC_RELAXED = 0.10

# ---------------------------------------------------------------------------
# Split prepare / screen execution
# ---------------------------------------------------------------------------
# Marker a --prepare-only run drops inside its output directory once ingest,
# tiling and filtering have finished, the candidate checkpoint is on disk and
# the pre-screen artefacts are written. It certifies exactly one thing: the
# memory-hungry front half of this run completed cleanly.
#
# Deliberately NOT `.flexiprobes_DONE`. That sentinel means the run is
# FINISHED — results.csv exists, the tarball is packaged — and every remaining
# chain link no-ops on sight of it. A prepared run has produced no probes at
# all, so writing DONE would strand the output directory looking complete and
# stop the screen job from ever starting.
#
# Written inside output_dir, not the bundle root, because it describes the
# state of THAT directory's checkpoint, travels with it in the results
# tarball, and stays correct if the bundle is moved. Dot-prefixed for the same
# reason `.candidate_cache/` is: `tar -C output .` picks it up while the GUI's
# results view never shows it.
PREPARED_MARKER = ".flexiprobes_PREPARED"
