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

"""Ranking engine: transparent composite scoring for probe candidates.

``composite_score`` is contractually in **[0, 1]** — including after the opt-in
primer3 adjustment, which is a bounded nudge (see ``PRIMER3_ADJUSTMENT_MIN`` /
``PRIMER3_ADJUSTMENT_MAX``) and not a free-floating additive term, and after the
TN-junction gate (see ``JUNCTION_BASE_MULTIPLIERS``), which is a multiplier in
[0, 1] and therefore cannot leave the range either.

Which molecule do the primer3 terms describe?
---------------------------------------------
The 50-bp mRNA window is never synthesised. The species in the tube are two
oligos per probe: ``lhs_5prime_handle + target_LHS`` (46 nt) and
``/5Phos/ + target_RHS + rhs_3prime_handle`` (37-63 nt depending on chemistry).
Hairpin / homodimer / heterodimer are therefore computed on the **target-binding
25-mers** ``target_LHS`` / ``target_RHS`` — real sub-molecules of the real
oligos, on the strand actually synthesised — rather than on the mRNA window.

They are deliberately *not* computed on the full assembled oligos, even though
that is the physically complete answer:

* The handles are constant within a chemistry, so their own self-structure is a
  constant offset that cannot re-order candidates. The handle-times-target cross
  structure is real and probe-specific, though: ranking by 25-mer structure
  and by full-oligo structure can disagree substantially. So this is a genuine
  limitation, not a no-op.
* primer3's thermodynamic aligner refuses any sequence longer than 60 nt for
  hairpin/homodimer, and the Flex v1 RHS oligo is 63 nt — it fails for 100% of
  Flex v1 probes, the flagship chemistry. Scoring handles for Flex v2/Visium but
  not Flex v1 would make scores non-comparable across chemistries in a pipeline
  that deliberately emits the same sites for all three.

Including the handles therefore needs a structure engine without the 60-nt
ceiling (ViennaRNA / NUPACK), which this module does not use. Chemistry
context itself is *not* the blocker — the pipeline knows the
chemistry at ranking time via ``config.modality``.
"""
from __future__ import annotations
import importlib.util
import logging
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Set
from .constants import (
    DTM_PER_MISMATCH,
    END_EXCLUSION_FRAC, GC_MAX_DEFAULT, GC_MIN_DEFAULT, OFF_TARGET_PIDENT_MIN,
    PROBE_ROLE_CONSENSUS, PROBE_ROLE_VARIANT,
    TM_DNA_CONC, TM_DNTP_CONC, TM_DV_CONC, TM_MV_CONC,
    mismatch_position_weight,
)
from .models import ProbeCandidate
from .tiling import reverse_complement

logger = logging.getLogger("flexprobe.ranking")


# ---------------------------------------------------------------------------
# primer3 heuristic-nudge constants (INTERNAL heuristics, NOT 10x-sourced).
# ---------------------------------------------------------------------------
#
# These are opt-in tie-breakers applied on top of the sourced composite score.
# They must stay *small* relative to that score — the whole point is that
# primer3 only separates probes the sourced axes rank similarly.
TM_TARGET_C = 60.0            # Tm the bonus peaks at
TM_BONUS_MAX = 0.05           # max Tm bonus
TM_BONUS_PER_C = 0.001        # bonus lost per °C away from TM_TARGET_C
TM_BALANCE_PER_C = 0.005      # penalty per °C of LHS/RHS Tm gap
TM_BALANCE_PENALTY_MAX = 0.05

# primer3 returns ΔG in **cal/mol**, not kcal/mol. Treating the raw value as
# kcal/mol, or scaling it without a cap, can turn a routine -13.6 kcal/mol
# self-dimer into a penalty larger than the entire [0, 1] composite score it is
# meant to nudge. We convert to kcal/mol and cap each penalty at the documented
# maximum.
DG_FLOOR_KCAL = -1.0          # weaker than this -> ignored (no penalty)
DG_SCALE_KCAL = -10.0         # at/below this -> full DG_PENALTY_MAX
DG_PENALTY_MAX = 0.05         # per-term cap (hairpin, homodimer, heterodimer)

# Worst-case total adjustment, used by the invariant tests and by the clamp
# below: +TM_BONUS_MAX at best, -(Tm-balance + 3 ΔG terms) at worst.
PRIMER3_ADJUSTMENT_MAX = TM_BONUS_MAX
PRIMER3_ADJUSTMENT_MIN = -(TM_BALANCE_PENALTY_MAX + 3 * DG_PENALTY_MAX)


def _dg_penalty(dg_cal_per_mol: float) -> float:
    """Map a primer3 ΔG (cal/mol) to a penalty in ``[0, DG_PENALTY_MAX]``.

    Linear ramp from ``DG_FLOOR_KCAL`` (no penalty) to ``DG_SCALE_KCAL`` (full
    penalty), then flat. The ramp is continuous at the floor, so there is no
    jump in penalty at -1 kcal/mol.
    """
    dg = float(dg_cal_per_mol) / 1000.0        # cal/mol -> kcal/mol
    if dg != dg or dg >= DG_FLOOR_KCAL:        # NaN or too weak to care about
        return 0.0
    frac = (dg - DG_FLOOR_KCAL) / (DG_SCALE_KCAL - DG_FLOOR_KCAL)
    return DG_PENALTY_MAX * min(1.0, max(0.0, frac))


@dataclass
class RankingWeights:
    """Weights for composite score components.

    WHERE THE FIVE MAIN WEIGHTS COME FROM. They are a fit against 10x's public
    human Flex probe set, rounded to two decimals. Candidates were generated
    for a sample of genes, the weights were fitted on a small subset of those
    genes, and they were evaluated on the held-out rest. The fit maximises how
    high 10x's shipped probes rank among the candidates of their own gene. It
    cannot validate the weights, because unshipped candidates are untested
    rather than known bad. The same weights transfer to 10x's mouse and rat
    sets without refitting, and label-permutation controls fall to roughly
    chance.

    The human fit is used rather than an average across species because the
    per-species fits disagree. ``specificity_margin`` is kept at 0.28 although
    the mouse and rat fits push it towards 0: 10x's shipped probes all passed
    10x's own off-target screen, so these labels cannot reward a specificity
    margin, and their silence is not evidence against it.
    """
    gc_proximity: float = 0.35
    gc_balance: float = 0.34
    specificity_margin: float = 0.28
    homopolymer_penalty: float = 0.03
    central_placement: float = 0.00
    # Reward probes that bind more of their gene's isoforms (group-by-gene
    # mode). 0 = disabled (the default), so ranking is unchanged on normal runs.
    isoform_coverage: float = 0.0
    # Where a probe's ON-target mismatches sit (multi-strain runs only; see
    # ``_on_target_mismatch_placement_score``). 0.10 against specificity's 0.28
    # (about a third of that axis), so a probe whose mismatches hug
    # the ligation nick is visibly demoted without the term ever being able to
    # overturn a specificity or GC verdict on its own. It is INAPPLICABLE, and
    # then drops out of both numerator and denominator, on any run with no
    # on-target mismatch data at all — which is every single-sequence run.
    mismatch_placement: float = 0.10

    def total(self) -> float:
        return (self.gc_proximity + self.gc_balance +
                self.specificity_margin + self.homopolymer_penalty +
                self.central_placement + self.isoform_coverage +
                self.mismatch_placement)


DEFAULT_WEIGHTS = RankingWeights()

# Weight given to isoform coverage when group-by-gene mode is enabled but the
# user didn't override it. Chosen to meaningfully reward pan-isoform probes
# without swamping GC/specificity (the pipeline applies it on top of the
# defaults via dataclasses.replace, so the other components keep their weights).
DEFAULT_ISOFORM_COVERAGE_WEIGHT = 0.20


# Ideal GC = the EMPIRICAL centre of 10x's own production probe sets, not the
# midpoint of the acceptance band.
#
# The band midpoint ((0.44 + 0.72) / 2 = 0.58) is not a good target: an
# acceptance band says which probes are *allowed*, not which are *typical*.
# Measured over all eight of 10x's public Chromium Flex probe-set CSVs
# (Human/Mouse/Rat, rows with included == TRUE):
#
#   439,656 probes / 879,312 25-mer halves
#   per-half GC   : median 0.52, mean 0.5065   (0.52 for every human and mouse
#                                               set, 0.48 for rat)
#   full 50-mer GC: median 0.50, mean 0.5065
#
# The two medians differ only because GC is quantised on a 1/25 grid per half
# and a 1/50 grid over the full window; both distributions centre on ~0.506.
# 0.52 (the per-half median 10x actually ships) is used here. A 0.58 target
# sits above the 75th percentile of 10x's own probes, so it would reward
# GC-rich windows that 10x itself does not preferentially pick.
_GC_IDEAL = 0.52
# Half-width of the band, so a probe at either accepted edge still scores >= 0.
# The band is asymmetric about the empirical ideal, so this normalises by the
# wider side (0.72 - 0.52 = 0.20): a probe at the GC-poor edge (0.44) keeps a
# partial score, one at the GC-rich edge (0.72) scores 0.
_GC_HALF_SPAN = max(_GC_IDEAL - GC_MIN_DEFAULT, GC_MAX_DEFAULT - _GC_IDEAL)


def _gc_proximity_score(gc: float, ideal: float = _GC_IDEAL) -> float:
    """Score 0-1: how close a single GC value is to the empirical 10x ideal
    (1 = perfect).

    This is the per-value primitive; the ranking terms below never call it on
    whole-window GC (``candidate.gc_content``) — see
    :func:`_gc_proximity_per_half_score`.

    NOTE: gc_proximity and gc_balance are correlated (a probe with both halves
    near the ideal scores high on both), so their combined weight is
    intentionally the tool's main compositional signal — they are not
    independent axes."""
    return max(0.0, 1.0 - abs(gc - ideal) / _GC_HALF_SPAN)


def _gc_proximity_per_half_score(
    gc_left: float, gc_right: float, ideal: float = _GC_IDEAL,
) -> float:
    """Score 0-1: GC proximity scored PER HALF, not on the whole window.

    The 50 nt window is never synthesised as one molecule — it becomes two
    separately-synthesised 25 nt oligos (see the module docstring) that must
    each independently hybridise. GC content is therefore a per-half property,
    not a per-window one, and this is what ``filters.py``'s hard GC-range gate
    already does: it checks ``gc_left`` and ``gc_right`` against
    ``[gc_min, gc_max]`` separately and rejects the window if EITHER half is
    out of range (``filters.py`` ~line 499-505) — the window is never allowed
    to pass on the strength of a compensating other half.

    This term is the continuous analogue of that same rule: the WORST of the
    two per-half proximity scores decides the term, matching both the filter's
    "either half fails -> window fails" gate above and the worst-of-N
    precedent already used elsewhere in this module for mismatch placement
    (``_on_target_mismatch_placement_score``: "The WORST strain decides ...
    averaging would let a probe that fails outright on one strain hide behind
    the others"). A probe is only as good as its worse half.

    Averaging the two halves instead (equivalently, scoring whole-window GC,
    since the halves are equal length) would let a GC-rich half and a GC-poor
    half cancel out: a 0.72/0.32 pair averages to the ideal 0.52 and would
    score a perfect 1.0 even though NEITHER half is anywhere near it. Scoring
    the worse half prevents that.
    """
    return min(
        _gc_proximity_score(gc_left, ideal),
        _gc_proximity_score(gc_right, ideal),
    )


def _gc_balance_score(gc_left: float, gc_right: float) -> float:
    """Score 0-1: how balanced GC is between halves."""
    return max(0.0, 1.0 - abs(gc_left - gc_right) / 0.5)


# ---------------------------------------------------------------------------
# TN ligation junction — mRNA base at index 25.
# ---------------------------------------------------------------------------
#
# The ligation junction wants 'A' at mRNA position 25 (so the LHS probe's 3'
# base is T). By default that is a HARD filter in filters.py, so every surviving
# candidate has A there and this term is a constant 1.0 — exactly inert, no
# effect on hard-mode ranking. It exists for the opt-in *soft* mode, where the
# floor relaxes to "not C" and A/G/T all survive filtering; the ranker must then
# put A-junction probes decisively above the rest.
#
# Why a MULTIPLIER rather than another weighted term in the sum:
#   * Inertness. An additive term needs its weight in the normalising
#     denominator, which in hard mode (where it is constant) compresses the
#     variable part of the score by that same factor — silently inflating the
#     bounded primer3 nudge relative to the sourced axes. A multiplier of
#     exactly 1.0 changes nothing at all.
#   * Dominance. Getting an additive term to guarantee "every A beats every
#     non-A" needs a weight larger than all other weights combined, which is
#     the same thing as a gate but harder to read and impossible to guarantee
#     once a user overrides the weights.
#   * It matches the biology: the junction base is a precondition for ligation,
#     not one quality axis among several. If the probe cannot ligate, how good
#     its GC balance is does not matter.
#
# Values. Measured over all eight of 10x's public Chromium Flex probe sets
# (n = 439,656;
# ``probe_seq`` there is the PROBE strand, so the mRNA junction base is
# complement(probe_seq[24])):
#
#     mRNA[25]:   A 51.35%   G 35.72%   T 12.05%   C 0.88%
#     background: A 25.91%   G 25.86%   T 23.44%   C 24.79%
#     enrichment: A 1.98x    G 1.38x    T 0.51x    C 0.04x
#
# So 10x's own production designer strongly prefers A, accepts G roughly 3x as
# often as T (2.96:1 by frequency, 2.70:1 by enrichment), and all but forbids C
# — which is what the soft mode's "not C" floor encodes. G is therefore ranked
# above T rather than treating the two as equal.
#
# Magnitudes are set so the *ceiling* of a non-A probe sits below the *floor* of
# any realistic A-junction probe: a G-junction probe scores at most 0.20, below
# the composite of any A-junction probe that passes the filters in practice.
# This is a FlexiProbes heuristic, not a 10x specification. The G:T ratio of
# 2.86 sits between the two empirical ratios above.
JUNCTION_BASE_MULTIPLIERS: Dict[str, float] = {
    "A": 1.00,   # the wanted junction base
    "G": 0.20,   # tolerated by 10x ~3x as often as T
    "T": 0.07,
    "C": 0.00,   # 0.9% of 10x's panel; excluded outright by soft mode's floor
}
# Ambiguity codes / anything unrecognised: cannot verify the junction, so no
# credit. A missing or too-short sequence is different — see below.
JUNCTION_UNKNOWN_MULTIPLIER = 0.0


def _junction_base_multiplier(sequence: str) -> float:
    """Gate in [0, 1] from the mRNA junction base of a target window.

    Returns 1.0 (fully inert) when there is no sequence to inspect, so callers
    that score bare feature vectors — tests, the GUI's what-if sliders — are not
    silently zeroed out. A present-but-unreadable base scores 0.
    """
    seq = (sequence or "").upper()
    if not seq:
        return 1.0
    half = len(seq) // 2
    if half >= len(seq):
        return 1.0
    return JUNCTION_BASE_MULTIPLIERS.get(seq[half], JUNCTION_UNKNOWN_MULTIPLIER)


def _specificity_margin_score(
    margin: float, max_margin: float = OFF_TARGET_PIDENT_MIN,
) -> float:
    """Score 0-1: normalised specificity margin.

    For passing probes ``specificity_margin`` is the identity gap between the
    off-target disqualification threshold and the closest near-miss off-target
    (0..off_target_pident_min), so we normalise by that threshold rather than by
    an arbitrary 100."""
    if max_margin <= 0:
        return 0.0
    return min(1.0, max(0.0, margin / max_margin))


# Runs at or below this length are treated as unpenalised. Above it the score
# decays linearly. Set to 2 because the chemistry filter's own ceiling is a run
# of 4, so a floor of 4 leaves the term constant on every candidate that
# survives filtering (see the docstring below).
HOMOPOLYMER_FREE_RUN = 2


def _homopolymer_penalty_score(max_run: int, threshold: int = 4) -> float:
    """Score 0-1: penalise long homopolymers (1 = no penalty).

    The score must be graded INSIDE the passing range.
    ``FilterThresholds.max_homopolymer_run`` defaults to 4 and rejects anything
    longer, so a score that is flat at 1.0 for every ``max_run <= threshold``
    would be constant on every candidate that reaches the ranker. It would then
    carry no ordering information and only dilute the other terms. Runs of 2,
    3 and 4 do all occur among passing candidates, and on 10x's public panels
    shorter runs are more likely to be shipped.

    The score therefore decays from ``HOMOPOLYMER_FREE_RUN`` and stays
    monotone across the threshold:
    run <=2 -> 1.0, 3 -> 0.75, 4 -> 0.5, 5 -> 0.25, 6+ -> 0.0.
    """
    del threshold  # retained for call compatibility; the floor drives the score
    if max_run <= HOMOPOLYMER_FREE_RUN:
        return 1.0
    return max(0.0, 1.0 - (max_run - HOMOPOLYMER_FREE_RUN) * 0.25)


def _central_placement_score(
    start: int, end: int, seq_length: int,
    end_frac: float = END_EXCLUSION_FRAC,
) -> float:
    """Score 0-1: end avoidance, flat across the middle ``1 - 2*end_frac``.

    A symmetric end-avoidance prior (``RankingWeights.central_placement``,
    0.00 by default): probes at the extreme 5'/3' ends of a cDNA are more
    exposed to incomplete/degraded ends and terminal secondary structure.
    "Central in the cDNA string" is arbitrary for spliced transcripts, and a
    3'-bias is the usual alternative, but Flex is a ligation-probe chemistry
    tiled across the whole transcript, not a 3'-capture assay, so a 3' bias is
    not justified here.

    PLATEAU, NOT A PEAK. A linear tent scoring 1.0 only at the exact midpoint
    would be an *attractor* on a long transcript: every window more than a few
    hundred bases off-centre would pay a visible score penalty, so the
    top-ranked candidates a gene offered would all come from a narrow band
    around the midpoint and the chosen probes would cluster there however
    roomy the transcript was. Signal from three probes 200 bp apart in the
    middle of a 12 kb mRNA is close to three measurements of one site.

    Being central is not itself the property worth scoring — being *away from
    the ends* is. So the term is flat at 1.0 across the middle
    ``1 - 2*end_frac`` of the transcript and expresses no preference within it;
    what spreads probes out is the spacing target (``MIN_PROBE_GAP_BP`` and the
    adaptive spread gap in ``auto_select``), which scales with the transcript.

    Outside the plateau the falloff is QUADRATIC in how far into the terminal
    ``end_frac`` the probe sits, reaching 0 at the very first/last base. A
    window 5% into the transcript scores ``(0.05/0.15)^2 = 0.11``, close to a
    linear tent's 0.10, while a window at 30% scores a full 1.0 where a tent
    would give 0.60. The ends are strongly discouraged, and the interior is not
    ranked by centredness. See ``END_EXCLUSION_FRAC`` for how the terminal band
    compares with the probes 10x ship.
    """
    if seq_length <= 0 or end <= start:
        return 0.0
    frac = ((start + end) / 2.0) / float(seq_length)
    frac = min(1.0, max(0.0, frac))
    if end_frac <= 0.0:
        return 1.0
    edge = min(frac, 1.0 - frac)          # distance to the nearer end, 0..0.5
    if edge >= end_frac:
        return 1.0
    return (edge / end_frac) ** 2


# ---------------------------------------------------------------------------
# ON-target mismatch placement (INTERNAL heuristic, NOT 10x-sourced).
# ---------------------------------------------------------------------------
#
# Multi-strain runs are the only ones that can have on-target mismatches at all:
# a single-sequence probe is a slice of its own target, so it matches perfectly
# by construction. In a pan-strain (``conserved``) design the probe is claimed
# to bind EVERY input record, and a window that carries a variant position pays
# for it against whichever strain that variant belongs to.
#
# ``estimated_dTm`` and ``likely_binds`` report that cost in results.csv, but
# reach no gate. Without a ranking term, two probes with identical GC and
# specificity but mismatches in completely different places — one at the outer
# end of a half, one straddling the ligation nick — would score as equals.
# They are not equals: the ligase
# will not seal across a nick-proximal mispair, so that probe simply fails on
# that strain, while an end-proximal mismatch sits where the duplex frays anyway.
#
# This term prices that difference using the SAME positional ramp as the ΔTm
# estimate itself (``constants.mismatch_position_weight``), so the CSV column
# and the rank tell one story instead of two.
#
# ``MISMATCH_PLACEMENT_DTM_FLOOR`` is where the sub-score bottoms out. -12.0 C
# is two nick-adjacent mismatches (2 x -6.0) or six probe-end ones — a probe
# that has stopped being a pan-strain probe on its worst strain by any reading.
# Anything worse is already unusable, so grading past it would only spend the
# term's resolution on candidates nobody will order. It is a heuristic scale,
# not a measured melting point.
MISMATCH_PLACEMENT_DTM_FLOOR = -12.0


def _on_target_mismatch_placement_score(
    candidate: ProbeCandidate,
) -> Optional[float]:
    """Sub-score in [0, 1] for WHERE a probe's on-target mismatches sit.

    1.0 = no on-target mismatches anywhere (or none that cost anything), 0.0 =
    at or past ``MISMATCH_PLACEMENT_DTM_FLOOR``. Between those it falls linearly
    in the position-weighted ΔTm, so for the same mismatch COUNT a probe whose
    mismatches sit near the probe ends outscores one whose mismatches sit near
    the ligation nick. On the standard 50-mer a single mismatch scores 0.833 at
    the very end of a half and 0.500 right at the nick.

    Returns **None** — not 1.0 — when the candidate carries no on-target
    mismatch data at all, i.e. ``strain_binding`` is empty. That is every
    single-sequence run, and the distinction matters: ``compute_composite_score``
    drops an inapplicable term from the normalising denominator as well as from
    the sum, exactly as :func:`prescreen_score` drops the weights of the terms
    BLAST has not supplied yet. Folding it in as a constant 1.0 instead would
    renormalise every single-sequence score by 1/1.10 and silently move numbers
    that this term must leave alone.

    The WORST strain decides, matching ``worst_strain_dTm`` in results.csv: a
    pan-strain probe is only as good as the strain it binds least well, and
    averaging would let a probe that fails outright on one strain hide behind
    the others it matches perfectly.
    """
    binding = getattr(candidate, "strain_binding", None) or []
    if not binding:
        return None
    half_length = len(candidate.sequence or "") // 2
    worst_weight = 0.0
    for entry in binding:
        positions = entry.get("mismatch_positions") or []
        worst_weight = max(worst_weight, sum(
            mismatch_position_weight(int(p), half_length) for p in positions
        ))
    dtm = DTM_PER_MISMATCH * worst_weight          # <= 0
    if MISMATCH_PLACEMENT_DTM_FLOOR >= 0:          # defensive; floor is negative
        return 1.0
    return 1.0 - min(1.0, max(0.0, dtm / MISMATCH_PLACEMENT_DTM_FLOOR))


def prescreen_score(
    candidate: ProbeCandidate,
    seq_length: int = 0,
    weights: RankingWeights = DEFAULT_WEIGHTS,
) -> float:
    """BLAST-free quality score in [0, 1] for *pre*-screening candidates.

    Uses only the features known before BLAST — GC proximity, GC balance,
    homopolymer penalty, central placement and the junction-base gate —
    dropping the specificity-margin and isoform-coverage terms (which need
    BLAST). Normalised by the sum of just those weights so it stays in [0, 1]
    and tracks the same quality notion the final composite score uses. Lets the
    pipeline BLAST only the best N windows per gene instead of every candidate.

    The junction gate is carried here too, so that in soft mode the per-gene
    pre-screen cap spends its budget on A-junction windows rather than filling
    up with G/T ones that the final ranking would then bury. In hard mode it is
    a constant 1.0 and the cap is unchanged.
    """
    gc_prox = _gc_proximity_per_half_score(candidate.gc_left, candidate.gc_right)
    gc_bal = _gc_balance_score(candidate.gc_left, candidate.gc_right)
    hp = _homopolymer_penalty_score(candidate.homopolymer_max)
    central = _central_placement_score(candidate.start, candidate.end, seq_length)
    denom = (weights.gc_proximity + weights.gc_balance +
             weights.homopolymer_penalty + weights.central_placement)
    if denom <= 0:
        return 0.0
    return _junction_base_multiplier(candidate.sequence) * (
        weights.gc_proximity * gc_prox +
        weights.gc_balance * gc_bal +
        weights.homopolymer_penalty * hp +
        weights.central_placement * central
    ) / denom


def select_top_per_gene(
    candidates: List[ProbeCandidate],
    *,
    gene_key: Callable[[ProbeCandidate], str],
    seq_length_of: Callable[[str], int],
    n: int,
    weights: RankingWeights = DEFAULT_WEIGHTS,
) -> List[ProbeCandidate]:
    """Keep the top-``n`` candidates per gene by :func:`prescreen_score`.

    ``gene_key`` maps a candidate to its biological-gene grouping key (e.g. the
    Ensembl gene for isoform-aware runs); ``seq_length_of`` maps a candidate's
    transcript id to its length (for the central-placement term). Genes with
    ``<= n`` candidates are kept whole. Deterministic: ties break by lowest
    ``start`` then ``probe_id``. The returned list preserves the input order of
    the kept candidates, so the downstream BLAST chunking stays stable.
    """
    if n <= 0:
        return list(candidates)
    by_gene: Dict[str, List[ProbeCandidate]] = {}
    for c in candidates:
        by_gene.setdefault(gene_key(c), []).append(c)
    keep_ids: Set[str] = set()
    for cands in by_gene.values():
        if len(cands) <= n:
            keep_ids.update(c.probe_id for c in cands)
            continue
        ordered = sorted(
            cands,
            key=lambda c: (
                -prescreen_score(c, seq_length_of(c.gene_id), weights),
                c.start, c.probe_id,
            ),
        )
        keep_ids.update(c.probe_id for c in ordered[:n])
    return [c for c in candidates if c.probe_id in keep_ids]


class TopPerGeneAccumulator:
    """Memory-bounded, streaming form of :func:`select_top_per_gene`.

    Feed candidates in any order with :meth:`add` / :meth:`add_many`, then call
    :meth:`result` for the kept set. Per gene it holds at most
    ``trim_factor * n`` candidates between trims and exactly ``n`` after the
    final trim, so peak memory is O(kept set) — not O(all candidates). That is
    what lets a whole-transcriptome run cap to a few million windows without ever
    materialising the hundreds of millions that pass filtering.

    The kept SET is identical to ``select_top_per_gene`` with the same ``n``,
    ``weights`` and tie-break (highest ``prescreen_score``, then lowest ``start``,
    then lowest ``probe_id``): an intermediate trim can never drop a global
    top-``n`` member, because such a member out-ranks all but < n candidates and
    so sits at position <= n in any subset that contains it. Only the *order* of
    the returned list differs from ``select_top_per_gene`` (here it is sorted by
    gene key, then start, then probe_id for determinism); order is irrelevant
    downstream — every probe is BLAST-scored independently and results are
    re-ranked per gene afterwards.
    """

    def __init__(
        self,
        *,
        gene_key: Callable[[ProbeCandidate], str],
        seq_length_of: Callable[[str], int],
        n: int,
        weights: RankingWeights = DEFAULT_WEIGHTS,
        trim_factor: int = 4,
    ):
        self._gene_key = gene_key
        self._seq_length_of = seq_length_of
        self._n = n
        self._weights = weights
        # Trim a gene's buffer once it grows past this; >= 2*n so trims amortise.
        self._limit = max(1, n) * max(2, trim_factor)
        self._buf: Dict[str, List[ProbeCandidate]] = {}

    def _sort_key(self, c: ProbeCandidate):
        return (
            -prescreen_score(c, self._seq_length_of(c.gene_id), self._weights),
            c.start, c.probe_id,
        )

    def add(self, c: ProbeCandidate) -> None:
        gk = self._gene_key(c)
        buf = self._buf.get(gk)
        if buf is None:
            buf = []
            self._buf[gk] = buf
        buf.append(c)
        if self._n > 0 and len(buf) > self._limit:
            buf.sort(key=self._sort_key)
            del buf[self._n:]

    def add_many(self, cands) -> None:
        for c in cands:
            self.add(c)

    def result(self) -> List[ProbeCandidate]:
        out: List[ProbeCandidate] = []
        for buf in self._buf.values():
            if self._n > 0 and len(buf) > self._n:
                buf.sort(key=self._sort_key)
                del buf[self._n:]
            out.extend(buf)
        # Deterministic order, independent of how candidates were fed in.
        out.sort(key=lambda c: (self._gene_key(c), c.start, c.probe_id))
        return out


def compute_composite_score(
    candidate: ProbeCandidate,
    seq_length: int = 0,
    weights: RankingWeights = DEFAULT_WEIGHTS,
) -> float:
    """Compute and assign composite score in [0, 1]; return the score."""
    gc_prox = _gc_proximity_per_half_score(candidate.gc_left, candidate.gc_right)
    gc_bal = _gc_balance_score(candidate.gc_left, candidate.gc_right)
    spec = _specificity_margin_score(candidate.specificity_margin)
    hp = _homopolymer_penalty_score(candidate.homopolymer_max)
    central = _central_placement_score(candidate.start, candidate.end, seq_length)
    # Already a 0-1 fraction; clamp defensively in case of bad input.
    iso = max(0.0, min(1.0, candidate.isoform_coverage))
    # None when this run has no on-target mismatch data (every single-sequence
    # run). The term then contributes to NEITHER the sum nor the denominator,
    # so such a run's score is bit-identical to a six-term score. See
    # _on_target_mismatch_placement_score.
    mm_raw = _on_target_mismatch_placement_score(candidate)
    mm = 1.0 if mm_raw is None else mm_raw
    mm_weight = 0.0 if mm_raw is None else weights.mismatch_placement

    score = (
        weights.gc_proximity * gc_prox +
        weights.gc_balance * gc_bal +
        weights.specificity_margin * spec +
        weights.homopolymer_penalty * hp +
        weights.central_placement * central +
        weights.isoform_coverage * iso +
        mm_weight * mm
    )
    # Normalise by the weight total so the composite score always stays in
    # [0, 1] no matter how the weights are set: they act as *relative*
    # importances and need not sum to 1. Dividing every candidate by the same
    # constant leaves ranking order untouched, and with the defaults (which sum
    # to 1.0) the score is identical to the un-normalised value — so enabling
    # the isoform-coverage weight cannot push scores above 1.0.
    #
    # The mismatch-placement weight participates only when the term does. Every
    # summand is a weight times a value in [0, 1] and the denominator is exactly
    # the sum of those same weights, so the quotient is a convex combination and
    # stays in [0, 1] whichever branch is taken — adding the seventh term cannot
    # break the invariant. Dropping an inapplicable term from the denominator is
    # what prescreen_score already does with the BLAST-dependent weights.
    #
    # Summed term by term in RankingWeights.total()'s own order rather than as
    # ``total() - mismatch_placement``: with the term applicable this IS total(),
    # and with it inapplicable it is bit-for-bit the six-weight total, with no
    # float subtraction that could perturb a single-sequence score in the last
    # ulp.
    total = (weights.gc_proximity + weights.gc_balance +
             weights.specificity_margin + weights.homopolymer_penalty +
             weights.central_placement + weights.isoform_coverage +
             mm_weight)
    score = score / total if total > 0 else 0.0
    # TN ligation junction gate. Exactly 1.0 (inert) in the default hard mode,
    # where the filter has already guaranteed 'A'; in soft mode it is what keeps
    # A-junction probes decisively above G/T ones. Multiplying by a value in
    # [0, 1] cannot leave the unit range. See JUNCTION_BASE_MULTIPLIERS.
    junction = _junction_base_multiplier(candidate.sequence)
    score *= junction
    candidate.composite_score = score
    # Expose the raw 0-1 sub-scores so a rank is auditable (which axis drove it).
    candidate.score_components = {
        "gc_proximity": gc_prox,
        "gc_balance": gc_bal,
        "specificity_margin": spec,
        "homopolymer_penalty": hp,
        "central_placement": central,
        "isoform_coverage": iso,
        # 1.0 here means "nothing to penalise" — either no on-target mismatches,
        # or (on a single-sequence run) no on-target mismatch data at all, in
        # which case the term was left out of the weighting entirely.
        "mismatch_placement": mm,
        # Multiplicative, not weighted-additive like the seven above.
        "junction_base": junction,
    }
    return score


def rank_candidates(
    candidates: List[ProbeCandidate],
    seq_length: int = 0,
    weights: RankingWeights = DEFAULT_WEIGHTS,
    primer3_scoring: bool = False,
) -> List[ProbeCandidate]:
    """Score and rank candidates deterministically.

    Parameters
    ----------
    candidates : list
        Must have passed filters and specificity.
    seq_length : int
        Full gene length (for central-placement component).
    weights : RankingWeights
        Component weights.
    primer3_scoring : bool
        If True and primer3 is available, add Tm/dimer adjustments.
        OFF by default so ranking is deterministic without primer3.

    Returns
    -------
    Sorted list (descending score). Each candidate gets .rank set.
    """
    for c in candidates:
        compute_composite_score(c, seq_length, weights)
        c.reduced_confidence = False

    if primer3_scoring:
        # Probe availability before per-item work, without an unused import
        # (find_spec keeps this lint-clean under pyflakes/ruff alike).
        if importlib.util.find_spec("primer3") is not None:
            _apply_primer3_adjustment(candidates)
        else:
            logger.warning(
                "primer3 not installed; skipping Tm/hairpin scoring. "
                "Composite scores are unaffected."
            )

    # A degenerate parent could not be scored on its own string above, so fill
    # its structure fields from the variants that WERE scored. Runs whatever
    # primer3_scoring was, because the Tm/GC ranges are populated at tile time
    # and are meaningful without primer3.
    aggregate_consensus_parents(candidates)

    # Deterministic sort: by score desc, then by probe_id asc (tie-break)
    candidates.sort(key=lambda c: (-c.composite_score, c.probe_id))
    for i, c in enumerate(candidates, 1):
        c.rank = i
    return candidates


def aggregate_consensus_parents(candidates: List[ProbeCandidate]) -> None:
    """Fill each consensus parent's aggregate fields from its own variants.

    ``_apply_primer3_adjustment`` skips any non-ACGT sequence and flags it
    ``reduced_confidence`` — correct, primer3 genuinely cannot score a mixed-base
    string, and that skip stays. But leaving the resulting fields blank would be
    the wrong answer to the question the review table is asked. Every variant of
    a site is synthesised into the SAME pool, so the pool behaves like its worst
    member: a parent whose second variant has a -9 kcal/mol hairpin is a risky
    probe, and a blank hairpin column sorts it to the same place as a probe with
    no structure at all — hiding exactly the probes carrying the most structural
    risk from anyone sorting by it.

    So a parent carries the RANGE where the question is "how much does the
    ambiguity cost" (Tm, GC — an 'R' has two resolutions and a narrow spread, an
    'N' four and a wide one) and the WORST CASE where the question is "will this
    pool work" (hairpin, heterodimer, Tm balance).

    ``reduced_confidence`` stays True on a parent, but on a parent it means
    "aggregated from variants", NOT "unscored".

    Tm/GC ranges are populated earlier, at tile time, because the filter funnel
    needs them to gate the parent at all (``pipeline._degenerate_stability``);
    they are only recomputed here if that never happened.
    """
    parents = [c for c in candidates if c.probe_role == PROBE_ROLE_CONSENSUS]
    if not parents:
        return
    kids: Dict[str, List[ProbeCandidate]] = {}
    for c in candidates:
        if c.probe_role == PROBE_ROLE_VARIANT and c.parent_probe_id:
            kids.setdefault(c.parent_probe_id, []).append(c)

    def _worst(vals: List[float], pick) -> float:
        real = [v for v in vals if v == v]          # drop NaN
        return pick(real) if real else float("nan")

    for p in parents:
        children = kids.get(p.probe_id, [])
        if not children:
            continue
        # Worst case = most negative ΔG (structure forms most readily), and the
        # largest LHS/RHS Tm gap (the half that melts off first decides).
        p.hairpin_dg_worst = _worst([c.hairpin_dg for c in children], min)
        p.hetero_dimer_dg_worst = _worst(
            [c.hetero_dimer_dg for c in children], min)
        p.tm_balance_max = _worst([c.tm_balance for c in children], max)
        # Ranges, if tile time did not already provide them (e.g. a candidate
        # set rebuilt outside the pipeline).
        if p.tm_lhs_min != p.tm_lhs_min:            # NaN
            half = len(p.sequence) // 2
            from ._utils import gc_fraction, safe_tm
            lhs = [safe_tm(c.sequence[:half]) for c in children]
            rhs = [safe_tm(c.sequence[half:]) for c in children]
            p.tm_lhs_min, p.tm_lhs_max = _worst(lhs, min), _worst(lhs, max)
            p.tm_rhs_min, p.tm_rhs_max = _worst(rhs, min), _worst(rhs, max)
            gl = [gc_fraction(c.sequence[:half]) for c in children]
            gr = [gc_fraction(c.sequence[half:]) for c in children]
            p.gc_lhs_min, p.gc_lhs_max = min(gl), max(gl)
            p.gc_rhs_min, p.gc_rhs_max = min(gr), max(gr)
        # Say out loud what reduced_confidence means on this row.
        p.reduced_confidence = True


def _target_oligo_halves(seq: str) -> tuple:
    """Return ``(target_LHS, target_RHS)`` for a sense-strand mRNA window.

    These are the two 25-mers that are actually *synthesised* — the
    target-binding halves of the LHS and RHS oligos — as built by
    :func:`flexprobe.chemistry.ChemistryConfig.split_tile`:

        target_LHS = RC(mRNA 3' half)      target_RHS = RC(mRNA 5' half)

    The mRNA window itself is never synthesised, so running structure
    prediction on it describes a molecule that does not exist in the tube.
    primer3's answers are not reverse-complement invariant in practice, so the
    strand matters: ranking the same probes by mRNA-window homodimer and by
    probe-strand homodimer can give noticeably different orders.
    """
    half = len(seq) // 2
    return reverse_complement(seq[half:]), reverse_complement(seq[:half])


def _apply_primer3_adjustment(candidates: List[ProbeCandidate]) -> None:
    """Optional primer3-based Tm and dimer scoring, per-candidate.

    Structure is computed on the two target-binding 25-mers (``target_LHS`` /
    ``target_RHS``), not on the mRNA window — see :func:`_target_oligo_halves`
    and the "which molecule?" note in the module header.

    The adjustment is bounded to ``[PRIMER3_ADJUSTMENT_MIN,
    PRIMER3_ADJUSTMENT_MAX]`` = [-0.20, +0.05] and the result is clamped into
    [0, 1], so ``composite_score`` keeps the range its docstring promises.

    A primer3 failure on one candidate (e.g. an ambiguous base it can't
    handle) only flags THAT candidate as reduced-confidence — it no longer
    taints the whole batch (C3). Candidates with non-ACGT bases are skipped
    cleanly rather than raising.
    """
    import primer3 as p3
    for c in candidates:
        seq = (c.sequence or "").upper()
        if not seq or any(b not in "ACGT" for b in seq):
            # primer3 rejects ambiguity codes — leave score unadjusted and
            # flag so the user knows Tm/hairpin weren't factored in.
            c.reduced_confidence = True
            continue
        try:
            salt = dict(mv_conc=TM_MV_CONC, dv_conc=TM_DV_CONC,
                        dntp_conc=TM_DNTP_CONC, dna_conc=TM_DNA_CONC)
            target_lhs, target_rhs = _target_oligo_halves(seq)
            tm = p3.calc_tm(seq, **salt)
            # Per-half Tm uniformity: the two 25-mer halves ligate
            # independently, so a big Tm gap means one half melts off first.
            # (Tm is RC-invariant under the nearest-neighbour model, so these
            # are the same numbers the mRNA halves gave.)
            tm_lhs = p3.calc_tm(target_lhs, **salt)
            tm_rhs = p3.calc_tm(target_rhs, **salt)
            c.tm_balance = abs(tm_lhs - tm_rhs)
            # Self / cross structure, on the synthesised oligos:
            #   hairpin / homodimer — worst (most negative) of the two halves
            #   heterodimer         — LHS half against RHS half, the two
            #                         species that share a tube
            # ΔG < 0 = structure forms.
            hairpin_dg = min(p3.calc_hairpin(target_lhs, **salt).dg,
                             p3.calc_hairpin(target_rhs, **salt).dg)
            # Persisted, not just consumed by the nudge below: a consensus
            # parent's ``hairpin_dg_worst`` is aggregated from its variants,
            # which is impossible if each variant's own hairpin is discarded
            # here.
            c.hairpin_dg = float(hairpin_dg)
            c.self_dimer_dg = float(min(
                p3.calc_homodimer(target_lhs, **salt).dg,
                p3.calc_homodimer(target_rhs, **salt).dg,
            ))
            c.hetero_dimer_dg = float(
                p3.calc_heterodimer(target_lhs, target_rhs, **salt).dg
            )
            # Internal heuristic nudges (NOT 10x-sourced), deliberately small so
            # primer3 only breaks ties among probes the sourced scores rank
            # similarly. Tunable via the constants at the top of this module.
            # Opt-in only.
            tm_bonus = max(
                0.0, TM_BONUS_MAX - abs(tm - TM_TARGET_C) * TM_BONUS_PER_C
            )
            tm_bal_penalty = min(TM_BALANCE_PENALTY_MAX,
                                 c.tm_balance * TM_BALANCE_PER_C)
            adjustment = (
                tm_bonus
                - tm_bal_penalty
                - _dg_penalty(hairpin_dg)
                - _dg_penalty(c.self_dimer_dg)
                - _dg_penalty(c.hetero_dimer_dg)
            )
            # Belt and braces: the arithmetic above cannot leave this window,
            # but the invariant is the published contract, so assert it here
            # rather than discovering a violation in a shipped CSV.
            adjustment = min(PRIMER3_ADJUSTMENT_MAX,
                             max(PRIMER3_ADJUSTMENT_MIN, adjustment))
            # The nudge perturbs the sourced score, so it scales with the same
            # junction gate that scaled it — otherwise a gated-down probe would
            # feel a proportionally huge nudge. 1.0 (no-op) in hard mode.
            adjustment *= (c.score_components or {}).get("junction_base", 1.0)
            c.composite_score = min(1.0, max(0.0, c.composite_score + adjustment))
            if c.score_components is not None:
                c.score_components["primer3_adjustment"] = adjustment
        except Exception as exc:
            logger.warning(
                "primer3 scoring failed for %s (%s); marking reduced-confidence",
                c.probe_id, exc,
            )
            c.reduced_confidence = True
