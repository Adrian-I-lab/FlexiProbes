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

"""Filter funnel for probe candidates."""
from __future__ import annotations
import importlib.util
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple
from ._utils import degenerate_count, is_nan, max_homopolymer, safe_tm
from .constants import (
    DEGENERATE_BASES, DEGENERATE_JUNCTION_CLEARANCE, PROBE_ROLE_CONSENSUS,
    PROBE_ROLE_PLAIN, TN_MODE_AUTO, TN_MODE_OFF, TN_MODE_SOFT, TN_MODE_STRICT_A,
)
from .models import ProbeCandidate, FilterThresholds, RunConfig

# Alignment gap characters. A tile spanning a gap ('-'/'.') — e.g. from an
# aligned / MSA-derived FASTA — can never be a valid probe, so it is dropped up
# front. NOTE this does NOT catch a consensus 'N' that came from a gap column:
# `consensus._consensus_from_aligned` turns a gap into 'N', which is a perfectly
# ordinary character here. Excluding indel windows is therefore the job of
# window selection (`consensus.degenerate_conserved_windows`, rule 1) and
# nothing downstream will catch it if that rule is relaxed.
_GAP_TILE_CHARS = frozenset("-.")

# Cached once: primer3 present? Used only to tell "we could not judge this
# parent's stability because primer3 is missing" (honest, pass) apart from "we
# could have, and the fan-out failed to give us the numbers" (a bug, reject).
_HAS_PRIMER3: Optional[bool] = None


def _primer3_available() -> bool:
    global _HAS_PRIMER3
    if _HAS_PRIMER3 is None:
        _HAS_PRIMER3 = importlib.util.find_spec("primer3") is not None
    return _HAS_PRIMER3


def _junction_band(
    thresholds: FilterThresholds, n: int, half: int,
) -> Tuple[int, int]:
    """Half-open offsets ``[lo, hi)`` of the ligation-junction clearance band.

    One definition, used by BOTH junction-clearance rules below (the
    degenerate-base rule and the concrete-mismatch rule) so they can never
    drift apart, and matching ``consensus.degenerate_conserved_windows``
    rule 3, which applies the same band at window-selection time.

    *** THE OFF-BY-ONE HERE IS AN OPEN INTERPRETATION, NOT A DERIVED RESULT. ***
    CG000621 Rev D says only:

      "Avoid common single nucleotide polymorphisms (SNPs) and potential
       mismatches at the ligation junction. Refer to the UCSC Genome Browser and
       the Single Nucleotide Polymorphism Database (dbSNP). If avoiding SNPs is
       not possible, SNPs and mismatches should be at least four bp away from
       the ligation junction."

    For a 50-mer at the default clearance of 4 this returns ``(21, 29)``, i.e.
    it excludes offsets 21-24 and 25-28 — the four bases immediately either side
    of the nick. Whether that is EXACTLY 10x's rule or one base stricter than it
    depends on something the technical note does not state: whether the base
    adjacent to the nick is "0 bp away" (making the fourth base out 3 bp away,
    so a four-bp exclusion is correct) or "1 bp away" (making the fourth base
    out 4 bp away, hence already permitted, so this band is one base too wide).
    Nothing in CG000621 or CG000839 resolves it, and neither reading is asserted
    here. The band is deliberately left on the stricter side, which is a
    flexiprobes choice, not a 10x specification — the same "slightly stricter
    than the letter, to improve reliability" stance the rest of this funnel
    takes. Do not narrow it without evidence that resolves the ambiguity.

    ``clearance`` 0 returns an empty band (``hi <= lo``), which every caller
    reads as "no junction rule"; use the ``reject_junction_*`` switches to turn
    a rule off deliberately rather than by zeroing the width.
    """
    clearance = max(0, int(getattr(
        thresholds, 'degenerate_junction_clearance',
        DEGENERATE_JUNCTION_CLEARANCE)))
    return max(0, half - clearance), min(n, half + clearance)


def on_target_record_ids(
    config: RunConfig,
    strain_ids: Sequence[str],
) -> Optional[FrozenSet[str]]:
    """Which input records a probe from ``config`` is REQUIRED to bind.

    This is the mode-awareness behind the ``junction_mismatch`` rule, and it is
    the part that has to be right: a junction-proximal mismatch is fatal against
    a record the probe must bind and DESIRABLE against a record it must
    discriminate from, so scoping it wrongly does not merely over-reject — in
    strain-unique mode it would reject exactly the windows the run exists to
    find. Returns ``None`` for "no scope; the rule does not apply", which is the
    safe direction (the gate stays inert rather than guessing).

    Per ``RunConfig.consensus_strategy``:

    ``conserved``
        Every record. A pan-strain probe's whole claim is that it binds all of
        them, so a mismatch to ANY of them inside the band disqualifies it.
    ``strain_unique``
        ONLY ``strain_unique_target``'s record. The probe is that strain's
        literal sequence and its mismatches against the others are the
        discrimination being designed for — 10x's own off-target rule wants
        them, and ``ProbeCandidate.offtarget_worst_junction_mm`` already notes
        that junction-proximal mismatches read backwards against a non-target.
    ``constitutive``
        Every record when the quorum is every record (``isoform_min_count`` 0 or
        >= the record count), which is the default and makes the mode's own
        rule "shared by all". Below that quorum the probe is designed for ONE
        isoform group and the records outside it are not targets, and WHICH
        group differs per window — so no run-level answer exists and this
        returns ``None``. That is not the end of the scoping: the pipeline
        carries each window's covered group onto the candidate as
        ``ProbeCandidate.covered_record_ids``, and
        :func:`candidate_on_target_records` prefers it over whatever this
        function returned. A partial-quorum run is therefore scoped
        per-candidate, not skipped.

    A single-sequence (non-multi-strain) run also returns ``None``: there is
    only one record, the probe is its literal subsequence, and no candidate on
    that path carries ``strain_binding`` for the rule to read anyway.
    """
    if not getattr(config, 'multi_strain', False) or not strain_ids:
        return None
    strategy = getattr(config, 'consensus_strategy', 'conserved')
    if strategy == 'strain_unique':
        idx = int(getattr(config, 'strain_unique_target', 0))
        if not (0 <= idx < len(strain_ids)):
            # Out of range is a config error the pipeline raises on separately;
            # do not silently fall back to "every record", which is the exact
            # inversion this function exists to prevent.
            return None
        return frozenset({strain_ids[idx]})
    if strategy == 'conserved':
        return frozenset(strain_ids)
    if strategy == 'constitutive':
        quorum = int(getattr(config, 'isoform_min_count', 0) or 0)
        if quorum <= 0 or quorum >= len(strain_ids):
            return frozenset(strain_ids)
        return None
    return None


def candidate_on_target_records(
    pc: ProbeCandidate,
    run_scope: Optional[FrozenSet[str]],
) -> Optional[FrozenSet[str]]:
    """The junction-mismatch scope for ONE candidate.

    :func:`on_target_record_ids` answers "which records must a probe from this
    RUN bind", which is the whole answer for every mode whose target set is the
    same for every window. Constitutive mode below its quorum is the exception:
    each window is the shared sequence of its own largest isoform group, so the
    target set is a property of the window and there is no run-level answer to
    give. The pipeline records that group on the candidate and this prefers it.

    Precedence is candidate-first and deliberate. ``covered_record_ids`` is only
    ever set where the covered group IS the set of records the probe must bind
    (see the field's own warning in ``models``), so when it is present it is the
    more specific and more correct scope; when it is absent — every other mode,
    every single-sequence run, every candidate built by a caller that knows
    nothing about this — the run-level scope is used unchanged and behaviour is
    exactly what it was before per-candidate scoping existed.

    Returns ``None`` only when neither is available, which leaves the rule
    inert. That remains the safe direction: an un-scoped junction-mismatch rule
    would reject exactly the windows a strain-unique run exists to find.
    """
    covered = getattr(pc, 'covered_record_ids', None)
    if covered:
        return frozenset(covered)
    return run_scope


def new_filter_counters() -> Dict[str, int]:
    """Fresh counter dict for the filter funnel."""
    return {
        'total': 0,
        'gap_or_invalid': 0,
        'junction_base': 0,
        # Deliberately NOT folded into 'junction_base'. Both rules police the
        # same band, but they are different evidence and have different blast
        # radii: 'junction_base' is a property of the probe string, readable
        # from the sequence in filter_report.json, while 'junction_mismatch' is
        # a statement about the input records and can only be checked on a
        # multi-strain run. Merging them would make a report unreadable and
        # would hide a mode-scoping bug (see on_target_record_ids) behind a
        # count that looks normal.
        'junction_mismatch': 0,
        'tn_junction': 0,
        'n_count': 0,
        # A consensus parent that arrived without a usable variant aggregate —
        # a fan-out bug, not a property of the window. Counted separately so it
        # can never be read as "the chemistry rejected these".
        'degenerate_unresolved': 0,
        'gc_range': 0,
        'tm_range': 0,
        'homopolymer': 0,
        'passed': 0,
    }


def evaluate_candidate(
    pc: ProbeCandidate,
    thresholds: FilterThresholds,
    counters: Dict[str, int],
    *,
    on_target_records: Optional[FrozenSet[str]] = None,
) -> bool:
    """Apply the filter funnel to ONE candidate in place.

    Updates ``pc.passed_filters`` / ``pc.filter_fail_reasons`` / GC + homopolymer
    fields and increments ``counters`` (including 'total' and 'passed').
    Returns ``pc.passed_filters``. This is the streaming-friendly core used by
    both ``apply_filters`` (list) and the pipeline's streaming loop.

    ``on_target_records`` is the RUN-level set of input-record ids this probe is
    required to bind, from :func:`on_target_record_ids`. A candidate that
    carries its own ``covered_record_ids`` overrides it for that candidate
    (:func:`candidate_on_target_records`); only constitutive mode sets that.
    Between them they are the ONLY thing that engages the ``junction_mismatch``
    rule — with neither present the rule stays inert, so a caller that knows
    nothing about the run's strain layout (a bare ``apply_filters`` in a test, a
    downstream tool) behaves exactly as it did before the rule existed.
    """
    counters['total'] += 1
    seq = pc.sequence
    seq_u = seq.upper()
    half = len(seq) // 2
    left = seq_u[:half]
    right = seq_u[half:]
    reasons: List[str] = []
    _filter_body(pc, thresholds, counters, seq_u, half, left, right, reasons,
                 on_target_records=on_target_records)
    if reasons:
        pc.passed_filters = False
        pc.filter_fail_reasons = reasons
    else:
        pc.passed_filters = True
        counters['passed'] += 1
    return pc.passed_filters


def apply_filters(
    candidates: List[ProbeCandidate],
    thresholds: FilterThresholds,
    *,
    on_target_records: Optional[FrozenSet[str]] = None,
) -> Dict[str, int]:
    """Apply the filter funnel to a list in-place; return counter dict.

    Thin wrapper over ``evaluate_candidate`` kept for callers/tests that work
    on a materialized list. The streaming pipeline calls ``evaluate_candidate``
    directly so it never holds every tile at once. ``on_target_records`` is
    forwarded unchanged; see :func:`evaluate_candidate`.
    """
    counters = new_filter_counters()
    for pc in candidates:
        evaluate_candidate(pc, thresholds, counters,
                           on_target_records=on_target_records)
    return counters


def _filter_body(
    pc: ProbeCandidate,
    thresholds: FilterThresholds,
    counters: Dict[str, int],
    seq_u: str,
    half: int,
    left: str,
    right: str,
    reasons: List[str],
    *,
    on_target_records: Optional[FrozenSet[str]] = None,
) -> None:
    """Run every filter check on one candidate, appending fail reasons.

    Mutates ``pc`` (gc + homopolymer fields), ``counters`` (per-reason), and
    ``reasons`` (the running fail list). Does not decide pass/fail — the
    caller does that from whether ``reasons`` is empty.
    """
    n = len(seq_u)

    # Gap filter — runs first. A tile that spans an alignment gap ('-'/'.')
    # can never be a valid probe; drop it before the chemistry/GC checks. This
    # lets aligned/consensus inputs be tiled in their original coordinates (so
    # start/end stay faithful) while gap windows are silently excluded.
    if _GAP_TILE_CHARS & set(seq_u):
        reasons.append('gap_or_invalid')
        counters['gap_or_invalid'] += 1
        return  # no point scoring GC/homopolymer on a gap tile

    # Junction-ambiguity filter. Supersedes the old two-base check, which
    # looked at exactly ``seq_u[half-1:half+1]`` — the pair straddling the nick
    # — and only for a literal 'N'. Both halves of that were too narrow:
    # CG000621 Rev D asks for mismatches at least four bp from the ligation
    # junction, and the ligase will not seal across a junction-proximal mispair,
    # so ambiguity anywhere in the clearance band makes the probe's behaviour on
    # the non-reference record unpredictable in the one place the chemistry
    # cannot tolerate it. An 'M' at offset 22 was as bad as an 'N' at 24 and the
    # old check saw neither. Band derived from the tile length, not a fixed 25.
    if thresholds.reject_junction_base_n and n >= 2 and half > 0:
        j_lo, j_hi = _junction_band(thresholds, n, half)
        if j_hi > j_lo and DEGENERATE_BASES & set(seq_u[j_lo:j_hi]):
            reasons.append('junction_base')
            counters['junction_base'] += 1

    # Junction-MISMATCH filter — the same band, for a probe that is already
    # concrete ACGT. The rule above can only see ambiguity written into the
    # probe string, and that is not the same thing as a mismatch:
    #
    #   "Avoid common single nucleotide polymorphisms (SNPs) and potential
    #    mismatches at the ligation junction. Refer to the UCSC Genome Browser
    #    and the Single Nucleotide Polymorphism Database (dbSNP). If avoiding
    #    SNPs is not possible, SNPs and mismatches should be at least four bp
    #    away from the ligation junction."
    #                                              — CG000621 Rev D
    #
    # Note what that sentence actually asks for. The PRIMARY instruction is to
    # avoid SNPs and mismatches at the junction outright; "at least four bp
    # away" is the CONCESSION for when avoiding them is impossible. This code
    # applies the concession as a hard gate, which is stricter than the letter
    # — a flexiprobes choice, not 10x's specification. See _junction_band for
    # the (unresolved) off-by-one in the band's width.
    #
    # Where the concrete mismatch comes from: on a multi-strain run every
    # candidate is scored against every input record at tile time
    # (pipeline._annotate_strain_binding -> consensus.score_probe_against_all_
    # strains), so ``pc.strain_binding[i]['mismatch_positions']`` holds the
    # exact window-relative offsets at which this probe mispairs record i. That
    # scoring reads an IUPAC code in the probe as "any of these bases", so a
    # degenerate parent contributes no mismatch at its own ambiguous positions
    # — those are the rule above's business, and the two do not double-count.
    #
    # *** THE SCOPE IS THE WHOLE SAFETY ARGUMENT. *** ``scope`` names the
    # records this probe is REQUIRED to bind, and nothing outside it is
    # consulted. A junction mismatch against a record the probe must bind is
    # fatal (the ligase will not seal across a junction-proximal mispair); the
    # identical mismatch against a record the probe is meant to DISCRIMINATE
    # from is the point of the design. Applying this rule un-scoped would not
    # over-reject at the margins, it would reject precisely the best
    # strain-unique windows. ``None`` (no scope) therefore means inert, never
    # "assume every record".
    #
    # It resolves per CANDIDATE, not per run: candidate_on_target_records takes
    # the window's own covered group when it has one and the run-level
    # ``on_target_records`` otherwise. Only constitutive mode sets a per-window
    # group; see on_target_record_ids for the per-mode run-level mapping.
    #
    # NOTE ON REACH: on most paths this is a backstop that does not fire.
    # Conserved windows are gap-free and either identical across all records or
    # carry their variant positions outside this band already
    # (consensus.degenerate_conserved_windows rules 1 and 3, applied with this
    # same clearance); a strain-unique probe is its own target's literal
    # sequence, so it cannot mismatch it; a default-quorum constitutive window
    # is shared verbatim by every record. On those paths this closes the case
    # the module docstring above already flags for the gap rule — "nothing
    # downstream will catch it if that rule is relaxed". Window selection is
    # not the only thing standing between a junction-proximal mismatch and an
    # ordered oligo.
    #
    # Constitutive mode below its quorum is no exception to that: its covered
    # group is grouped BY exact window sequence, so every member matches the
    # probe byte for byte and cannot contribute a mismatch of any kind. The
    # check still runs there, through the per-candidate scope, because it is
    # the one mode whose window selection applies no junction rule of its own
    # (constitutive_start_positions picks the LARGEST group; it never asks
    # where the isoforms disagree). All four modes therefore fail the same way
    # if window selection ever stops guaranteeing what it currently
    # guarantees.
    #
    # A candidate with no ``strain_binding`` at all passes: there is no evidence
    # either way, and the "cannot judge -> pass" convention used for Tm below is
    # the right one here too. On a multi-strain run the pipeline annotates every
    # candidate before it reaches this funnel, so an unannotated candidate means
    # the caller is not on that path — not that the probe is clean.
    # ``strain_binding`` is tested before the scope is resolved so that a
    # single-sequence run — which is every whole-transcriptome run, millions of
    # tiles — short-circuits on an empty list and never makes the call.
    if (thresholds.reject_junction_mismatch and pc.strain_binding
            and n >= 2 and half > 0):
        scope = candidate_on_target_records(pc, on_target_records)
        j_lo, j_hi = _junction_band(thresholds, n, half)
        if scope and j_hi > j_lo:
            for entry in pc.strain_binding:
                if entry.get('strain_id', '') not in scope:
                    continue
                positions = entry.get('mismatch_positions') or ()
                if any(j_lo <= int(p) < j_hi for p in positions):
                    reasons.append('junction_mismatch')
                    counters['junction_mismatch'] += 1
                    break   # one disqualifying record is enough; count once

    # TN-junction filter. ``seq_u[half]`` is the mRNA base opposite the LHS
    # probe's 3' (ligation) end; when it is 'A' the LHS oligo ends in T, giving
    # the TN motif. Which bases *fail* is set by ``thresholds.tn_junction_mode``,
    # resolved upstream by chemistry.resolve_tn_mode (citations in
    # constants.TN_MODE_*):
    #
    #   strict_a — HARD, and the default wherever the gate is engaged: reject
    #              anything but 'A'. CG000621 states it as a requirement, and
    #              CG000839 still recommends it "for a simplified and robust
    #              design process".
    #   soft     — explicit opt-in, Flex v1/v2 only: reject only 'C'. CG000839
    #              says the TN motif "is not always strictly required, and
    #              other nucleotide motifs can also function effectively", and
    #              10x's own panels put A there only ~51% of the time while
    #              avoiding C ~99% of the time. Under this mode 'A' is a
    #              preference expressed in ranking, not a gate.
    #   off      — no junction-identity gate (the long-standing CLI default).
    #
    # ``junction_base`` and ``tn_junction_ok`` are recorded on EVERY candidate
    # whatever the mode, so the ranking layer can always prefer a true TN
    # junction — that is what keeps A-junction probes at the top under `soft`.
    if n > half:
        junction_base = seq_u[half]
        pc.junction_base = junction_base
        pc.tn_junction_ok = junction_base == 'A'
        mode = getattr(thresholds, 'tn_junction_mode', TN_MODE_AUTO)
        if mode == TN_MODE_AUTO:
            # Unresolved (a bare FilterThresholds, e.g. a direct apply_filters
            # call in a test or a downstream tool). Fall back to the legacy
            # switch, which is exactly what this code did before the mode
            # existed — never to a gate the caller did not ask for.
            mode = (TN_MODE_STRICT_A
                    if getattr(thresholds, 'require_tn_junction', True)
                    else TN_MODE_OFF)
        failed = (
            (mode == TN_MODE_STRICT_A and junction_base != 'A')
            or (mode == TN_MODE_SOFT and junction_base == 'C')
        )
        if failed:
            reasons.append('tn_junction')
            counters['tn_junction'] += 1

    # Degenerate-base budget per half. Counts EVERY non-ACGT base, not just a
    # literal 'N': the consensus builder emits an IUPAC code wherever the
    # disagreeing bases are all in ACGT and only falls back to 'N' for a gap,
    # so an 'N'-only count would let every two-way code (the ones that dominate
    # a real strain alignment) through unmeasured. The threshold keeps the
    # name `max_n_per_half` for continuity with the CLI flag and with existing
    # run_config.json files; at its default of 1 the semantics are
    # <=1 degenerate position per 25 nt half, so <=2 per 50-mer.
    if degenerate_count(left) > thresholds.max_n_per_half:
        reasons.append('n_count_left')
        counters['n_count'] += 1
    elif degenerate_count(right) > thresholds.max_n_per_half:
        reasons.append('n_count_right')
        counters['n_count'] += 1

    # GC per half — count G/C once per half. ALWAYS computed (ranking and
    # results.csv use it), but only a *gate* when explicitly asked for.
    #
    # 10x publishes GC 44-72% per half (CG000839 p.2 / CG000621 p.2, p.12). As
    # a hard filter it is off by default because 10x does not honour it
    # themselves: 8,940 of the 106,918 halves in their public human
    # v1.1.0 panel (8.36%) fall outside the band, down to GC 0.28 — while not
    # one of those halves falls below 56.32 C. The band is a GC proxy for duplex
    # stability written for a mammalian transcriptome (median half GC 0.52); on
    # AT-rich targets it rejects thermodynamically fine windows wholesale. The
    # Tm gate below is the real constraint; GC survives as a ranking term.
    gc_left_count = left.count('G') + left.count('C')
    gc_right_count = right.count('G') + right.count('C')
    gc_l = gc_left_count / len(left) if left else 0.0
    gc_r = gc_right_count / len(right) if right else 0.0
    pc.gc_content = (gc_left_count + gc_right_count) / n if n else 0.0
    # A consensus parent's own string has no meaningful GC — an 'S' is G+C with
    # certainty and counts 0, an 'N' dilutes the half by 4% — so it reports the
    # midpoint of its variants' GC range instead, and is GATED on the range's
    # ends. Every variant goes into the same synthesis pool, so one variant
    # outside the band puts the pool outside the band.
    is_parent = getattr(pc, 'probe_role', PROBE_ROLE_PLAIN) == PROBE_ROLE_CONSENSUS
    gc_l_lo, gc_l_hi, gc_r_lo, gc_r_hi = gc_l, gc_l, gc_r, gc_r
    if is_parent and not is_nan(pc.gc_lhs_min) and not is_nan(pc.gc_rhs_min):
        gc_l_lo, gc_l_hi = pc.gc_lhs_min, pc.gc_lhs_max
        gc_r_lo, gc_r_hi = pc.gc_rhs_min, pc.gc_rhs_max
        gc_l = (gc_l_lo + gc_l_hi) / 2.0
        gc_r = (gc_r_lo + gc_r_hi) / 2.0
        pc.gc_content = (gc_l + gc_r) / 2.0
    pc.gc_left = gc_l
    pc.gc_right = gc_r
    if getattr(thresholds, 'enforce_gc_range', True):
        if not (thresholds.gc_min <= gc_l_lo and gc_l_hi <= thresholds.gc_max):
            reasons.append('gc_range_left')
            counters['gc_range'] += 1
        elif not (thresholds.gc_min <= gc_r_lo and gc_r_hi <= thresholds.gc_max):
            reasons.append('gc_range_right')
            counters['gc_range'] += 1

    # Per-half Tm gate — the duplex-stability filter that replaces the GC band.
    #
    # *** RELATIVE, NOT AN IN-ASSAY MELTING TEMPERATURE. *** primer3 calc_tm
    # under the PCR-primer salt conditions in constants.py (10x does not publish
    # the Flex hybridisation buffer). The floor is derived under those same
    # conditions from 10x's own shipped panels, so gate and value are in the
    # same units and the comparison is meaningful — but the number must not be
    # read as the temperature at which these probes melt in the Flex assay.
    #
    # Computed after the cheap character-level checks so a window already
    # rejected for gaps/N/homopolymer costs no primer3 call. NaN (primer3
    # missing, or an ambiguity code it will not score) means "cannot judge", and
    # an unjudgeable window is passed rather than silently rejected — the same
    # honest-sentinel convention as _utils.safe_tm and results.csv.
    # *** A DEGENERATE PARENT MUST NEVER REACH safe_tm. *** primer3 raises on
    # any non-ACGT input, so safe_tm returns NaN for a mixed-base half, and the
    # "cannot judge -> pass" convention below would wave it straight through.
    # Since the GC band was demoted to a ranking term, Tm is THE duplex-
    # stability gate — an unjudged parent is an ungated parent. So a parent is
    # gated on the RANGE its concrete variants span (populated at tile time by
    # `pipeline._degenerate_stability`): if any variant's half falls outside the
    # band, the whole pool does, because they are all synthesised into it.
    tm_min = getattr(thresholds, 'tm_min', None)
    tm_max = getattr(thresholds, 'tm_max', None)
    tm_gated = tm_min is not None or tm_max is not None
    if is_parent:
        if pc.n_variants <= 0:
            # A consensus parent with no variants cannot be judged by anything
            # and cannot be resolved into an orderable oligo either. That is a
            # fan-out bug, so it fails closed rather than riding through on the
            # NaN convention.
            reasons.append('degenerate_unresolved')
            counters['degenerate_unresolved'] += 1
        elif tm_gated:
            ranges = (
                ('left', pc.tm_lhs_min, pc.tm_lhs_max),
                ('right', pc.tm_rhs_min, pc.tm_rhs_max),
            )
            if all(is_nan(lo) for _s, lo, _hi in ranges) and _primer3_available():
                # primer3 could have scored these variants — it scores plain
                # ACGT — so a missing aggregate means the fan-out never ran.
                reasons.append('degenerate_unresolved')
                counters['degenerate_unresolved'] += 1
            else:
                for side, lo, hi in ranges:
                    if is_nan(lo) or is_nan(hi):
                        continue
                    if (tm_min is not None and lo < tm_min) or \
                            (tm_max is not None and hi > tm_max):
                        reasons.append(f'tm_range_{side}')
                        counters['tm_range'] += 1
                        break
    elif tm_gated and left and right:
        for side, seq_half in (('left', left), ('right', right)):
            tm = safe_tm(seq_half)
            if is_nan(tm):
                continue
            if (tm_min is not None and tm < tm_min) or \
                    (tm_max is not None and tm > tm_max):
                reasons.append(f'tm_range_{side}')
                counters['tm_range'] += 1
                break

    # Homopolymer run
    hp = max_homopolymer(seq_u)
    pc.homopolymer_max = hp
    if hp > thresholds.max_homopolymer_run:
        reasons.append('homopolymer')
        counters['homopolymer'] += 1
