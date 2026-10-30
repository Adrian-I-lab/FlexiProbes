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

"""Greedy spaced top-N probe selection per gene.

Used by the GUI's whole-transcriptome / auto-mode and by the bulk-export
project flow when a run has no user-supplied selections.

Spacing model
-------------
Three gaps, in decreasing order of ambition:

``spread`` (adaptive, see :func:`spread_gap_for`)
    What an unconstrained gene should give you: the probes divided evenly over
    the middle ``1 - 2*END_EXCLUSION_FRAC`` of the transcript. Scales with the
    gene, so three probes on a 12 kb mRNA land kilobases apart instead of merely
    clearing a fixed rule. Never costs a probe — see the search below.
``MIN_PROBE_GAP_BP`` (150)
    The spacing the design targets. Anything below it is reported.
``RELAXED_PROBE_GAP_BP`` (40)
    The floor used only when a gene cannot yield ``PROBE_GAP_RELAX_MIN_PROBES``
    probes at 150. Genes that land here are flagged, not silently accepted.

The greedy walk is *monotone in the gap*: it visits candidates in a fixed
score order and rejects one only for being too close to an already-kept probe,
so lowering the gap can never lose a pick that a higher gap kept. That is what
makes both the binary search for the widest workable spread and the "relax only
if it actually buys a probe" fallback sound, rather than heuristics that might
trade one probe for another.

Compared with a plain fixed-gap greedy walk, the adaptive spread raises the
median closest gap and the share of the middle 70% that probes cover, most of
all on long transcripts. The cost is a small loss of probes, almost all on
short transcripts that cannot fit three probes at 40 bp. Spreading probes out
puts somewhat more of them in the terminal 15% of the transcript, but the share
stays below the one a placement-blind design gives, and below the share in
10x's own published human panel (about 26%).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import (
    Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple,
)

from .constants import (
    END_EXCLUSION_FRAC, END_EXCLUSION_FRAC_RELAXED, GAP_TIER_GREEN,
    GAP_TIER_NAMES, MIN_PROBE_GAP_BP, MIN_PROBE_GAP_FLOOR_BP,
    PROBE_GAP_RELAX_MIN_PROBES, PROBE_ROLE_CONSENSUS, PROBE_ROLE_VARIANT,
    RELAXED_PROBE_GAP_BP, gap_tier,
)


#: The spacing ladder, best rung first: ``(gap floor bp, terminal fraction
#: avoided at each end)``. Selection stops at the first rung that yields
#: ``PROBE_GAP_RELAX_MIN_PROBES`` probes, so a gene pays only for the tightness
#: it actually needs. Gap and terminal band drop TOGETHER on the last rung
#: because a gene reaches it only when its transcript cannot hold three
#: separated probes — on such a transcript the middle 70% is usually the binding
#: constraint, and relaxing the gap alone would buy nothing.
SPACING_LADDER: Tuple[Tuple[int, float], ...] = (
    (MIN_PROBE_GAP_BP, END_EXCLUSION_FRAC),              # green  150 / 15%
    (RELAXED_PROBE_GAP_BP, END_EXCLUSION_FRAC),          # orange   40 / 15%
    (MIN_PROBE_GAP_FLOOR_BP, END_EXCLUSION_FRAC_RELAXED),  # red    15 / 10%
)


def ladder_from(min_gap: int = MIN_PROBE_GAP_BP,
                relaxed_gap: int = RELAXED_PROBE_GAP_BP,
                floor_gap: int = MIN_PROBE_GAP_FLOOR_BP,
                end_frac: float = END_EXCLUSION_FRAC,
                end_frac_relaxed: float = END_EXCLUSION_FRAC_RELAXED,
                ) -> Tuple[Tuple[int, float], ...]:
    """Build a ladder from explicit gaps, dropping rungs that cannot help.

    A rung whose gap is not strictly tighter than the one above it is removed
    rather than run: it would return the identical pick (the walk is
    deterministic in the gap) and then report a relaxation that never happened.
    So ``--relaxed-gap 150`` (equal to the default min) drops that rung only:
    the floor rung, 15 bp by default, is still tighter than 150 and still runs.
    No fallback at all needs both relaxed and floor >= min. Dropping a rung is
    also the mechanism behind the GUI's spacing filter: excluding the red tier
    is just dropping its rung.
    """
    rungs = [(min_gap, end_frac)]
    for gap, frac in ((relaxed_gap, end_frac), (floor_gap, end_frac_relaxed)):
        if gap < rungs[-1][0]:
            rungs.append((gap, frac))
    return tuple(rungs)


def _gap(a: tuple, b: tuple) -> int:
    """Gap in bp between two ``[start, end)`` intervals (negative = overlap)."""
    (s1, e1), (s2, e2) = a, b
    if e1 <= s2:
        return s2 - e1
    if e2 <= s1:
        return s1 - e2
    return -1


def spread_gap_for(gene_length: int, n_per_gene: int,
                   probe_length: int = 50,
                   min_gap: int = MIN_PROBE_GAP_BP,
                   end_frac: float = END_EXCLUSION_FRAC) -> int:
    """The gap ``n_per_gene`` probes would sit at if spread evenly over the
    usable middle of a ``gene_length`` transcript.

    "Usable" is the middle ``1 - 2*end_frac`` — the same span
    ``ranking._central_placement_score`` scores flat — because a spread target
    that counted the terminal 15% would aim probes at positions the ranking is
    simultaneously pushing them away from, and the two would just fight.

    This is a CEILING derived from the transcript, floored only at ``min_gap``.
    It deliberately does NOT floor at ``PREFERRED_PROBE_GAP_BP``: on a 1 kb
    transcript, three 50-mers at 300 bp gaps need 750 nt, i.e. 75% of it, so a
    300 bp floor would *force* probes into the terminal 15% that the placement
    score is pushing them out of. The share of chosen probes in the terminal
    15% would then be about the uniform rate, as if the end-avoidance term were
    switched off. Deriving the target from the transcript instead keeps a
    3-probe layout inside the middle 70% for every transcript long enough to
    hold one.

    ``min_gap`` still wins when the even-spread figure falls below it: the
    design target outranks the aspiration, and a 3-probe layout at 150 bp needs
    only 450 nt, which fits the middle 70% of anything from ~650 nt up.

    This is only ever an *aspiration* — :func:`auto_pick_with_report` lowers it
    as far as needed rather than dropping a probe for it.
    """
    if gene_length <= 0 or n_per_gene <= 1:
        return min_gap
    usable = gene_length * (1.0 - 2.0 * end_frac)
    # n probes occupy n*probe_length; the n-1 gaps share what is left.
    slack = usable - n_per_gene * probe_length
    if slack <= 0:
        return min_gap
    return max(min_gap, int(slack // (n_per_gene - 1)))


@dataclass
class GeneSpacing:
    """What spacing a single gene actually ended up with."""
    gene_id: str
    n_probes: int
    min_gap: Optional[int]       # closest chosen pair in bp; None if <2 probes
    gap_used: int                # the rung's gap floor
    end_frac_used: float         # the rung's terminal fraction
    target_gap: int              # the gap the top rung asked for
    relaxed: bool                # a rung below the top one won
    gene_length: int = 0

    @property
    def tier(self) -> int:
        """GREEN / ORANGE / RED, from the spacing actually achieved.

        Read off ``min_gap`` rather than off the rung that produced it: a gene
        can settle on the orange rung and still end up with 400 bp between its
        probes, because the rung is a floor and the search reaches above it.
        Reporting such a gene as orange would overstate how crowded the panel
        is. The tier a user sees is always the one their probes earn.
        """
        return gap_tier(self.min_gap)

    @property
    def tight(self) -> bool:
        """True when some chosen pair sits closer than the design target."""
        return self.min_gap is not None and self.min_gap < self.target_gap


@dataclass
class AutoPick:
    """Selections plus why each gene's spacing is what it is."""
    selections: Dict[str, List[str]] = field(default_factory=dict)
    spacing: Dict[str, GeneSpacing] = field(default_factory=dict)

    @property
    def relaxed_genes(self) -> List[str]:
        """Genes whose probes had to be packed below ``MIN_PROBE_GAP_BP``."""
        return sorted(g for g, s in self.spacing.items() if s.tight)

    def genes_by_tier(self) -> Dict[int, List[str]]:
        out: Dict[int, List[str]] = {}
        for g, sp in self.spacing.items():
            out.setdefault(sp.tier, []).append(g)
        return {k: sorted(v) for k, v in out.items()}

    def tier_counts(self) -> Dict[str, int]:
        """Probes per spacing tier, by name. Counts PROBES, not genes.

        Probes because that is the unit a user orders and pays for, and because
        gene counts hide the shape of the problem: 200 crowded genes carrying
        two probes each is a very different panel from 200 carrying five. The
        three counts always sum to the probe total.
        """
        counts = {name: 0 for name in GAP_TIER_NAMES.values()}
        for gid, ids in self.selections.items():
            sp = self.spacing.get(gid)
            tier = sp.tier if sp else GAP_TIER_GREEN
            counts[GAP_TIER_NAMES[tier]] += len(ids)
        return counts

    def summary(self) -> Dict[str, int]:
        tight = self.relaxed_genes
        out: Dict[str, int] = {
            "genes": len(self.selections),
            "probes": sum(len(v) for v in self.selections.values()),
            "genes_below_min_gap": len(tight),
            "min_gap_bp": MIN_PROBE_GAP_BP,
            "relaxed_gap_bp": RELAXED_PROBE_GAP_BP,
            "floor_gap_bp": MIN_PROBE_GAP_FLOOR_BP,
        }
        for name, n in self.tier_counts().items():
            out[f"probes_{name}"] = n
        by_tier = self.genes_by_tier()
        for tier, name in GAP_TIER_NAMES.items():
            out[f"genes_{name}"] = len(by_tier.get(tier, ()))
        return out


def _walk(ranked: List[Tuple[Mapping[str, str], Tuple[int, int]]],
          n_per_gene: int, gap: int,
          variants_of: Mapping[str, List[str]]) -> Tuple[List[str], List[tuple]]:
    """One greedy pass at a fixed ``gap``. Returns (probe ids, footprints).

    Candidates are walked in the caller's score order; each is kept only if it
    sits at least ``gap`` bp from every already-kept probe (which also rules out
    overlap). Monotone in ``gap`` — see the module docstring.
    """
    picked_ids: List[str] = []
    picked_iv: List[tuple] = []
    for r, (s, e) in ranked:
        if any(_gap((s, e), iv) < gap for iv in picked_iv):
            continue
        pid = r.get("probe_id", "")
        if r.get("probe_role") == PROBE_ROLE_CONSENSUS:
            kids = variants_of.get(pid, [])
            if not kids:
                continue        # unorderable and nothing to order for it
            picked_ids.extend(kids)
        else:
            picked_ids.append(pid)
        picked_iv.append((s, e))
        # One family = one pick, so count SITES, not emitted ids.
        if len(picked_iv) >= n_per_gene:
            break
    return picked_ids, picked_iv


def _best_spread(pool, n_per_gene: int, min_gap: int, length: int,
                 probe_len: int, variants_of, end_frac: float
                 ) -> Tuple[List[str], List[tuple], int]:
    """Widest gap in ``[min_gap, spread_gap_for(...)]`` that still fills the quota.

    Binary search, exact because the walk is monotone in the gap (module
    docstring): if a gap fills the quota then every smaller gap does too, so the
    predicate "fills the quota" is a step function and bisection finds its edge.
    Falls back to the ``min_gap`` result when the quota is unreachable at any
    gap, which is the answer callers want anyway — as many probes as exist.
    """
    if not pool:
        return [], [], min_gap
    ids, ivs = _walk(pool, n_per_gene, min_gap, variants_of)
    gap_used = min_gap
    if len(ivs) < n_per_gene:
        return ids, ivs, gap_used
    lo = min_gap
    hi = spread_gap_for(length, n_per_gene, probe_len, min_gap, end_frac)
    while lo < hi:                       # invariant: lo fills the quota
        mid = (lo + hi + 1) // 2
        t_ids, t_ivs = _walk(pool, n_per_gene, mid, variants_of)
        if len(t_ivs) >= n_per_gene:
            lo, ids, ivs, gap_used = mid, t_ids, t_ivs, mid
        else:
            hi = mid - 1
    return ids, ivs, gap_used


def _closest_gap(ivs: List[tuple]) -> Optional[int]:
    """Smallest gap between any two footprints, or None when there is no pair.

    None rather than a sentinel number: a one-probe gene has no spacing, which
    is a different statement from "its spacing is very large", and only the
    first is true.
    """
    if len(ivs) < 2:
        return None
    ivs = sorted(ivs)
    return min(_gap(ivs[i], ivs[i + 1]) for i in range(len(ivs) - 1))


def auto_pick_with_report(
    rows: Iterable[Mapping[str, str]],
    *,
    n_per_gene: int,
    score_key: str = "composite_score",
    min_gap: int = MIN_PROBE_GAP_BP,
    relaxed_gap: int = RELAXED_PROBE_GAP_BP,
    floor_gap: int = MIN_PROBE_GAP_FLOOR_BP,
    relax_below: int = PROBE_GAP_RELAX_MIN_PROBES,
    end_frac: float = END_EXCLUSION_FRAC,
    end_frac_relaxed: float = END_EXCLUSION_FRAC_RELAXED,
    ladder: Optional[Sequence[Tuple[int, float]]] = None,
    gene_length_of: Optional[Callable[[str], int]] = None,
) -> AutoPick:
    """Spaced top-N pick per gene, with the spacing each gene settled at.

    Each gene walks :data:`SPACING_LADDER` — ``(150 bp, 15%)``, ``(40 bp, 15%)``,
    ``(15 bp, 10%)`` — and STOPS at the first rung yielding
    ``min(n_per_gene, relax_below)`` probes. A gene therefore only ever pays for
    the tightness it needs, and the rung it stopped on is recorded.

    Within one rung, in order:

    1. **Spread, interior first.** Restrict candidates to the middle
       ``1 - 2*end_frac`` of the transcript and binary-search the largest gap in
       ``[gap, spread_gap_for(...)]`` that still yields ``n_per_gene`` probes. On
       a roomy transcript that gap is far above the rung's floor and is what
       stops the top-scoring windows from bunching; on a tight one it collapses
       back to the floor at no cost. Monotonicity (module docstring) is what
       makes the search exact rather than approximate.

       The interior restriction is not decoration: maximising a minimum gap, by
       itself, pushes the outer probes to the extreme ends of the transcript.
       Without it the same run put 29% of chosen probes in the terminal 15% —
       the rate you get from ignoring placement entirely.
    2. **Widen.** If the interior cannot fill the quota, repeat over every
       candidate at the same gap. Ends are strongly discouraged, not forbidden;
       a gene whose only spaced probes sit near an end still gets them.

    A rung is kept only if it produced MORE probes than the best rung so far, so
    descending can add probes but never trade one away, and a gene never reports
    a relaxation that bought it nothing.

    Pass ``ladder`` to override the rungs — this is how the GUI's spacing filter
    works: excluding the red tier is passing a ladder without its rung, so those
    probes are never candidates rather than being picked and then filtered out.
    :func:`ladder_from` builds one from explicit gap values.

    ``gene_length_of`` maps a gene id to its transcript length; without it the
    length is taken as the furthest probe end seen for that gene, which
    under-reads a transcript whose 3' end produced no passing probe and so only
    ever makes the spread target more conservative, never less.

    **A degenerate family counts once.** A consensus parent and its variants are
    one probe at one site, not N+1 probes: they occupy the same 50 bp, they are
    ordered together and they are reviewed together. Counting each row against
    the cap would spend a "top 5 per gene" budget on two or three sites. So the
    parent row is the cap unit — it carries the spacing and the score — and the
    ids returned for it are its VARIANTS', because the parent itself is not
    orderable. Variant rows are never walked on their own.
    """
    if n_per_gene <= 0:
        return AutoPick()
    rungs = tuple(ladder) if ladder else ladder_from(
        min_gap, relaxed_gap, floor_gap, end_frac, end_frac_relaxed)
    if not rungs:
        return AutoPick()
    target_gap = rungs[0][0]

    rows = list(rows)
    # parent probe_id -> its variants' probe_ids, in row order.
    variants_of: Dict[str, List[str]] = {}
    for r in rows:
        parent = r.get("parent_probe_id") or ""
        if parent and r.get("probe_role") == PROBE_ROLE_VARIANT:
            variants_of.setdefault(parent, []).append(r.get("probe_id", ""))

    by_gene: Dict[str, List[Mapping[str, str]]] = {}
    for r in rows:
        if r.get("probe_role") == PROBE_ROLE_VARIANT:
            continue        # represented by its parent
        by_gene.setdefault(r.get("gene_id", ""), []).append(r)

    def _score(r: Mapping[str, str]) -> float:
        try:
            return float(r.get(score_key, 0) or 0)
        except (TypeError, ValueError):
            return 0.0

    def _interval(r: Mapping[str, str]) -> tuple:
        try:
            return int(r.get("start", 0)), int(r.get("end", 0))
        except (TypeError, ValueError):
            return 0, 0

    out = AutoPick()
    for gid, cands in by_gene.items():
        ranked = []
        for r in sorted(cands, key=lambda r: -_score(r)):
            s, e = _interval(r)
            if e > s:
                ranked.append((r, (s, e)))
        if not ranked:
            continue

        if gene_length_of is not None:
            length = int(gene_length_of(gid) or 0)
        else:
            length = max(e for _, (_s, e) in ranked)
        probe_len = min(e - s for _, (s, e) in ranked)

        quota = min(n_per_gene, relax_below)
        ids: List[str] = []
        ivs: List[tuple] = []
        gap_used, frac_used = rungs[0]
        for rung_gap, rung_frac in rungs:
            # Interior first. Maximising the minimum gap on its own drives the
            # outer probes to the very ends of the transcript — that is what
            # maximising a minimum does — so the search has to be told where it
            # may not go. Without this step the share of chosen probes in the
            # terminal 15% is about what you get by ignoring placement
            # altogether. Widening only when the interior
            # cannot fill the quota is the "strongly discouraged, not forbidden"
            # reading: a gene whose only spaced probes sit near an end gets them.
            interior = [
                (r, iv) for r, iv in ranked
                if rung_frac * length <= (iv[0] + iv[1]) / 2.0
                <= (1 - rung_frac) * length
            ] if length > 0 else []
            r_ids, r_ivs, r_gap = _best_spread(
                interior, n_per_gene, rung_gap, length, probe_len, variants_of,
                rung_frac)
            if len(r_ivs) < n_per_gene:
                w_ids, w_ivs, w_gap = _best_spread(
                    ranked, n_per_gene, rung_gap, length, probe_len, variants_of,
                    rung_frac)
                if len(w_ivs) > len(r_ivs):
                    r_ids, r_ivs, r_gap = w_ids, w_ivs, w_gap
            # Descend only for probes we did not already have. Without this a
            # looser rung could swap a well-spaced pick for an equally-sized but
            # tighter one, and the gene would be reported as relaxed for nothing.
            if len(r_ivs) > len(ivs):
                ids, ivs, gap_used, frac_used = r_ids, r_ivs, r_gap, rung_frac
            if len(ivs) >= quota:
                break

        if not ids:
            continue
        out.selections[gid] = ids
        out.spacing[gid] = GeneSpacing(
            gene_id=gid, n_probes=len(ivs), min_gap=_closest_gap(ivs),
            gap_used=gap_used, end_frac_used=frac_used, target_gap=target_gap,
            relaxed=gap_used < target_gap, gene_length=length,
        )
    return out


def auto_pick_per_gene(
    rows: Iterable[Mapping[str, str]],
    *,
    n_per_gene: int,
    score_key: str = "composite_score",
    min_gap: int = MIN_PROBE_GAP_BP,
    relaxed_gap: int = RELAXED_PROBE_GAP_BP,
    floor_gap: int = MIN_PROBE_GAP_FLOOR_BP,
    relax_below: int = PROBE_GAP_RELAX_MIN_PROBES,
    end_frac: float = END_EXCLUSION_FRAC,
    end_frac_relaxed: float = END_EXCLUSION_FRAC_RELAXED,
    ladder: Optional[Sequence[Tuple[int, float]]] = None,
    gene_length_of: Optional[Callable[[str], int]] = None,
) -> Dict[str, List[str]]:
    """Return ``{gene_id: [probe_id, ...]}``.

    Thin wrapper over :func:`auto_pick_with_report` for callers that only want
    the picks. Use the report form if you need to know which genes had to be
    packed below ``min_gap`` — that is the fact worth surfacing to a user, and
    it is invisible in this return type.
    """
    return auto_pick_with_report(
        rows, n_per_gene=n_per_gene, score_key=score_key, min_gap=min_gap,
        relaxed_gap=relaxed_gap, floor_gap=floor_gap, relax_below=relax_below,
        end_frac=end_frac, end_frac_relaxed=end_frac_relaxed, ladder=ladder,
        gene_length_of=gene_length_of,
    ).selections
