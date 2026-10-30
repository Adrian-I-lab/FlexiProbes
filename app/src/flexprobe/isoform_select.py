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

"""Coverage-maximising probe selection per gene.

The default :func:`flexprobe.auto_select.auto_pick_per_gene` picks the top-N
highest-scoring probes per gene, spaced apart. This module instead picks up to
N spaced probes per gene that, *between them*, bind as many of the gene's
isoforms as physically possible — a greedy weighted set-cover — and reports the
genes where full coverage (or N spaced probes) cannot be achieved.

It is used by the desktop review dialog ("Maximise isoform coverage" action)
and by the headless ``run_manager`` auto-pick path.

Why it needs the input FASTA
----------------------------
``results.csv`` keys probes by the *transcript* they were tiled on (the
``gene_id`` column is an ENST id) and records only the *count* of sibling
isoforms each probe binds — not which ones. To group probes by their real gene
(ENSG) and run a genuine "cover every isoform" optimisation we recover, for each
candidate, the exact set of the gene's transcripts whose sequence contains the
probe's 50-bp target. That comes from the run's input cDNA FASTA, whose Ensembl
headers carry ``gene:ENSG...``.

Coverage model: a probe binds an isoform when the isoform's sequence contains
the probe's ``mrna_target_50bp`` (exact match by default; ``max_mismatches`` >0
loosens this toward the engine's >=90%-identity BLAST rule). On sibling isoforms
that share identical exonic sequence this reproduces the engine's BLAST-based
``isoform_hits`` exactly; on divergent isoforms it is deliberately conservative.

Spacing is enforced in *sequence space*: two probes (which may come from
different source transcripts) conflict only if, on an isoform they both bind,
their footprints overlap or sit closer than ``min_gap`` bp.
"""
from __future__ import annotations

import heapq
import re
from dataclasses import dataclass, field
from typing import (
    Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple,
)

from .constants import (
    GAP_TIER_GREEN, GAP_TIER_NAMES, MIN_PROBE_GAP_BP, MIN_PROBE_GAP_FLOOR_BP,
    PREFERRED_PROBE_GAP_BP, PROBE_GAP_RELAX_MIN_PROBES, PROBE_ROLE_CONSENSUS,
    PROBE_ROLE_VARIANT,
    RELAXED_PROBE_GAP_BP, gap_tier,
)

_GENE_RE = re.compile(r"gene:(\S+)")
_SYMBOL_RE = re.compile(r"gene_symbol:(\S+)")


# --------------------------------------------------------------------------- #
# Gene index (parsed from the run's input FASTA)
# --------------------------------------------------------------------------- #
@dataclass
class GeneIndex:
    """ENST->ENSG mapping, each gene's transcript set, sequences, symbols."""
    gene_of: Dict[str, str]               # transcript id -> gene key
    group_members: Dict[str, Set[str]]    # gene key -> {transcript ids}
    seqs: Dict[str, str]                  # transcript id -> uppercased sequence
    symbol_of: Dict[str, str]             # gene key -> gene symbol (if known)

    def gene_key(self, transcript: str) -> str:
        return self.gene_of.get(transcript, transcript)


def build_gene_index(fasta_path: str,
                     needed_transcripts: Optional[Set[str]] = None) -> GeneIndex:
    """Parse the cDNA FASTA into a :class:`GeneIndex`.

    Headers are read for every record (cheap) to build the ENST->ENSG map and
    each gene's full isoform set. Sequences are loaded only for transcripts
    whose *gene* is referenced by ``needed_transcripts`` — i.e. the genes that
    produced candidate probes — so peak RAM tracks the selected genes, not the
    whole transcriptome. Pass ``needed_transcripts=None`` to load everything.
    """
    gene_of: Dict[str, str] = {}
    symbol_of: Dict[str, str] = {}
    group_members: Dict[str, Set[str]] = {}

    # Pass 1: headers -> gene_of / group_members / symbol_of.
    with open(fasta_path) as fh:
        for line in fh:
            if not line or line[0] != ">":
                continue
            tid = line[1:].split(None, 1)[0]
            m = _GENE_RE.search(line)
            gkey = m.group(1) if m else tid
            gene_of[tid] = gkey
            group_members.setdefault(gkey, set()).add(tid)
            sm = _SYMBOL_RE.search(line)
            if sm:
                symbol_of[gkey] = sm.group(1)

    # Which transcripts' sequences do we actually need? Every sibling of a gene
    # that produced a candidate (so coverage and spacing see the whole gene).
    if needed_transcripts is None:
        wanted: Optional[Set[str]] = None
    else:
        wanted = set()
        for t in needed_transcripts:
            wanted.update(group_members.get(gene_of.get(t, t), {t}))

    seqs = _load_sequences(fasta_path, wanted)
    return GeneIndex(gene_of, group_members, seqs, symbol_of)


def _load_sequences(path: str, wanted: Optional[Set[str]]) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    cur: Optional[str] = None
    keep = True
    buf: List[str] = []
    with open(path) as fh:
        for line in fh:
            if line and line[0] == ">":
                if cur is not None and keep:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split(None, 1)[0]
                keep = wanted is None or cur in wanted
                buf = []
            elif keep:
                buf.append(line.strip().upper())
        if cur is not None and keep:
            seqs[cur] = "".join(buf)
    return seqs


# --------------------------------------------------------------------------- #
# Candidate / result records
# --------------------------------------------------------------------------- #
@dataclass
class _Cand:
    probe_id: str
    transcript: str
    start: int
    end: int
    target: str
    score: float
    hit_set: Set[str] = field(default_factory=set)
    footprint: Dict[str, int] = field(default_factory=dict)
    # The SITE this candidate occupies, which is not its id. Every variant child
    # of a degenerate consensus window is built with the PARENT's own start/end
    # (see the child construction in ``flexprobe.pipeline``), so a family's
    # siblings are one probe at one place. Defaults to the probe's own id, which
    # makes every ordinary candidate its own site and leaves the non-degenerate
    # path byte-identical.
    site_key: str = ""
    # True when ``hit_set`` came from the persisted ``isoform_hit_ids`` column
    # (no FASTA needed); False means it must be recomputed from sequence.
    from_column: bool = False


def gene_index_from_manifest(gene_isoforms: Mapping[str, Iterable[str]],
                             symbol_of: Optional[Mapping[str, str]] = None
                             ) -> "GeneIndex":
    """Build a sequence-free :class:`GeneIndex` from a manifest's
    ``gene_isoforms`` map (gene -> its transcript ids). Used for the FASTA-free
    selection path, where membership comes from the ``isoform_hit_ids`` column
    and only the gene->isoform universe is needed for flagging."""
    gene_of: Dict[str, str] = {}
    group_members: Dict[str, Set[str]] = {}
    for gene, ts in gene_isoforms.items():
        members = set(ts)
        group_members[gene] = members
        for t in members:
            gene_of[t] = gene
    return GeneIndex(gene_of, group_members, {}, dict(symbol_of or {}))


@dataclass
class GeneFlag:
    """One gene that could not be fully satisfied."""
    gene: str
    isoforms_total: int
    probes_selected: int
    isoforms_covered: int
    missing_isoforms: List[str]
    reasons: List[str]
    # Where the isoform universe (the coverage_fraction denominator) came from:
    #   "authoritative" — the gene's full isoform set is known (manifest/FASTA
    #                      map), so the fraction is trustworthy.
    #   "floor"         — the set was inferred from the candidates only, so the
    #                      true isoform count may be higher and the fraction is
    #                      optimistic (an upper bound).
    coverage_basis: str = "authoritative"
    # Spacing actually achieved. ``min_gap_bp`` is the closest chosen pair on any
    # isoform both probes bind (None when the gene got fewer than two probes);
    # ``gap_used_bp`` is the floor the winning pass ran at. A gene whose
    # ``reasons`` contain ``spacing_relaxed`` only has its probes because the
    # selector dropped from MIN_PROBE_GAP_BP to RELAXED_PROBE_GAP_BP — the
    # probes are real, but they are packed tighter than the design targets and
    # that has to travel with them rather than being inferable from the count.
    min_gap_bp: Optional[int] = None
    gap_used_bp: int = MIN_PROBE_GAP_BP
    # Spacing tier the gene's probes actually earned — read off ``min_gap_bp``,
    # not off the rung that produced it, because a gene can settle on the orange
    # rung and still come out with 400 bp between its probes.
    tier: int = GAP_TIER_GREEN

    @property
    def tier_name(self) -> str:
        return GAP_TIER_NAMES[self.tier]


@dataclass
class CoveragePick:
    """Result of :func:`coverage_pick`.

    ``selections`` is keyed by *source transcript* (the probe's results.csv
    ``gene_id``) so it drops straight into the review dialog / ``export_selected
    _probes`` plumbing. ``flags`` lists genes that fell short of full coverage or
    N spaced probes. ``summary`` carries run-level tallies.
    """
    selections: Dict[str, List[str]]
    flags: List[GeneFlag]
    summary: Dict[str, int]


# --------------------------------------------------------------------------- #
# Matching / spacing
# --------------------------------------------------------------------------- #
def mismatches_for_identity(seq_len: int, min_identity: float = 0.90) -> int:
    """Mismatch budget over ``seq_len`` bp equivalent to a ``min_identity`` rule.

    The FASTA-membership path (``coverage_pick(max_mismatches=...)``) counts
    isoform coverage by exact/near-exact substring match, while the engine's
    BLAST isoform_coverage uses >= 90% identity. They diverge unless the caller
    aligns them: pass ``max_mismatches=mismatches_for_identity(tile_len)`` (5 for
    a 50-mer at 90%) so the FASTA path matches the BLAST rule."""
    return max(0, int(round(seq_len * (1.0 - min_identity))))


def _find_mm(hay: str, needle: str, max_mm: int) -> int:
    """First index where ``needle`` matches ``hay`` with <= max_mm mismatches
    (no indels), else -1. Only used when max_mm > 0."""
    n, m = len(hay), len(needle)
    for i in range(0, n - m + 1):
        mm = 0
        ok = True
        win = hay[i:i + m]
        for a, b in zip(needle, win):
            if a != b:
                mm += 1
                if mm > max_mm:
                    ok = False
                    break
        if ok:
            return i
    return -1


def _membership(cands: List[_Cand], siblings: Sequence[str],
                seqs: Mapping[str, str], max_mm: int) -> None:
    sib_seqs = [(t, seqs[t]) for t in siblings if t in seqs]
    for c in cands:
        for tid, seq in sib_seqs:
            if tid == c.transcript:
                # Use the probe's true tiled position on its own transcript
                # (results.csv ``start``) — ``find`` would return the first
                # occurrence, which is wrong for repetitive targets and would
                # corrupt the spacing check between probes on this transcript.
                c.hit_set.add(tid)
                c.footprint[tid] = c.start
                continue
            pos = seq.find(c.target)
            if pos < 0 and max_mm > 0:
                pos = _find_mm(seq, c.target, max_mm)
            if pos >= 0:
                c.hit_set.add(tid)
                c.footprint[tid] = pos
        # A probe always binds the transcript it was tiled from, even if that
        # record is missing from the FASTA.
        if c.transcript not in c.hit_set:
            c.hit_set.add(c.transcript)
            c.footprint.setdefault(c.transcript, c.start)


def _conflicts(a: _Cand, b: _Cand, length: int, min_gap: int) -> bool:
    # One site cannot crowd itself. Measuring a family's siblings against each
    # other is a category error: they share a footprint by construction, so the
    # answer is always -1 bp and a correctly spaced panel reports itself as
    # entirely crowded. Two DIFFERENT families sitting close still conflict, and
    # the suite carries that as a control.
    if a.site_key and a.site_key == b.site_key:
        return False
    for tid in a.footprint.keys() & b.footprint.keys():
        s1, s2 = a.footprint[tid], b.footprint[tid]
        e1, e2 = s1 + length, s2 + length
        if e1 <= s2:
            gap = s2 - e1
        elif e2 <= s1:
            gap = s1 - e2
        else:
            gap = -1
        if gap < min_gap:
            return True
    return False


def _greedy(cands: List[_Cand], n: int, length: int, min_gap: int,
            preferred_gap: int) -> List[_Cand]:
    """Greedy set-cover with a soft spacing preference.

    ``min_gap`` is a hard constraint (candidates closer than it to a picked
    probe are rejected). ``preferred_gap`` (>= min_gap) is a soft target: at each
    step we take the candidate that adds the most new isoform coverage, breaking
    ties first toward probes that also clear the *preferred* spacing from every
    pick, then toward higher composite score. So selection reaches for the
    roomier spacing but won't sacrifice coverage (or a needed probe) to get it.
    Once coverage saturates, remaining slots fill with the best spaced probes."""
    selected: List[_Cand] = []
    selected_ids: Set[str] = set()   # membership by probe_id, not field-value ==
    selected_sites: Set[str] = set()  # the quota is spent per SITE, not per id
    covered: Set[str] = set()
    pool = sorted(cands, key=lambda c: -c.score)
    while len(selected_sites) < n:
        best: Optional[_Cand] = None
        best_key: Tuple[int, int, float] = (-1, -1, float("-inf"))
        for c in pool:
            if c.probe_id in selected_ids:
                continue
            if any(_conflicts(c, s, length, min_gap) for s in selected):
                continue
            meets_pref = 0 if any(
                _conflicts(c, s, length, preferred_gap) for s in selected) else 1
            key = (len(c.hit_set - covered), meets_pref, c.score)
            if key > best_key:
                best_key, best = key, c
        if best is None:
            break
        selected.append(best)
        selected_ids.add(best.probe_id)
        selected_sites.add(best.site_key or best.probe_id)
        covered |= best.hit_set
        # A degenerate family is ONE site, and only a child is orderable, so
        # choosing one child takes its siblings with it and the family spends one
        # quota slot. ``auto_select._walk`` already does exactly this, under the
        # comment "One family = one pick". Without it the two selectors disagree
        # about what a family is, and coverage maximisation silently drops the
        # other resolutions of a degenerate window.
        if best.site_key:
            for sib in pool:
                if sib.probe_id in selected_ids:
                    continue
                if sib.site_key == best.site_key:
                    selected.append(sib)
                    selected_ids.add(sib.probe_id)
                    covered |= sib.hit_set
    return selected


def _closest_gap(chosen: List[_Cand], length: int) -> Optional[int]:
    """Tightest spacing among the chosen probes, over isoforms they share.

    None when no two chosen probes land on a common isoform — there is then no
    pair whose spacing is defined, which is not the same as a comfortable gap.
    """
    best: Optional[int] = None
    for i, a in enumerate(chosen):
        for b in chosen[i + 1:]:
            if a.site_key and a.site_key == b.site_key:
                continue          # one site, not a pair. See _conflicts.
            for tid in a.footprint.keys() & b.footprint.keys():
                s1, s2 = a.footprint[tid], b.footprint[tid]
                e1, e2 = s1 + length, s2 + length
                gap = s2 - e1 if e1 <= s2 else (s1 - e2 if e2 <= s1 else -1)
                best = gap if best is None else min(best, gap)
    return best


def _greedy_ladder(cands: List[_Cand], n: int, length: int,
                   preferred_gap: int, gaps: Sequence[int],
                   relax_below: int) -> Tuple[List[_Cand], int]:
    """:func:`_greedy` down a ladder of gap floors, stopping when the quota is met.

    Returns ``(chosen, gap_used)``. A looser rung runs only when the one above
    could not reach ``min(n, relax_below)`` probes, and its result is kept only
    when it genuinely yields MORE probes — so a gene that was already fine is
    never re-packed, and a gene never reports a relaxation that bought it
    nothing.

    Unlike ``auto_select``'s ladder this one carries no terminal-band rung: see
    :func:`coverage_pick` for why position is not used to prune candidates here.
    """
    quota = min(n, relax_below)
    best: List[_Cand] = []
    gap_used = gaps[0] if gaps else 0
    for gap in gaps:
        chosen = _greedy(cands, n, length, gap, max(preferred_gap, gap))
        if len(chosen) > len(best):
            best, gap_used = chosen, gap
        if len(best) >= quota:
            break
    return best, gap_used


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #
def coverage_pick(
    rows: Iterable[Mapping[str, str]],
    index: GeneIndex,
    *,
    n_per_gene: int = 3,
    min_gap: int = MIN_PROBE_GAP_BP,
    preferred_gap: int = PREFERRED_PROBE_GAP_BP,
    relaxed_gap: int = RELAXED_PROBE_GAP_BP,
    floor_gap: int = MIN_PROBE_GAP_FLOOR_BP,
    relax_below: int = PROBE_GAP_RELAX_MIN_PROBES,
    max_mismatches: int = 0,
    require_pass: bool = True,
    keep_reduced: bool = True,
    top_k: int = 600,
) -> CoveragePick:
    """Pick up to ``n_per_gene`` spaced probes per gene maximising isoform
    coverage. See module docstring for the model.

    Membership source (per probe), in order of preference:
      1. the persisted ``isoform_hit_ids`` column (no FASTA needed); else
      2. exact/fuzzy sequence match against ``index.seqs`` (the FASTA path).
    A run produced by a current pipeline carries the column, so ``index`` can be
    a sequence-free :func:`gene_index_from_manifest` and no FASTA is required.

    ``min_gap`` is the spacing target; ``preferred_gap`` (>= min_gap) is the
    soft target selection reaches for (see :func:`_greedy`). A gene that cannot
    yield ``min(n_per_gene, relax_below)`` probes at ``min_gap`` drops to
    ``relaxed_gap`` and then to ``floor_gap``, stopping at the first rung that
    fills the quota (see :func:`_greedy_ladder`), and is flagged
    ``spacing_relaxed``. Every gene's achieved spacing and tier is reported in
    ``GeneFlag.min_gap_bp`` / ``GeneFlag.tier`` whether it relaxed or not, and
    the run-level probe counts per tier are in ``summary``.

    Pass ``relaxed_gap``/``floor_gap`` >= the rung above to drop that rung — the
    same mechanism the GUI spacing filter uses to make orange or red probes
    un-selectable rather than selected-then-warned-about.

    Deliberately UNLIKE :func:`flexprobe.auto_select.auto_pick_with_report`,
    this selector runs no adaptive spread search and no interior-first pass.
    Its objective is isoform coverage, and both of those work by shrinking the
    candidate pool — which here would trade away the very thing it exists to
    maximise. End avoidance still reaches it, through the placement term inside
    ``composite_score`` (see ``ranking._central_placement_score``), which is the
    tie-break once coverage saturates."""
    if n_per_gene <= 0:
        return CoveragePick({}, [], {"genes": 0, "probes": 0, "flagged": 0})

    # Group candidate rows by gene (ENSG), keeping top-K by score per gene.
    heaps: Dict[str, List[Tuple[float, int, _Cand]]] = {}
    counter = 0
    for row in rows:
        if require_pass and not _truthy(row.get("specificity_pass", "")):
            continue
        if not keep_reduced and _truthy(row.get("reduced_confidence", "")):
            continue
        transcript = row.get("gene_id", "")
        hit_ids = _parse_hit_ids(row.get("isoform_hit_ids"))
        target = (row.get("mrna_target_50bp") or "").strip().upper()
        # A degenerate consensus parent is not something isoform coverage can be
        # picked from: it is never orderable, and its concrete variants are in
        # these same rows. Excluded on ROLE, on BOTH membership paths.
        #
        # This replaces an `"N" in target` test that ran only on the FASTA-
        # membership branch — a parent carrying a persisted isoform_hit_ids
        # column took the other branch and was silently included, so whether a
        # parent reached the selection depended on which membership source the
        # run happened to have. Role is the property actually being tested, and
        # it holds identically on both paths.
        if row.get("probe_role") == PROBE_ROLE_CONSENSUS:
            continue
        if not hit_ids and not target:
            continue
        try:
            score = float(row.get("composite_score") or 0.0)
            start = int(row["start"]); end = int(row["end"])
        except (TypeError, ValueError, KeyError):
            continue
        _pid = row.get("probe_id", "")
        _role = str(row.get("probe_role") or "").strip()
        _parent = str(row.get("parent_probe_id") or "").strip()
        cand = _Cand(probe_id=_pid, transcript=transcript,
                     start=start, end=end, target=target, score=score,
                     site_key=(_parent if _role == PROBE_ROLE_VARIANT and _parent
                               else _pid))
        if hit_ids:
            cand.from_column = True
            cand.hit_set = set(hit_ids)
            cand.hit_set.add(transcript)       # always binds its own transcript
            # Column path: we know WHICH isoforms a probe binds (hit_ids) but not
            # its coordinate on each, so the footprint is source-transcript-
            # local. Spacing is therefore enforced per source transcript; two
            # probes from different transcripts that overlap on a shared isoform
            # are not de-conflicted here. The FASTA-membership path reconstructs
            # true cross-isoform footprints; pass a FASTA if you need that.
            cand.footprint = {transcript: start}
        gene = index.gene_key(transcript)
        counter += 1
        h = heaps.setdefault(gene, [])
        if len(h) < top_k:
            heapq.heappush(h, (score, counter, cand))
        elif score > h[0][0]:
            heapq.heapreplace(h, (score, counter, cand))

    # Rungs, tightest-first, with any that cannot loosen further removed: a rung
    # equal to the one above returns the identical pick and would then report a
    # relaxation that never happened.
    gap_rungs: List[int] = [min_gap]
    for g in (relaxed_gap, floor_gap):
        if g < gap_rungs[-1]:
            gap_rungs.append(g)

    selections: Dict[str, List[str]] = {}
    flags: List[GeneFlag] = []
    n_probes = 0
    n_relaxed = 0
    tier_probes: Dict[int, int] = {t: 0 for t in GAP_TIER_NAMES}
    tier_genes: Dict[int, int] = {t: 0 for t in GAP_TIER_NAMES}

    for gene, heap in heaps.items():
        cands = [c for _, _, c in heap]
        length = (cands[0].end - cands[0].start) if cands else 50

        # Recompute membership from sequence only for candidates that lack the
        # persisted column (older runs). With no FASTA either, fall back to the
        # probe binding just its own source transcript.
        need = [c for c in cands if not c.from_column]
        if need:
            if index.seqs:
                sibs = sorted(index.group_members.get(gene)
                              or {c.transcript for c in need})
                _membership(need, sibs, index.seqs, max_mismatches)
            else:
                for c in need:
                    c.hit_set = {c.transcript}
                    c.footprint = {c.transcript: c.start}

        # The gene's full isoform universe: the manifest/FASTA map if known, else
        # everything the candidates actually touch (a floor, not the true total).
        universe = set(index.group_members.get(gene, set()))
        coverage_basis = "authoritative"
        if not universe:
            coverage_basis = "floor"   # inferred from candidates → optimistic
            universe = {c.transcript for c in cands}
            for c in cands:
                universe |= c.hit_set

        chosen, gap_used = _greedy_ladder(
            cands, n_per_gene, length, preferred_gap, gap_rungs, relax_below)
        closest = _closest_gap(chosen, length)

        covered: Set[str] = set()
        for c in chosen:
            covered |= c.hit_set
            # Key by the probe's source transcript so the existing per-gene
            # (per-transcript) review UI and export keep working unchanged.
            selections.setdefault(c.transcript, []).append(c.probe_id)
            n_probes += 1
        covered &= universe

        tier = gap_tier(closest)
        tier_probes[tier] += len(chosen)
        if chosen:
            tier_genes[tier] += 1

        reasons: List[str] = []
        if not chosen:
            reasons.append("no_passing_probes")
        if len(chosen) < n_per_gene:
            reasons.append(f"only_{len(chosen)}_spaced_probes")
        if covered != universe:
            reasons.append("incomplete_isoform_coverage")
        if closest is not None and closest < min_gap:
            reasons.append("spacing_relaxed")
            n_relaxed += 1
        if reasons:
            flags.append(GeneFlag(
                gene=gene, isoforms_total=len(universe),
                probes_selected=len(chosen), isoforms_covered=len(covered),
                missing_isoforms=sorted(universe - covered), reasons=reasons,
                coverage_basis=coverage_basis,
                min_gap_bp=closest, gap_used_bp=gap_used, tier=tier,
            ))

    summary = {"genes": len(heaps), "probes": n_probes, "flagged": len(flags),
               "genes_below_min_gap": n_relaxed, "min_gap_bp": min_gap,
               "relaxed_gap_bp": relaxed_gap, "floor_gap_bp": floor_gap}
    for t, name in GAP_TIER_NAMES.items():
        summary[f"probes_{name}"] = tier_probes[t]
        summary[f"genes_{name}"] = tier_genes[t]
    return CoveragePick(selections, flags, summary)


def _truthy(v) -> bool:
    return str(v).strip().lower() in {"true", "1", "yes"}


def _parse_hit_ids(raw) -> List[str]:
    """Parse the ';'-joined ``isoform_hit_ids`` cell into a list of ids."""
    if not raw:
        return []
    return [t for t in (s.strip() for s in str(raw).split(";")) if t]


def results_transcripts(rows: Iterable[Mapping[str, str]]) -> Set[str]:
    """The set of source transcripts (``gene_id`` column) present in rows —
    used to scope :func:`build_gene_index` sequence loading."""
    return {r.get("gene_id", "") for r in rows if r.get("gene_id")}


# --------------------------------------------------------------------------- #
# CSV output (shared by the CLI, the HPC job, and the standalone script)
# --------------------------------------------------------------------------- #
SELECTION_FIELDS = [
    # spacing_tier is green/orange/red for the GENE this probe belongs to —
    # spacing is a property of a pair, so it cannot be attributed to one probe
    # alone, and the gene is the smallest unit it is true of. Carried here so a
    # tier count can be taken straight off the order list without joining
    # against flagged_genes.csv.
    "gene", "gene_symbol", "spacing_tier", "probe_rank_in_gene", "probe_id",
    "source_transcript", "start", "end", "isoforms_total",
    "engine_isoform_hits", "isoform_hit_ids", "composite_score",
    "specificity_status", "reduced_confidence", "mrna_target_50bp",
    "lhs_oligo", "rhs_oligo",
]
FLAGGED_FIELDS = [
    "gene", "gene_symbol", "isoforms_total", "probes_selected",
    "isoforms_covered", "coverage_fraction", "coverage_basis",
    "flags", "missing_isoforms",
    # APPENDED, never inserted, so a script reading an older flagged_genes.csv
    # keeps working. min_gap_bp is the tightest chosen pair (blank when the gene
    # got fewer than two probes); gap_used_bp is the floor the winning pass ran
    # at, so gap_used_bp < MIN_PROBE_GAP_BP is exactly the relaxation pass.
    # spacing_tier is the green/orange/red the user sees everywhere else.
    "min_gap_bp", "gap_used_bp", "spacing_tier",
]


def write_selection_csvs(pick: CoveragePick, index: GeneIndex,
                         rows: Iterable[Mapping[str, str]],
                         out_path: str, flagged_path: str) -> None:
    """Write the per-probe selection CSV and the flagged-genes CSV.

    ``out_path`` lists the chosen probes (gene-ordered, ranked by score within
    each gene). ``flagged_path`` lists genes that couldn't be fully covered or
    couldn't fit the requested number of spaced probes.
    """
    import csv as _csv

    rows_by_id = {r.get("probe_id", ""): r for r in rows}
    by_gene: Dict[str, List[str]] = {}
    for transcript, ids in pick.selections.items():
        by_gene.setdefault(index.gene_key(transcript), []).extend(ids)
    # Every gene's tier, not just the flagged ones: a green gene raises no flag,
    # so without this the column would be blank exactly where it says "fine".
    tier_of = {f.gene: f.tier_name for f in pick.flags}

    def _score(pid: str) -> float:
        try:
            return float(rows_by_id[pid].get("composite_score") or 0.0)
        except (KeyError, TypeError, ValueError):
            return 0.0

    with open(out_path, "w", newline="") as fo:
        w = _csv.DictWriter(fo, fieldnames=SELECTION_FIELDS)
        w.writeheader()
        for gene in sorted(by_gene):
            ids = sorted(by_gene[gene], key=lambda pid: -_score(pid))
            total = len(index.group_members.get(gene, ())) or 1
            for rank, pid in enumerate(ids, 1):
                r = rows_by_id.get(pid, {})
                w.writerow({
                    "gene": gene,
                    "gene_symbol": index.symbol_of.get(gene, ""),
                    "spacing_tier": tier_of.get(gene, "green"),
                    "probe_rank_in_gene": rank,
                    "probe_id": pid,
                    "source_transcript": r.get("gene_id", ""),
                    "start": r.get("start", ""), "end": r.get("end", ""),
                    "isoforms_total": total,
                    "engine_isoform_hits": r.get("isoform_hits", ""),
                    "isoform_hit_ids": r.get("isoform_hit_ids", ""),
                    "composite_score": r.get("composite_score", ""),
                    "specificity_status": r.get("specificity_status", ""),
                    "reduced_confidence": r.get("reduced_confidence", ""),
                    "mrna_target_50bp": r.get("mrna_target_50bp", ""),
                    "lhs_oligo": r.get("lhs_oligo", ""),
                    "rhs_oligo": r.get("rhs_oligo", ""),
                })

    with open(flagged_path, "w", newline="") as ff:
        w = _csv.DictWriter(ff, fieldnames=FLAGGED_FIELDS)
        w.writeheader()
        for f in sorted(pick.flags, key=lambda x: x.gene):
            frac = f.isoforms_covered / max(f.isoforms_total, 1)
            w.writerow({
                "gene": f.gene,
                "gene_symbol": index.symbol_of.get(f.gene, ""),
                "isoforms_total": f.isoforms_total,
                "probes_selected": f.probes_selected,
                "isoforms_covered": f.isoforms_covered,
                "coverage_fraction": f"{frac:.4f}",
                "coverage_basis": f.coverage_basis,
                "flags": ";".join(f.reasons),
                "missing_isoforms": ";".join(f.missing_isoforms),
                "min_gap_bp": "" if f.min_gap_bp is None else f.min_gap_bp,
                "gap_used_bp": f.gap_used_bp,
                "spacing_tier": f.tier_name,
            })
