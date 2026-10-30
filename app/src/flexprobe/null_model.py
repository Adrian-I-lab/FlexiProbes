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
# flexiprobes. If not, see <https://www.gnu.org/licenses/>.
"""Null-model calibration of the off-target specificity gate.

A specificity screen only means something if it rejects *specific* sequence more
often than it rejects noise. This module measures the second number, so a run
can state its own gate's false-positive rate instead of leaving it unknown.

The method:

1. Take a bounded random sample of the windows that cleared the filter funnel.
2. Generate ``shuffles`` dinucleotide-preserving shuffles of each
   (Altschul-Erikson Eulerian-path shuffle). Dinucleotide shuffling is the right
   null: BLAST word seeding and duplex stability both track dinucleotide
   composition, so a mononucleotide shuffle understates how often noise finds a
   match and a uniform-random control understates it badly.
3. Screen the shuffled windows through *exactly* the same path as the real ones
   — same 25 nt half split, same BLAST parameters, same databases, same
   ``flexprobe.specificity`` classifier.
4. Report the fraction of pure-noise windows the gate rejects: the gate's
   empirical false-positive rate on composition-matched sequence, with a
   conservative 95% Wilson interval alongside it (see ``summarise``) --
   ``shuffles`` copies of one window are correlated, not independent draws,
   so the interval is computed on the sampled-WINDOW count, not the larger
   shuffle-inflated total.

A gate that rejects a third of composition-matched random sequence is not
measuring specificity; it is measuring how big the database is. The rate can
vary widely between runs that differ only in target and database, which is why
it is reported per run rather than assumed.

**Cost.** Bounded by construction: ``max_windows`` caps the sample regardless of
run size, so the extra work is at most ``max_windows * shuffles`` BLAST queries
(default 200 * 20 = 4,000) no matter whether the run screened 40 windows or 35
million. That is negligible beside a whole-transcriptome screen and small even
beside a viral segment run.
"""
from __future__ import annotations

import logging
import math
import os
import random
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)

_NUCS = "ACGT"

# Defaults for a run that asks for calibration without tuning it. 20 shuffles
# per window is the floor at which a per-window rate is worth quoting; 200
# windows keeps the whole calibration under ~4,000 extra queries.
DEFAULT_SHUFFLES = 20
DEFAULT_MAX_WINDOWS = 200
REPORT_FILE = "null_model.json"


# ---------------------------------------------------------------------------
# Dinucleotide-preserving shuffle (Altschul & Erikson 1985; Eulerian-path form)
# ---------------------------------------------------------------------------
def _count_and_lists(s: str) -> Tuple[Dict[str, Dict[str, int]],
                                      Dict[str, List[str]]]:
    lists: Dict[str, List[str]] = {x: [] for x in _NUCS}
    dinucl = {x: {y: 0 for y in _NUCS} for x in _NUCS}
    for i in range(len(s) - 1):
        x, y = s[i], s[i + 1]
        lists[x].append(y)
        dinucl[x][y] += 1
    return dinucl, lists


def _choose_edge(x: str, dinucl: Dict[str, Dict[str, int]],
                 rng: random.Random) -> Optional[str]:
    denom = sum(dinucl[x][y] for y in _NUCS)
    if denom == 0:
        return None
    z = rng.random()
    acc = 0
    for y in _NUCS:
        acc += dinucl[x][y]
        if z < acc / denom:
            dinucl[x][y] -= 1
            return y
    for y in _NUCS:                      # float-rounding fallback
        if dinucl[x][y] > 0:
            dinucl[x][y] -= 1
            return y
    return None


def _connected_to_last(edge: Dict[str, str], last: str) -> bool:
    """Every vertex with a chosen out-edge must reach ``last`` by following them."""
    for x in edge:
        seen = set()
        v: Optional[str] = x
        while v is not None and v != last and v not in seen:
            seen.add(v)
            v = edge.get(v)
        if v != last:
            return False
    return True


def _eulerian(s: str, rng: random.Random, tries: int = 200):
    last = s[-1]
    for _ in range(tries):
        dinucl, lists = _count_and_lists(s)
        edge: Dict[str, str] = {}
        ok = True
        for x in _NUCS:
            if x == last or not lists[x]:
                continue
            y = _choose_edge(x, dinucl, rng)
            if y is None:
                ok = False
                break
            edge[x] = y
        if ok and _connected_to_last(edge, last):
            return lists, edge
    return None


def dinucleotide_shuffle(s: str, rng: random.Random) -> str:
    """Shuffle ``s`` preserving its exact dinucleotide counts (and both ends).

    Falls back to a mononucleotide shuffle of the interior for sequences that
    are too short, contain non-ACGT characters, or hit the pathological case
    where no Eulerian path is found within ``tries``. The fallback preserves
    length and base composition but not dinucleotide counts, so it is a weaker
    null — it is counted separately in the report rather than hidden.
    """
    s = s.upper()
    if len(s) < 3 or any(c not in _NUCS for c in s):
        body = list(s[1:-1])
        rng.shuffle(body)
        return s[0] + "".join(body) + s[-1]
    res = _eulerian(s, rng)
    if res is None:                       # pathological; degrade gracefully
        body = list(s[1:-1])
        rng.shuffle(body)
        return s[0] + "".join(body) + s[-1]
    lists, edge = res
    for x, y in edge.items():
        lists[x].remove(y)
    for x in _NUCS:
        rng.shuffle(lists[x])
    for x, y in edge.items():
        lists[x].append(y)               # last edge out of each vertex
    out = [s[0]]
    ptr = {x: 0 for x in _NUCS}
    cur = s[0]
    for _ in range(len(s) - 1):
        nxt = lists[cur][ptr[cur]]
        ptr[cur] += 1
        out.append(nxt)
        cur = nxt
    return "".join(out)


def preserves_dinucleotides(original: str, shuffled: str) -> bool:
    """True when ``shuffled`` has the exact dinucleotide counts of ``original``."""
    if len(original) != len(shuffled):
        return False
    return _count_and_lists(original.upper())[0] == \
        _count_and_lists(shuffled.upper())[0]


# ---------------------------------------------------------------------------
# Building the null population
# ---------------------------------------------------------------------------
def sample_windows(candidates: List[Any], *, max_windows: int,
                   rng: random.Random) -> List[Any]:
    """Bounded, reproducible sample of the windows to calibrate against.

    Sampled rather than truncated: taking the first N would take them all from
    the first gene / the 5' end of the first segment, whose composition is not
    the run's composition.
    """
    if max_windows <= 0 or len(candidates) <= max_windows:
        return list(candidates)
    return rng.sample(list(candidates), max_windows)


def build_null_candidates(candidates: List[Any], *, shuffles: int,
                          rng: random.Random) -> Tuple[List[Any], int]:
    """Shuffled stand-ins for ``candidates``, ready for the real screen.

    Returns ``(null_candidates, n_degraded)`` where ``n_degraded`` counts the
    shuffles that fell back to mononucleotide shuffling (see
    ``dinucleotide_shuffle``).
    """
    from .models import ProbeCandidate

    out: List[Any] = []
    degraded = 0
    for c in candidates:
        for k in range(shuffles):
            seq = dinucleotide_shuffle(c.sequence, rng)
            if not preserves_dinucleotides(c.sequence, seq):
                degraded += 1
            out.append(ProbeCandidate(
                # Same gene_id, so on-target classification treats the shuffled
                # window exactly as it treated the real one. The probe_id is
                # namespaced so a null row can never be mistaken for a design.
                gene_id=c.gene_id,
                probe_id=f"{c.probe_id}__null{k:03d}",
                sequence=seq,
                start=c.start,
                end=c.end,
            ))
    return out, degraded


# z for a 95% two-sided normal interval (scipy.stats.norm.ppf(0.975)).
# Hardcoded rather than pulling in scipy for one constant -- this module has
# no numeric dependency today and a 95% CI is what a Methods section expects
# by default.
_Z_95 = 1.959963984540054


def _wilson_interval(phat: float, n: float, z: float = _Z_95) -> Tuple[float, float]:
    """Wilson score interval for a binomial proportion ``phat`` over ``n`` trials.

    Preferred over the naive normal (Wald) interval ``phat +/- z*sqrt(phat(1-
    phat)/n)``: Wald can extend below 0 or above 1 and collapses to a
    zero-width interval at phat=0, which is routine here: a run can
    reject 0/N of the null population. Wilson stays inside [0, 1] and degrades sensibly at
    small n or extreme phat instead.
    """
    if n <= 0:
        return (0.0, 1.0)          # no evidence at all -- the widest honest answer
    denom = 1.0 + z * z / n
    center = (phat + z * z / (2.0 * n)) / denom
    half = (z / denom) * math.sqrt(phat * (1.0 - phat) / n + z * z / (4.0 * n * n))
    return (max(0.0, center - half), min(1.0, center + half))


def summarise(null_candidates: List[Any], *, n_real_windows: int,
              n_sampled: int, shuffles: int, degraded: int,
              seed: int) -> Dict[str, Any]:
    """Turn a screened null population into the reported statistics.

    Reports a rate only when every null window was recounted EXACTLY. A
    conservative or legacy recount means positions were reconstructed rather
    than compared, and both over-reject rather than over-pass -- so they bias
    the false-positive rate upward, which is precisely the direction that
    corrupts the one conclusion this number is drawn for ("a high value means
    the gate is rejecting on database size"). A number that can only be wrong
    one way is not a calibration, so an unmeasurable run says so and omits the
    rate rather than shipping an inflated one.
    """
    from .specificity import RECOUNT_EXACT

    n = len(null_candidates)
    rejected = sum(1 for c in null_candidates if not c.specificity_pass)
    by_status: Dict[str, int] = {}
    for c in null_candidates:
        by_status[c.specificity_status] = by_status.get(c.specificity_status, 0) + 1
    by_basis: Dict[str, int] = {}
    for c in null_candidates:
        b = c.offtarget_recount_basis
        by_basis[b] = by_basis.get(b, 0) + 1
    inexact = n - by_basis.get(RECOUNT_EXACT, 0)
    fp_rate = (rejected / n) if n else 0.0
    if inexact:
        return {
            "windows_in_run": n_real_windows,
            "windows_sampled": n_sampled,
            "shuffles_per_window": shuffles,
            "null_windows_screened": n,
            "seed": seed,
            "degraded_shuffles": degraded,
            # No `false_positive_rate` key AT ALL, deliberately. A reader doing
            # payload.get("false_positive_rate") gets None and cannot mistake it
            # for a measurement; a sentinel value like -1 or 0.0 could be plotted.
            "measurable": False,
            "null_rejected": rejected,
            "null_status_counts": by_status,
            "null_recount_basis_counts": by_basis,
            "unmeasurable_reason": (
                f"{inexact} of {n} null windows could not be recounted exactly. "
                "Their off-target mismatch counts were reconstructed rather than "
                "compared, which over-rejects and never over-passes, so any rate "
                "computed from them would be an unknown amount too HIGH. Usual "
                "cause: the off-target reference cannot serve sequence retrieval "
                "(a BLAST database built without -parse_seqids). No rate is "
                "reported."
            ),
            "interpretation": (
                "NOT MEASURED. See unmeasurable_reason. Do not substitute zero "
                "or treat the absence of a rate as a low rate."
            ),
        }
    # CONFIDENCE INTERVAL. The naive binomial interval on `n` (=
    # windows_sampled * shuffles) treats all `n` outcomes as independent, but
    # they are not: `shuffles` dinucleotide-shuffled copies of the SAME parent
    # window share that window's composition and its position relative to
    # whatever off-target sequence it happens to be near, so a window that is
    # "easy" (or "hard") for the gate to reject tends to stay easy (or hard)
    # across its own shuffles. That is pseudo-replication -- the effective
    # number of INDEPENDENT trials is closer to windows_sampled than to
    # windows_sampled * shuffles, and a flat-n interval computed on the larger
    # number is too narrow, understating the true uncertainty in the exact
    # direction that would make a borderline rate look more precise than it is.
    #
    # This uses windows_sampled itself as the effective sample size, which is
    # the WORST-CASE (most conservative) bound: it is what the effective n
    # would be if every shuffle of a window were perfectly correlated with its
    # siblings (intraclass correlation = 1), i.e. one window contributes
    # exactly as much independent evidence as one of its shuffles alone would.
    # A full cluster-robust estimate (e.g. fitting the actual intra-cluster
    # correlation via a one-way ANOVA over per-window rejection rates) would
    # recover a tighter interval when the true correlation is less than 1, but
    # most runs here sample only tens of windows -- too few to estimate an ICC
    # reliably -- and this number exists to warn against OVER-trusting the
    # point estimate, not to produce the narrowest defensible interval. Taking
    # the conservative bound errs in the direction that keeps the number
    # honest: it can only ever be too wide, never too narrow.
    eff_n = n_sampled if n_sampled > 0 else n
    ci_lo, ci_hi = _wilson_interval(fp_rate, eff_n)
    return {
        "measurable": True,
        "null_recount_basis_counts": by_basis,
        # What was calibrated.
        "windows_in_run": n_real_windows,
        "windows_sampled": n_sampled,
        "shuffles_per_window": shuffles,
        "null_windows_screened": n,
        "seed": seed,
        "degraded_shuffles": degraded,
        # The number the whole module exists to produce.
        "null_rejected": rejected,
        "false_positive_rate": round(fp_rate, 6),
        # 95% Wilson score interval around false_positive_rate, using
        # windows_sampled (not null_windows_screened) as the effective n to
        # account for shuffles-of-the-same-window pseudo-replication -- see
        # the comment above. Present only alongside a reported rate: there is
        # nothing to bound when the point estimate itself is suppressed as
        # unmeasurable (the branch above never reaches this far).
        "false_positive_rate_ci95": {
            "lower": round(ci_lo, 6),
            "upper": round(ci_hi, 6),
            "effective_n": eff_n,
            "method": "wilson_score_conservative_cluster_n",
        },
        "null_status_counts": by_status,
        "interpretation": (
            "Fraction of composition-matched random 50-mers the off-target "
            "gate rejects. This is the gate's false-positive rate, not an "
            "error rate of the run: a high value means the gate is rejecting "
            "on database size rather than on specificity, and the run's "
            "rejections should be read with that in mind. "
            "false_positive_rate_ci95 is a conservative 95% interval that "
            "treats each sampled window (not each shuffle) as one independent "
            "observation, since shuffles of the same window are correlated."
        ),
    }


def write_report(output_dir: str, payload: Dict[str, Any]) -> str:
    """Write ``null_model.json`` next to the run's other artefacts."""
    import json

    path = os.path.join(output_dir, REPORT_FILE)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=2)
    os.replace(tmp, path)
    return path
