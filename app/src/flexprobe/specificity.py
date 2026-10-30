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

"""Specificity engine: strain-rule evaluation and off-target rejection."""
from __future__ import annotations
import logging
from enum import Enum
from typing import List, Dict, Any, Sequence, Set, Optional, Tuple
from .blast import reverse_complement as _reverse_complement
from .models import ProbeCandidate
from .constants import (
    OFF_TARGET_PIDENT_MIN, OFF_TARGET_ALN_MIN,
    ON_TARGET_PIDENT_MIN, ON_TARGET_ALN_MIN,
    OFF_TARGET_MIN_MISMATCH, OFF_TARGET_FULL_HALF_DEFAULT,
    OFF_TARGET_REQUIRE_BOTH_HALVES_DEFAULT, OFF_TARGET_SAME_SUBJECT_DEFAULT,
    OFF_TARGET_JUNCTION_DUPLEX_BP,
    HALF_LENGTH_DEFAULT, JUNCTION_WINDOW,
)

logger = logging.getLogger("flexprobe.specificity")


class StrainRule(Enum):
    """Supported strain-matching modes."""
    EXACT_SET = "exact_set"
    PAN_STRAIN = "pan_strain"
    STRAIN_UNIQUE = "strain_unique"


# Characters that separate the sequence id from trailing fields in a BLAST
# subject id (FASTA headers like ``OROV_BeAn19991|Large|ACC123`` or
# ``KR013237.1``). Deliberately EXCLUDES '_' because strain/gene ids legitimately
# contain underscores — so a target ``OROV`` must NOT match ``OROV_BeAn19991``.
_ID_DELIMS = "|.: /\t"


def id_matches(sseqid: str, tid: str) -> bool:
    """True if BLAST subject ``sseqid`` refers to target ``tid``.

    A match is an exact equality OR ``tid`` followed by a header delimiter
    (so ``OROV_BeAn19991`` matches ``OROV_BeAn19991|Large|ACC`` but a short
    prefix like ``OROV`` does NOT spuriously match ``OROV_BeAn19991``). A bare
    ``startswith`` would mis-classify off-targets whose accession merely begins
    with a target id.
    """
    if not tid:
        return False
    if sseqid == tid:
        return True
    if sseqid.startswith(tid):
        return sseqid[len(tid):len(tid) + 1] in _ID_DELIMS
    return False


def hit_query_coverage(hit: Dict[str, Any], half_len: int) -> Optional[float]:
    """Fraction of the probe half this HSP actually aligned, or None.

    BLAST's own ``qcovhsp`` field is exactly ``(qend - qstart + 1) / qlen``,
    so it is derived here rather than requested as an extra output column. The column
    would add a few bytes to every row of a file that already runs to hundreds
    of megabytes per chunk under the exact recount.

    Returns None when the hit carries no usable query range (a hand-built
    fixture, a truncated row), which callers treat as "unknown" rather than
    folding a fabricated 0 or 1 into a minimum.
    """
    if half_len <= 0:
        return None
    try:
        qstart = int(hit["qstart"])
        qend = int(hit["qend"])
    except (KeyError, TypeError, ValueError):
        return None
    # blastn reports qstart <= qend for the query even on a minus-strand hit
    # (the strand shows in sstart/send), but abs() costs nothing and makes this
    # independent of that guarantee.
    span = abs(qend - qstart) + 1
    return max(0.0, min(1.0, span / half_len))


def hit_half(hit: Dict[str, Any]) -> str:
    """Return which probe half a BLAST hit came from: ``'L'``, ``'R'`` or ``''``.

    The BLAST screen queries each 25 bp half separately under the query ids
    ``<probe_id>__L`` / ``<probe_id>__R`` (see pipeline._run_screen), and
    the raw ``qseqid`` survives into the parsed hit dict. That suffix is the
    only record of *which* half a hit belongs to, and the 10x off-target rule
    is stated per half ("in at least one of the LHS or RHS probes"), so it has
    to be recovered here rather than discarded with the id.

    Returns ``''`` for a hit with no half suffix (e.g. a hand-built fixture or
    a cached hit file from a whole-tile query); callers treat that as
    "half unknown", which is handled conservatively — see
    ``evaluate_specificity``.
    """
    qid = str(hit.get("qseqid", ""))
    if qid.endswith("__L"):
        return "L"
    if qid.endswith("__R"):
        return "R"
    return ""


# How a probe's off-target mismatch counts were arrived at. Recorded per
# candidate (ProbeCandidate.offtarget_recount_basis) and surfaced in results.csv
# + the manifest, so a run can never be read as exactly-recounted when it was
# only bounded.
#   exact         — every off-target hit was recounted against the real subject
#                   bases; the count IS the number of mismatches.
#   conservative  — at least one hit's subject window could not be retrieved, so
#                   its unaligned query positions were counted as MATCHES. The
#                   count is a LOWER BOUND, so the screen may over-reject.
#   legacy_local  — off_target_full_half=False; the local-alignment rule ran and
#                   no full-half count was made at all.
#   not_screened  — no BLAST hits were scored for this probe (BLAST skipped or
#                   failed).
RECOUNT_EXACT = "exact"
RECOUNT_CONSERVATIVE = "conservative"
RECOUNT_LEGACY_LOCAL = "legacy_local"
RECOUNT_NOT_SCREENED = "not_screened"


def hit_junction_mismatches(
    hit: Dict[str, Any],
    half_len: int,
    junction_window: int = JUNCTION_WINDOW,
) -> Optional[int]:
    """Mismatches of this off-target hit that sit NEXT TO the ligation junction.

    THE SIGN IS INVERTED HERE, and that is the whole point of the field.

    For a strain the probe is meant to bind, a junction-proximal mismatch is the
    worst kind: the ligase will not seal across it, so the probe fails. See
    ``consensus._binding_score_for_window``, which weights each mismatch by how
    far it sits from the nick — a linear ramp from 3x at the junction down to 1x
    at the end of the probe, so -6.0 C for a mismatch abutting the nick against
    -2.0 C for one at the outer end. It is a ramp rather than a bonus applied
    inside a fixed window: the destabilisation does not step off a cliff at any
    particular distance.

    For an OFF-target, that same mismatch is PROTECTIVE. A probe whose four
    mismatches against some other gene all sit at the junction still cannot
    ligate there and generates no false signal. A probe whose four mismatches sit
    at the outer ends of its halves is barely destabilised, hybridises, ligates,
    and puts reads on the wrong gene. The flat ">= 5 mismatches in one half"
    count treats those two as identical. They are not remotely equivalent, and
    the second is the one that costs you an experiment.

    The junction-proximal zone of a half is the ``junction_window`` bases
    ADJACENT TO THE NICK — the LAST bases of the LHS half and the FIRST bases of
    the RHS half — because the nick sits between them. It is not the middle of
    the half.

    Returns None when the hit cannot be attributed to a half, or when the
    aligned segments are unavailable, or when the HSP does not reach the
    junction zone at all: in each case the answer is unknown, and an unknown
    must not be reported as zero (zero here means "nothing protective", the
    dangerous reading).
    """
    which = hit_half(hit)
    if which not in ("L", "R") or half_len <= 0 or junction_window <= 0:
        return None
    qseq = str(hit.get("qseq", "") or "").upper()
    sseq = str(hit.get("sseq", "") or "").upper()
    if not qseq or not sseq or len(qseq) != len(sseq):
        return None
    try:
        qstart = int(hit["qstart"])
    except (KeyError, TypeError, ValueError):
        return None
    # The zone, as 0-based query offsets within this half.
    if which == "L":
        lo, hi = half_len - junction_window, half_len - 1
    else:
        lo, hi = 0, junction_window - 1
    qpos = qstart - 1                       # 0-based cursor into the half
    covered = False
    n = 0
    for q, sbase in zip(qseq, sseq):
        if q == "-":
            continue                        # consumes no query position
        if lo <= qpos <= hi:
            covered = True
            if q != sbase:
                n += 1
        qpos += 1
    return n if covered else None


#: Key that marks a junction-duplex record inside a probe's hit list. Such a
#: record is NOT an alignment and must never reach the mismatch-count rule; see
#: :func:`split_junction_records`.
JUNCTION_DUPLEX_KEY = "junction_duplex"

#: ``specificity_status`` for a pair the 10x count rule would have shipped and
#: the junction-duplex gate rejected. Deliberately distinct from
#: ``failed_off_target``: that status means the VENDOR's sentence rejected the
#: probe, this one means ours did, and a reader of results.csv has to be able to
#: separate the two populations — the whole claim of this gate is that it only
#: ever rejects probes the count rule passes.
STATUS_FAILED_JUNCTION_DUPLEX = "failed_junction_duplex"

#: Key that marks an UNSCREENABLE-HALF record inside a probe's hit list. Like a
#: junction record it is not an alignment; unlike one it is not evidence either,
#: it is the screen saying "I could not search this half at all". See
#: :func:`split_screen_records` and ``screen.expand_degenerate_half``.
UNSCREENABLE_HALF_KEY = "unscreenable_half"

#: ``specificity_status`` for a pair one of whose halves the screen refused to
#: search — in practice, a half carrying more degenerate bases than the engine will
#: enumerate concrete resolutions for.
#:
#: It is a REJECTION, not a warning, and that is the whole point of the status.
#: If such a half were silently dropped from the search, the pair would keep
#: whatever verdict its other half earned, and a probe with a perfect
#: off-target under the unsearched half would be certified clean. A screen that
#: cannot answer must say so in the column the pipeline already reads, not
#: return the answer it would have given had the question been asked.
STATUS_UNSCREENABLE_HALF = "unscreenable_half"

# Bases that can actually pair. Anything else — an N, an IUPAC ambiguity code, a
# record-separator sentinel — is not a matched base pair and terminates a clean
# run. See junction_clean_duplex for which way that errs.
_PAIRING_BASES = frozenset("ACGT")


def split_screen_records(
    hits: List[Dict[str, Any]],
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    """``(alignments, junction records, unscreenable-half records)``.

    All three travel in the same per-probe list because that list is the only
    channel between the screen and this module (and the only thing the resume
    cache persists). Only the first kind is an alignment: the other two have no
    ``qseq``/``sseq``, no half attribution and no mismatch count, so letting one
    reach ``classify_hit`` would score it as a 0-mismatch off-target on BOTH
    halves and reject every probe. Every consumer of a hit list in this module
    therefore filters first.
    """
    aln: List[Dict[str, Any]] = []
    jd: List[Dict[str, Any]] = []
    uns: List[Dict[str, Any]] = []
    for h in hits:
        if h.get(JUNCTION_DUPLEX_KEY):
            jd.append(h)
        elif h.get(UNSCREENABLE_HALF_KEY):
            uns.append(h)
        else:
            aln.append(h)
    return aln, jd, uns


def split_junction_records(
    hits: List[Dict[str, Any]],
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """``(alignments, junction records)`` — the two-way form of the split.

    Kept because most callers only ever want the alignments. An
    unscreenable-half record is dropped by BOTH lists here: it is not an
    alignment and it is not junction evidence, and a caller that does not know
    about it must not be handed it as either. ``evaluate_specificity`` — the one
    place the refusal has to change a verdict — calls
    :func:`split_screen_records` instead.
    """
    aln, jd, _uns = split_screen_records(hits)
    return aln, jd


def junction_clean_duplex(
    query: str,
    subject: str,
    nick: Optional[int] = None,
) -> Tuple[int, int]:
    """``(left, right)``: matched bases either side of the ligation nick.

    THE SIGN IS INVERTED HERE TOO, and in the same direction as
    :func:`hit_junction_mismatches` — but this measures the thing that function
    can only infer, so read the two together.

    ``hit_junction_mismatches`` argues that mismatches next to the nick are
    PROTECTIVE against an off-target: the ligase will not seal across a
    junction-proximal mispair, so a probe whose mismatches bunch there makes no
    library molecule and no false signal. This function states the same
    biology from the other end. A long unbroken MATCHED run spanning the nick is
    the exact condition under which the ligase *will* seal, so it is the direct
    measure of danger rather than the absence of a protective feature. A pair
    can satisfy 10x's ">= 5 mismatches in one half" comfortably — with every one
    of those mismatches out at the oligo ends — and still present 20 perfectly
    matched bases across the junction. That pair hybridises, ligates, and puts
    reads on the wrong gene.

    ``query`` is the probe's target window (the mRNA sense sequence it is tiled
    from) and ``subject`` is the CONTIGUOUS reference sequence opposite it, in
    the same frame and of the same length. Contiguity is load-bearing and is
    what makes this stricter than the per-half rule rather than merely
    different: two runs measured in one unbroken stretch of reference are
    necessarily adjacent on ONE molecule in ONE gene, which is the condition
    ``off_target_same_subject`` approximates by grouping ids. Here it is not an
    approximation.

    ``nick`` is the index in ``query`` of the first base 3' of the ligation
    junction — i.e. the start of the RHS half — and defaults to the midpoint.
    ``left`` counts backwards from ``nick - 1`` and ``right`` forwards from
    ``nick``; each stops at the first position that is not an identical A/C/G/T
    pair. They are reported separately, never summed: the ligase needs a stable
    duplex on EACH side, so 15+3 is not 9+9 (see
    :data:`constants.OFF_TARGET_JUNCTION_DUPLEX_BP`).

    UNKNOWNS. Anything that is not an identical ACGT pair terminates the run:
    an N or IUPAC code on either side, a record-separator sentinel, or the
    subject simply running out (a site near the end of a transcript). All of
    those are positions where no base pair forms, so none of them extends a
    duplex. Unlike ``hit_junction_mismatches`` there is no ``None`` case, and
    the reason is that the answer is never partial: this walks contiguous
    sequence outward from the nick, so the first unknown IS the end of the
    duplex and nothing beyond it can matter. What the caller must still
    distinguish is a measurement of 0 from a scan that never ran — see
    ``ProbeCandidate.offtarget_clean_left_bp``, which reserves -1 for the
    second.

    The error direction, stated plainly: an ambiguous reference SHORTENS a run,
    which makes this gate fire LESS often. That is the permissive direction, and
    it is deliberate — a run of Ns is not evidence of a ligatable duplex, and
    10x's count rule still sits underneath this one either way.
    """
    n = len(query)
    if nick is None:
        nick = n // 2
    if n == 0 or not 0 <= nick <= n:
        return (0, 0)
    q = query.upper()
    sub = subject.upper()
    left = 0
    for j in range(nick - 1, -1, -1):
        if j >= len(sub) or q[j] != sub[j] or q[j] not in _PAIRING_BASES:
            break
        left += 1
    right = 0
    for j in range(nick, n):
        if j >= len(sub) or q[j] != sub[j] or q[j] not in _PAIRING_BASES:
            break
        right += 1
    return (left, right)


def _hsp_internal_mismatches(hit: Dict[str, Any]) -> Tuple[int, bool]:
    """Mismatched query positions INSIDE the HSP, and whether that is exact.

    Counted column by column from the aligned segments (``qseq``/``sseq``, added
    to :data:`blast.BLASTN_OUTFMT_COLS` for exactly this) rather than derived
    from ``pident``, which is rounded to three decimals and says nothing about
    where the mismatches are. A column whose subject base is a gap is a query
    position that does not base-pair, so it counts; a column whose *query* base
    is a gap consumes no query position at all, so it does not.

    Two fallbacks, for cached hits that carry no aligned segments and for
    hand-built fixtures: the reported ``mismatch`` column, then ``length * (100 - pident)``.
    Both ignore gap columns and so can only UNDER-count — the safe direction —
    but neither is exact, hence the flag.

    Note what the second fallback is NOT. Using ``pident`` to count mismatches
    WITHIN the HSP is a faithful reconstruction of the column BLAST would have
    reported. What must never happen is different and worse: extrapolating that
    identity across the WHOLE half, so query positions BLAST never looked at are
    scored as mismatches. Nothing here extends past ``length``.
    """
    qseq = str(hit.get("qseq", "") or "").upper()
    sseq = str(hit.get("sseq", "") or "").upper()
    if qseq and sseq and len(qseq) == len(sseq):
        mm = 0
        for q, s in zip(qseq, sseq):
            if q == "-":
                continue          # gap in the query: no query position here
            if q != s:
                mm += 1           # mismatch, or a subject gap (unpaired)
        return mm, True
    try:
        if "mismatch" in hit:
            return max(0, int(hit["mismatch"])), False
        aln_len = int(hit.get("length", 0))
        pident = float(hit.get("pident", 100.0))
        return max(0, int(round(aln_len * (100.0 - pident) / 100.0))), False
    except (TypeError, ValueError):
        return 0, False


def _flank_mismatches(
    hit: Dict[str, Any],
    half_len: int,
    query_seq: str,
    subjects: Any,
) -> Optional[int]:
    """Mismatches in the query positions the HSP never covered, or None.

    ``None`` means "could not be determined" and the caller must treat those
    positions as matches (the conservative fallback). Returns 0 without touching
    the database when the HSP already spans the whole half — the common case, and
    the reason an on-target-heavy run costs no retrievals.

    Coordinates: BLAST always reports the query on the plus strand, so
    ``qstart <= qend``; a minus-strand hit is signalled by ``sstart > send``.
    Walking one base LEFT along the query therefore walks left along the subject
    on a plus-strand hit and RIGHT on a minus-strand one, and the retrieved
    window is reverse-complemented in the latter case so it reads in the probe's
    frame.
    """
    try:
        qstart = int(hit.get("qstart", 0))
        qend = int(hit.get("qend", 0))
        sstart = int(hit.get("sstart", 0))
        send = int(hit.get("send", 0))
    except (TypeError, ValueError):
        return None
    if qstart < 1 or qend < qstart or qend > half_len:
        return None
    left_n = qstart - 1
    right_n = half_len - qend
    if left_n == 0 and right_n == 0:
        return 0                  # HSP spans the half; nothing was skipped
    if subjects is None or len(query_seq) < half_len:
        return None
    sseqid = str(hit.get("sseqid", "") or "")
    if not sseqid or sstart < 1 or send < 1:
        return None
    minus = send < sstart

    total = 0
    for n, q_lo, s_lo, s_hi in (
        # (positions, query slice start (0-based), subject range 1-based)
        (left_n, 0,
         (sstart - left_n) if not minus else (sstart + 1),
         (sstart - 1) if not minus else (sstart + left_n)),
        (right_n, qend,
         (send + 1) if not minus else (send - right_n),
         (send + right_n) if not minus else (send - 1)),
    ):
        if n <= 0:
            continue
        _wc = getattr(subjects, "window_clamped", None)
        if _wc is None:
            # Duck-typed provider with only the all-or-nothing window().
            # Use its semantics as they are rather than assume a shape.
            _w = subjects.window(sseqid, s_lo, s_hi)
            got = None if _w is None else (_w, n - len(_w))
        else:
            got = _wc(sseqid, s_lo, s_hi)
        if got is None:
            # The SUBJECT could not be retrieved. Genuinely unknown, so the whole
            # hit degrades to the conservative count rather than guessing.
            return None
        win, missing = got
        if missing:
            # The flank runs off the end of the subject record. This is NOT
            # unknown: there is no sequence there, so the probe cannot base-pair
            # at those positions and they are mismatches by definition. Counting
            # them as matches (the conservative path) would make an off-target
            # look closer than it is and reject a probe that is fine. It is also the COMMON case, not an edge case:
            # any hit within 25 nt of a transcript end has one.
            total += missing
        if len(win) + missing != n:
            return None                     # arithmetic disagrees; do not guess
        if minus:
            win = _reverse_complement(win)
        # Align the comparison to the part that actually exists. On a left flank
        # the missing positions are the leftmost ones, so the retrieved bases
        # correspond to the RIGHTMOST query positions of that flank.
        q_part = query_seq[q_lo:q_lo + n].upper()
        if missing and n - len(win) > 0:
            q_part = q_part[n - len(win):] if q_lo == 0 else q_part[:len(win)]
        total += sum(1 for a, b in zip(q_part, win) if a != b)
    return total


def full_half_mismatch_detail(
    hit: Dict[str, Any],
    half_len: int,
    *,
    query_seq: str = "",
    subjects: Any = None,
) -> Tuple[int, str]:
    """``(mismatches over the full half, basis)`` — the recount.

    See :func:`full_half_mismatches` for the rule and its rationale. ``basis`` is
    :data:`RECOUNT_EXACT` when every one of the half's positions was actually
    compared, :data:`RECOUNT_CONSERVATIVE` when some could not be and were
    counted as matches.
    """
    internal, internal_exact = _hsp_internal_mismatches(hit)
    flank = _flank_mismatches(hit, half_len, query_seq, subjects)
    if flank is None or not internal_exact:
        return min(internal, half_len), RECOUNT_CONSERVATIVE
    return min(internal + flank, half_len), RECOUNT_EXACT


def full_half_mismatches(
    hit: Dict[str, Any],
    half_len: int,
    *,
    query_seq: str = "",
    subjects: Any = None,
) -> int:
    """Mismatches between a probe half and a BLAST subject over the FULL half.

    10x counts "five mismatches" over the whole 25 nt probe half, but BLAST
    reports a *local* alignment, which may cover only part of the query. Any
    query position that is not identically paired — a reported mismatch, a gap
    column, or a position outside the HSP — is a position where the probe does
    not base-pair, so all three count the same way here.

    This is what makes ``off_target_full_half`` differ from the legacy local
    form: a 25-mer with 6 mismatches clustered at its ends still contains a
    clean 20 nt core, which the legacy >=20 bp / >=80% rule scores as a
    disqualifying off-target even though the full half is 6 mismatches away
    and cannot hybridise.

    The count is made in two parts:

    * INSIDE the HSP, column by column from ``qseq``/``sseq``;
    * OUTSIDE it, by retrieving the subject bases opposite the unaligned query
      positions from the BLAST DB (``blast.SubjectWindows``, i.e. ``blastdbcmd
      -entry <sseqid> -range a-b``) and comparing them to the probe's own bases.

    ``query_seq`` is the probe half itself (the BLAST query), which the hit dict
    does not carry — only the aligned slice of it. ``subjects`` is a
    :class:`blast.SubjectWindows`. Omit either and the count degrades to the
    conservative form below.

    WHY THE COUNT IS NOT INFERRED FROM HSP GEOMETRY (pinned by
    ``test_off_target_full_half_underestimates_mismatches``)
    ------------------------------------------------------------------
    The tempting shortcut is::

        identities = round(length * pident / 100)
        mismatches = half_len - identities        # WRONG

    which assumes the reported HSP spans the half. Often it does not. Under
    blastn-short's default -3 mismatch penalty, a 25-mer sitting 4 mismatches
    from a subject scores 21 - 12 = 9 at full length but 11 as a bare 11 nt
    exact core, and BLAST reports the core. The shortcut then counts the
    unaligned query positions as mismatches, so the true count of 4 reads as
    25 - 11 = 14 and the probe is certified specific. The error is always in
    the unsafe direction: the estimate can only overstate how specific a probe
    is. run_blastn passes -penalty -1, which makes the same hit score
    21 - 4 = 17 at full length, so the full alignment usually wins. Trimming
    is rarer but not impossible, so the recount below stays.

    The blind spot lives in the 3-5 mismatch band. That is where a closely
    related strain sits, and where an unrelated host transcript rarely does,
    so it matters most when screening one viral strain against another.

    Softening ``-penalty`` alone does NOT make the geometry reliable. Trimming
    one edge mismatch stays profitable at -2 (lose 1 match, drop 1 mismatch)
    and score-neutral at -1, so no scoring choice guarantees a full-length HSP,
    and lower penalties start to cost genuinely clean probes.

    THE CONSERVATIVE FALLBACK
    -------------------------
    When the subject window cannot be retrieved — no ``blastdbcmd``, a DB built
    without ``-parse_seqids``, a subject missing from the DB, a flank running off
    the end of the record, a malformed hit — the unaligned query positions are
    counted as MATCHES. That makes the result a LOWER BOUND on the mismatch
    count, so a hit can be called disqualifying when it is not, and the screen
    over-rejects. That is the exact inversion of the geometry error and the
    only safe direction for a specificity screen. The geometry formula is never used as a
    fallback: being wrong in the unsafe direction is worse than being blunt.
    """
    return full_half_mismatch_detail(
        hit, half_len, query_seq=query_seq, subjects=subjects,
    )[0]


# Sentinel subject group for an off-target hit whose subject cannot be
# identified (empty/whitespace sseqid). Its halves are unioned into EVERY other
# subject group as well as standing alone, so an unidentifiable hit can only
# ever make the per-subject screen stricter, never weaker.
UNKNOWN_SUBJECT = "\x00unknown"


def off_target_subject_key(sseqid: str) -> str:
    """Normalise a BLAST subject id to the *gene* the 10x rule talks about.

    10x's sentence is about "matches to off-target **genes**", not about
    individual accessions, so the same gene arriving under several ids must
    collapse to one group. The id is cut at the first header delimiter and
    lower-cased, which folds together:

        ``KR013237.1`` / ``KR013237.2``          -> ``kr013237``  (version)
        ``ENSG0000001|ENST0000009``              -> ``ensg0000001`` (isoform)
        ``OROV_BeAn19991|Large|ACC``             -> ``orov_bean19991``

    Deliberately errs on the side of over-merging. Merging two subjects that
    are really distinct can only move halves INTO a shared group, making
    rejection more likely; splitting one subject across two groups would let a
    genuinely non-specific pair through. Only the safe direction is acceptable
    for a specificity screen, so the aggressive normalisation is the correct
    trade-off here even though it occasionally over-groups.
    """
    sid = (sseqid or "").strip()
    if not sid:
        return UNKNOWN_SUBJECT
    for i, ch in enumerate(sid):
        if ch in _ID_DELIMS:
            sid = sid[:i]
            break
    return sid.lower() or UNKNOWN_SUBJECT


# An Ensembl cDNA/ncRNA header carries the gene it belongs to in a
# ``gene:ENSG...`` token, but the RECORD ID is the transcript — so
# ``off_target_subject_key`` can only ever return a transcript key for it. See
# gene_scoping_degraded.
_ENSEMBL_GENE_TOKEN = "gene:"

# Above this many records a reference is a transcriptome, not a hand-built panel
# of a handful of targets. Only there does "every record is its own subject
# group" mean the grouping failed rather than that the panel really has one
# record per gene. Well below any real transcriptome (GRCh38 cdna.all has
# ~190,000 records) and well above any viral or custom panel this tool is used
# with.
_TRANSCRIPTOME_RECORDS = 1000


def gene_scoping_degraded(headers: Sequence[str]) -> Optional[str]:
    """Why per-gene off-target scoping will not work on this reference, or None.

    ``off_target_same_subject`` (on by default) reads 10x's sentence as being
    about off-target **genes**: a pair is rejected only when BOTH halves have a
    disqualifying match to the SAME gene. The gene is derived from the subject
    id by :func:`off_target_subject_key`, which is all the screen has — it never
    sees an annotation.

    On an Ensembl ``cdna.all`` FASTA, a common whole-transcriptome reference,
    that derivation cannot succeed. The record
    id is ``ENST00000631435.1``; the gene appears only later in the description,
    as ``gene:ENSG...``. So every key is a TRANSCRIPT and "same gene" silently
    becomes "same transcript".

    THE DIRECTION IS PERMISSIVE — it under-rejects. Finer keys split one gene's
    isoforms into separate groups, so a pair whose LHS is disqualifyingly close
    to one isoform and whose RHS is disqualifyingly close to another isoform of
    the SAME gene no longer has both halves in one group, and passes. Those two
    halves are adjacent on a real molecule; the pair ligates and puts reads on a
    gene it was not designed for. (``off_target_subject_key``'s own docstring
    states the rule this is the bad half of: over-merging can only make the
    screen stricter, splitting can let a non-specific pair through.)

    Returns the warning text, or None when the reference is fine or when there
    is nothing to go on. Two signals, both conservative:

    * Ensembl-style headers whose keys are all transcript-level. Precise, and
      it is the common case in practice.
    * A transcriptome-scale reference in which no two records share a key at
      all, so no grouping is happening whatever the id convention is.

    A small custom or viral panel is deliberately NOT flagged: one record per
    gene is what those references genuinely are, and warning on every such run
    would train the reader to ignore the message.
    """
    n = len(headers)
    if n == 0:
        return None
    ensembl = 0
    gene_keyed = 0
    keys = set()
    for h in headers:
        keys.add(off_target_subject_key(h))
        i = h.find(_ENSEMBL_GENE_TOKEN)
        if i < 0:
            continue
        ensembl += 1
        rest = h[i + len(_ENSEMBL_GENE_TOKEN):].split(None, 1)
        gene = rest[0] if rest else ""
        if gene and off_target_subject_key(h) == off_target_subject_key(gene):
            gene_keyed += 1
    if ensembl and gene_keyed == 0:
        return (
            "Off-target scoping is per-TRANSCRIPT, not per-gene, on this "
            f"reference: {ensembl} of its {n} records name their gene only in "
            "an Ensembl 'gene:ENSG...' description field, and the screen groups "
            "off-target hits by RECORD ID (an ENST...), which no record shares. "
            "off_target_same_subject then requires both probe halves to hit the "
            "same TRANSCRIPT, so a pair whose two halves hit two isoforms of one "
            "gene is NOT rejected. That is the permissive direction — the screen "
            "under-rejects. To scope per gene, rebuild the reference with the "
            "gene id first in each header (e.g. '>ENSG...|ENST...')."
        )
    if ensembl == 0 and n >= _TRANSCRIPTOME_RECORDS and len(keys) == n:
        return (
            f"Off-target scoping is per-RECORD on this reference: all {n} of "
            "its record ids reduce to distinct subject keys, so no two records "
            "are ever treated as the same gene. off_target_same_subject then "
            "requires both probe halves to hit the same record, which is the "
            "permissive direction — a pair whose halves hit two records of one "
            "gene is not rejected."
        )
    return None


def off_target_reject_pident(
    off_target_pident_min: float = OFF_TARGET_PIDENT_MIN,
    *,
    full_half: bool = OFF_TARGET_FULL_HALF_DEFAULT,
    min_mismatch: int = OFF_TARGET_MIN_MISMATCH,
    half_len: int = HALF_LENGTH_DEFAULT,
) -> float:
    """The identity (%) at which an off-target match becomes disqualifying.

    In full-half mode this is the identity of a hit sitting exactly on the
    mismatch boundary — ``min_mismatch - 1`` mismatches across the whole half,
    e.g. 4/25 => 84.0%. In legacy local-alignment mode it is simply
    ``off_target_pident_min``. Used to express the reported specificity margin
    in the same units as the rule that is actually in force.
    """
    if not full_half or half_len <= 0:
        return float(off_target_pident_min)
    keep = max(0, half_len - max(0, min_mismatch - 1))
    return 100.0 * keep / half_len


def classify_hit(
    hit: Dict[str, Any],
    target_ids: Set[str],
    *,
    off_target_pident_min: float = OFF_TARGET_PIDENT_MIN,
    off_target_aln_min: int = OFF_TARGET_ALN_MIN,
    on_target_pident_min: float = ON_TARGET_PIDENT_MIN,
    on_target_aln_min: int = ON_TARGET_ALN_MIN,
    off_target_full_half: bool = OFF_TARGET_FULL_HALF_DEFAULT,
    off_target_min_mismatch: int = OFF_TARGET_MIN_MISMATCH,
    off_target_half_len: int = HALF_LENGTH_DEFAULT,
    full_half_mismatch_count: Optional[int] = None,
) -> str:
    """Classify a BLAST hit as 'on_target', 'off_target', or 'ignore'.

    The gate is deliberately the 10x-sourced mismatch rule (CG000621 /
    CG000839: ">= 5 mismatches in at least one of the LHS or RHS probes"), NOT
    bitscore/e-value. Switching the *decision* to bitscore would sever
    traceability to the vendor spec and change which probes pass. Instead the
    bitscore/e-value of the closest near-miss are *reported* per probe (see
    evaluate_specificity / results.csv) so a user can apply their own
    DB-size-aware judgement: the pass/fail rule stays, and the length- and
    DB-robust statistics are surfaced beside it.

    ``off_target_full_half`` (default True) applies the rule as written:
    mismatches are counted over the whole ``off_target_half_len``-nt half and
    the hit is off-target when there are fewer than ``off_target_min_mismatch``
    of them. Set it False for the legacy local-alignment translation
    (``off_target_pident_min`` over ``off_target_aln_min`` bp), which is
    retained for reproducing result sets made with that rule.

    NOTE: this classifies a single (half, subject) alignment. Whether the
    *pair* is rejected is decided in ``evaluate_specificity``, because the 10x
    rule is a disjunction over the two halves.

    ``full_half_mismatch_count`` lets the caller pass in a count it has already
    made (``evaluate_specificity`` does, because it also needs to record whether
    that count was exact). When omitted the count is made here — without the
    probe sequence or a subject-retrieval handle, so it takes the CONSERVATIVE
    path. Callers that want the exact recount must supply the count.
    """
    pident = float(hit.get("pident", 0))
    aln_len = int(hit.get("length", 0))
    sseqid = hit.get("sseqid", "")

    is_target = any(id_matches(sseqid, tid) for tid in target_ids)

    if is_target:
        if pident >= on_target_pident_min and aln_len >= on_target_aln_min:
            return "on_target"
        return "ignore"
    if off_target_full_half:
        n_mm = (full_half_mismatch_count if full_half_mismatch_count is not None
                else full_half_mismatches(hit, off_target_half_len))
        if n_mm < off_target_min_mismatch:
            return "off_target"
        return "ignore"
    if pident >= off_target_pident_min and aln_len >= off_target_aln_min:
        return "off_target"
    return "ignore"


def _matched_strains(
    hits: List[Dict[str, Any]],
    target_ids: Set[str],
    *,
    on_target_pident_min: float = ON_TARGET_PIDENT_MIN,
    on_target_aln_min: int = ON_TARGET_ALN_MIN,
) -> Set[str]:
    """Return the set of target IDs that have an on-target hit."""
    matched: Set[str] = set()
    for h in hits:
        sid = h.get("sseqid", "")
        for tid in target_ids:
            if id_matches(sid, tid):
                pident = float(h.get("pident", 0))
                aln_len = int(h.get("length", 0))
                if pident >= on_target_pident_min and aln_len >= on_target_aln_min:
                    matched.add(tid)
    return matched


def _worst_junction_duplex(
    junction_records: List[Dict[str, Any]],
    target_ids: Set[str],
    excused: Any,
) -> Tuple[Optional[int], Optional[int]]:
    """Longest junction-spanning clean duplex against a NON-target gene.

    ``(None, None)`` when no junction record reached us at all — the scan did
    not run, which is not the same claim as "it ran and found nothing" and must
    not be flattened into ``(0, 0)``. ``(0, 0)`` is what a scan that found no
    non-target site returns, and the record list is the only evidence that
    separates the two.

    Target and excused subjects are dropped HERE rather than in the screen
    because the screen has no idea which gene a probe belongs to: under
    gene-aware mode the on-target set is per probe (its own gene's transcripts),
    and the probe's own transcript is by construction a perfect 25+25 site. A
    scan that filtered by itself would either reject every probe against its own
    gene or have to be told the target set twice.

    "Worst" is the largest ``min(left, right)``, ties broken on the longer total
    run. It is the minimum that decides, because the gate is per side.
    """
    if not junction_records:
        return (None, None)
    best: Tuple[int, int] = (0, 0)
    best_key = (-1, -1)
    for rec in junction_records:
        for site in rec.get("sites") or ():
            try:
                sseqid, left, right = str(site[0]), int(site[1]), int(site[2])
            except (IndexError, TypeError, ValueError):
                continue
            if any(id_matches(sseqid, tid) for tid in target_ids):
                continue
            if excused({"sseqid": sseqid}):
                continue
            key = (min(left, right), left + right)
            if key > best_key:
                best_key = key
                best = (left, right)
    return best


def evaluate_specificity(
    candidate: ProbeCandidate,
    hits: List[Dict[str, Any]],
    target_ids: Set[str],
    strain_rule: StrainRule = StrainRule.EXACT_SET,
    required_strains: Optional[Set[str]] = None,
    off_target_pident_min: float = OFF_TARGET_PIDENT_MIN,
    off_target_aln_min: int = OFF_TARGET_ALN_MIN,
    on_target_pident_min: float = ON_TARGET_PIDENT_MIN,
    on_target_aln_min: int = ON_TARGET_ALN_MIN,
    extra_on_target_ids: Optional[Set[str]] = None,
    require_on_target: bool = False,
    off_target_full_half: bool = OFF_TARGET_FULL_HALF_DEFAULT,
    off_target_min_mismatch: int = OFF_TARGET_MIN_MISMATCH,
    off_target_half_len: int = HALF_LENGTH_DEFAULT,
    off_target_require_both_halves: bool = OFF_TARGET_REQUIRE_BOTH_HALVES_DEFAULT,
    off_target_same_subject: bool = OFF_TARGET_SAME_SUBJECT_DEFAULT,
    off_target_junction_duplex_bp: int = OFF_TARGET_JUNCTION_DUPLEX_BP,
    # Bases each side of the nick that count as junction-proximal when the
    # closest off-target is profiled (FilterThresholds.junction_window,
    # --junction-window).
    junction_window: int = JUNCTION_WINDOW,
    subjects: Any = None,
) -> None:
    """Evaluate specificity for a single probe candidate (in-place).

    BLAST is primarily an **off-target screen**: a probe is REJECTED when its
    off-target hits satisfy the 10x rejection rule below. An on-target BLAST
    hit is NOT required by default — a probe is tiled from its target, so its
    on-target binding is given; the DBs the user supplies are things the probe
    must NOT hit. Requiring an on-target hit would produce ZERO probes whenever
    the target is not itself in a supplied DB, or when the DB names the target
    differently from the input FASTA record.

    Off-target rejection rule
    -------------------------
    CG000621 Rev D / CG000839 Rev A, verbatim: *"Matches to off-target genes
    should have at least five mismatches in at least one of the LHS or RHS
    probes to prevent efficient hybridization."*

    With ``off_target_require_both_halves=True`` (the default, and what the
    sentence says) the pair is rejected only when **both** halves have an
    off-target match closer than the mismatch floor. This is the biologically
    correct rule for a ligation chemistry: RTL only produces a library molecule
    when both halves hybridise and the ligase seals the nick between them, so a
    pair whose LHS is >= 5 mismatches from an off-target generates no ligation
    product there and therefore no false signal — rejecting it discards a
    perfectly specific probe.

    Set it False for the either-half reading, in which a hit on EITHER half
    disqualifies the pair. That is strictly more conservative but is not the
    published rule; it is kept for reproducing result sets made with it.

    ``off_target_same_subject`` (default True) scopes that
    disjunction per off-target GENE, which is what 10x's sentence says: "matches
    to off-target **genes**". An LHS that binds gene X and an RHS that binds gene
    Y are never adjacent on one molecule, so the ligase has nothing to seal and
    no product forms on either. A pair is therefore rejected only when both
    halves have a disqualifying match to the SAME subject (normalised by
    ``off_target_subject_key``). Set it False to score the disjunction over the
    whole probe instead, kept for reproducing result sets made with that
    reading. It has no effect when ``off_target_require_both_halves`` is
    False, since that mode never consults the second half at all.

    The default is also consistent with 10x's own choices. On 10x's public Flex
    v1 human probe panel, screened against GENCODE, the per-gene reading
    rejects far fewer of the probes 10x ships than the per-probe reading does,
    while rejecting almost as many of the probes 10x marks ``included=FALSE``.
    See :data:`constants.OFF_TARGET_SAME_SUBJECT_DEFAULT` for why per-LOCUS
    scoping is NOT a further improvement.

    Junction-duplex gate (OURS, not 10x's)
    --------------------------------------
    ``off_target_junction_duplex_bp`` (default 9, 0 disables) adds a second,
    STRICTER rejection on top of the count rule above: a pair is rejected when
    some non-target gene presents at least that many perfectly matched bases on
    BOTH sides of the ligation nick, in contiguous reference sequence. See
    :func:`junction_clean_duplex` for what is measured and
    :data:`constants.OFF_TARGET_JUNCTION_DUPLEX_BP` for why 9.

    It is applied AFTER the count rule, never instead of it, and only to pairs
    the count rule PASSED — so a probe carrying
    :data:`STATUS_FAILED_JUNCTION_DUPLEX` is by construction one that 10x's own
    sentence would have shipped. This gate can only ever reject more than the
    published rule, never less, which is why it can be on by default without
    contradicting the spec it is built on.

    The evidence is recorded on the candidate whether or not the gate is
    enabled (``offtarget_clean_left_bp`` / ``offtarget_clean_right_bp``), and
    -1 there means the screen never measured it — the BLAST engine emits no
    junction records, so a BLAST run reports -1 and this gate never fires.

    Half attribution comes from the ``__L`` / ``__R`` query-id suffix (see
    ``hit_half``). A hit whose half cannot be determined — a hand-built fixture,
    or a cached hit file from a whole-tile query — is counted against BOTH
    halves, so an un-attributable hit can never make the screen *weaker* than
    the either-half reading.

    ``require_on_target=True`` restores the strict mode (the probe must have an
    on-target hit and satisfy the strain rule) — for diagnostics where the
    target IS in the DB under a matching id.

    ``extra_on_target_ids`` are accessions excused from off-target rejection
    (e.g. the target gene's other isoforms).

    ``subjects`` is a :class:`blast.SubjectWindows` used to retrieve the subject
    bases opposite the *unaligned* part of each probe half, so the full-half
    mismatch count is a real count rather than an inference from HSP geometry
    (see ``full_half_mismatches``). Pass None — or omit it — and every hit that
    BLAST trimmed falls back to the conservative count. Which path each probe
    took is recorded on ``candidate.offtarget_recount_basis``; a run is only
    exactly recounted if that says so.
    """
    # The junction and unscreenable-half records travel in the same list and are
    # not alignments; take them out BEFORE anything counts a mismatch. See
    # split_screen_records.
    hits, junction_records, unscreenable = split_screen_records(hits)
    # FAIL CLOSED. The screen refused to search one of this pair's halves, so
    # the 10x rule below cannot be evaluated: its whole shape is "at least one
    # half is far from every off-target", and a half that was never searched is
    # not evidence of distance. Returning "passed" here is the exact failure
    # this record exists to prevent, so the refusal becomes the verdict.
    if unscreenable:
        candidate.specificity_pass = False
        candidate.specificity_margin = 0.0
        candidate.specificity_status = STATUS_UNSCREENABLE_HALF
        return
    on_target_scores: List[float] = []
    off_target_scores: List[float] = []
    # Off-target identities split by which probe half produced the hit. The 10x
    # rule is a disjunction over these two sets, so they cannot be merged.
    off_by_half: Dict[str, List[float]] = {"L": [], "R": []}
    # subject group -> set of halves with a disqualifying match to it. Only
    # consulted under off_target_same_subject; see off_target_subject_key.
    off_halves_by_subject: Dict[str, Set[str]] = {}
    # Near-miss off-targets: non-target hits that clear the alignment-length
    # floor but fall *just short* on identity (so they're "ignore", not
    # disqualifying). The best one measures how close the probe came to being
    # non-specific. It sets a discriminating specificity margin and is
    # reported as off-target evidence.
    nearmiss_off_hits: List[Dict[str, Any]] = []

    def _excused(hit: Dict[str, Any]) -> bool:
        if not extra_on_target_ids:
            return False
        _sid = hit.get("sseqid", "")
        return any(id_matches(_sid, iid) for iid in extra_on_target_ids)

    # The BLAST query for a hit is one 25 nt HALF of the candidate's 50 bp
    # window (pipeline._run_screen writes them as ``<probe_id>__L`` /
    # ``__R``). The hit dict carries only the ALIGNED slice of that query, so the
    # recount needs the half itself to compare the positions BLAST skipped.
    _seq = candidate.sequence or ""
    _mid = len(_seq) // 2
    _query_halves = {"L": _seq[:_mid], "R": _seq[_mid:]}

    def _query_seq(hit: Dict[str, Any]) -> str:
        # A hit with no __L/__R suffix cannot be attributed to a half, so there
        # is no query to compare against and the recount degrades to
        # conservative — consistent with how ``hit_half`` says such hits are
        # handled everywhere else.
        return _query_halves.get(hit_half(hit), "")

    # Basis of the counts made for THIS probe. Starts optimistic and is demoted
    # by the first hit that could not be recounted exactly: a probe is only
    # "exact" if every off-target count behind its verdict was.
    recount_basis = (RECOUNT_EXACT if off_target_full_half
                     else RECOUNT_LEGACY_LOCAL)
    n_fallback_hits = 0
    # Smallest fraction of a probe half that BLAST actually ALIGNED, across the
    # non-target hits behind this verdict. `offtarget_recount_basis` says whether
    # unexamined positions had to be reconstructed; this says how much of the
    # half was unexamined in the first place, which is the other half of the same
    # question. A probe that passed with every hit covering all 25 nt is a
    # stronger pass than one that passed on a 9 nt core plus a database lookup,
    # and this field lets results.csv tell them apart.
    #
    # Deliberately NOT read from a BLAST `qcovhsp` column: it is exactly
    # (qend - qstart + 1) / half_len by BLAST's own definition, so asking for
    # the column would buy a
    # derivable number at the cost of a wider row on every hit — at precisely the
    # volumes where hit I/O is the binding constraint.
    min_qcov: Optional[float] = None
    # (full-half mismatches, junction-proximal mismatches) of the closest
    # disqualifying off-target. See hit_junction_mismatches for why the second
    # number is protective rather than harmful on an off-target.
    worst_off: Optional[Tuple[int, Optional[int]]] = None

    for hit in hits:
        n_mm: Optional[int] = None
        if off_target_full_half:
            n_mm, basis = full_half_mismatch_detail(
                hit, off_target_half_len,
                query_seq=_query_seq(hit), subjects=subjects,
            )
            if basis != RECOUNT_EXACT:
                # Only demote on hits that are NOT on-target: an on-target hit
                # never goes through the full-half rule, so an inexact count for
                # it says nothing about the off-target verdict.
                if not any(id_matches(hit.get("sseqid", ""), tid)
                           for tid in target_ids):
                    recount_basis = RECOUNT_CONSERVATIVE
                    n_fallback_hits += 1
        cls = classify_hit(
            hit, target_ids,
            off_target_pident_min=off_target_pident_min,
            off_target_aln_min=off_target_aln_min,
            on_target_pident_min=on_target_pident_min,
            on_target_aln_min=on_target_aln_min,
            off_target_full_half=off_target_full_half,
            off_target_min_mismatch=off_target_min_mismatch,
            off_target_half_len=off_target_half_len,
            full_half_mismatch_count=n_mm,
        )
        if cls == "on_target":
            on_target_scores.append(float(hit.get("pident", 0)))
        elif cls == "off_target":
            if _excused(hit):
                continue  # a known isoform / accepted cross-hit — not off-target
            pid = float(hit.get("pident", 0))
            _cov = hit_query_coverage(hit, off_target_half_len)
            if _cov is not None and (min_qcov is None or _cov < min_qcov):
                min_qcov = _cov
            # Keep the profile of the CLOSEST off-target — the one that decides
            # how much trouble this probe is in. Fewer mismatches is closer; ties
            # break toward the more dangerous (fewer protective junction ones).
            _jm = hit_junction_mismatches(hit, off_target_half_len, junction_window)
            _n = n_mm if n_mm is not None else off_target_min_mismatch
            if (worst_off is None or _n < worst_off[0]
                    or (_n == worst_off[0] and (_jm or 0) < (worst_off[1] or 0))):
                worst_off = (_n, _jm)
            off_target_scores.append(pid)
            which = hit_half(hit)
            # No suffix => can't attribute; charge it to both halves so an
            # un-attributable hit is never treated more leniently.
            halves = (which,) if which else ("L", "R")
            for k in halves:
                off_by_half[k].append(pid)
            off_halves_by_subject.setdefault(
                off_target_subject_key(hit.get("sseqid", "")), set()
            ).update(halves)
        else:  # "ignore" — capture non-target near-misses over the length floor
            sid = hit.get("sseqid", "")
            is_target = any(id_matches(sid, tid) for tid in target_ids)
            if is_target or _excused(hit):
                continue
            # Tracked BEFORE the near-miss length floor. A trimmed core sits
            # below that floor by definition, and it is exactly the hit whose
            # coverage matters — it is the one whose "distant enough" verdict
            # rested on reconstructed bases rather than on aligned ones.
            _cov = hit_query_coverage(hit, off_target_half_len)
            if _cov is not None and (min_qcov is None or _cov < min_qcov):
                min_qcov = _cov
            if int(hit.get("length", 0)) >= off_target_aln_min:
                nearmiss_off_hits.append(hit)

    # Set BEFORE any of the early returns below, so every exit path — rejected,
    # no-on-target, passed — carries an honest account of how its counts were
    # made. A verdict whose provenance depends on which branch it left by is
    # exactly the kind of gap this field exists to close.
    candidate.offtarget_recount_basis = recount_basis
    # Same reasoning as the basis above: set before every early return, so a
    # rejected probe reports the coverage of the evidence that rejected it and a
    # passing one reports the coverage of the evidence it survived. NaN means
    # there were no non-target hits at all — the strongest kind of pass, and not
    # the same thing as coverage of zero.
    candidate.offtarget_min_qcov = (
        float("nan") if min_qcov is None else float(min_qcov)
    )
    # The junction-duplex evidence, and the verdict it implies. Measured here —
    # beside the recount basis and for the same reason — so that a probe
    # rejected by the COUNT rule below still reports what its junction looked
    # like. Left as the -1 model default when no junction record reached us,
    # which is how a BLAST run (or a hit cache without junction records) says
    # "not measured".
    jd_left, jd_right = _worst_junction_duplex(
        junction_records, target_ids, _excused,
    )
    if jd_left is not None and jd_right is not None:
        candidate.offtarget_clean_left_bp = int(jd_left)
        candidate.offtarget_clean_right_bp = int(jd_right)
    if worst_off is not None:
        candidate.offtarget_worst_mismatches = int(worst_off[0])
        candidate.offtarget_worst_junction_mm = (
            -1 if worst_off[1] is None else int(worst_off[1])
        )
    if subjects is not None and off_target_full_half:
        # Aggregate only; the caller logs ONE line per run. A whole-
        # transcriptome screen has tens of millions of hits, so per-hit logging
        # is itself an outage.
        subjects.n_fallback += n_fallback_hits
        subjects.n_exact += max(0, len(hits) - n_fallback_hits)

    _strain_kw: Dict[str, Any] = dict(
        on_target_pident_min=on_target_pident_min,
        on_target_aln_min=on_target_aln_min,
    )

    # 1) Core off-target screen. Per CG000621/CG000839 the pair survives as long
    #    as at least ONE half is >= off_target_min_mismatch away from every
    #    off-target — that half fails to hybridise, the nick never forms, and no
    #    ligation product (hence no signal) is generated. So reject only when
    #    BOTH halves carry a disqualifying off-target match.
    if off_target_require_both_halves and off_target_same_subject:
        # Stricter-literal reading: "matches to off-target GENES" scoped per
        # gene. Both halves must be close to the SAME subject before a ligation
        # product is possible. Halves whose subject could not be identified are
        # unioned into every group (and stand alone), so an attribution gap can
        # only make this stricter — never let a non-specific pair through.
        unknown = off_halves_by_subject.get(UNKNOWN_SUBJECT, set())
        rejected = any(
            len(halves | unknown) == 2
            for subject, halves in off_halves_by_subject.items()
            if subject != UNKNOWN_SUBJECT
        ) or len(unknown) == 2
    elif off_target_require_both_halves:
        rejected = bool(off_by_half["L"]) and bool(off_by_half["R"])
    else:
        # Legacy: a hit on EITHER half counts against the probe.
        rejected = bool(off_target_scores)
    if rejected:
        candidate.specificity_pass = False
        best_off = max(off_target_scores)
        best_on = max(on_target_scores) if on_target_scores else 0.0
        candidate.specificity_margin = best_on - best_off
        candidate.specificity_status = "failed_off_target"
        return

    # 1b) Junction-spanning clean duplex. OURS, not 10x's, and applied only to
    #     pairs the count rule above has already passed — so this can only ever
    #     reject MORE than the published sentence, never less, and the status it
    #     writes identifies exactly the population that distinction creates: a
    #     probe 10x's own rule would have shipped.
    #
    #     No `off_target_require_both_halves` / `off_target_same_subject`
    #     equivalent is needed here, and adding one would be wrong. Both of those
    #     knobs exist to approximate "are these two halves adjacent on one
    #     molecule?" from a pair of independent per-half alignments; this
    #     measurement is taken in one unbroken stretch of reference on one
    #     record, so adjacency is not inferred, it is the coordinate system.
    if (off_target_junction_duplex_bp > 0
            and jd_left is not None and jd_right is not None
            and min(jd_left, jd_right) >= off_target_junction_duplex_bp):
        candidate.specificity_pass = False
        # Not the count rule's margin: nothing here rests on percent identity,
        # and reporting best_on - best_off would attribute this rejection to a
        # number that did not make it. 0.0 is the honest floor for a probe that
        # can ligate on a gene it was not designed for.
        candidate.specificity_margin = 0.0
        candidate.specificity_status = STATUS_FAILED_JUNCTION_DUPLEX
        return
    if off_target_scores:
        # The pair passed despite a disqualifying off-target on at least one
        # half, so it must not keep the full-credit margin of a genuinely clean
        # probe. Report the non-target hits of whichever half(s) the pair's
        # specificity now rests on.
        #
        # Normally exactly one half is dirty and the other carries the pair, so
        # the clean half's own worst off-target is the honest near-miss. Under
        # off_target_same_subject BOTH halves can be dirty (to different genes)
        # and the pair still survives; then neither half is a fallback, so every
        # non-target hit counts and the margin collapses toward 0 — which is the
        # correct signal that this probe survives only on the per-subject
        # reading.
        if off_by_half["L"] and off_by_half["R"]:
            carrying = {"L", "R"}
        else:
            carrying = {"R" if off_by_half["L"] else "L"}
        nearmiss_off_hits = [
            h for h in hits
            if hit_half(h) in carrying and not _excused(h)
            and not any(id_matches(h.get("sseqid", ""), tid) for tid in target_ids)
        ]

    # 2) Off-target clean. Only enforce on-target / strain coverage when the
    #    caller asked for it (strict diagnostic mode). Otherwise the probe
    #    passes the off-target screen regardless of whether the target was in a
    #    supplied DB.
    if require_on_target:
        if not on_target_scores:
            candidate.specificity_pass = False
            candidate.specificity_margin = 0.0
            candidate.specificity_status = "no_on_target_hit"
            return
        matched = _matched_strains(hits, target_ids, **_strain_kw)
        ok = True
        if strain_rule == StrainRule.EXACT_SET:
            ok = (required_strains or target_ids).issubset(matched)
        elif strain_rule == StrainRule.PAN_STRAIN:
            ok = bool(matched)
        elif strain_rule == StrainRule.STRAIN_UNIQUE:
            ok = len(matched) == 1
        if not ok:
            candidate.specificity_pass = False
            candidate.specificity_margin = 0.0
            candidate.specificity_status = "no_on_target_hit"
            return

    candidate.specificity_pass = True
    candidate.specificity_status = "passed"
    # Record quantitative off-target evidence: the closest near-miss (highest
    # identity) and how many near-misses there were.
    candidate.offtarget_nearmiss_count = len(nearmiss_off_hits)
    best_hit = (max(nearmiss_off_hits, key=lambda h: float(h.get("pident", 0)))
                if nearmiss_off_hits else None)
    best_nearmiss = float(best_hit.get("pident", 0)) if best_hit else 0.0
    if best_hit is not None:
        candidate.best_offtarget_pident = best_nearmiss
        candidate.best_offtarget_bitscore = float(best_hit.get("bitscore", 0))
        candidate.best_offtarget_evalue = float(best_hit.get("evalue", 0))
    # Margin = how far the closest near-miss off-target sits BELOW the
    # disqualification threshold. No near-miss at all → full margin. A near-miss
    # right at the threshold → ~0. This discriminates among passing probes
    # (best_on - best_off would be ~0 for every passing probe in the default
    # off-target-only mode). The threshold has to follow whichever
    # off-target rule is in force, else the margin is measured against a number
    # that no longer decides anything.
    candidate.specificity_margin = max(
        0.0,
        off_target_reject_pident(
            off_target_pident_min,
            full_half=off_target_full_half,
            min_mismatch=off_target_min_mismatch,
            half_len=off_target_half_len,
        ) - best_nearmiss,
    )


def batch_evaluate(
    candidates: List[ProbeCandidate],
    hits_by_probe: Dict[str, List[Dict[str, Any]]],
    target_ids: Set[str],
    strain_rule: StrainRule = StrainRule.EXACT_SET,
    required_strains: Optional[Set[str]] = None,
    *,
    off_target_pident_min: float = OFF_TARGET_PIDENT_MIN,
    off_target_aln_min: int = OFF_TARGET_ALN_MIN,
    on_target_pident_min: float = ON_TARGET_PIDENT_MIN,
    on_target_aln_min: int = ON_TARGET_ALN_MIN,
    extra_on_target_ids: Optional[Set[str]] = None,
    require_on_target: bool = False,
    off_target_full_half: bool = OFF_TARGET_FULL_HALF_DEFAULT,
    off_target_min_mismatch: int = OFF_TARGET_MIN_MISMATCH,
    off_target_half_len: int = HALF_LENGTH_DEFAULT,
    off_target_require_both_halves: bool = OFF_TARGET_REQUIRE_BOTH_HALVES_DEFAULT,
    off_target_same_subject: bool = OFF_TARGET_SAME_SUBJECT_DEFAULT,
    off_target_junction_duplex_bp: int = OFF_TARGET_JUNCTION_DUPLEX_BP,
    junction_window: int = JUNCTION_WINDOW,
    subjects: Any = None,
    gene_of: Optional[Dict[str, str]] = None,
    group_members: Optional[Dict[str, Set[str]]] = None,
) -> int:
    """Evaluate specificity for many candidates. Returns count passing.

    Gene-aware mode (``gene_of`` + ``group_members`` both given): instead of one
    global ``target_ids`` set, each probe's on-target set is the transcripts of
    *its own gene* (``group_members[gene_of[probe.gene_id]]``). Sibling isoforms
    therefore never count as off-target, and only cross-gene hits disqualify a
    probe. Each candidate is also annotated with its isoform coverage
    (``isoform_hits`` / ``isoform_total`` / ``isoform_coverage``) — how many of
    its gene's transcripts it binds — for downstream scoring / gating.
    """
    gene_aware = gene_of is not None and group_members is not None
    _off_kw: Dict[str, Any] = dict(
        off_target_pident_min=off_target_pident_min,
        off_target_aln_min=off_target_aln_min,
        on_target_pident_min=on_target_pident_min,
        on_target_aln_min=on_target_aln_min,
        extra_on_target_ids=extra_on_target_ids,
        require_on_target=require_on_target,
        off_target_full_half=off_target_full_half,
        off_target_min_mismatch=off_target_min_mismatch,
        off_target_half_len=off_target_half_len,
        off_target_require_both_halves=off_target_require_both_halves,
        off_target_same_subject=off_target_same_subject,
        off_target_junction_duplex_bp=off_target_junction_duplex_bp,
        junction_window=junction_window,
        subjects=subjects,
    )
    passed = 0
    for c in candidates:
        probe_hits = hits_by_probe.get(c.probe_id, [])
        if gene_aware:
            # gene_aware already implies both are non-None; assert so the type
            # checker narrows the Optionals (no behaviour change).
            assert gene_of is not None and group_members is not None
            gkey = gene_of.get(c.gene_id, c.gene_id)
            siblings = group_members.get(gkey) or {c.gene_id}
            evaluate_specificity(
                c, probe_hits, siblings, strain_rule, required_strains,
                **_off_kw,
            )
            matched = _matched_strains(
                # Alignments only: a junction record carries no identity or
                # length and would never match, but it must not be walked as
                # though it might.
                split_junction_records(probe_hits)[0], siblings,
                on_target_pident_min=on_target_pident_min,
                on_target_aln_min=on_target_aln_min,
            )
            # The probe is tiled from its own transcript, so it binds that one
            # by construction even if the DB lacks a self-record to BLAST against.
            matched.add(c.gene_id)
            total = len(siblings) or 1
            bound = matched & siblings
            c.isoform_total = total
            c.isoform_hits = len(bound)
            c.isoform_coverage = c.isoform_hits / total
            # Persist *which* isoforms so downstream selection needs no FASTA.
            c.isoform_hit_ids = sorted(bound)
        else:
            evaluate_specificity(
                c, probe_hits, target_ids, strain_rule, required_strains,
                **_off_kw,
            )
        if c.specificity_pass:
            passed += 1
    return passed
