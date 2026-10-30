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

"""Results table writer (CSV + optional XLSX).

One row per ranked tile. Columns include the original mRNA target window
(sense strand, as it appears in the input FASTA), the 50-bp probe-binding
sequence (reverse complement), the two 25-bp target halves, the assembled
LHS/RHS oligos, GC stats, melting temperatures, and the ranking outputs.
"""
from __future__ import annotations
import csv
import os
import tempfile
from typing import List

from ..constants import SCREEN_ENGINE_DEFAULT
from ..models import ProbeCandidate
from ..errors import WriterError
from .._utils import fmt_tm as _tm, is_nan as _is_nan


def _fmt_qcov(value) -> str:
    """Render offtarget_min_qcov: ``n/a`` when there were no non-target hits.

    Same convention as fmt_tm — the sentinel must stay distinguishable from a
    real 0.0, because "nothing to compare against" is the strongest pass and
    "BLAST aligned none of the half" would be the weakest.
    """
    if value is None or _is_nan(value):
        return "n/a"
    try:
        return f"{float(value):.4f}"
    except (TypeError, ValueError):
        return "n/a"


BASE_HEADERS = [
    "rank", "probe_id", "gene_id", "start", "end",
    # --- the user's requested "original gene/genome sequence" + "matching probe seq (RC)" ---
    "mrna_target_50bp", "probe_rc_50bp",
    "target_LHS_25bp", "target_RHS_25bp",
    "lhs_oligo", "rhs_oligo",
    # --- per-half stats ---
    "gc_content_full", "gc_lhs_half", "gc_rhs_half",
    "tm_target_lhs", "tm_target_rhs", "tm_full_50",
    # junction_base = the mRNA base opposite the LHS probe's 3' end (index 25);
    # tn_junction_ok = that base is 'A' (a true TN motif). Reported for every
    # probe whatever TN gate ran, so a reviewer can see the junction identity
    # of probes admitted under the soft ('not C') gate.
    "homopolymer_max", "junction_base", "tn_junction_ok",
    # --- specificity / ranking ---
    # Which search engine produced the hits behind every specificity column
    # below. Constant across a run, and repeated on every row on purpose: a
    # results.csv is read on its own, detached from the manifest that would
    # otherwise carry it, and the engine changes what the neighbouring columns
    # MEAN. Under 'kmer', offtarget_recount_basis is 'exact' by construction
    # rather than by successful retrieval, offtarget_min_qcov is always 1.0,
    # best_offtarget_bitscore/_evalue are undefined (nan), and
    # offtarget_nearmiss_count covers a narrower mismatch band. Two runs cannot
    # be compared on those columns without first comparing this one.
    "screen_engine",
    "specificity_pass", "specificity_status", "specificity_margin",
    # How the off-target mismatch counts behind that verdict were made:
    # exact / conservative / legacy_local / not_screened. "conservative" means
    # the counts are LOWER BOUNDS (a subject window could not be retrieved, so
    # unaligned probe positions were counted as matches) — the screen may
    # over-reject. Emitted so a results.csv can never be read as exactly
    # recounted when it was only bounded.
    "offtarget_recount_basis",
    "offtarget_min_qcov",
    "offtarget_worst_mismatches",
    "offtarget_worst_junction_mm",
    # The junction-spanning clean duplex against the worst NON-target gene:
    # matched bases immediately 5' of the ligation nick and immediately 3' of
    # it, in contiguous reference sequence. This is the evidence behind a
    # `failed_junction_duplex` verdict, and it is emitted for every probe,
    # including ones the gate was not asked to judge.
    #
    # Same -1 convention as offtarget_worst_junction_mm above, and for the same
    # reason: -1 is "never measured" (the BLAST engine emits no junction scan,
    # and neither does a hit cache written before 0.6.1), which must stay
    # distinguishable from a measured 0. It is NOT interchangeable with a blank:
    # a blank in this file means an inapplicable aggregate (see _num), and this
    # column applies to every row.
    "offtarget_clean_left_bp", "offtarget_clean_right_bp",
    # Quantitative off-target evidence for a passing probe:
    # closest near-miss off-target identity / bitscore / e-value + how many.
    "offtarget_nearmiss_count", "best_offtarget_pident",
    "best_offtarget_bitscore", "best_offtarget_evalue",
    # --- isoform coverage (group-by-gene mode; 1/1/1.0 when not grouped) ---
    "isoform_coverage", "isoform_hits", "isoform_total",
    # The specific transcript ids the probe binds (';'-joined). Empty when not
    # gene-aware. Lets the per-gene selection run without the input FASTA.
    "isoform_hit_ids",
    # composite_score is a NON-CALIBRATED preference score (weighted sum of the
    # score_* sub-scores below), not a probability of probe success.
    "composite_score", "reduced_confidence",
    # Optional primer3 thermodynamics (0 unless primer3_scoring was on).
    # Relative indicators, not in-assay values.
    "tm_balance", "self_dimer_dg", "hetero_dimer_dg",
    # Per-component 0-1 sub-scores that make up composite_score.
    "score_gc_proximity", "score_gc_balance", "score_specificity_margin",
    "score_homopolymer_penalty", "score_central_placement",
    "score_isoform_coverage",
    # TN-junction gate. MULTIPLIES the weighted sum above rather than adding to
    # it; 1.0 (inert) whenever the junction filter is in its default hard mode.
    "score_junction_base",
    # Bounded primer3 nudge added on top (blank unless primer3_scoring was on).
    # In [-0.20, +0.05] before the junction gate scales it; see ranking.py.
    "score_primer3_adjustment",
    # --- Degenerate (consensus) probes: the parent / variant model ---
    # APPENDED, never inserted: every column above keeps its position so a
    # spreadsheet, script or notebook built against an earlier results.csv
    # still reads. Blank/0 on every row of a run that used no tolerance.
    #   probe_role           — "" (plain), "consensus" (never orderable) or
    #                          "variant" (the concrete ACGT oligo)
    #   parent_probe_id      — the consensus row a variant resolves from
    #   resolves_to          — ';'-joined record ids this variant matches
    #   degenerate_positions — ','-joined 0-based offsets into the 50-mer
    #   iupac_codes          — the code at each of those offsets, same order
    #   n_variants           — distinct variants under the parent
    "probe_role", "parent_probe_id", "resolves_to",
    "degenerate_positions", "iupac_codes", "n_variants",
    # Parent aggregates over its variants — ranges where the question is how
    # much the ambiguity costs, worst case where it is whether the pool works.
    # Blank on plain and variant rows; hairpin_dg is per-row, not an aggregate.
    "tm_lhs_min", "tm_lhs_max", "tm_rhs_min", "tm_rhs_max",
    "gc_lhs_min", "gc_lhs_max", "gc_rhs_min", "gc_rhs_max",
    "hairpin_dg", "hairpin_dg_worst", "hetero_dimer_dg_worst", "tm_balance_max",
    # ';'-joined input record ids this window was designed to COVER. Populated
    # only by constitutive (isoform-collapse) mode, where each probe is the
    # shared sequence of its largest isoform group: below the quorum that group
    # is a strict subset and this is the only place a reader learns which
    # isoforms the probe is actually for. Blank on every other row.
    #
    # Deliberately NOT folded into resolves_to, which sits under the parent /
    # variant model above and is documented as the record ids a VARIANT row
    # matches. These probes are plain rows (probe_role ""), so reusing that
    # column would put a non-blank resolves_to on a row no reader expects one
    # on, and would narrow the junction-mismatch scope of conserved-mode
    # variants as a side effect. Separate field, separate column.
    "covered_record_ids",
    # --- sub-scores added after the column order was frozen ---
    # Out of order on purpose: it belongs with the other ``score_*`` columns
    # above but cannot be inserted there without shifting every column after it.
    # Where a probe's ON-target mismatches sit relative to the ligation nick,
    # 1.0 = nothing to penalise. Always 1.0 on a single-sequence run, which has
    # no on-target mismatches and where the term is not weighted at all.
    "score_mismatch_placement",
]


def _num(value, fmt: str = "{:.2f}") -> str:
    """Render an aggregate: blank for NaN/None, formatted otherwise.

    Blank rather than ``nan`` because these columns are empty on the majority
    of rows (every plain and variant row) and a column of ``nan`` reads as a
    failed computation rather than an inapplicable one.
    """
    if value is None:
        return ""
    try:
        f = float(value)
    except (TypeError, ValueError):
        return ""
    if f != f:                      # NaN
        return ""
    return fmt.format(f)

# Sub-score component keys (from ProbeCandidate.score_components) in header order.
_SCORE_COMPONENTS = [
    "gc_proximity", "gc_balance", "specificity_margin",
    "homopolymer_penalty", "central_placement", "isoform_coverage",
    "junction_base", "primer3_adjustment",
]
# Sub-scores added after the column order was frozen. They render at the END of
# the base columns, not next to their siblings above, because every column in
# BASE_HEADERS keeps the position it has always had — see the note there.
_SCORE_COMPONENTS_APPENDED = [
    "mismatch_placement",
]

# Per-strain columns: emitted once per strain in the alignment when the
# candidate list comes from a multi-strain run. Header names use the strain
# IDs as a suffix so users can sort / filter by strain in Excel.
_PER_STRAIN_COLS = [
    "mismatches", "lhs_mm", "rhs_mm", "junction_mm",
    "mismatch_positions", "est_dTm", "likely_binds",
]


def _strain_headers(strain_ids):
    """Return the per-strain columns + the cross-strain summary scalars."""
    out = []
    for sid in strain_ids:
        for col in _PER_STRAIN_COLS:
            out.append(f"{col}__{sid}")
    if strain_ids:
        out.extend([
            "min_strain_mismatches",
            "max_strain_mismatches",
            "worst_strain_dTm",
        ])
    return out


def _row(c: ProbeCandidate, strain_ids=None,
         screen_engine: str = SCREEN_ENGINE_DEFAULT):
    base = [
        c.rank, c.probe_id, c.gene_id, c.start, c.end,
        c.sequence, c.probe_rc,
        c.target_lhs, c.target_rhs,
        c.lhs_oligo, c.rhs_oligo,
        f"{c.gc_content:.4f}", f"{c.gc_left:.4f}", f"{c.gc_right:.4f}",
        _tm(c.tm_target_lhs), _tm(c.tm_target_rhs), _tm(c.tm_full_50),
        c.homopolymer_max, c.junction_base, c.tn_junction_ok,
        screen_engine,
        c.specificity_pass, c.specificity_status, f"{c.specificity_margin:.2f}",
        c.offtarget_recount_basis,
        _fmt_qcov(c.offtarget_min_qcov),
        c.offtarget_worst_mismatches, c.offtarget_worst_junction_mm,
        c.offtarget_clean_left_bp, c.offtarget_clean_right_bp,
        c.offtarget_nearmiss_count, f"{c.best_offtarget_pident:.2f}",
        f"{c.best_offtarget_bitscore:.1f}", f"{c.best_offtarget_evalue:.2e}",
        f"{c.isoform_coverage:.4f}", c.isoform_hits, c.isoform_total,
        ";".join(c.isoform_hit_ids),
        f"{c.composite_score:.6f}", c.reduced_confidence,
        f"{c.tm_balance:.2f}", f"{c.self_dimer_dg:.1f}", f"{c.hetero_dimer_dg:.1f}",
    ]
    comps = c.score_components or {}
    base.extend(
        f"{comps[k]:.4f}" if k in comps else "" for k in _SCORE_COMPONENTS
    )
    base.extend([
        c.probe_role, c.parent_probe_id, ";".join(c.resolves_to),
        ",".join(str(p) for p in c.degenerate_positions),
        ",".join(c.iupac_codes), c.n_variants,
        _num(c.tm_lhs_min), _num(c.tm_lhs_max),
        _num(c.tm_rhs_min), _num(c.tm_rhs_max),
        _num(c.gc_lhs_min, "{:.4f}"), _num(c.gc_lhs_max, "{:.4f}"),
        _num(c.gc_rhs_min, "{:.4f}"), _num(c.gc_rhs_max, "{:.4f}"),
        _num(c.hairpin_dg, "{:.1f}"), _num(c.hairpin_dg_worst, "{:.1f}"),
        _num(c.hetero_dimer_dg_worst, "{:.1f}"), _num(c.tm_balance_max),
        ";".join(c.covered_record_ids),
    ])
    base.extend(
        f"{comps[k]:.4f}" if k in comps else ""
        for k in _SCORE_COMPONENTS_APPENDED
    )
    if not strain_ids:
        return base
    # Pull the strain dicts keyed by strain_id so column order matches the
    # header order regardless of how the pipeline ordered the list.
    by_sid = {d.get("strain_id", ""): d for d in (c.strain_binding or [])}
    for sid in strain_ids:
        d = by_sid.get(sid, {})
        positions = d.get("mismatch_positions") or []
        base.extend([
            d.get("total_mismatches", ""),
            d.get("lhs_mismatches", ""),
            d.get("rhs_mismatches", ""),
            d.get("junction_mismatches", ""),
            ",".join(str(p) for p in positions),
            f"{d.get('estimated_dTm', 0.0):.1f}" if d else "",
            d.get("likely_binds", ""),
        ])
    base.extend([
        c.min_strain_mismatches, c.max_strain_mismatches,
        f"{c.worst_strain_dTm:.1f}",
    ])
    return base


def write_results_table(
    candidates: List[ProbeCandidate],
    output_path: str,
    write_xlsx: bool = False,
    overwrite: bool = False,
    screen_engine: str = SCREEN_ENGINE_DEFAULT,
) -> str:
    if os.path.exists(output_path) and not overwrite:
        raise WriterError(f"Output exists and overwrite=False: {output_path}")

    # Detect a multi-strain run from the candidates themselves: if any
    # candidate carries a non-empty strain_binding list, expand the
    # header with one column-group per strain.
    strain_ids = []
    for c in candidates:
        if c.strain_binding:
            strain_ids = [d.get("strain_id", "") for d in c.strain_binding]
            break
    headers = BASE_HEADERS + _strain_headers(strain_ids)

    dir_name = os.path.dirname(output_path) or '.'
    os.makedirs(dir_name, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=dir_name, suffix='.tmp')
    try:
        with os.fdopen(fd, 'w', newline='') as fh:
            writer = csv.writer(fh)
            writer.writerow(headers)
            for c in candidates:
                writer.writerow(_row(c, strain_ids, screen_engine))
        os.replace(tmp_path, output_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise

    if write_xlsx:
        try:
            import openpyxl
            # write_only streams rows straight to the archive instead of holding
            # the whole cell grid in memory. No wb.active in this
            # mode — create the sheet explicitly.
            wb = openpyxl.Workbook(write_only=True)
            ws = wb.create_sheet('Probes')
            ws.append(headers)
            for c in candidates:
                ws.append(_row(c, strain_ids, screen_engine))
            xlsx_path = output_path.rsplit('.', 1)[0] + '.xlsx'
            wb.save(xlsx_path)
        except ImportError:
            pass  # openpyxl optional

    return output_path
