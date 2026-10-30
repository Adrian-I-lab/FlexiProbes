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

"""Discovery-mode review dialog.

Opens for runs in ``ready_for_review`` status. For each gene the user sees:

  • a horizontal map of the mRNA with selected probes as coloured bars
    (probes closer than the minimum spacing turn red) and variant positions
    marked as ticks,
  • below that map, when the run has more than one sequence, one lane per
    sequence — every isoform of the gene (isoform-aware runs) or every input
    record (multi-strain / homology runs) — with each selected probe drawn
    solid on the lanes it binds and hatched-hollow on the lanes it does not,
    so a probe missing from an isoform is visible at a glance,
  • a filterable, sortable table of every ranked candidate with checkboxes,
  • per-sequence mismatch / ΔTm columns when the run was multi-sequence
    (homology mode),
  • a "Pick best N (spaced)" shortcut,
  • a "Maximise isoform coverage (all genes)" action that picks, for every
    gene, up to max_per_gene spaced probes covering as many isoforms as
    possible (gene-aware runs; needs the input cDNA FASTA),
  • a counter showing "Selected: N / max_per_gene".

"Confirm & export" writes the order CSVs restricted to the user's
selections and flips the run status to ``complete``.
"""
from __future__ import annotations

import os
import threading
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from PySide6.QtCore import QEventLoop, Qt, QThread, Signal
from PySide6.QtGui import QBrush, QColor, QFont, QFontMetrics, QPainter, QPen
from PySide6.QtWidgets import (
    QAbstractItemView, QCheckBox, QComboBox, QDialog, QDialogButtonBox,
    QDoubleSpinBox,
    QFrame, QGridLayout, QHBoxLayout, QHeaderView, QLabel, QLineEdit,
    QMessageBox, QProgressDialog, QPushButton, QScrollArea, QSizePolicy,
    QSpinBox, QTableWidget, QTableWidgetItem, QToolTip, QVBoxLayout, QWidget,
)

from flexprobe.auto_select import auto_pick_per_gene, ladder_from
from flexprobe.constants import (
    GAP_TIER_GREEN, GAP_TIER_NAMES, GAP_TIER_ORANGE, GAP_TIER_RED,
    MIN_PROBE_GAP_BP, MIN_PROBE_GAP_FLOOR_BP, PREFERRED_PROBE_GAP_BP,
    PROBE_GAP_RELAX_MIN_PROBES, PROBE_ROLE_CONSENSUS, PROBE_ROLE_VARIANT,
    RELAXED_PROBE_GAP_BP,
)
from flexprobe.isoform_select import (
    build_gene_index, coverage_pick, gene_index_from_manifest,
    results_transcripts,
)
from flexprobe.select_export import (
    ResultsLoadCancelled, export_selected_probes, load_gene_lengths,
    load_modality_info, load_results_rows, load_selections,
    load_variant_positions, save_selections,
)

from .reveal import reveal_in_finder  # re-exported for existing importers
from .run_history import RunRecord
from .run_manager import load_run_config
from .theme import muted_hex

__all__ = ["ReviewDialog", "reveal_in_finder"]


@dataclass
class RunData:
    """Everything the review dialog needs off disk for one run.

    A plain value with no Qt in it, so the whole load — including the grouping
    pass over every row, which on a whole-transcriptome run is millions of
    iterations — can happen on a worker thread and arrive as one object.
    """
    rows_by_gene: Dict[str, List[Dict[str, str]]] = field(default_factory=dict)
    gene_lengths: Dict[str, int] = field(default_factory=dict)
    strain_ids: List[str] = field(default_factory=list)
    gene_isoforms: Dict[str, Any] = field(default_factory=dict)
    isoform_universe: Dict[str, Tuple[str, List[str]]] = field(default_factory=dict)
    has_iso_col: bool = False
    has_strain_cols: bool = False
    variants: List[int] = field(default_factory=list)
    selections: Dict[str, List[str]] = field(default_factory=dict)


def read_run_data(output_dir: str, *, progress=None,
                  should_cancel=None) -> RunData:
    """Read and shape one run's review data. No Qt, no widgets, no globals.

    Raises ``ResultsLoadCancelled`` when ``should_cancel()`` goes true partway
    through the CSV.
    """
    data = RunData()
    try:
        rows = load_results_rows(output_dir, progress=progress,
                                 should_cancel=should_cancel)
    except FileNotFoundError:
        rows = []
    for r in rows:
        data.rows_by_gene.setdefault(r.get("gene_id", ""), []).append(r)
    data.gene_lengths = load_gene_lengths(output_dir)
    for gid, rs in data.rows_by_gene.items():
        if data.gene_lengths.get(gid):
            continue
        try:
            data.gene_lengths[gid] = max(int(r.get("end", 0) or 0) for r in rs)
        except ValueError:
            data.gene_lengths[gid] = 0
    info = load_modality_info(output_dir)
    data.strain_ids = [str(s) for s in (info.get("strain_ids") or [])]
    # Gene -> isoform-ids map (manifest) + whether results carry the
    # isoform_hit_ids column. Together they let the coverage pick run without
    # the input FASTA.
    gi = info.get("gene_isoforms")
    data.gene_isoforms = gi if isinstance(gi, dict) else {}
    data.has_iso_col = any(
        (r.get("isoform_hit_ids") or "").strip()
        for rs in data.rows_by_gene.values() for r in rs
    )
    # Transcript -> its gene's full isoform list, inverted once here so the
    # lane block never has to walk the manifest map per repaint.
    for gene, transcripts in data.gene_isoforms.items():
        if not isinstance(transcripts, (list, tuple)):
            continue
        members = [str(t) for t in transcripts]
        for t in members:
            data.isoform_universe[t] = (str(gene), members)
    # Do the per-record columns actually exist? A homology run that predates
    # them lists strain_ids in the manifest but carries no mismatches__<id>
    # columns, and lanes built from absent columns would be fiction.
    first = next((r for rs in data.rows_by_gene.values() for r in rs), {})
    data.has_strain_cols = any(
        f"mismatches__{sid}" in first for sid in data.strain_ids
    )
    data.variants = load_variant_positions(output_dir)
    existing = load_selections(output_dir)
    for gid, ids in existing.items():
        valid_ids = {r.get("probe_id") for r in data.rows_by_gene.get(gid, [])}
        data.selections[gid] = [pid for pid in ids if pid in valid_ids]
    return data


class _RunDataWorker(QThread):
    """Reads one run's results.csv + manifest off the UI thread.

    The largest results.csv on this project is 1.17 GB; parsing it inline froze
    the dialog (and the whole app) for as long as the read took, with no
    progress and no way out. Shaped after ``_CoverageWorker`` below — the
    sibling worker in this same dialog: done/failed signals, parented to the
    dialog, deleteLater on finish — with the ``threading.Event`` cancel flag
    taken from ``run_manager._DesignWorker.request_cancel``.
    """
    progress = Signal(int)          # 0-100 through the CSV
    done = Signal(object)           # RunData
    failed = Signal(str)            # error message
    cancelled = Signal()

    def __init__(self, output_dir: str, parent=None):
        super().__init__(parent)
        self._output_dir = output_dir
        self._cancel = threading.Event()

    def request_cancel(self) -> None:
        self._cancel.set()

    def run(self) -> None:
        try:
            data = read_run_data(
                self._output_dir,
                progress=self.progress.emit,
                should_cancel=self._cancel.is_set,
            )
        except ResultsLoadCancelled:
            self.cancelled.emit()
            return
        except Exception as exc:    # surfaced in a dialog on the UI thread
            self.failed.emit(str(exc))
            return
        self.done.emit(data)


class _CoverageWorker(QThread):
    """Runs the per-gene coverage pick off the UI thread.

    Preferred path (current runs): membership comes from the persisted
    ``isoform_hit_ids`` column + the manifest's gene->isoform map, so the index
    is built in memory and no FASTA is read. Fallback (older runs): build the
    index from the input cDNA FASTA (hundreds of MB), cached across clicks."""
    done = Signal(object, object)   # (CoveragePick, GeneIndex)
    failed = Signal(str)            # error message

    def __init__(self, *, rows, n_per_gene, min_gap, preferred_gap,
                 relaxed_gap=RELAXED_PROBE_GAP_BP,
                 floor_gap=MIN_PROBE_GAP_FLOOR_BP,
                 gene_isoforms=None, fasta=None, cached_index=None, parent=None):
        super().__init__(parent)
        self._rows = rows
        self._n = n_per_gene
        self._min_gap = min_gap
        self._preferred_gap = preferred_gap
        self._relaxed_gap = relaxed_gap
        self._floor_gap = floor_gap
        self._gene_isoforms = gene_isoforms
        self._fasta = fasta
        self._index = cached_index

    def run(self) -> None:
        try:
            index = self._index
            if index is None:
                if self._gene_isoforms is not None:
                    index = gene_index_from_manifest(self._gene_isoforms)
                else:
                    index = build_gene_index(
                        self._fasta, results_transcripts(self._rows))
            pick = coverage_pick(
                self._rows, index, n_per_gene=self._n, min_gap=self._min_gap,
                preferred_gap=self._preferred_gap,
                relaxed_gap=self._relaxed_gap, floor_gap=self._floor_gap)
            self.done.emit(pick, index)
        except Exception as exc:  # surfaced in a dialog on the UI thread
            self.failed.emit(str(exc))


# Probe-bar colours by closest-neighbour spacing. Keyed by the shared
# GAP_TIER_* vocabulary so the map, the legend, the exports and the CLI can
# never disagree about what "orange" means.
_GAP_TIER_COLORS = {GAP_TIER_RED: "#d62728", GAP_TIER_ORANGE: "#e8820c",
                    GAP_TIER_GREEN: "#2ca02c"}

# The tier encoded WITHOUT hue. Red/green fill alone excludes roughly one man in
# twelve, every screen reader, and any greyscale print of the figure — so each
# bar is stamped with its tier letter and, for the two non-ideal tiers, overlaid
# with a fill pattern that survives desaturation.
_GAP_TIER_GLYPHS = {GAP_TIER_RED: "R", GAP_TIER_ORANGE: "O",
                    GAP_TIER_GREEN: "G"}
_GAP_TIER_PATTERNS = {
    GAP_TIER_RED: Qt.BDiagPattern,      # dense backward hatch: crowded
    GAP_TIER_ORANGE: Qt.FDiagPattern,   # lighter forward hatch: tight
    GAP_TIER_GREEN: None,               # ideal: plain fill, nothing to flag
}
# The tier in WORDS, for tooltips and accessible descriptions. GAP_TIER_NAMES
# gives the colour word; this gives what the colour MEANS, which is the part a
# reader who cannot see the colour actually needs.
_GAP_TIER_WORDS = {
    GAP_TIER_RED: "crowded — closer than the relaxed rung",
    GAP_TIER_ORANGE: "tight — between the relaxed rung and the design target",
    GAP_TIER_GREEN: "well spaced — at or above the design target",
}


def tier_phrase(tier: int) -> str:
    """e.g. ``"R (red): crowded — closer than the relaxed rung"``."""
    return (f"{_GAP_TIER_GLYPHS.get(tier, '?')} "
            f"({GAP_TIER_NAMES.get(tier, '?')}): "
            f"{_GAP_TIER_WORDS.get(tier, '')}")

# ---------------------------------------------------------------------------
# The run's OWN spacing configuration
#
# MIN_PROBE_GAP_BP and friends are only the application defaults. A run
# designed with `--min-gap 300` records 300 in its run_config.json sidecar, and
# every number this dialog shows, filters on, colours by and re-picks with has
# to come from THAT, not from the module constants — otherwise a 300 bp run is
# displayed and re-picked as if it were a 150 bp one, and the labels lie about
# it. Reachable through the Export-HPC → Import → Review workflow.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class RunSpacing:
    """The three spacing rungs this run was actually configured with.

    ``from_run_config`` is False when no ``run_config.json`` could be read — a
    legacy run, or a hand-made output folder. In that case these ARE the module
    defaults, and the dialog says so out loud rather than presenting them as if
    they were the run's own settings.
    """
    min_gap: int = MIN_PROBE_GAP_BP
    relaxed_gap: int = RELAXED_PROBE_GAP_BP
    floor_gap: int = MIN_PROBE_GAP_FLOOR_BP
    from_run_config: bool = False

    def tier(self, gap: Optional[int]) -> int:
        """This run's tier for a closest-neighbour gap (None = no neighbour).

        Same shape as ``flexprobe.constants.gap_tier``, but against the run's
        own rungs instead of the defaults.
        """
        if gap is None or gap >= self.min_gap:
            return GAP_TIER_GREEN
        if gap >= self.relaxed_gap:
            return GAP_TIER_ORANGE
        return GAP_TIER_RED


def run_spacing_for(output_dir: str) -> RunSpacing:
    """Read a run's configured spacing, falling back to the module defaults."""
    cfg = load_run_config(output_dir) if output_dir else None
    if cfg is None:
        return RunSpacing()
    try:
        return RunSpacing(
            min_gap=int(cfg.min_probe_gap_bp),
            relaxed_gap=int(cfg.relaxed_probe_gap_bp),
            floor_gap=int(cfg.floor_probe_gap_bp),
            from_run_config=True,
        )
    except (AttributeError, TypeError, ValueError):
        # A sidecar from a version that predates these fields. Defaults, and
        # flagged as defaults.
        return RunSpacing()


def _spacing_filters(sp: RunSpacing):
    """The spacing filter's three settings, loosest first, in THIS run's units.

    Each is the ladder truncated: refusing a tier means its rung is never
    walked, so those probes are never candidates — as opposed to being picked
    and then warned about.
    """
    return (
        ("Allow all (green / orange / red)", sp.floor_gap, GAP_TIER_RED),
        (f"No red (>= {sp.relaxed_gap} bp apart)", sp.relaxed_gap,
         GAP_TIER_ORANGE),
        (f"Green only (>= {sp.min_gap} bp apart)", sp.min_gap, GAP_TIER_GREEN),
    )


def _interval_gap(a, b) -> int:
    """Gap in bp between two ``[start, end)`` intervals (negative = overlap)."""
    (s1, e1), (s2, e2) = a, b
    if e1 <= s2:
        return s2 - e1
    if e2 <= s1:
        return s1 - e2
    return -1


def site_key(row: Mapping[str, Any]) -> str:
    """The id of the SITE a row occupies, which is not the id of the row.

    A degenerate consensus parent and every one of its variant children are
    built on one and the same ``[start, end)``. See the child construction in
    ``flexprobe.pipeline``, which passes the parent's own ``start``/``end``
    straight through. They are one probe at one place: ``auto_select`` walks the
    parent, counts the family once against the gene's quota ("one family = one
    pick") and emits the children because only a child is orderable.

    So anything that measures SPACING has to work in sites. Measuring it per
    emitted id compares a family's siblings against each other, which is a
    category error: the answer is always -1 bp, and a correctly spaced panel
    reports itself as entirely crowded.
    """
    role = str(row.get("probe_role") or "").strip()
    parent = str(row.get("parent_probe_id") or "").strip()
    if role == PROBE_ROLE_VARIANT and parent:
        return parent
    return str(row.get("probe_id", ""))


def _probe_gap_tiers(probes, min_gap: int, preferred_gap: int):
    """Per-probe spacing tier from its closest neighbour on the track:
    0 = closer than min_gap (red), 1 = min..preferred (orange), 2 = >= preferred
    or no neighbour (green).

    Two probes of the same SITE are not neighbours. See :func:`site_key`. A
    probe dict may carry a ``"site"``; without one it is its own site, which is
    what every non-degenerate run gives and is the behaviour this had before.
    """
    ivs = [(int(p["start"]), int(p["end"])) for p in probes]
    # The index fallback keeps two id-less rows from collapsing into one site.
    sites = [str(p.get("site") or p.get("probe_id") or f"#{i}")
             for i, p in enumerate(probes)]
    tiers = []
    for i, (s1, e1) in enumerate(ivs):
        # None until the first neighbour is seen; -1 encodes an overlap, so the
        # annotation has to be explicit or the initialiser types it as None.
        best: Optional[int] = None
        for j, (s2, e2) in enumerate(ivs):
            if i == j or sites[i] == sites[j]:
                continue
            g = _interval_gap((s1, e1), (s2, e2))
            best = g if best is None else min(best, g)
        if best is None or best >= preferred_gap:
            tiers.append(2)
        elif best >= min_gap:
            tiers.append(1)
        else:
            tiers.append(0)
    return tiers


# ---------------------------------------------------------------------------
# Multi-lane model
# ---------------------------------------------------------------------------

# How a probe relates to one lane's sequence. Deliberately NOT carried by hue
# alone — hue already encodes the spacing tier — so each state also has its own
# fill pattern and outline and survives a greyscale screenshot in a paper.
MARK_ABSENT = 0      # does not bind this sequence: hollow, hatched, dashed edge
MARK_MISMATCH = 1    # binds, but imperfectly: solid + diagonal hatch overlay
MARK_BINDS = 2       # binds cleanly: solid, exactly like the summary lane

# A lane block of 40 isoforms would otherwise produce a ~1400 px widget; the
# dialog scrolls instead, and past this many lanes the tail is summarised in
# the caption rather than drawn (nobody reads lane 200).
MAX_LANES = 60

# Vertical geometry of the lane block, in px. LANE_PITCH must leave room for
# the small label drawn ABOVE each track: labels sit above rather than in a
# left-hand gutter precisely so every lane track spans the same x range as the
# summary track and the shared x-axis really is shared.
LANE_LABEL_H = 12
LANE_TRACK_H = 14
LANE_PITCH = LANE_LABEL_H + LANE_TRACK_H + 6
LANE_BLOCK_TOP_PAD = 6

# The summary block keeps exactly the height the single-track widget always had,
# so a single-sequence run is pixel-identical to the pre-lane dialog.
SUMMARY_BLOCK_H = 110

# The most vertical space the map is ever allowed to take in the dialog. Past
# this the lane block scrolls: 40 isoforms must not push the table off-screen.
MAX_MAP_VIEWPORT_H = 340


@dataclass(frozen=True)
class LaneMark:
    """One selected probe as it appears on ONE lane.

    ``start``/``end`` are in the lane block's shared coordinate system (see
    :class:`LaneSet.coord_caption`) — not necessarily an offset into this lane's
    own sequence, which older runs simply do not record.
    """
    probe_id: str
    rank: int
    start: int
    end: int
    state: int          # MARK_ABSENT | MARK_MISMATCH | MARK_BINDS
    detail: str = ""    # tooltip line: mismatches / why it is absent


@dataclass(frozen=True)
class SequenceLane:
    """One sequence — an isoform or an input record — and its probe marks."""
    label: str
    marks: Tuple[LaneMark, ...] = ()
    length: int = 0     # this sequence's OWN length in bp, 0 when unknown
    note: str = ""      # right-aligned annotation, e.g. "2/3 probes bind"


@dataclass(frozen=True)
class LaneSet:
    """The lane block for one gene: the lanes plus the caption that says which
    coordinate system they are drawn in.

    The caption is not decoration. Alignment coordinates and per-record
    coordinates are different things, and the column-based isoform path records
    *which* isoforms a probe binds without recording *where* on each — so the
    caption is the only place the user learns which of those they are looking
    at. An empty ``lanes`` means "nothing extra to show": the widget then
    renders the single summary track and nothing else.
    """
    lanes: Tuple[SequenceLane, ...] = ()
    coord_caption: str = ""
    kind: str = ""      # 'isoform' | 'record' | ''


def _lane_state_brushes(tier_color: QColor, state: int):
    """(brush, pen) for one lane mark. See MARK_* for why pattern, not hue."""
    if state == MARK_ABSENT:
        # Hollow + backward-diagonal hatch in neutral grey: unmistakably "not
        # here" even when the page is printed in black and white.
        return (QBrush(QColor("#9e9e9e"), Qt.BDiagPattern),
                QPen(QColor("#9e9e9e"), 1, Qt.DashLine))
    if state == MARK_MISMATCH:
        return QBrush(tier_color), QPen(tier_color.darker(160), 1, Qt.DotLine)
    return QBrush(tier_color), QPen(tier_color.darker(140), 1)


# ---------------------------------------------------------------------------
# Visualization widget
# ---------------------------------------------------------------------------

class ProbeMapWidget(QWidget):
    """Horizontal mRNA track with selected-probe rectangles + variant ticks,
    optionally followed by one lane per isoform / input record.

    With no lanes the widget is exactly the single-track map it has always been
    — that path is untouched, because single-sequence runs are the common case
    and must not shift by a pixel.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._gene_id = ""
        self._gene_length = 0
        self._probes: List[Dict] = []
        self._gap_tiers: List[int] = []
        self._tier_by_pid: Dict[str, int] = {}
        self._variants: List[int] = []
        self._lane_set = LaneSet()
        self._overflow = 0
        # Defaults until the dialog hands us the run's own rungs; a bare widget
        # (tests, a preview) still colours sensibly.
        self._spacing = RunSpacing()
        self.setMinimumHeight(SUMMARY_BLOCK_H)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.setMouseTracking(True)
        self.setAccessibleName("Probe map")

    def _describe(self) -> str:
        """What a screen reader is told about the map, in words only.

        The map is a QPainter drawing: it exposes no child widgets, so without
        this it is an unlabelled rectangle to assistive technology."""
        if not self._probes:
            return (f"Probe map for {self._gene_id or 'no gene'}: "
                    "no probes selected.")
        by_tier: Dict[int, int] = {}
        for t in self._gap_tiers:
            by_tier[t] = by_tier.get(t, 0) + 1
        parts = ", ".join(
            f"{by_tier[t]} {GAP_TIER_NAMES[t]}"
            for t in (GAP_TIER_GREEN, GAP_TIER_ORANGE, GAP_TIER_RED)
            if by_tier.get(t)
        )
        lanes = (f", {len(self._lane_set.lanes)} sequence lane(s)"
                 if self._lane_set.lanes else "")
        return (f"Probe map for {self._gene_id}, {self._gene_length} bp: "
                f"{len(self._probes)} selected probe(s) by spacing tier — "
                f"{parts}{lanes}.")

    def set_spacing(self, spacing: RunSpacing) -> None:
        """Adopt the run's configured rungs. Call BEFORE set_data — the tiers
        are computed there and are only correct against the right thresholds."""
        self._spacing = spacing

    def set_data(self, gene_id, gene_length, probes, variants=None,
                 lane_set: Optional[LaneSet] = None) -> None:
        self._gene_id = gene_id
        self._gene_length = max(0, int(gene_length))
        self._probes = list(probes)
        # Tiers follow the SELECTION policy, not two arbitrary numbers:
        # red = under the relaxation floor (the selector would never emit this,
        # so it can only be a manual pick), orange = between the floor and the
        # design target (what the relaxation pass produces — legal but flagged),
        # green = at or above the target.
        self._gap_tiers = (
            _probe_gap_tiers(self._probes, self._spacing.relaxed_gap,
                             self._spacing.min_gap)
            if self._probes else []
        )
        self._tier_by_pid = {
            str(pr.get("probe_id", "")): (
                self._gap_tiers[i] if i < len(self._gap_tiers) else 2)
            for i, pr in enumerate(self._probes)
        }
        self._variants = list(variants or [])
        lane_set = lane_set or LaneSet()
        # Truncate here rather than in the caller so no caller can hand us a
        # 500-isoform gene and blow the widget up; the caller reports the
        # overflow in its caption.
        self._overflow = max(0, len(lane_set.lanes) - MAX_LANES)
        if self._overflow:
            lane_set = LaneSet(lanes=tuple(lane_set.lanes[:MAX_LANES]),
                               coord_caption=lane_set.coord_caption,
                               kind=lane_set.kind)
        self._lane_set = lane_set
        # Fixed height, so the enclosing scroll area knows how much there is to
        # scroll; the cap on what is *visible* is the scroll area's job.
        self.setFixedHeight(self.content_height())
        self.setAccessibleDescription(self._describe())
        self.update()

    @property
    def lanes(self) -> Tuple[SequenceLane, ...]:
        return self._lane_set.lanes

    def content_height(self) -> int:
        """Full height needed to draw everything — summary block plus lanes."""
        if not self._lane_set.lanes:
            return SUMMARY_BLOCK_H
        return (SUMMARY_BLOCK_H + LANE_BLOCK_TOP_PAD
                + LANE_PITCH * len(self._lane_set.lanes))

    @property
    def has_overlap(self) -> bool:
        """Any chosen pair closer than the hard minimum (red bars)."""
        return any(t == 0 for t in self._gap_tiers)

    @property
    def has_tight(self) -> bool:
        """Any chosen pair in the acceptable-but-tight band (orange bars)."""
        return any(t == 1 for t in self._gap_tiers)

    # -- geometry ------------------------------------------------------
    # Shared by paintEvent and the hover handler so a tooltip can never point
    # at a lane the painter drew somewhere else.

    def _track_span(self) -> Tuple[int, int]:
        """(x, width) of the track area — identical for the summary and every
        lane, which is what makes the x-axis genuinely common."""
        margin_x = 28
        return margin_x, max(1, self.width() - 2 * margin_x)

    def _summary_block_h(self) -> int:
        # Without lanes the summary owns the whole widget, exactly as before.
        return SUMMARY_BLOCK_H if self._lane_set.lanes else self.height()

    def _lane_rows(self) -> List[Tuple[int, int, int]]:
        """(lane_index, label_baseline_y, track_top_y) for each drawn lane."""
        top = self._summary_block_h() + LANE_BLOCK_TOP_PAD
        rows = []
        for i in range(len(self._lane_set.lanes)):
            y = top + i * LANE_PITCH
            rows.append((i, y + LANE_LABEL_H - 2, y + LANE_LABEL_H))
        return rows

    def paintEvent(self, _event) -> None:
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        w = self.width()
        h = self._summary_block_h()
        margin_x = 28
        track_y = h // 2 + 6
        track_h = 22
        track_x = margin_x
        track_w = max(1, w - 2 * margin_x)

        def x_for(pos: int) -> int:
            if not self._gene_length:
                return track_x
            return track_x + int(round(
                track_w * max(0, min(pos, self._gene_length)) / self._gene_length
            ))

        # Track background.
        p.setPen(QPen(QColor("#bbbbbb"), 1))
        p.setBrush(QBrush(QColor("#ececec")))
        p.drawRoundedRect(track_x, track_y - track_h // 2, track_w, track_h, 4, 4)

        # Gene label.
        p.setPen(QColor("#444"))
        p.drawText(track_x, track_y - track_h // 2 - 18,
                   f"{self._gene_id}  (length: {self._gene_length} bp"
                   + (f", {len(self._variants)} variant positions)" if self._variants else ")"))

        # Endpoint labels.
        p.setPen(QColor("#666"))
        p.drawText(track_x, track_y + track_h // 2 + 16, "0")
        if self._gene_length:
            lbl = str(self._gene_length)
            fm = QFontMetrics(p.font())
            p.drawText(track_x + track_w - fm.horizontalAdvance(lbl),
                       track_y + track_h // 2 + 16, lbl)

        # Variant ticks (small marks below the track) — these show the user
        # which positions differ between strains, so they can see when a
        # probe spans a variant.
        if self._gene_length and self._variants:
            p.setPen(QPen(QColor("#c0392b"), 1))
            ty = track_y + track_h // 2 + 1
            for vp in self._variants:
                vx = x_for(vp)
                p.drawLine(vx, ty, vx, ty + 5)

        if self._gene_length and self._probes:
            for idx, probe in enumerate(self._probes):
                x1, x2 = x_for(int(probe["start"])), x_for(int(probe["end"]))
                band_w = max(2, x2 - x1)
                band_y = track_y - track_h // 2 - 6
                band_h = track_h + 12
                # Colour by spacing tier: green (ideal) / orange (tight) / red
                # (too close), so the user sees spacing quality at a glance.
                tier = self._gap_tiers[idx] if idx < len(self._gap_tiers) else 2
                fill = QColor(_GAP_TIER_COLORS[tier])
                p.setPen(QPen(fill.darker(140), 1))
                p.setBrush(QBrush(fill))
                p.drawRoundedRect(x1, band_y, band_w, band_h, 3, 3)
                # Non-colour cue #1: a hatch over the two tiers that need
                # attention, so a desaturated screenshot still separates them.
                pattern = _GAP_TIER_PATTERNS[tier]
                if pattern is not None:
                    p.setBrush(QBrush(QColor(255, 255, 255, 130), pattern))
                    p.setPen(Qt.NoPen)
                    p.drawRoundedRect(x1, band_y, band_w, band_h, 3, 3)
                # Non-colour cue #2: the tier LETTER beside the rank. The
                # legend, the tooltip and this glyph all use the same G/O/R.
                label = f"#{probe.get('rank', idx + 1)} {_GAP_TIER_GLYPHS[tier]}"
                fm = QFontMetrics(p.font())
                if band_w < fm.horizontalAdvance(label) + 6:
                    # Too narrow for both — the tier is the part that cannot be
                    # recovered from anywhere else on the bar, so it wins.
                    label = _GAP_TIER_GLYPHS[tier]
                if band_w >= fm.horizontalAdvance(label) + 6:
                    p.setPen(QColor("white"))
                    p.drawText(x1 + 3, band_y + band_h - 5, label)

        if self._lane_set.lanes:
            self._paint_lanes(p, x_for, track_x, track_w)
        p.end()

    def _paint_lanes(self, p: QPainter, x_for, track_x: int, track_w: int) -> None:
        """Draw the per-sequence lane block under the summary track."""
        small = QFont(p.font())
        small.setPointSizeF(max(6.5, small.pointSizeF() - 2.0))
        p.setFont(small)
        fm = QFontMetrics(small)

        for i, label_y, top_y in self._lane_rows():
            lane = self._lane_set.lanes[i]
            # Label above the track (never a left gutter — a gutter would shift
            # the lane's x origin away from the summary track's).
            caption = lane.label
            if lane.length:
                caption += f"  ({lane.length:,} bp)"
            p.setPen(QColor("#555"))
            p.drawText(track_x, label_y,
                       fm.elidedText(caption, Qt.ElideMiddle, track_w - 120))
            if lane.note:
                p.setPen(QColor("#888"))
                p.drawText(track_x + track_w - fm.horizontalAdvance(lane.note),
                           label_y, lane.note)

            p.setPen(QPen(QColor("#cccccc"), 1))
            p.setBrush(QBrush(QColor("#f4f4f4")))
            p.drawRoundedRect(track_x, top_y, track_w, LANE_TRACK_H, 3, 3)

            for mark in lane.marks:
                x1, x2 = x_for(mark.start), x_for(mark.end)
                band_w = max(2, x2 - x1)
                tier = self._tier_for(mark.probe_id)
                brush, pen = _lane_state_brushes(QColor(_GAP_TIER_COLORS[tier]),
                                                 mark.state)
                p.setBrush(brush)
                p.setPen(pen)
                p.drawRoundedRect(x1, top_y + 1, band_w, LANE_TRACK_H - 2, 2, 2)
                if mark.state == MARK_MISMATCH:
                    # Solid fill, then a light hatch on top: still obviously
                    # "filled" (it binds) but visibly textured (imperfectly).
                    p.setBrush(QBrush(QColor(255, 255, 255, 150), Qt.FDiagPattern))
                    p.setPen(Qt.NoPen)
                    p.drawRoundedRect(x1, top_y + 1, band_w,
                                      LANE_TRACK_H - 2, 2, 2)

        if self._overflow:
            last_top = self._lane_rows()[-1][2]
            p.setPen(QColor("#888"))
            p.drawText(track_x, last_top + LANE_TRACK_H + 11,
                       f"… {self._overflow:,} further sequence(s) not drawn")

    def _tier_for(self, probe_id: str) -> int:
        # Precomputed in set_data: paintEvent runs on every selection change and
        # every resize, so it must not scan the probe list per lane mark.
        return self._tier_by_pid.get(probe_id, 2)

    # -- hover ---------------------------------------------------------

    def _pos_at(self, x: float) -> int:
        track_x, track_w = self._track_span()
        frac = (x - track_x) / track_w if track_w else 0.0
        return int(round(self._gene_length * max(0.0, min(1.0, frac))))

    def mouseMoveEvent(self, event) -> None:
        if not self._gene_length:
            return
        pos = self._pos_at(event.position().x())
        y = event.position().y()

        # Over the lane block? Report binding + mismatch state for THAT lane.
        for i, _label_y, top_y in self._lane_rows():
            if top_y <= y <= top_y + LANE_TRACK_H:
                lane = self._lane_set.lanes[i]
                hits = [m for m in lane.marks if m.start <= pos < m.end]
                if hits:
                    QToolTip.showText(
                        event.globalPosition().toPoint(),
                        f"<b>{lane.label}</b><br>" + "<br>".join(
                            f"#{m.rank} {m.probe_id} ({m.start}-{m.end}) — "
                            f"{m.detail}" for m in hits
                        ),
                        self,
                    )
                return

        if not self._probes:
            return
        hits_p = [pr for pr in self._probes if pr["start"] <= pos < pr["end"]]
        if hits_p:
            QToolTip.showText(
                event.globalPosition().toPoint(),
                "<br>".join(
                    f"#{p.get('rank', '?')} {p['probe_id']} ({p['start']}-{p['end']})"
                    # The tier IN WORDS. Without this the spacing verdict was
                    # only ever available as a fill colour, so it did not exist
                    # for a colour-blind reader or a screen reader at all.
                    + f" — spacing {tier_phrase(self._tier_for(p['probe_id']))}"
                    + self._binding_summary(p["probe_id"])
                    for p in hits_p
                ),
                self,
            )

    def _binding_summary(self, probe_id: str) -> str:
        """"binds 2 / 3 sequences" for the summary-lane tooltip. Empty when
        there are no lanes, so single-sequence runs keep today's exact text."""
        lanes = self._lane_set.lanes
        if not lanes:
            return ""
        bound = sum(
            1 for lane in lanes
            for m in lane.marks
            if m.probe_id == probe_id and m.state != MARK_ABSENT
        )
        return f" — binds {bound} / {len(lanes)} sequences"


# ---------------------------------------------------------------------------
# Review dialog
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class _Col:
    """One review-table column. Drives BOTH the header and the per-row cell, so
    column order/content live in a single registry — no hardcoded cell indices
    to drift out of sync when a column is added or moved."""
    label: str
    width: int
    kind: str           # 'select' | 'num' | 'gc' | 'tm' | 'spec' | 'text'
    source: str = ""    # results.csv column the cell reads (empty for 'select')


# Base (always-present) columns. Strain columns are appended dynamically.
# GC is shown three ways: the full 50-bp target and each 25-bp half (LHS/RHS),
# all as percentages (the underlying results.csv stores them as fractions).
_BASE_COLS = [
    _Col("Select", 56, "select"),
    _Col("Rank", 52, "num", "rank"),
    _Col("Start", 64, "num", "start"),
    _Col("End", 64, "num", "end"),
    _Col("GC% 50bp", 72, "gc", "gc_content_full"),
    _Col("GC% LHS", 68, "gc", "gc_lhs_half"),
    _Col("GC% RHS", 68, "gc", "gc_rhs_half"),
    _Col("Tm LHS", 70, "tm", "tm_target_lhs"),
    _Col("Tm RHS", 70, "tm", "tm_target_rhs"),
    _Col("Score", 80, "num", "composite_score"),
    _Col("Spec", 110, "spec", "specificity_status"),
    # The EVIDENCE behind that verdict, not just the verdict. Four separate
    # defects in the specificity path shipped a plausible number for an unknown
    # one, and every one of them was invisible in a results table that showed
    # only pass/fail. These are the columns that make the difference legible:
    #   Recount   — exact, or a lower bound because a subject could not be read
    #   Min qcov  — the least of the half BLAST actually aligned, worst case
    #   OT mm     — how far the CLOSEST off-target sits over the full half
    #   OT junc   — how many of those mismatches sit beside the nick, which on
    #               an OFF-target is protective: the ligase will not seal there
    #   Dup L/R   — the other side of that same coin, and the evidence behind a
    #               `failed_junction_duplex` verdict: matched bases either side
    #               of the nick against the worst non-target gene. -1 = never
    #               measured (only the k-mer engine records these; a
    #               BLAST-engine run leaves them -1), which is not 0.
    _Col("Recount", 92, "text", "offtarget_recount_basis"),
    _Col("Min qcov", 78, "num", "offtarget_min_qcov"),
    _Col("OT mm", 62, "num", "offtarget_worst_mismatches"),
    _Col("OT junc", 66, "num", "offtarget_worst_junction_mm"),
    _Col("Dup L", 60, "num", "offtarget_clean_left_bp"),
    _Col("Dup R", 60, "num", "offtarget_clean_right_bp"),
    # Blank on a normal run; "consensus" marks a degenerate parent, which is
    # never orderable, and "variant" its concrete per-record resolution.
    _Col("Role", 80, "text", "probe_role"),
    _Col("Probe ID", 200, "text", "probe_id"),
]

# The single checkbox column — derived from the registry so it can never drift.
SELECT_COL = next(i for i, c in enumerate(_BASE_COLS) if c.kind == "select")


class ReviewDialog(QDialog):
    """Modal review-and-select dialog for one discovery-mode run."""

    confirmed = Signal(str)   # run_id — emitted after a successful export

    def __init__(self, record: RunRecord, parent=None):
        super().__init__(parent)
        self._record = record
        # FIRST, before any label, colour, filter or re-pick is built: every one
        # of them reads these numbers, and reading the module defaults instead
        # is the bug this exists to close.
        self._spacing = run_spacing_for(record.output_dir)
        self.setWindowTitle(f"Review: {record.name}")
        self.setModal(True)
        self.resize(1120, 760)
        self.setMinimumSize(760, 500)
        self.setSizeGripEnabled(True)

        self._rows_by_gene: Dict[str, List[Dict[str, str]]] = {}
        self._gene_lengths: Dict[str, int] = {}
        self._variants: List[int] = []
        self._selections: Dict[str, List[str]] = {}
        self._max_per_gene = 3
        self._strain_ids: List[str] = []
        self._gene_index = None        # cached GeneIndex for coverage mode
        self._cov_worker: Optional[_CoverageWorker] = None
        # Set when the user cancels the load, or it fails: exec() then returns
        # Rejected without ever showing a half-built dialog.
        self._load_aborted = False
        self._load_error = ""
        # Overwritten by _adopt_run_data; needed as a default because an
        # aborted load never reaches it and the filter bar is built regardless.
        self._has_degenerate = False

        self._load_run_data()
        self._cols = list(_BASE_COLS) + self._strain_columns()

        outer = QVBoxLayout(self)
        outer.setContentsMargins(12, 12, 12, 12)
        outer.setSpacing(8)

        self._header_label = QLabel()
        self._header_label.setWordWrap(True)
        outer.addWidget(self._header_label)
        self._refresh_header()

        # Gene picker + per-gene limit + auto-pick.
        gene_row = QHBoxLayout()
        # Every label/field pair below is joined with setBuddy: that is what
        # gives the field its accessible name from its own visible label, and
        # what makes the label's Alt-mnemonic move focus to the field.
        gene_lbl = QLabel("&Gene:")
        gene_row.addWidget(gene_lbl)
        self.gene_combo = QComboBox()
        gene_lbl.setBuddy(self.gene_combo)
        self.gene_combo.setAccessibleName("Gene")
        self.gene_combo.setAccessibleDescription(
            "Which gene's candidate probes the map and table below show.")
        for gid in sorted(self._rows_by_gene):
            self.gene_combo.addItem(
                f"{gid}  ({len(self._rows_by_gene[gid])} candidates)", userData=gid,
            )
        gene_row.addWidget(self.gene_combo, 1)
        gene_row.addSpacing(16)
        max_lbl = QLabel("Ma&x per gene:")
        gene_row.addWidget(max_lbl)
        self.max_spin = QSpinBox(); self.max_spin.setRange(1, 50)
        self.max_spin.setValue(self._max_per_gene)
        max_lbl.setBuddy(self.max_spin)
        self.max_spin.setAccessibleName("Maximum probes per gene")
        self.max_spin.setAccessibleDescription(
            "The most probes that may be selected for any one gene.")
        gene_row.addWidget(self.max_spin)
        gene_row.addSpacing(16)
        spacing_lbl = QLabel("&Spacing:")
        gene_row.addWidget(spacing_lbl)
        self.spacing_filter = QComboBox()
        spacing_lbl.setBuddy(self.spacing_filter)
        self.spacing_filter.setAccessibleName("Spacing filter")
        for label, floor, min_tier in _spacing_filters(self._spacing):
            self.spacing_filter.addItem(label, userData=(floor, min_tier))
        self.spacing_filter.setToolTip(
            "How crowded a pair of chosen probes is allowed to get.\n\n"
            f"Auto-selection walks a ladder — {self._spacing.min_gap} bp apart, "
            f"then {self._spacing.relaxed_gap} bp, then "
            f"{self._spacing.floor_gap} bp — and stops at the first rung that "
            f"finds {PROBE_GAP_RELAX_MIN_PROBES} probes (or Max per gene, if "
            "smaller). This setting removes "
            "the lower rungs, so tighter probes are never candidates rather "
            "than being picked and then flagged.\n\n"
            "Tightening it costs probes on genes whose transcript is too short "
            "to hold separated ones — those genes come back with fewer probes, "
            "not with worse spacing.\n\n"
            f"{self._spacing_source_sentence()}"
        )
        self.spacing_filter.setAccessibleDescription(
            "How close two selected probes are allowed to be. "
            + self._spacing_source_sentence())
        gene_row.addWidget(self.spacing_filter)
        self.autopick_btn = QPushButton("&Pick best N (spaced)")
        self.autopick_btn.setToolTip(
            "Auto-select the top-N highest-scoring probes for this gene, spread "
            "as widely as the transcript allows and at least "
            f"{self._spacing.min_gap} bp apart. If fewer than "
            f"{PROBE_GAP_RELAX_MIN_PROBES} probes fit at that spacing the gap "
            f"drops to {self._spacing.relaxed_gap} bp, then to "
            f"{self._spacing.floor_gap} bp, and the tighter pairs are flagged "
            "orange or red. Use the Spacing box to forbid those rungs outright. "
            "You can still adjust the selection after.\n\n"
            f"{self._spacing_source_sentence()}"
        )
        gene_row.addWidget(self.autopick_btn)
        self.autopick_btn.setAccessibleName("Pick best N spaced probes")
        self.coverage_btn = QPushButton("&Maximise isoform coverage (all genes)")
        self.coverage_btn.setAccessibleName(
            "Maximise isoform coverage across all genes")
        self.coverage_btn.setToolTip(
            "Across every gene, auto-pick up to 'Max per gene' spaced probes that "
            "between them bind as many of the gene's isoforms as possible. Genes "
            "that can't be fully covered are reported afterwards."
        )
        gene_row.addWidget(self.coverage_btn)
        if self._coverage_source()[0] is None:
            self.coverage_btn.setEnabled(False)
            self.coverage_btn.setToolTip(
                "Isoform-coverage selection needs either the persisted "
                "isoform_hit_ids column or the run's input cDNA FASTA, neither "
                "of which is available for this run."
            )
        outer.addLayout(gene_row)

        # Visualization.
        viz_frame = QFrame(); viz_frame.setFrameShape(QFrame.StyledPanel)
        viz_l = QVBoxLayout(viz_frame); viz_l.setContentsMargins(8, 8, 8, 8)
        self.map_widget = ProbeMapWidget()
        self.map_widget.set_spacing(self._spacing)
        # The map lives in a height-capped scroll area so a gene with dozens of
        # isoforms scrolls instead of growing the dialog to 2000 px. Frameless
        # and sized to the content when there is only the summary track, so a
        # single-sequence run looks exactly as it did before lanes existed.
        self.map_scroll = QScrollArea()
        self.map_scroll.setWidget(self.map_widget)
        self.map_scroll.setWidgetResizable(True)
        self.map_scroll.setFrameShape(QFrame.NoFrame)
        self.map_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.map_scroll.setFixedHeight(SUMMARY_BLOCK_H)
        viz_l.addWidget(self.map_scroll)
        # Which coordinate system the lanes are drawn in. Hidden entirely when
        # there are no lanes; never silently omitted when there are.
        self.coord_label = QLabel()
        self.coord_label.setWordWrap(True)
        self.coord_label.setStyleSheet(f"QLabel {{ color:{muted_hex()}; }}")
        self.coord_label.setVisible(False)
        viz_l.addWidget(self.coord_label)
        # Lane legend — pattern, not hue, so it reads in a greyscale figure.
        self.lane_legend = QLabel(
            "Lanes: <b>solid</b> = probe binds this sequence&nbsp;&nbsp;"
            "<b>textured</b> = binds with mismatches&nbsp;&nbsp;"
            "<b>hollow / hatched, dashed edge</b> = does not bind"
        )
        self.lane_legend.setTextFormat(Qt.RichText)
        self.lane_legend.setStyleSheet(f"QLabel {{ color:{muted_hex()}; }}")
        self.lane_legend.setVisible(False)
        viz_l.addWidget(self.lane_legend)
        # Spacing legend — each tier carries its LETTER as well as its hue, and
        # the letter is what the probe bars are stamped with, so the map is
        # readable without colour vision and in a greyscale figure.
        sp = self._spacing
        legend = QLabel(
            f"Spacing: "
            f"<span style='color:{_GAP_TIER_COLORS[GAP_TIER_GREEN]}'>●</span> "
            f"<b>G</b> green, ≥{sp.min_gap} bp on target&nbsp;&nbsp;"
            f"<span style='color:{_GAP_TIER_COLORS[GAP_TIER_ORANGE]}'>●</span> "
            f"<b>O</b> orange, {sp.relaxed_gap}–{sp.min_gap} bp tight&nbsp;&nbsp;"
            f"<span style='color:{_GAP_TIER_COLORS[GAP_TIER_RED]}'>●</span> "
            f"<b>R</b> red, &lt;{sp.relaxed_gap} bp crowded "
            f"(floor {sp.floor_gap} bp)"
        )
        legend.setTextFormat(Qt.RichText)
        legend.setStyleSheet(f"QLabel {{ color:{muted_hex()}; }}")
        legend.setAccessibleName("Probe spacing legend")
        viz_l.addWidget(legend)
        # WHERE those three numbers came from. Never omitted: a dialog that
        # silently substitutes the application defaults for a run's own settings
        # is exactly the failure this line exists to make visible.
        self.spacing_source_label = QLabel(self._spacing_source_sentence())
        self.spacing_source_label.setWordWrap(True)
        self.spacing_source_label.setAccessibleName("Spacing threshold source")
        self.spacing_source_label.setStyleSheet(
            f"QLabel {{ color:{muted_hex()}; }}" if sp.from_run_config
            else "QLabel { color:#8a5a00; }")
        viz_l.addWidget(self.spacing_source_label)
        # Running tier tally for the whole selection — the number a user needs
        # before ordering, and the one the CLI prints at the end of a headless
        # run, so the two agree.
        self.tier_label = QLabel()
        self.tier_label.setTextFormat(Qt.RichText)
        viz_l.addWidget(self.tier_label)
        self.overlap_warning = QLabel()
        self.overlap_warning.setStyleSheet("color:#b00020; font-weight:bold;")
        viz_l.addWidget(self.overlap_warning)
        outer.addWidget(viz_frame)

        # Filter bar (U3).
        outer.addWidget(self._build_filter_bar())

        # Candidates table.
        self.table = QTableWidget(0, len(self._cols))
        self.table.setHorizontalHeaderLabels([c.label for c in self._cols])
        for i, c in enumerate(self._cols):
            self.table.setColumnWidth(i, c.width)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSortingEnabled(True)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setAccessibleName("Candidate probes")
        self.table.setAccessibleDescription(
            "Every ranked candidate for the selected gene. Tick the Select "
            "column to include a probe in the export.")
        outer.addWidget(self.table, 1)

        self.counter_label = QLabel(); self.counter_label.setStyleSheet(f"color:{muted_hex()};")
        self.counter_label.setAccessibleName("Selection counter")
        outer.addWidget(self.counter_label)
        self.tier_label.setAccessibleName("Spacing tally of the selection")
        self.overlap_warning.setAccessibleName("Spacing warning")

        btns = QDialogButtonBox()
        self.save_btn = btns.addButton("Save &draft", QDialogButtonBox.ActionRole)
        self.save_btn.setAccessibleName("Save selections as a draft")
        # "&&" is a literal ampersand in Qt; "e" carries the mnemonic.
        self.confirm_btn = btns.addButton("Confirm && &export",
                                          QDialogButtonBox.AcceptRole)
        self.confirm_btn.setAccessibleName("Confirm and export order CSVs")
        cancel_btn = btns.addButton(QDialogButtonBox.Cancel)
        outer.addWidget(btns)

        # An explicit, reading-order tab chain. Qt's default follows widget
        # CREATION order, which here interleaves the map's decorations with the
        # controls; this is the order the dialog is actually used in.
        self._set_tab_order([
            self.gene_combo, self.max_spin, self.spacing_filter,
            self.autopick_btn, self.coverage_btn,
            self.f_gc_min, self.f_gc_max, self.f_tm_min, self.f_score_min,
            self.f_text, self.f_clear, self.show_parents_check,
            self.table,
            self.save_btn, self.confirm_btn, cancel_btn,
        ])

        self.gene_combo.currentIndexChanged.connect(self._on_gene_changed)
        self.max_spin.valueChanged.connect(self._on_max_changed)
        self.spacing_filter.currentIndexChanged.connect(
            self._on_spacing_filter_changed)
        self.autopick_btn.clicked.connect(self._on_autopick)
        self.coverage_btn.clicked.connect(self._on_coverage_pick)
        self.table.itemChanged.connect(self._on_item_changed)
        self.save_btn.clicked.connect(self._on_save_draft)
        self.confirm_btn.clicked.connect(self._on_confirm)
        cancel_btn.clicked.connect(self.reject)

        if self._load_aborted:
            # The load was cancelled or failed. exec() reports it and returns
            # Rejected; do NOT also claim the run has no candidates, which is a
            # statement about the data rather than about the load.
            self.confirm_btn.setEnabled(False)
        elif self.gene_combo.count() == 0:
            QMessageBox.warning(self, "No candidates",
                                "results.csv is empty or unreadable; nothing to review.")
            self.confirm_btn.setEnabled(False)
        else:
            self._on_gene_changed(0)

    # ------------------------------------------------------------------
    # Accessibility
    # ------------------------------------------------------------------

    def _set_tab_order(self, widgets) -> None:
        """Chain ``widgets`` into an explicit tab order, in the given order."""
        for first, second in zip(widgets, widgets[1:]):
            self.setTabOrder(first, second)

    # ------------------------------------------------------------------
    # Spacing configuration
    # ------------------------------------------------------------------

    def _spacing_source_sentence(self) -> str:
        """One sentence naming where the three rungs came from.

        Shown in the dialog and repeated in the two spacing tooltips. When there
        is no run_config.json the fallback is stated, not hidden: the numbers on
        screen would otherwise be presented as this run's settings when they are
        the application's."""
        sp = self._spacing
        if sp.from_run_config:
            return (f"Spacing thresholds {sp.min_gap} / {sp.relaxed_gap} / "
                    f"{sp.floor_gap} bp are this run's own "
                    f"(read from run_config.json).")
        return (f"\u26a0 No run_config.json in this run's output folder, so the "
                f"spacing thresholds shown ({sp.min_gap} / {sp.relaxed_gap} / "
                f"{sp.floor_gap} bp) are the application defaults, not "
                f"necessarily what this run was designed at. Colours, the "
                f"Spacing filter and re-picks all use these defaults.")

    # ------------------------------------------------------------------
    # Strain columns
    # ------------------------------------------------------------------

    def _strain_columns(self):
        """One (mm, ΔTm) column pair per strain, for multi-strain runs (U1)."""
        cols = []
        for sid in self._strain_ids:
            cols.append(_Col(f"mm:{sid}", 70, "num", f"mismatches__{sid}"))
            cols.append(_Col(f"ΔTm:{sid}", 74, "num", f"est_dTm__{sid}"))
        return cols

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_run_data(self) -> None:
        """Read the run off a worker thread, with progress and a live Cancel.

        The parse is genuinely off the UI thread; this call blocks on a NESTED
        event loop rather than on the read, so the progress bar animates, Cancel
        is clickable, and the window manager keeps getting served. Cancel sets
        the worker's flag, which the streaming CSV reader polls, so the read
        actually stops instead of merely being hidden.
        """
        worker = _RunDataWorker(self._record.output_dir, self)
        loop = QEventLoop(self)

        # 0-100 determinate: the reader knows the file size, and "47% of 1.17
        # GB" is the only honest thing to show for a read that takes minutes.
        progress = QProgressDialog(self)
        progress.setLabelText(
            f"Loading candidates from {os.path.basename(self._record.output_dir)}…")
        progress.setWindowTitle("Loading run")
        progress.setRange(0, 100)
        progress.setWindowModality(Qt.WindowModal)
        # Nothing at all for the small runs that finish in a blink; the bar
        # appears by itself once the read is visibly slow.
        progress.setMinimumDuration(400)
        progress.setAutoClose(False)
        progress.setAutoReset(False)

        result: Dict[str, Any] = {}

        def _done(data):
            result["data"] = data
            loop.quit()

        def _failed(msg):
            result["error"] = msg
            loop.quit()

        def _cancelled():
            result["cancelled"] = True
            loop.quit()

        worker.progress.connect(progress.setValue)
        worker.done.connect(_done)
        worker.failed.connect(_failed)
        worker.cancelled.connect(_cancelled)
        progress.canceled.connect(worker.request_cancel)

        worker.start()
        progress.setValue(0)
        loop.exec()
        # QProgressDialog.close() emits canceled() from its own closeEvent;
        # take the handler off before closing so a finished load is never
        # retro-actively cancelled.
        progress.canceled.disconnect(worker.request_cancel)
        # Never leave the QThread running past this point — the dialog owns it
        # and Qt would destroy it out from under the thread.
        worker.wait()
        progress.close()
        worker.deleteLater()

        if result.get("cancelled"):
            self._load_aborted = True
            return
        if "error" in result:
            self._load_aborted = True
            self._load_error = str(result["error"])
            return
        self._adopt_run_data(result.get("data") or RunData())

    def _adopt_run_data(self, data: RunData) -> None:
        self._rows_by_gene = data.rows_by_gene
        self._gene_lengths = data.gene_lengths
        self._strain_ids = data.strain_ids
        self._gene_isoforms = data.gene_isoforms
        self._isoform_universe = data.isoform_universe
        self._has_iso_col = data.has_iso_col
        self._has_strain_cols = data.has_strain_cols
        self._variants = data.variants
        self._selections = data.selections
        # Computed once, here, rather than per repaint: _populate_table runs on
        # every filter keystroke and this would otherwise walk the whole run.
        self._has_degenerate = any(
            str(r.get("probe_role") or "").strip() == PROBE_ROLE_CONSENSUS
            for rows in self._rows_by_gene.values() for r in rows
        )

    def exec(self):
        """Refuse to show a dialog whose data never arrived.

        A cancelled or failed load leaves every table, map and counter empty;
        showing that would look like "this run has no candidates", which is a
        different and much worse claim than "you stopped the load".
        """
        if self._load_aborted:
            if self._load_error:
                QMessageBox.critical(self.parent(), "Could not open run",
                                     self._load_error)
            return QDialog.Rejected
        return super().exec()

    # ------------------------------------------------------------------
    # Filter bar
    # ------------------------------------------------------------------

    def _build_filter_bar(self) -> QWidget:
        box = QFrame(); box.setFrameShape(QFrame.StyledPanel)
        g = QGridLayout(box); g.setContentsMargins(8, 6, 8, 6); g.setHorizontalSpacing(10)
        box.setAccessibleName("Candidate filters")
        g.addWidget(QLabel("Filter:"), 0, 0)

        # A placeholder is not a label: it vanishes on focus and is not exposed
        # as one to assistive technology. Every control here gets a real
        # accessible name, and a buddy wherever it has a visible label.
        self.f_gc_min = QDoubleSpinBox(); self.f_gc_min.setRange(0, 1); self.f_gc_min.setSingleStep(0.01); self.f_gc_min.setValue(0.0)
        self.f_gc_max = QDoubleSpinBox(); self.f_gc_max.setRange(0, 1); self.f_gc_max.setSingleStep(0.01); self.f_gc_max.setValue(1.0)
        gc_lbl = QLabel("G&C"); gc_lbl.setBuddy(self.f_gc_min)
        self.f_gc_min.setAccessibleName("Minimum GC fraction")
        self.f_gc_max.setAccessibleName("Maximum GC fraction")
        g.addWidget(gc_lbl, 0, 1); g.addWidget(self.f_gc_min, 0, 2); g.addWidget(QLabel("–"), 0, 3); g.addWidget(self.f_gc_max, 0, 4)

        self.f_tm_min = QDoubleSpinBox(); self.f_tm_min.setRange(0, 120); self.f_tm_min.setValue(0.0)
        tm_lbl = QLabel("&Tm ≥"); tm_lbl.setBuddy(self.f_tm_min)
        self.f_tm_min.setAccessibleName("Minimum Tm of both probe halves")
        self.f_tm_min.setAccessibleDescription(
            "Rows are kept only when BOTH half melting temperatures are at or "
            "above this. Zero disables the filter.")
        g.addWidget(tm_lbl, 0, 5); g.addWidget(self.f_tm_min, 0, 6)

        self.f_score_min = QDoubleSpinBox(); self.f_score_min.setRange(0, 1); self.f_score_min.setSingleStep(0.01); self.f_score_min.setValue(0.0)
        score_lbl = QLabel("Sc&ore ≥"); score_lbl.setBuddy(self.f_score_min)
        self.f_score_min.setAccessibleName("Minimum composite score")
        g.addWidget(score_lbl, 0, 7); g.addWidget(self.f_score_min, 0, 8)

        self.f_text = QLineEdit(); self.f_text.setPlaceholderText("probe_id contains…")
        self.f_text.setAccessibleName("Probe ID contains")
        self.f_text.setAccessibleDescription(
            "Keep only candidates whose probe ID contains this text.")
        g.addWidget(self.f_text, 0, 9)

        self.f_clear = QPushButton("C&lear")
        self.f_clear.setAccessibleName("Clear all candidate filters")
        g.addWidget(self.f_clear, 0, 10)
        g.setColumnStretch(9, 1)

        # Degenerate consensus parents are OFF by default. A parent is a
        # mixed-base IUPAC string: it is never ordered and never exported, and
        # it sits on the same coordinates as its own children, so leaving it in
        # the table asks the user to discount a row on every single family. It
        # is still worth being able to see, because the parent is the only row
        # that shows the degenerate positions and the IUPAC codes, hence a toggle
        # rather than a deletion. Absent on a run with no families at all: a
        # control for something the data does not contain is noise.
        self.show_parents_check = QCheckBox("Show consensus &parents")
        self.show_parents_check.setChecked(False)
        self.show_parents_check.setAccessibleName(
            "Show degenerate consensus parent rows")
        self.show_parents_check.setAccessibleDescription(
            "Off by default. A consensus parent is a mixed-base sequence that "
            "is never ordered and never exported; its variant children are the "
            "concrete oligonucleotides. Turn this on to see the parent row, "
            "which carries the IUPAC codes and the degenerate positions.")
        self.show_parents_check.setVisible(self._has_degenerate)
        g.addWidget(self.show_parents_check, 0, 11)

        for w in (self.f_gc_min, self.f_gc_max, self.f_tm_min, self.f_score_min):
            w.valueChanged.connect(self._populate_table)
        self.f_text.textChanged.connect(self._populate_table)
        self.f_clear.clicked.connect(self._clear_filters)
        self.show_parents_check.toggled.connect(self._populate_table)
        return box

    def _clear_filters(self) -> None:
        self.f_gc_min.setValue(0.0); self.f_gc_max.setValue(1.0)
        self.f_tm_min.setValue(0.0); self.f_score_min.setValue(0.0)
        self.f_text.clear()
        self._populate_table()

    def _passes_filter(self, r: Dict[str, str]) -> bool:
        def _f(key):
            try:
                return float(r.get(key, "") or "nan")
            except ValueError:
                return float("nan")
        gc = _f("gc_content_full")
        if gc == gc:  # not NaN
            if not (self.f_gc_min.value() <= gc <= self.f_gc_max.value()):
                return False
        tm_min = self.f_tm_min.value()
        if tm_min > 0:
            # Filter on the same per-half Tm columns shown in the table
            # (Tm LHS / Tm RHS). The floor means "both halves melt at or above
            # this" — the weaker half limits hybridisation. n/a Tm rows are
            # excluded only when the user sets a Tm floor.
            tm_lhs = _f("tm_target_lhs")
            tm_rhs = _f("tm_target_rhs")
            if not (tm_lhs == tm_lhs and tm_rhs == tm_rhs
                    and tm_lhs >= tm_min and tm_rhs >= tm_min):
                return False
        score = _f("composite_score")
        if score == score and score < self.f_score_min.value():
            return False
        txt = self.f_text.text().strip().lower()
        if txt and txt not in r.get("probe_id", "").lower():
            return False
        return True

    # ------------------------------------------------------------------
    # Header / counters
    # ------------------------------------------------------------------

    def _refresh_header(self) -> None:
        n_genes = len(self._rows_by_gene)
        n_cand = sum(len(v) for v in self._rows_by_gene.values())
        strain_note = (
            f"  ·  homology mode ({len(self._strain_ids)} sequences: "
            f"{', '.join(self._strain_ids)})" if self._strain_ids else ""
        )
        screened = (self._record.stats or {}).get("specificity_screened")
        warn = ("  ·  <span style='color:#b00020'>⚠ probes were NOT "
                "off-target screened</span>") if screened is False else ""
        self._header_label.setText(
            f"<b>{self._record.name}</b> — {self._record.modality_label} · "
            f"{n_genes} gene(s), {n_cand} candidates{strain_note}{warn}<br>"
            f"Pick up to N probes per gene, then <b>Confirm &amp; export</b>."
        )

    def _current_gene(self) -> Optional[str]:
        idx = self.gene_combo.currentIndex()
        return self.gene_combo.itemData(idx) if idx >= 0 else None

    def selection_tier_counts(self) -> Dict[str, int]:
        """Selected SITES per spacing tier, across the WHOLE selection.

        Each site is tiered by its OWN nearest other site in the same gene,
        exactly as the map colours its bars, and counted once however many
        variant ids it carries. So the three tiers add up to the same site
        count the "Selected" counter shows.

        Until 2026-09-22 every id in a gene took the tier of that gene's single
        closest pair and was counted per id: one crowded pair turned a gene's
        well-spaced probes red here while the map drew them green, and 3
        degenerate sites read "green 6" beside "Selected: 3 / 3".

        Recomputed from the selected intervals rather than carried from whatever
        made the selection, so hand-edited picks are counted on the same footing
        as auto-picked ones and the number cannot drift out of date.
        """
        counts = {name: 0 for name in GAP_TIER_NAMES.values()}
        for gid, ids in self._selections.items():
            by_id = {r.get("probe_id"): r for r in self._rows_by_gene.get(gid, [])}
            # One interval per SITE, not per id: a degenerate family's variants
            # share a footprint. See site_key.
            by_site: Dict[str, Tuple[int, int]] = {}
            for pid in ids:
                r = by_id.get(pid)
                if r is None:
                    continue
                try:
                    by_site.setdefault(
                        site_key(r),
                        (int(r.get("start") or 0), int(r.get("end") or 0)))
                except (TypeError, ValueError):
                    pass
            ivs = list(by_site.values())
            for i, iv in enumerate(ivs):
                closest: Optional[int] = None
                for j, other in enumerate(ivs):
                    if i == j:
                        continue
                    g = _interval_gap(iv, other)
                    closest = g if closest is None else min(closest, g)
                counts[GAP_TIER_NAMES[self._spacing.tier(closest)]] += 1
        return counts

    def _refresh_counter(self) -> None:
        gid = self._current_gene()
        if gid is None:
            self.counter_label.setText(""); return
        n = len(self._sites(gid, self._selections.get(gid, [])))
        total = sum(len(self._sites(g, v)) for g, v in self._selections.items())
        shown = self.table.rowCount()
        all_n = len(self._rows_by_gene.get(gid, []))
        tiers = self.selection_tier_counts()
        self.counter_label.setText(
            f"Selected for {gid}: {n} / {self._max_per_gene}  ·  "
            f"Total across genes: {total}  ·  showing {shown}/{all_n} candidates"
        )
        # Always all three tiers, even at zero: a count that disappears when the
        # news is good is a count the user has to go and verify by hand.
        self.tier_label.setText(
            "Spacing of selected probe sites:  "
            f"<span style='color:{_GAP_TIER_COLORS[GAP_TIER_GREEN]}'>●</span> "
            f"green {tiers['green']}&nbsp;&nbsp;"
            f"<span style='color:{_GAP_TIER_COLORS[GAP_TIER_ORANGE]}'>●</span> "
            f"orange {tiers['orange']}&nbsp;&nbsp;"
            f"<span style='color:{_GAP_TIER_COLORS[GAP_TIER_RED]}'>●</span> "
            f"red {tiers['red']}"
        )

    # ------------------------------------------------------------------
    # Table population
    # ------------------------------------------------------------------

    @staticmethod
    def _make_cell(col: "_Col", r: Dict[str, str], pid: str, selected_set) -> QTableWidgetItem:
        """Build one table cell from a column spec + a results row, dispatching
        on ``col.kind``. Keeps cell construction declarative so adding/moving a
        column is a single registry edit."""
        if col.kind == "select":
            sel = QTableWidgetItem()
            sel.setData(Qt.UserRole, pid)
            if str(r.get("probe_role") or "").strip() == PROBE_ROLE_CONSENSUS:
                # Visible (with the toggle on) but never orderable: a parent is
                # a mixed-base string, and FlexiProbes only exports
                # oligonucleotides that are concrete observed sequences. Offering a
                # tick box would be offering something the exporter refuses.
                sel.setFlags(Qt.ItemIsEnabled | Qt.ItemIsSelectable)
                sel.setToolTip(
                    "Degenerate consensus parent, not orderable. Its variant "
                    "children carry the concrete oligonucleotides; select those.")
                return sel
            sel.setFlags(Qt.ItemIsEnabled | Qt.ItemIsUserCheckable | Qt.ItemIsSelectable)
            sel.setCheckState(Qt.Checked if pid in selected_set else Qt.Unchecked)
            return sel
        value = r.get(col.source, "")
        if col.kind == "gc":
            return _gc_item(value)
        if col.kind == "tm":
            return _tm_item(value)
        if col.kind == "spec":
            return _spec_item(value)
        if col.kind == "num":
            return _num_item(value)
        return QTableWidgetItem(str(value))

    def _row_is_visible(self, r: Dict[str, str]) -> bool:
        """Whether a row belongs in the table at all, before any user filter.

        The only rows this hides are degenerate consensus parents, and only
        while the toggle is off. This is a VIEW decision and nothing else:
        hiding a parent does not change what is selected, spaced or exported,
        because a parent was never any of those things.
        """
        if str(r.get("probe_role") or "").strip() != PROBE_ROLE_CONSENSUS:
            return True
        return self.show_parents_check.isChecked()

    def _populate_table(self) -> None:
        gid = self._current_gene()
        if gid is None:
            return
        rows = [r for r in self._rows_by_gene.get(gid, [])
                if self._passes_filter(r) and self._row_is_visible(r)]
        selected_set = set(self._selections.get(gid, []))

        self.table.blockSignals(True)
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)
        for r in rows:
            row = self.table.rowCount()
            self.table.insertRow(row)
            pid = r.get("probe_id", "")
            for col_idx, col in enumerate(self._cols):
                self.table.setItem(row, col_idx, self._make_cell(col, r, pid, selected_set))
        self.table.setSortingEnabled(True)
        self.table.blockSignals(False)
        self._refresh_counter()

    # ------------------------------------------------------------------
    # Lane block (one lane per isoform / per input record)
    # ------------------------------------------------------------------

    @staticmethod
    def _split_ids(value) -> List[str]:
        """Split a ';'-joined id column. Tolerates the column being absent —
        every one of these is optional and older results.csv files lack them."""
        return [s for s in (str(value or "")).split(";") if s.strip()]

    def _record_lane_mark(self, row, probe, sid: str) -> LaneMark:
        """How ``probe`` relates to input record ``sid``.

        Reads the degenerate-consensus contract (probe_role / resolves_to)
        defensively: those columns do not exist in older results.csv files, and
        when they are missing or blank this falls straight back to the
        per-record likely_binds / mismatches columns.
        """
        role = str(row.get("probe_role") or "").strip()
        resolves = self._split_ids(row.get("resolves_to"))
        note = ""
        if role == PROBE_ROLE_VARIANT and resolves:
            # A variant is the concrete oligo for its own record(s) only; on any
            # other lane the sibling variant is the one that binds.
            if sid not in resolves:
                return LaneMark(
                    probe_id=probe["probe_id"], rank=probe["rank"],
                    start=probe["start"], end=probe["end"], state=MARK_ABSENT,
                    detail=f"variant resolves to {', '.join(resolves)} — not "
                           "this record")
            note = "variant for this record"
        elif role == PROBE_ROLE_CONSENSUS:
            note = "degenerate consensus parent"

        mm_raw = row.get(f"mismatches__{sid}")
        binds_raw = row.get(f"likely_binds__{sid}")
        try:
            mm = int(float(mm_raw)) if str(mm_raw or "").strip() else None
        except (TypeError, ValueError):
            mm = None
        binds_s = str(binds_raw or "").strip().lower()
        if binds_s in ("false", "0", "no"):
            state = MARK_ABSENT
        elif mm:
            state = MARK_MISMATCH
        else:
            state = MARK_BINDS
        bits = [note] if note else []
        bits.append("does not bind" if state == MARK_ABSENT else "binds")
        if mm is not None:
            bits.append(f"{mm} mismatch(es)")
        dtm = str(row.get(f"est_dTm__{sid}") or "").strip()
        if dtm:
            bits.append(f"ΔTm {dtm}")
        return LaneMark(probe_id=probe["probe_id"], rank=probe["rank"],
                        start=probe["start"], end=probe["end"], state=state,
                        detail=", ".join(b for b in bits if b))

    def _build_lane_set(self, gid: str, probes: Sequence[Dict[str, Any]]) -> LaneSet:
        """Lanes for the current gene, or an empty LaneSet when the run has
        nothing extra to show (single sequence, no isoform data, old columns
        missing) — in which case the map falls back to today's single track."""
        rows_by_pid = {r.get("probe_id", ""): r
                       for r in self._rows_by_gene.get(gid, [])}

        # (a) Multi-sequence / multi-strain: results.csv coordinates ARE
        # alignment coordinates (every record was aligned to the consensus the
        # probes were tiled on), so one shared x-axis is genuinely correct here.
        if self._strain_ids and self._has_strain_cols:
            lanes: List[SequenceLane] = []
            for sid in self._strain_ids:
                rec_marks = tuple(
                    self._record_lane_mark(rows_by_pid.get(p["probe_id"], {}), p, sid)
                    for p in probes
                )
                bound = sum(1 for m in rec_marks if m.state != MARK_ABSENT)
                lanes.append(SequenceLane(
                    label=sid, marks=rec_marks,
                    note=(f"{bound}/{len(rec_marks)} bind" if rec_marks else ""),
                ))
            return LaneSet(
                lanes=tuple(lanes), kind="record",
                coord_caption=(
                    f"Lanes: {len(lanes)} input sequence(s). "
                    "<b>x-axis = alignment (consensus) coordinates</b>, shared by "
                    "every lane and by the summary track above — the probes were "
                    "tiled on the consensus, so one position means the same base "
                    "on every record."),
            )

        # (b) Isoform-aware: membership comes from isoform_hit_ids, which says
        # WHICH isoforms a probe binds but not WHERE on each of them.
        gene_entry = self._isoform_universe.get(gid)
        if self._has_iso_col and gene_entry and len(gene_entry[1]) > 1:
            gene, members = gene_entry
            # Source transcript first: its coordinates are the axis everything
            # else is drawn against, so it belongs at the top of the block.
            ordered = [gid] + sorted(t for t in members if t != gid)
            lanes = []
            for tid in ordered:
                marks: List[LaneMark] = []
                for p in probes:
                    row = rows_by_pid.get(p["probe_id"], {})
                    hits = set(self._split_ids(row.get("isoform_hit_ids")))
                    hits.add(gid)          # a probe always binds its own source
                    binds = tid in hits
                    marks.append(LaneMark(
                        probe_id=p["probe_id"], rank=p["rank"],
                        start=p["start"], end=p["end"],
                        state=MARK_BINDS if binds else MARK_ABSENT,
                        detail=("binds (listed in isoform_hit_ids)" if binds
                                else "not in isoform_hit_ids — does not bind "
                                     "this isoform"),
                    ))
                bound = sum(1 for m in marks if m.state != MARK_ABSENT)
                lanes.append(SequenceLane(
                    label=tid + ("  [source]" if tid == gid else ""),
                    marks=tuple(marks),
                    length=self._gene_lengths.get(tid, 0),
                    note=(f"{bound}/{len(marks)} bind" if marks else ""),
                ))
            extra = max(0, len(lanes) - MAX_LANES)
            return LaneSet(
                lanes=tuple(lanes), kind="isoform",
                coord_caption=(
                    f"Lanes: {len(lanes)} isoform(s) of {gene}"
                    + (f" ({extra} not drawn)" if extra else "") + ". "
                    "<b>x-axis = coordinates on the source transcript "
                    f"{gid}</b> — this run records which isoforms each probe "
                    "binds (isoform_hit_ids) but not the probe's offset on each, "
                    "so a probe is drawn at its source-transcript position on "
                    "every lane it binds. Each isoform's own length is shown "
                    "beside its label."),
            )

        return LaneSet()

    # ------------------------------------------------------------------
    # Visualization sync
    # ------------------------------------------------------------------

    def _sync_map_height(self) -> None:
        """Show the whole map when it fits, scroll it when it does not."""
        self.map_scroll.setFixedHeight(
            min(self.map_widget.content_height(), MAX_MAP_VIEWPORT_H))

    def _refresh_viz(self) -> None:
        gid = self._current_gene()
        if gid is None:
            self.map_widget.set_data("", 0, [])
            self._sync_map_height()
            self.overlap_warning.setText(""); return
        sel_ids = set(self._selections.get(gid, []))
        probes: List[Dict[str, Any]] = []
        for r in self._rows_by_gene.get(gid, []):
            if r.get("probe_id") in sel_ids:
                try:
                    probes.append({
                        "probe_id": r.get("probe_id", ""),
                        # Carried so the map's tiering can tell a family's
                        # siblings apart from genuine neighbours. See site_key.
                        "site": site_key(r),
                        "rank": int(r.get("rank") or 0),
                        "start": int(r.get("start") or 0),
                        "end": int(r.get("end") or 0),
                    })
                except ValueError:
                    pass
        probes.sort(key=lambda p: p["start"])
        # Built once per selection change, never in paintEvent.
        lane_set = self._build_lane_set(gid, probes)
        self.map_widget.set_data(gid, self._gene_lengths.get(gid, 0), probes,
                                 self._variants, lane_set)
        self.coord_label.setText(lane_set.coord_caption)
        self.coord_label.setVisible(bool(lane_set.lanes))
        self.lane_legend.setVisible(bool(lane_set.lanes))
        self._sync_map_height()
        if self.map_widget.has_overlap:
            self.overlap_warning.setStyleSheet("color:#b00020; font-weight:bold;")
            self.overlap_warning.setText(
                f"⚠ Selected probes are closer than "
                f"{self._spacing.relaxed_gap} bp "
                "(red bars, stamped R). Auto-selection only reaches this on "
                f"the ladder's last rung ({self._spacing.floor_gap} bp). "
                "Their "
                "transcript may simply be too short to hold this many separated "
                "probes; consider fewer probes for this gene.")
        elif self.map_widget.has_tight:
            self.overlap_warning.setStyleSheet(f"color:{_GAP_TIER_COLORS[1]};")
            self.overlap_warning.setText(
                f"Some selected probes are {self._spacing.relaxed_gap}–"
                f"{self._spacing.min_gap} bp apart (orange, stamped O) — closer "
                f"than the {self._spacing.min_gap} bp target. Usable, and "
                f"carried into the exports as such; ≥{self._spacing.min_gap} bp "
                f"(green, G) is the aim.")
        else:
            self.overlap_warning.setText("")

    # ------------------------------------------------------------------
    # Event handlers
    # ------------------------------------------------------------------

    def _on_gene_changed(self, _idx) -> None:
        self._populate_table(); self._refresh_viz(); self._refresh_counter()

    def _spacing_floor(self) -> int:
        """Tightest gap the user currently permits, in bp."""
        data = self.spacing_filter.currentData()
        return int(data[0]) if data else self._spacing.floor_gap

    def _min_allowed_tier(self) -> int:
        data = self.spacing_filter.currentData()
        return int(data[1]) if data else GAP_TIER_RED

    def _pick_ladder(self):
        """The spacing ladder truncated to what the filter permits.

        Built from the floor rather than by filtering results, so a forbidden
        tier is never reachable: the selector simply has fewer rungs to walk and
        returns fewer probes for a gene that needed one.
        """
        floor = self._spacing_floor()
        return ladder_from(self._spacing.min_gap,
                           max(self._spacing.relaxed_gap, floor),
                           max(self._spacing.floor_gap, floor))

    def _on_spacing_filter_changed(self, _idx) -> None:
        """Drop existing selections the new filter no longer permits.

        A filter that only governed future picks would leave the user looking at
        red bars they just asked never to see, and would let those probes reach
        the order CSV. Manual picks are dropped by the same rule as automatic
        ones — the filter is about what may be ordered, not about who chose it.
        """
        floor = self._spacing_floor()
        dropped = 0
        for gid, ids in list(self._selections.items()):
            keep = self._prune_to_floor(gid, ids, floor)
            if len(keep) != len(ids):
                dropped += len(ids) - len(keep)
                if keep:
                    self._selections[gid] = keep
                else:
                    self._selections.pop(gid, None)
        self._populate_table(); self._refresh_viz(); self._refresh_counter()
        if dropped:
            self.overlap_warning.setStyleSheet(f"color:{_GAP_TIER_COLORS[1]};")
            self.overlap_warning.setText(
                f"Spacing filter dropped {dropped} already-selected probe(s) "
                f"that sat closer than {floor} bp.")

    def _prune_to_floor(self, gid, ids, floor: int):
        """Keep the highest-ranked subset of ``ids`` no two closer than ``floor``.

        Walks in the table's own rank order so the survivor set is the one the
        user would expect to keep — the best probes — rather than whichever
        happened to come first in the file.
        """
        by_id = {r.get("probe_id"): r for r in self._rows_by_gene.get(gid, [])}

        def _rank(pid):
            try:
                return int(by_id.get(pid, {}).get("rank") or 1 << 30)
            except (TypeError, ValueError):
                return 1 << 30

        # Measured between SITES (site_key). A family's variants share one
        # footprint, so comparing ids read siblings as -1 bp apart and dropped
        # every variant after the first, even at the default 15 bp floor.
        keep: List[str] = []
        kept_sites: Dict[str, Tuple[int, int]] = {}
        for pid in sorted(ids, key=_rank):
            r = by_id.get(pid)
            if r is None:
                keep.append(pid)         # not a row we can place; leave it alone
                continue
            site = site_key(r)
            if site in kept_sites:
                keep.append(pid)         # a sibling of a site already kept
                continue
            try:
                s0, e0 = int(r.get("start") or 0), int(r.get("end") or 0)
            except (TypeError, ValueError):
                keep.append(pid)
                continue
            if any(_interval_gap((s0, e0), iv) < floor
                   for iv in kept_sites.values()):
                continue
            keep.append(pid)
            kept_sites[site] = (s0, e0)
        # Preserve the caller's original ordering among the survivors.
        survivors = set(keep)
        return [pid for pid in ids if pid in survivors]

    # ---- the per-gene limit counts SITES ----
    #
    # auto_select counts a degenerate family once ("one family = one pick") and
    # emits its children, so a correct pick of N sites can be more than N ids.
    # Every limit check and count below therefore works in sites (site_key),
    # never in ids. Counting ids showed "6 / 3" on a correct pick, cut families
    # in half when the limit was lowered, and refused ticks under the limit.

    def _site_by_id(self, gid) -> Dict[str, str]:
        return {str(r.get("probe_id", "")): site_key(r)
                for r in self._rows_by_gene.get(gid, [])}

    def _sites(self, gid, ids) -> List[str]:
        """Distinct sites of ``ids``, in first-seen order."""
        of = self._site_by_id(gid)
        out: List[str] = []
        for pid in ids:
            site = of.get(pid, pid)
            if site not in out:
                out.append(site)
        return out

    def _within_site_limit(self, gid, ids, limit) -> List[str]:
        """``ids`` restricted to their first ``limit`` sites. A family is kept
        or dropped whole, never split."""
        of = self._site_by_id(gid)
        keep = set(self._sites(gid, ids)[:limit])
        return [pid for pid in ids if of.get(pid, pid) in keep]

    def _on_max_changed(self, value) -> None:
        self._max_per_gene = int(value)
        for gid, ids in list(self._selections.items()):
            if len(self._sites(gid, ids)) > self._max_per_gene:
                self._selections[gid] = self._within_site_limit(
                    gid, ids, self._max_per_gene)
        self._populate_table(); self._refresh_viz(); self._refresh_counter()

    def _on_autopick(self) -> None:
        gid = self._current_gene()
        if gid is None:
            return
        rows = self._rows_by_gene.get(gid, [])
        picked = auto_pick_per_gene(
            rows, n_per_gene=self._max_per_gene, ladder=self._pick_ladder(),
            gene_length_of=lambda g: self._gene_lengths.get(g, 0),
        )
        self._selections[gid] = picked.get(gid, [])
        self._populate_table(); self._refresh_viz(); self._refresh_counter()

    # ------------------------------------------------------------------
    # Isoform-coverage auto-pick (whole run, gene-aware)
    # ------------------------------------------------------------------

    def _coverage_source(self):
        """How to get isoform membership for the coverage pick, preferring the
        FASTA-free path. Returns ``("column", gene_isoforms)`` when the run
        persisted ``isoform_hit_ids`` + a gene map, ``("fasta", path)`` when the
        input FASTA is on disk, else ``(None, None)``."""
        if self._has_iso_col and self._gene_isoforms:
            return "column", self._gene_isoforms
        path = getattr(self._record, "fasta_path", "") or ""
        if path and os.path.isfile(path):
            return "fasta", path
        return None, None

    def _on_coverage_pick(self) -> None:
        # Re-entrancy guard: a coverage pass is already running. The progress
        # dialog is modal, but guard anyway so a stray double-trigger can't spawn
        # a second worker and orphan the first.
        if self._cov_worker is not None and self._cov_worker.isRunning():
            return
        mode, payload = self._coverage_source()
        if mode is None:
            QMessageBox.warning(
                self, "Can't maximise isoform coverage",
                "This run doesn't carry the isoform_hit_ids column and its input "
                "cDNA FASTA couldn't be found on disk, so isoform membership is "
                "unavailable. Re-run the design, or re-import with the FASTA.")
            return
        if any(self._selections.values()):
            if QMessageBox.question(
                self, "Replace current selections?",
                "This picks probes for every gene to maximise isoform coverage "
                "and will replace your current selections. Continue?",
                QMessageBox.Yes | QMessageBox.No,
            ) != QMessageBox.Yes:
                return

        # Flatten every candidate row back out for the whole-run pass.
        rows = [r for rs in self._rows_by_gene.values() for r in rs]

        msg = ("Maximising isoform coverage across all genes…"
               if mode == "column" else
               "Reading the cDNA FASTA and maximising isoform coverage…\n"
               "(first run loads the FASTA; this can take a few seconds)")
        # Built via the parent-only overload rather than
        # QProgressDialog(msg, None, 0, 0, self): the 5-argument form types its
        # second parameter as the cancel-button *label*, so passing None to mean
        # "no cancel button" is not a valid call even though Qt tolerates it.
        # setCancelButton(None) below is the supported way to say that, and
        # min == max == 0 is what makes the bar indeterminate.
        progress = QProgressDialog(self)
        progress.setLabelText(msg)
        progress.setMinimum(0)
        progress.setMaximum(0)
        progress.setWindowTitle("Maximise isoform coverage")
        progress.setWindowModality(Qt.WindowModal)
        progress.setMinimumDuration(0)
        progress.setCancelButton(None)
        progress.show()

        # Parent the worker to the dialog so it's owned by the Qt object tree,
        # and free it once it finishes. Disable the trigger while it runs.
        worker = _CoverageWorker(
            rows=rows, n_per_gene=self._max_per_gene,
            min_gap=self._spacing.min_gap,
            # PREFERRED_PROBE_GAP_BP is the isoform-coverage selector's OWN
            # spreading target and has no RunConfig field (see models.RunConfig)
            # — it is genuinely a constant, unlike the three rungs above it.
            preferred_gap=PREFERRED_PROBE_GAP_BP,
            relaxed_gap=max(self._spacing.relaxed_gap, self._spacing_floor()),
            floor_gap=max(self._spacing.floor_gap, self._spacing_floor()),
            gene_isoforms=(payload if mode == "column" else None),
            fasta=(payload if mode == "fasta" else None),
            cached_index=self._gene_index,
            parent=self,
        )
        self._cov_worker = worker
        self.coverage_btn.setEnabled(False)

        def _done(pick, index):
            self._gene_index = index            # cache for repeat runs
            progress.close()
            self._apply_coverage_pick(pick)
            self._cov_worker = None
            self.coverage_btn.setEnabled(True)

        def _fail(msg):
            progress.close()
            QMessageBox.critical(self, "Coverage selection failed", msg)
            self._cov_worker = None
            self.coverage_btn.setEnabled(True)

        worker.done.connect(_done)
        worker.failed.connect(_fail)
        worker.finished.connect(worker.deleteLater)
        worker.start()

    def closeEvent(self, event) -> None:
        # Don't let Qt destroy a running QThread out from under us — that's the
        # classic "QThread: Destroyed while thread is still running" crash.
        worker = self._cov_worker
        if worker is not None and worker.isRunning():
            worker.wait()
        super().closeEvent(event)

    def _apply_coverage_pick(self, pick) -> None:
        # Replace selections with the coverage result (keyed by source
        # transcript, which matches our per-gene rows), capped to max-per-gene.
        valid_by_gene = {
            gid: {r.get("probe_id") for r in rs}
            for gid, rs in self._rows_by_gene.items()
        }
        self._selections = {}
        for gid, ids in pick.selections.items():
            keep = [p for p in ids if p in valid_by_gene.get(gid, set())]
            if keep:
                self._selections[gid] = keep[:self._max_per_gene]
        self._populate_table(); self._refresh_viz(); self._refresh_counter()

        s = pick.summary
        msg = (
            f"Selected {s['probes']:,} probes across {s['genes']:,} genes.\n\n"
            f"No flags (≤{self._max_per_gene} spaced probes reach every "
            f"isoform): {s['genes'] - s['flagged']:,}\n"
            f"Flagged (incomplete coverage, fewer than "
            f"{self._max_per_gene} spaceable probes, or spacing relaxed below "
            f"{self._spacing.min_gap} bp): {s['flagged']:,}"
        )
        msg += (
            f"\n\nProbe spacing:"
            f"\n   green  (>= {s.get('min_gap_bp', self._spacing.min_gap)} bp): "
            f"{s.get('probes_green', 0):,}"
            f"\n   orange ({s.get('relaxed_gap_bp', self._spacing.relaxed_gap)}-"
            f"{s.get('min_gap_bp', self._spacing.min_gap)} bp): "
            f"{s.get('probes_orange', 0):,}"
            f"\n   red    ({s.get('floor_gap_bp', self._spacing.floor_gap)}-"
            f"{s.get('relaxed_gap_bp', self._spacing.relaxed_gap)} bp): "
            f"{s.get('probes_red', 0):,}"
        )
        if s.get("probes_red"):
            msg += (
                f"\n\nThe {s.get('genes_red', 0):,} red gene(s) reached the last "
                f"rung of the spacing ladder: their transcripts cannot hold "
                f"{PROBE_GAP_RELAX_MIN_PROBES} separated probes. Consider fewer "
                f"probes for them, or set the Spacing box to exclude red."
            )
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Information)
        box.setWindowTitle("Isoform coverage selection")
        box.setText(msg)
        save_btn = None
        if pick.flags:
            save_btn = box.addButton("Save flagged genes…", QMessageBox.ActionRole)
        box.addButton(QMessageBox.Ok)
        box.exec()
        if save_btn is not None and box.clickedButton() is save_btn:
            self._save_flagged_report(pick)

    def _save_flagged_report(self, pick) -> None:
        import csv
        path = os.path.join(self._record.output_dir, "flagged_genes.csv")
        idx = self._gene_index
        try:
            with open(path, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["gene", "gene_symbol", "isoforms_total",
                            "probes_selected", "isoforms_covered",
                            "coverage_fraction", "coverage_basis",
                            "flags", "missing_isoforms",
                            "min_gap_bp", "gap_used_bp", "spacing_tier"])
                for f in sorted(pick.flags, key=lambda x: x.gene):
                    sym = idx.symbol_of.get(f.gene, "") if idx else ""
                    frac = f.isoforms_covered / max(f.isoforms_total, 1)
                    w.writerow([f.gene, sym, f.isoforms_total, f.probes_selected,
                                f.isoforms_covered, f"{frac:.4f}",
                                f.coverage_basis,
                                ";".join(f.reasons),
                                ";".join(f.missing_isoforms),
                                "" if f.min_gap_bp is None else f.min_gap_bp,
                                f.gap_used_bp, f.tier_name])
        except OSError as exc:
            QMessageBox.critical(self, "Could not write report", str(exc))
            return
        QMessageBox.information(self, "Flagged genes saved",
                                f"Wrote:\n{path}")

    def _blocked_by_spacing_filter(self, gid, pid) -> Optional[int]:
        """Gap to the nearest already-selected probe when the filter forbids it.

        Returns the offending gap, or None when the pick is allowed. Checked
        before the tick is honoured rather than warned about after: the filter's
        contract is that forbidden spacing never reaches the order CSV, and a
        warning the user can ignore is not that.
        """
        floor = self._spacing_floor()
        by_id = {r.get("probe_id"): r for r in self._rows_by_gene.get(gid, [])}
        row = by_id.get(pid)
        if row is None:
            return None
        try:
            iv = (int(row.get("start") or 0), int(row.get("end") or 0))
        except (TypeError, ValueError):
            return None
        worst: Optional[int] = None
        # Siblings of the same site are not neighbours: they share a footprint
        # by construction and read as -1 bp. See site_key.
        site = site_key(row)
        for other in self._selections.get(gid, []):
            if other == pid:
                continue
            o = by_id.get(other)
            if o is None or site_key(o) == site:
                continue
            try:
                oiv = (int(o.get("start") or 0), int(o.get("end") or 0))
            except (TypeError, ValueError):
                continue
            g = _interval_gap(iv, oiv)
            if g < floor and (worst is None or g < worst):
                worst = g
        return worst

    def _on_item_changed(self, item) -> None:
        if item.column() != SELECT_COL:
            return
        gid = self._current_gene()
        if gid is None:
            return
        pid = item.data(Qt.UserRole)
        current = self._selections.setdefault(gid, [])
        if item.checkState() == Qt.Checked:
            if pid in current:
                pass
            elif (self._site_by_id(gid).get(pid, pid)
                  not in self._sites(gid, current)
                  and len(self._sites(gid, current)) >= self._max_per_gene):
                self.table.blockSignals(True)
                item.setCheckState(Qt.Unchecked)
                self.table.blockSignals(False)
                QMessageBox.information(
                    self, "Limit reached",
                    f"You've already chosen {self._max_per_gene} probe sites for "
                    f"{gid}. Uncheck one or raise the per-gene limit.",
                )
                return
            elif (bad := self._blocked_by_spacing_filter(gid, pid)) is not None:
                self.table.blockSignals(True)
                item.setCheckState(Qt.Unchecked)
                self.table.blockSignals(False)
                floor = self._spacing_floor()
                QMessageBox.information(
                    self, "Too close under the current spacing filter",
                    f"This probe sits {bad} bp from one you have already chosen "
                    f"for {gid}, and the Spacing filter is set to allow no "
                    f"closer than {floor} bp.\n\n"
                    "Loosen the Spacing box if you want it anyway.",
                )
                return
            else:
                current.append(pid)
        elif pid in current:
            current.remove(pid)
        self._refresh_viz(); self._refresh_counter()

    def _on_save_draft(self) -> None:
        save_selections(self._record.output_dir, self._selections, note="draft")
        QMessageBox.information(self, "Draft saved",
                                f"Selections written to:\n{self._record.output_dir}/selections.json")

    def _on_confirm(self) -> None:
        sel = {gid: list(ids) for gid, ids in self._selections.items() if ids}
        if not sel:
            QMessageBox.warning(self, "Nothing selected",
                                "Select at least one probe before confirming.")
            return
        if self.map_widget.has_overlap:
            if QMessageBox.question(
                self, "Probes too close",
                f"Some selected probes are closer than "
                f"{self._spacing.relaxed_gap} bp "
                "on the gene track (red). Auto-selection only goes this tight "
                "on its last rung "
                f"({self._spacing.floor_gap} bp), for a gene too short to hold "
                "separated probes. Export anyway?",
                QMessageBox.Yes | QMessageBox.No,
            ) != QMessageBox.Yes:
                return
        try:
            save_selections(self._record.output_dir, sel, note="confirmed")
            written = export_selected_probes(self._record.output_dir, sel, overwrite=True)
        except Exception as exc:
            QMessageBox.critical(self, "Export failed", str(exc))
            return
        QMessageBox.information(self, "Done", "Order CSVs written:\n" + "\n".join(written.values()))
        self.confirmed.emit(self._record.run_id)
        self.accept()


# ---------------------------------------------------------------------------
# Cell-item helpers
# ---------------------------------------------------------------------------

def _num_item(value) -> QTableWidgetItem:
    text = "" if value is None else str(value)
    item = QTableWidgetItem()
    try:
        item.setData(Qt.DisplayRole, float(text))
    except (TypeError, ValueError):
        item.setText(text)
    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
    return item


def _gc_item(value) -> QTableWidgetItem:
    """GC cell. results.csv stores GC as a fraction (0–1); show it as a
    percentage (e.g. 0.44 → 44.0). Stored as a number so the column still
    sorts numerically; the '%' is conveyed by the 'GC% …' header."""
    item = QTableWidgetItem()
    try:
        item.setData(Qt.DisplayRole, round(float(value) * 100.0, 1))
    except (TypeError, ValueError):
        item.setText("" if value is None else str(value))
    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
    return item


def _tm_item(value) -> QTableWidgetItem:
    """Tm cell that renders ``n/a`` in muted italic (C5)."""
    text = "" if value is None else str(value)
    item = QTableWidgetItem()
    if text.strip().lower() in ("n/a", "nan", ""):
        item.setText("n/a")
        item.setForeground(QBrush(QColor("#999")))
        f = item.font(); f.setItalic(True); item.setFont(f)
        item.setData(Qt.UserRole, -1.0)   # sort n/a to the bottom
    else:
        try:
            item.setData(Qt.DisplayRole, float(text))
        except (TypeError, ValueError):
            item.setText(text)
    item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
    return item


def _spec_item(status: str) -> QTableWidgetItem:
    """Specificity-status cell; non-'passed' values are highlighted."""
    item = QTableWidgetItem(status or "")
    if status and status != "passed":
        item.setForeground(QBrush(QColor("#b00020")))
        if status in ("not_screened", "blast_failed"):
            item.setToolTip("This probe was NOT off-target screened — no "
                            "reference was given, or the screen failed.")
    return item


