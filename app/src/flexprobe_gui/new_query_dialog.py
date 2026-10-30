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

"""Modal dialog for configuring a probe-design query.

Layout (top to bottom):
  • Inputs: the target FASTA (or the bulk FASTA list) and the off-target
            references, i.e. what the run is designed from and screened against
  • Output: output base folder, overwrite, run name
  • Chemistry — modality + multiplex/singleplex toggle + Probe Barcode
                (Flex v1 only). Selecting a modality auto-fills everything else.
  • Homology mode — single sequence vs. a consensus over several related sequences
  • Workflow: run mode + project
  • Off-target screening: engine and screen options, on by default (k-mer
                          engine, BLAST optional)
  • Advanced options (foldable, collapsed by default):
        – Filters & tile  — pre-filled from chemistry; user-editable
        – Ranking weights — composite-score weight overrides

The actual handle sequences (LHS 5', RHS 3', Probe Barcode tables) live in
``flexprobe.chemistry`` and are looked up by the pipeline at run time.
"""
from __future__ import annotations

import logging
import os
from typing import Optional

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox, QComboBox, QDialog, QDialogButtonBox, QDoubleSpinBox,
    QFileDialog, QFormLayout, QFrame, QGroupBox, QHBoxLayout, QLabel,
    QLineEdit, QListWidget, QMessageBox, QPlainTextEdit,
    QPushButton, QScrollArea, QSizePolicy, QSpinBox, QToolButton,
    QVBoxLayout, QWidget,
)

from flexprobe.chemistry import FLEX_V1_PROBE_BARCODES, get_chemistry
from flexprobe.constants import (
    OFF_TARGET_JUNCTION_DUPLEX_MIN_BP, REFERENCE_STRAND_BOTH,
    REFERENCE_STRAND_PLUS, SCREEN_ENGINE_DEFAULT, TM_MIN_DEFAULT,
    TN_MODE_OFF, TN_MODE_SOFT, TN_MODE_STRICT_A,
)
from flexprobe.models import (
    FilterThresholds, RunConfig, ScreenResources, SpecificityThresholds,
)
from flexprobe.ranking import DEFAULT_WEIGHTS
from flexprobe.screen import ENGINE_BLAST, ENGINE_KMER

from .modality import DEFAULT_MODALITY, MODALITY_DEFAULTS, MODALITY_PRESETS
from .theme import link_hex, muted_hex
from .project_history import ProjectHistory
from .query_drafts import QueryDraftStore
from ._nq_draft import DraftMixin
from ._nq_build import (
    ALL_TEMPLATES, BuildMixin, OT_COUNT_FULL_HALF, OT_COUNT_LOCAL,
    OT_RULE_CROSS_SUBJECT, OT_RULE_EITHER_HALF, OT_RULE_SAME_SUBJECT,
    ot_rule_code,
)

#: A probe half is 25 nt, so a junction duplex needs at most 25 bp per side.
_JUNCTION_DUPLEX_MAX_BP = 25
#: The k-mer screen reports hits with at most 5 mismatches per half
#: (screen.KmerScreen.MAX_MISMATCH_REPORTED). A pair is rejected at
#: min_mismatch - 1 mismatches, so 6 is the highest threshold it can check
#: completely; above that it would silently miss qualifying sites.
_OT_MIN_MISMATCH_MAX = 6
from .settings import AppSettings

logger = logging.getLogger("flexprobe_gui.new_query_dialog")


def _wrap(layout) -> QWidget:
    w = QWidget()
    w.setLayout(layout)
    return w


class _CollapsibleSection(QWidget):
    """A foldable section: a clickable header with a disclosure arrow that
    shows / hides its content. Collapsed by default so the dialog opens compact.
    """

    def __init__(self, title: str, parent=None):
        super().__init__(parent)
        self.toggle = QToolButton()
        self.toggle.setText(title)
        self.toggle.setCheckable(True)
        self.toggle.setChecked(False)
        self.toggle.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.toggle.setArrowType(Qt.RightArrow)
        self.toggle.setStyleSheet(
            "QToolButton { border: none; font-weight: bold; padding: 4px 0; }"
        )
        self.toggle.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.toggle.toggled.connect(self._on_toggled)

        self._content = QWidget()
        self._content_layout = QVBoxLayout(self._content)
        self._content_layout.setContentsMargins(0, 0, 0, 0)
        self._content.setVisible(False)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        lay.addWidget(self.toggle)
        lay.addWidget(self._content)

    def _on_toggled(self, checked: bool) -> None:
        self.toggle.setArrowType(Qt.DownArrow if checked else Qt.RightArrow)
        self._content.setVisible(checked)

    def add_widget(self, w: QWidget) -> None:
        self._content_layout.addWidget(w)

    def set_expanded(self, expanded: bool) -> None:
        self.toggle.setChecked(bool(expanded))

    def is_expanded(self) -> bool:
        return self.toggle.isChecked()


class NewQueryDialog(DraftMixin, BuildMixin, QDialog):
    def __init__(
        self,
        parent=None,
        *,
        projects: Optional[ProjectHistory] = None,
        bulk_fasta_paths: Optional[list] = None,
        settings: Optional[AppSettings] = None,
        drafts: Optional[QueryDraftStore] = None,
    ):
        super().__init__(parent)
        # Optional ProjectHistory injected by the home window so projects
        # created in this dialog show up in the home page filter too. If
        # missing (e.g. tests), the dialog still works but the project
        # combo holds only "(none)" + "+ New project…".
        self._projects = projects or ProjectHistory()
        # App preferences (default output dir) + the saved-draft store. Both
        # fall back to on-disk defaults so the dialog is usable standalone.
        self._settings = settings or AppSettings()
        self._drafts = drafts or QueryDraftStore()
        # Bulk mode: every FASTA in this list becomes its own run sharing
        # the dialog's chemistry / filters / BLAST settings. Empty list =
        # ordinary single-FASTA mode.
        self._bulk_paths = list(bulk_fasta_paths or [])
        self._is_bulk = bool(self._bulk_paths)
        # Run name provenance (see DraftMixin._restore_run_name): True while the
        # box holds a name restored from a draft rather than one typed here,
        # and the FASTA it was restored against. A restored name is cleared
        # when the FASTA changes; a typed one is kept.
        self._run_name_restored = False
        self._run_name_fasta = ""
        # Record pickers (template / target) are filled from the FASTA, which a
        # draft sets AFTER the selection. A selection that finds an empty combo
        # waits here, keyed by widget attribute, until the combo is filled.
        self._pending_combo_data: dict = {}
        self.setWindowTitle(
            f"New bulk query — {len(self._bulk_paths)} FASTAs"
            if self._is_bulk else "New query"
        )
        self.setModal(True)
        # Resizable window with a sensible default and a true minimum so users
        # can shrink it below the Mac dock. Content scrolls if the user makes
        # the window short.
        # 680 px was narrower than the form itself (694 px), so the dialog
        # opened with a sideways scrollbar and clipped Browse buttons.
        self.resize(780, 700)
        self.setMinimumSize(460, 320)
        self.setSizeGripEnabled(True)

        # Outer = scroll area (content) + sticky button row at the bottom so
        # Start/Cancel are always reachable regardless of window height.
        page = QVBoxLayout(self)
        page.setContentsMargins(0, 0, 0, 0)
        page.setSpacing(0)

        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        page.addWidget(scroll, 1)

        content = QWidget()
        scroll.setWidget(content)
        outer = QVBoxLayout(content)
        outer.setContentsMargins(12, 12, 12, 12)

        # Inputs first: what the run is designed from and screened against is
        # the first thing a user has to supply, and everything below it has a
        # usable default. Built in display order; the only ordering constraint
        # is that the screening group precedes the ranking group (it reads
        # self.group_by_gene), and all of these precede the signal wiring.
        outer.addWidget(self._build_inputs_group())
        outer.addWidget(self._build_output_group())
        outer.addWidget(self._build_chemistry_group())
        outer.addWidget(self._build_strain_group())
        outer.addWidget(self._build_workflow_group())
        outer.addWidget(self._build_blast_group())

        # Advanced options — the fiddly knobs (filters/tile + ranking weights)
        # are pre-filled from chemistry and rarely need touching, so they live
        # in a foldable section that opens collapsed.
        self.advanced_section = _CollapsibleSection("Advanced options")
        self.advanced_section.add_widget(self._build_filters_group())
        self.advanced_section.add_widget(self._build_specificity_group())
        self.advanced_section.add_widget(self._build_spacing_group())
        self.advanced_section.add_widget(self._build_input_calibration_group())
        self.advanced_section.add_widget(self._build_ranking_group())
        outer.addWidget(self.advanced_section)
        outer.addStretch(1)

        # --- Buttons: pinned to the bottom of the dialog, OUTSIDE the scroll
        # area, so they remain reachable when the window is short (e.g. when
        # docked against the macOS dock).
        btn_strip = QWidget(self)
        btn_strip.setStyleSheet("background: palette(window);")
        btn_layout = QHBoxLayout(btn_strip)
        btn_layout.setContentsMargins(12, 8, 12, 12)
        # Draft save/load on the left so the user can stash a half-configured
        # query and pick it back up after restarting.
        self.save_draft_btn = QPushButton("Save draft…")
        self.save_draft_btn.setToolTip(
            "Save these settings under a name so you can reload and finish the "
            "query later."
        )
        self.save_draft_btn.clicked.connect(self._on_save_draft)
        btn_layout.addWidget(self.save_draft_btn)
        self.load_draft_btn = QPushButton("Load draft…")
        self.load_draft_btn.clicked.connect(self._on_load_draft)
        btn_layout.addWidget(self.load_draft_btn)
        # Export an HPC job bundle instead of running locally — for memory-heavy
        # whole-transcriptome runs that don't fit on a laptop.
        self.export_hpc_btn = QPushButton("Export for HPC…")
        self.export_hpc_btn.setToolTip(
            "Write a PBS/SLURM job bundle to run this design on a cluster, then "
            "import the results back into the app."
        )
        self.export_hpc_btn.clicked.connect(self._on_export_hpc)
        btn_layout.addWidget(self.export_hpc_btn)
        btn_layout.addStretch(1)
        btns = QDialogButtonBox(QDialogButtonBox.Cancel)
        start_label = (
            f"Start {len(self._bulk_paths)} runs" if self._is_bulk else "Start run"
        )
        self.start_btn = btns.addButton(start_label, QDialogButtonBox.AcceptRole)
        btns.accepted.connect(self._validate_and_accept)
        btns.rejected.connect(self.reject)
        btn_layout.addWidget(btns)
        page.addWidget(btn_strip, 0)

        # Wire signals AFTER all widgets exist so initial fan-out works.
        self.modality.currentTextChanged.connect(self._on_modality_changed)
        self.mode.currentTextChanged.connect(self._on_mode_changed)
        self.barcode_combo.currentTextChanged.connect(self._refresh_handle_preview)
        self.run_mode.currentTextChanged.connect(self._on_run_mode_changed)
        self.new_project_btn.clicked.connect(self._on_new_project)
        self.strain_mode_combo.currentTextChanged.connect(self._on_strain_mode_changed)
        self.strategy_combo.currentTextChanged.connect(self._on_strain_mode_changed)
        self.fasta_edit.textChanged.connect(self._on_fasta_path_changed)
        self.fasta_edit.textChanged.connect(self._on_fasta_changed_for_run_name)
        self.run_name.textEdited.connect(self._on_run_name_edited)
        self._on_modality_changed(self.modality.currentText())
        self._on_run_mode_changed(self.run_mode.currentText())
        self._on_strain_mode_changed("")

    def _build_chemistry_group(self):
        # --- Chemistry / modality (drives auto-populate + visibility) ---
        chem_box = QGroupBox("Chemistry")
        chem_outer = QVBoxLayout(chem_box)
        chem_l = QFormLayout()
        chem_l.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)
        chem_outer.addLayout(chem_l)

        self.modality = QComboBox()
        self.modality.addItems(list(MODALITY_PRESETS.keys()))
        self.modality.setCurrentText(DEFAULT_MODALITY)
        chem_l.addRow("Modality:", self.modality)

        self.mode = QComboBox()
        self.mode.addItems(["Multiplex", "Singleplex"])
        chem_l.addRow("Experiment mode:", self.mode)

        self.barcode_label = QLabel("Probe Barcode:")
        self.barcode_combo = QComboBox()
        # Populated dynamically per modality
        chem_l.addRow(self.barcode_label, self.barcode_combo)

        # Multiplex sample count — visible only for Flex v1 multiplex. The user
        # states how many samples they are multiplexing; the run uses Probe
        # Barcodes BC001..BC00N (N samples). results.csv keeps a single generic
        # {BC}-placeholder probe per candidate, and every selected probe is
        # repeated once per barcode at export time (one pool per sample, named
        # after its barcode).
        self.multi_bc_label = QLabel("Number of multiplex samples:")
        multi_bc_wrap = QHBoxLayout(); multi_bc_wrap.setContentsMargins(0, 0, 0, 0)
        self.multi_count = QSpinBox()
        self.multi_count.setRange(2, len(FLEX_V1_PROBE_BARCODES))
        self.multi_count.setValue(4)
        self.multi_count.setToolTip(
            "How many samples you are multiplexing. The run uses Probe Barcodes "
            "BC001 upwards, one per sample: 4 samples → BC001-BC004, 8 → BC001-BC008, "
            "16 → BC001-BC016."
        )
        multi_bc_wrap.addWidget(self.multi_count)
        self.multi_count_hint = QLabel()
        self.multi_count_hint.setStyleSheet(f"color: {muted_hex()};")
        multi_bc_wrap.addWidget(self.multi_count_hint)
        multi_bc_wrap.addStretch(1)
        self.multi_bc_widget = _wrap(multi_bc_wrap)
        chem_l.addRow(self.multi_bc_label, self.multi_bc_widget)
        self.multi_count.valueChanged.connect(self._on_multi_count_changed)
        self._on_multi_count_changed(self.multi_count.value())
        # Hidden by default — visibility is set per modality + mode.
        self.multi_bc_label.setVisible(False)
        self.multi_bc_widget.setVisible(False)

        # Read-only handle preview lives OUTSIDE the form layout so the
        # monospace string can use the full width of the group box and word-
        # wrap cleanly instead of being clipped to the form's right column.
        handle_caption = QLabel("Handle preview:")
        handle_caption.setStyleSheet("margin-top: 6px;")
        chem_outer.addWidget(handle_caption)

        self.handle_preview = QLabel()
        self.handle_preview.setWordWrap(True)
        self.handle_preview.setTextFormat(Qt.RichText)
        self.handle_preview.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.handle_preview.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
        self.handle_preview.setStyleSheet(
            "background: rgba(0,0,0,0.04); padding: 8px; border-radius: 4px;"
        )
        chem_outer.addWidget(self.handle_preview)

        return chem_box

    def _build_strain_group(self):
        # --- Homology mode: build a consensus across several related sequences? ---
        strain_box = QGroupBox("Homology mode")
        strain_l = QFormLayout(strain_box)
        self.strain_mode_combo = QComboBox()
        self.strain_mode_combo.addItems([
            "Single sequence per gene (default)",
            "Homology consensus (one FASTA = N related sequences of one fragment)",
        ])
        strain_l.addRow("Mode:", self.strain_mode_combo)

        self.strategy_combo = QComboBox()
        self.strategy_combo.addItems([
            "Conserved — probes hit every sequence",
            "Sequence-unique — probes distinguish one sequence",
            "Constitutive (isoform-collapse) — probes shared by ≥N records",
        ])
        self.strategy_label = QLabel("Strategy:")
        strain_l.addRow(self.strategy_label, self.strategy_combo)

        # Which record supplies the base at every variant position of the
        # designed probe. A single run only ever emits that record's flavour of
        # each window, so "every sequence" fans the run out and the union is the
        # complete pan-strain set. ALL_TEMPLATES is the sentinel for the fan-out.
        self.template_combo = QComboBox()
        self.template_combo.setEditable(False)
        self.template_label = QLabel("Consensus template:")
        strain_l.addRow(self.template_label, self.template_combo)

        # Constitutive / isoform-collapse: minimum records a window must be
        # shared by (0 = all). Shown only for the constitutive strategy.
        isoform_row = QHBoxLayout()
        isoform_row.setContentsMargins(0, 0, 0, 0)
        self.isoform_min = QSpinBox()
        self.isoform_min.setRange(0, 999)
        self.isoform_min.setValue(0)
        isoform_row.addWidget(self.isoform_min)
        isoform_row.addWidget(QLabel("isoforms (0 = all)"))
        isoform_row.addStretch(1)
        self.isoform_min_label = QLabel("Min isoforms shared:")
        # The whole row is hidden with its label, unit text included. Hiding only
        # the spin box left "isoforms (0 = all)" floating in the form.
        self.isoform_min_row = _wrap(isoform_row)
        strain_l.addRow(self.isoform_min_label, self.isoform_min_row)

        # Conserved strategy: how much disagreement a pan-record window may
        # carry. 0 = identical everywhere, the pre-0.5.0 rule. Opens at 1 to
        # match RunConfig.conserved_max_degenerate_per_half, which is what an
        # unflagged CLI run gets; the dialog used to open at 0, so the same
        # untouched form produced exact-identity windows in the app and
        # degenerate ones on the command line.
        degen_row = QHBoxLayout()
        degen_row.setContentsMargins(0, 0, 0, 0)
        self.conserved_degen = QSpinBox()
        self.conserved_degen.setRange(0, 5)
        self.conserved_degen.setValue(1)
        degen_row.addWidget(self.conserved_degen)
        degen_row.addWidget(QLabel("per 25 nt half (0 = identical everywhere)"))
        degen_row.addStretch(1)
        self.conserved_degen_label = QLabel("Variant positions allowed:")
        self.conserved_degen.setToolTip(
            "0: a window is used only if every record carries exactly the "
            "same 50 bases.\n\n"
            "1 (default, and what an unflagged command-line run uses) or more: "
            "a window may differ at up to this many positions in each "
            "25 nt half, provided none falls within the junction clearance — "
            "CG000621 Rev D asks for mismatches at least 4 bp from the ligation "
            "junction, and the ligase will not seal across one. Such a window is "
            "designed as a degenerate consensus probe, which is never ordered, "
            "plus one concrete probe per input record, which is.\n\n"
            "Exact identity is severe on real strains: two strains of one "
            "virus can leave few or no usable windows on a segment at 0, and "
            "many times more at 1 per half.\n\n"
            "Windows spanning an insertion or deletion are excluded at every "
            "setting: one record has a base the other does not, so a probe there "
            "bulges rather than mispairs and the downstream half is "
            "frame-shifted."
        )
        self.conserved_degen_row = _wrap(degen_row)
        strain_l.addRow(self.conserved_degen_label, self.conserved_degen_row)

        self.unique_target_combo = QComboBox()
        self.unique_target_combo.setEditable(False)
        self.unique_target_label = QLabel("Target sequence:")
        strain_l.addRow(self.unique_target_label, self.unique_target_combo)

        unique_minmm_row = QHBoxLayout()
        unique_minmm_row.setContentsMargins(0, 0, 0, 0)
        self.unique_min_mm = QSpinBox()
        self.unique_min_mm.setRange(1, 25)
        self.unique_min_mm.setValue(5)
        unique_minmm_row.addWidget(self.unique_min_mm)
        unique_minmm_row.addWidget(QLabel(
            "bp (≥5 matches the 10x off-target threshold)"
        ))
        unique_minmm_row.addStretch(1)
        self.unique_min_mm_label = QLabel("Min mismatches vs others:")
        self.unique_min_mm_row = _wrap(unique_minmm_row)
        strain_l.addRow(self.unique_min_mm_label, self.unique_min_mm_row)

        # Strain-unique junction duplex: the strain-unique twin of the
        # off-target junction-duplex gate. Same 0-or-6..25 bp rule.
        _sujd = RunConfig(fasta_path="", output_dir="").strain_unique_junction_duplex_bp
        ujd_row = QHBoxLayout(); ujd_row.setContentsMargins(0, 0, 0, 0)
        self.unique_jd = QCheckBox("Reject pairs matching another sequence "
                                   "across the junction by")
        self.unique_jd.setChecked(_sujd > 0)
        self.unique_jd_bp = QSpinBox()
        self.unique_jd_bp.setRange(OFF_TARGET_JUNCTION_DUPLEX_MIN_BP,
                                   _JUNCTION_DUPLEX_MAX_BP)
        self.unique_jd_bp.setValue(_sujd if _sujd > 0
                                   else OFF_TARGET_JUNCTION_DUPLEX_MIN_BP)
        self.unique_jd_bp.setSuffix(" bp per side")
        self.unique_jd.toggled.connect(self.unique_jd_bp.setEnabled)
        self.unique_jd.setToolTip(
            "A discriminating probe must not ligate on the OTHER sequences "
            "either. This rejects a pair that matches another input sequence "
            "perfectly across the ligation nick for at least this many bases "
            f"on each side (default {_sujd}). "
            f"{OFF_TARGET_JUNCTION_DUPLEX_MIN_BP} bp is the shortest the scan "
            "can see.")
        self.unique_jd_bp.setToolTip(self.unique_jd.toolTip())
        ujd_row.addWidget(self.unique_jd)
        ujd_row.addWidget(self.unique_jd_bp)
        ujd_row.addStretch(1)
        self.unique_jd_label = QLabel("Junction duplex vs others:")
        self.unique_jd_row = _wrap(ujd_row)
        strain_l.addRow(self.unique_jd_label, self.unique_jd_row)

        self.unique_per_half = QCheckBox(
            "Require the mismatches within a single probe half"
        )
        self.unique_per_half.setChecked(True)
        self.unique_per_half.setToolTip(
            "On (recommended): the target must differ from every other record "
            "by the threshold above within its LHS half or its RHS half. This "
            "is the unit 10x's ≥5 mismatch rule is written in and the unit "
            "the off-target screen uses.\n\n"
            "Off: the mismatches are counted across the whole 50 bp window. "
            "Mismatches split across the junction (e.g. 3 LHS + 2 RHS) then "
            "count as 5, but leave both halves able to bind the non-target "
            "record — so the probe still hybridises and ligates there and "
            "does not actually discriminate."
        )
        strain_l.addRow("", self.unique_per_half)

        self.strain_hint = QLabel()
        self.strain_hint.setStyleSheet(f"color: {muted_hex()};")
        self.strain_hint.setWordWrap(True)
        strain_l.addRow(self.strain_hint)
        return strain_box

    def _build_workflow_group(self):
        # --- Workflow: pick how the run handles selection + which project ---
        wf_box = QGroupBox("Workflow")
        wf_l = QFormLayout(wf_box)
        self.run_mode = QComboBox()
        self.run_mode.addItems([
            "Manual review (default)",
            "Auto-pick top N per gene (whole-transcriptome)",
        ])
        wf_l.addRow("Run mode:", self.run_mode)
        self.auto_n = QSpinBox()
        self.auto_n.setRange(1, 50)
        self.auto_n.setValue(3)
        self.auto_n_label = QLabel("Probes per gene:")
        wf_l.addRow(self.auto_n_label, self.auto_n)

        self.project_combo = QComboBox()
        self.new_project_btn = QPushButton("+ New…")
        proj_row = QHBoxLayout(); proj_row.setContentsMargins(0, 0, 0, 0)
        proj_row.addWidget(self.project_combo, 1)
        proj_row.addWidget(self.new_project_btn)
        wf_l.addRow("Project:", _wrap(proj_row))
        self._refresh_project_combo()
        return wf_box

    def _build_inputs_group(self):
        # --- Inputs: target FASTA(s) + off-target references ---
        io_box = QGroupBox("Inputs")
        io_l = QFormLayout(io_box)
        io_l.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)
        self.fasta_edit = QLineEdit()
        fh = QHBoxLayout(); fh.setContentsMargins(0, 0, 0, 0)
        fh.addWidget(self.fasta_edit)
        b1 = QPushButton("Browse…"); b1.clicked.connect(self._browse_fasta); fh.addWidget(b1)
        self.fasta_row_label = QLabel("Target FASTA:")
        self._fasta_row_widget = _wrap(fh)
        io_l.addRow(self.fasta_row_label, self._fasta_row_widget)

        # Bulk mode shows the dropped/picked files as a read-only list and
        # hides the single-FASTA picker.
        if self._is_bulk:
            self.bulk_list = QListWidget()
            self.bulk_list.setFixedHeight(140)
            for p in self._bulk_paths:
                self.bulk_list.addItem(os.path.basename(p))
            self.bulk_list.setToolTip(
                "\n".join(self._bulk_paths)
            )
            io_l.addRow(QLabel(f"FASTAs ({len(self._bulk_paths)}):"), self.bulk_list)
            self.fasta_row_label.setVisible(False)
            self._fasta_row_widget.setVisible(False)
            # Pre-fill fasta_edit so any latent dialog code that reads it
            # doesn't choke (it's hidden so the user can't change it).
            self.fasta_edit.setText(self._bulk_paths[0])

        # Each line is either an existing BLAST DB prefix OR a raw FASTA path.
        # FASTAs are auto-indexed (makeblastdb -parse_seqids) right before the
        # run, cached by SHA-256 next to the FASTA, so re-runs are instant.
        # blastn searches across every line at once. The accession index is not
        # optional: the off-target recount retrieves subject bases through
        # `blastdbcmd -entry`, and a DB without one is rejected up front.
        self.db_edit = QPlainTextEdit()
        self.db_edit.setPlaceholderText(
            "One reference FASTA per line (the k-mer engine indexes it and "
            "caches the index beside it). Under the BLAST engine a line may "
            "instead be a BLAST DB prefix:\n"
            "/Users/me/blast_dbs/human_tx\n"
            "/Users/me/refs/other_viruses.fasta"
        )
        self.db_edit.setFixedHeight(72)
        db_btn_row = QHBoxLayout(); db_btn_row.setContentsMargins(0, 0, 0, 0)
        add_fasta_btn = QPushButton("Add FASTA (auto-index)…")
        add_fasta_btn.setToolTip(
            "Pick a FASTA of off-target sequences. The default k-mer engine "
            "searches it directly and caches its index beside it after the "
            "first run. Under the BLAST engine a database is built from it "
            "automatically (with -parse_seqids), so you never need to run "
            "makeblastdb yourself."
        )
        add_fasta_btn.clicked.connect(self._browse_blast_fasta)
        db_btn_row.addWidget(add_fasta_btn)
        add_db_btn = QPushButton("Add existing DB…")
        add_db_btn.setToolTip(
            "BLAST engine only. Add a BLAST database prefix. It must have been "
            "built with -parse_seqids, or the run stops with instructions.\n\n"
            "The default k-mer engine needs the FASTA the database was built "
            "from. It finds that FASTA automatically when it sits beside the "
            "database with a .fa or .fasta extension, and stops otherwise.")
        add_db_btn.clicked.connect(self._browse_db)
        db_btn_row.addWidget(add_db_btn)
        db_btn_row.addStretch(1)
        io_l.addRow("Off-target references:", self.db_edit)
        io_l.addRow("", _wrap(db_btn_row))
        return io_box

    def _build_output_group(self):
        # --- Output folder + run identification ---
        out_box = QGroupBox("Output")
        out_l = QFormLayout(out_box)
        self.out_edit = QLineEdit(self._settings.default_output_dir())
        oh = QHBoxLayout(); oh.setContentsMargins(0, 0, 0, 0)
        oh.addWidget(self.out_edit)
        b2 = QPushButton("Browse…"); b2.clicked.connect(self._browse_out); oh.addWidget(b2)
        out_l.addRow("Output base dir:", _wrap(oh))
        self.overwrite = QCheckBox("Overwrite existing output files")
        self.overwrite.setChecked(True)
        out_l.addRow("", self.overwrite)

        self.run_name = QLineEdit()
        self.run_name.setPlaceholderText("(auto from FASTA filename)")
        self.run_name_label = QLabel("Run name:")
        out_l.addRow(self.run_name_label, self.run_name)
        # In bulk mode the run name is auto-derived per-file; show a hint
        # instead of an editable field.
        if self._is_bulk:
            self.run_name.setVisible(False)
            self.run_name_label.setVisible(False)
            bulk_hint = QLabel(
                f"Bulk mode: {len(self._bulk_paths)} FASTAs → {len(self._bulk_paths)} runs. "
                f"Run names are taken from each FASTA filename."
            )
            bulk_hint.setStyleSheet(f"color: {muted_hex()};")
            bulk_hint.setWordWrap(True)
            out_l.addRow(bulk_hint)
        return out_box

    def _build_filters_group(self):
        # --- Filters & tile ---
        self.filt_box = QGroupBox("Filters and tiling (auto-filled by chemistry)")
        filt_l = QFormLayout(self.filt_box)
        self.tile_size = QSpinBox(); self.tile_size.setRange(15, 200)
        self.tile_size_label = QLabel("Tile size (bp):")
        filt_l.addRow(self.tile_size_label, self.tile_size)
        self.tile_step = QSpinBox(); self.tile_step.setRange(1, 100)
        self.tile_step_label = QLabel("Tile step (bp):")
        filt_l.addRow(self.tile_step_label, self.tile_step)
        self.gc_min = QDoubleSpinBox(); self.gc_min.setRange(0.0, 1.0)
        self.gc_min.setSingleStep(0.01); self.gc_min.setDecimals(2)
        self.gc_min_label = QLabel("GC min (per half):")
        filt_l.addRow(self.gc_min_label, self.gc_min)
        self.gc_max = QDoubleSpinBox(); self.gc_max.setRange(0.0, 1.0)
        self.gc_max.setSingleStep(0.01); self.gc_max.setDecimals(2)
        self.gc_max_label = QLabel("GC max (per half):")
        filt_l.addRow(self.gc_max_label, self.gc_max)
        self.enforce_gc = QCheckBox("Reject probes outside the GC range")
        self.enforce_gc.setToolTip(
            "Off by default. 10x publishes GC 44-72% per half, but 8.4% of the "
            "halves in their own human panel fall outside it (down to 0.28) "
            "while none falls below 56.3 C — so the Tm floor below is the real "
            "duplex-stability gate and GC only contributes to ranking. Turn "
            "this on to enforce the published GC range as a hard filter."
        )
        filt_l.addRow("", self.enforce_gc)
        self.tm_min = QDoubleSpinBox(); self.tm_min.setRange(0.0, 100.0)
        self.tm_min.setSingleStep(0.5); self.tm_min.setDecimals(1)
        self.tm_min_label = QLabel("Tm min (per half, C):")
        self.tm_min.setToolTip(
            "Per-half melting-temperature floor. RELATIVE, not an in-assay Tm: "
            "primer3 under PCR-primer salt conditions, because 10x does not "
            "publish the Flex hybridisation buffer. The default is derived "
            "under those same conditions from the minimum per-half Tm of 10x's "
            "shipped panels (56.32 C human / 56.35 C mouse). 0 disables it."
        )
        filt_l.addRow(self.tm_min_label, self.tm_min)
        _tm_max_default = FilterThresholds().tm_max
        self.tm_max = QDoubleSpinBox(); self.tm_max.setRange(0.0, 100.0)
        self.tm_max.setSingleStep(0.5); self.tm_max.setDecimals(1)
        self.tm_max.setValue(0.0 if _tm_max_default is None
                             else float(_tm_max_default))
        self.tm_max.setToolTip(
            "Per-half melting-temperature ceiling, same relative primer3 scale "
            "as Tm min. 0 = no ceiling (the default: 10x's shipped panels show "
            "no upper bound to calibrate against). Must be at least Tm min.")
        # 0 is "no ceiling"; say so rather than show a 0 C ceiling.
        self.tm_max.setSpecialValueText("off")
        filt_l.addRow("Tm max (per half, C):", self.tm_max)
        self.max_n = QSpinBox(); self.max_n.setRange(0, 25)
        self.max_n_label = QLabel("Max degenerate bases per half:")
        self.max_n.setToolTip(
            "How many ambiguous positions a 25 nt half may carry. Counts EVERY "
            "non-ACGT base — 'N' and every IUPAC code (M R W S Y K V H D B) — "
            "not just a literal N. The setting is named max_n_per_half in "
            "config.json and on the --max-n-per-half flag.\n\n"
            "This matters in homology mode: the consensus builder emits a "
            "two-way code wherever the records disagree but all agree on a "
            "base, and falls back to 'N' only for a gap. At 0 every such "
            "window is discarded. Default 1 (i.e. up to 2 per 50-mer), which "
            "is what an unflagged command-line run uses.\n\n"
            "Separately, and at every setting, no degenerate base is allowed "
            "within the junction clearance — the ligase will not seal across "
            "an ambiguous position beside the nick."
        )
        filt_l.addRow(self.max_n_label, self.max_n)
        # Up to 50: a value at or above the tile length switches the filter
        # off, which some paper runs used (max_homopolymer_run=50).
        self.max_hp = QSpinBox(); self.max_hp.setRange(1, 50)
        filt_l.addRow("Max homopolymer run:", self.max_hp)
        # Three-way, not a tick box. The engine has had strict_a / soft / off
        # since the mode was added, but the dialog only had a boolean — so the
        # soft gate, which is the whole point of CG000839's relaxation, was
        # reachable from the CLI and from a hand-written config and from
        # nowhere in the app. The filter report even explained a mode the app
        # could not produce.
        self.tn_junction_mode = QComboBox()
        self.tn_junction_mode.addItem(
            "Hard — reject unless A (default)", TN_MODE_STRICT_A)
        self.tn_junction_mode.addItem(
            "Soft — reject only C (Flex v1 / v2 only)", TN_MODE_SOFT)
        self.tn_junction_mode.addItem(
            "Off — no junction-identity check", TN_MODE_OFF)
        self.tn_junction_mode.setToolTip(
            "The base in the target mRNA opposite the LHS probe's 3' "
            "(ligation) end — window position 26, index 25. When it is A the "
            "LHS oligo ends in T, which is the TN motif.\n\n"
            "HARD (default): reject anything but A. CG000621 Rev D states it "
            "as a requirement, and CG000839 Rev A, having relaxed it, still "
            "recommends it \"for a simplified and robust design process\".\n\n"
            "SOFT: reject only C. CG000839 says the TN motif \"is not always "
            "strictly required, and other nucleotide motifs can also function "
            "effectively\"; 10x's own panels put A there only ~51% of the time "
            "while avoiding C ~99% of the time. Under soft, A becomes a "
            "ranking preference rather than a gate, so A-junction probes still "
            "sort to the top. Flex only — CG000621's Visium section states the "
            "requirement flatly, so Visium refuses this and the option is "
            "disabled there.\n\n"
            "OFF: no identity check at all. junction_base and tn_junction_ok are still "
            "recorded on every candidate whichever you pick."
        )
        self.tn_junction_label = QLabel("Ligation junction (TN):")
        filt_l.addRow(self.tn_junction_label, self.tn_junction_mode)
        # Explains, on the form, why "Soft" just greyed out or reverted.
        self.tn_junction_note = QLabel()
        self.tn_junction_note.setWordWrap(True)
        self.tn_junction_note.setStyleSheet(f"color: {muted_hex()};")
        self.tn_junction_note.setVisible(False)
        filt_l.addRow("", self.tn_junction_note)
        self.primer3 = QCheckBox("Enable primer3 Tm/hairpin scoring")
        filt_l.addRow("", self.primer3)

        # Parallelism for multi-gene (whole-transcriptome) runs.
        self.workers = QSpinBox()
        self.workers.setRange(1, 64)
        self.workers.setValue(1)
        self.workers.setToolTip(
            "Worker processes for TILING + FILTERING on multi-gene inputs. "
            "1 = serial. Has no effect on single-gene or homology-mode runs.\n\n"
            "This is not the screen. On a whole transcriptome the off-target "
            "screen is where nearly all the time goes — set 'Screen worker "
            "processes' in the Off-target screening section for that."
        )
        filt_l.addRow("Tiling worker processes:", self.workers)
        return self.filt_box

    def _build_specificity_group(self):
        """Off-target specificity rules and reference handling.

        Exposes the same specificity models as the command line
        (cross-subject pooling, either-half rejection, the junction-duplex
        threshold), so any run the CLI can make can be rebuilt here. Every
        default below is read from the engine's own dataclasses.
        """
        sp = SpecificityThresholds()
        sr = ScreenResources()
        box = QGroupBox("Off-target specificity rules")
        form = QFormLayout(box)

        self.ot_rule = QComboBox()
        self.ot_rule.addItem("Both halves hit the same off-target gene "
                             "(default, 10x rule)", OT_RULE_SAME_SUBJECT)
        self.ot_rule.addItem("Both halves hit any off-target genes "
                             "(cross-subject)", OT_RULE_CROSS_SUBJECT)
        self.ot_rule.addItem("Either half hits an off-target (strictest)",
                             OT_RULE_EITHER_HALF)
        self.ot_rule.setCurrentIndex(self.ot_rule.findData(ot_rule_code(
            sp.off_target_require_both_halves, sp.off_target_same_subject)))
        self.ot_rule.setToolTip(
            "When a probe pair counts as non-specific.\n\n"
            "Same gene (default): rejected only when BOTH halves match the SAME "
            "off-target gene too closely, which is what 10x's rule says. "
            "Cross-subject: both halves must fail, but they may fail on "
            "different genes, which rejects more. Either half: one failing "
            "half is enough; the strictest, and it ignores that both halves "
            "must bind for ligation.")
        form.addRow("A pair is rejected when:", self.ot_rule)

        self.ot_min_mm = QSpinBox()
        self.ot_min_mm.setRange(1, _OT_MIN_MISMATCH_MAX)
        self.ot_min_mm.setValue(int(sp.off_target_min_mismatch))
        self.ot_min_mm.setToolTip(
            f"An off-target match counts only if it has FEWER than this many "
            f"mismatches in a 25 nt half (default "
            f"{sp.off_target_min_mismatch}, the 10x rule). Capped at "
            f"{_OT_MIN_MISMATCH_MAX}: the k-mer screen reports matches with up "
            "to 5 mismatches, so a higher value could not be checked "
            "completely.")
        form.addRow("Off-target mismatches needed per half:", self.ot_min_mm)

        self.ot_counting = QComboBox()
        self.ot_counting.addItem("Across the full 25 nt half (default)",
                                 OT_COUNT_FULL_HALF)
        self.ot_counting.addItem("Local alignment, >= 80% over >= 20 bp "
                                 "(legacy, BLAST engine only)", OT_COUNT_LOCAL)
        self.ot_counting.setCurrentIndex(self.ot_counting.findData(
            OT_COUNT_FULL_HALF if sp.off_target_full_half else OT_COUNT_LOCAL))
        self.ot_counting.setToolTip(
            "How mismatches are counted. The legacy local-alignment rule is "
            "offered only under BLAST: the k-mer engine always compares the "
            "full 25 nt half, so under it the local rule would not be a local "
            "alignment at all.")
        form.addRow("Mismatch counting:", self.ot_counting)

        dup_row = QHBoxLayout(); dup_row.setContentsMargins(0, 0, 0, 0)
        self.junction_duplex = QCheckBox("Reject junction-spanning duplexes of")
        self.junction_duplex.setChecked(sp.off_target_junction_duplex_bp > 0)
        self.junction_duplex_bp = QSpinBox()
        self.junction_duplex_bp.setRange(OFF_TARGET_JUNCTION_DUPLEX_MIN_BP,
                                         _JUNCTION_DUPLEX_MAX_BP)
        self.junction_duplex_bp.setValue(
            int(sp.off_target_junction_duplex_bp)
            if sp.off_target_junction_duplex_bp > 0
            else OFF_TARGET_JUNCTION_DUPLEX_MIN_BP)
        self.junction_duplex_bp.setSuffix(" bp per side")
        self.junction_duplex.toggled.connect(self.junction_duplex_bp.setEnabled)
        dup_row.addWidget(self.junction_duplex)
        dup_row.addWidget(self.junction_duplex_bp)
        dup_row.addStretch(1)
        self.junction_duplex.setToolTip(
            "Rejects a pair that matches an off-target gene perfectly across "
            "the ligation nick for at least this many bases on each side, "
            f"even when its mismatch count passes (default "
            f"{sp.off_target_junction_duplex_bp} bp). Our criterion, not "
            f"10x's; it only ever rejects more. {OFF_TARGET_JUNCTION_DUPLEX_MIN_BP} "
            "bp is the shortest the scan can see. k-mer engine only: under "
            "BLAST it is switched off.")
        self.junction_duplex_bp.setToolTip(self.junction_duplex.toolTip())
        form.addRow("Junction duplex:", _wrap(dup_row))

        self.exclude_alt = QCheckBox(
            "Drop ALT haplotype and PATCH scaffold records from the reference")
        self.exclude_alt.setChecked(bool(sr.exclude_alt_scaffolds))
        self.exclude_alt.setToolTip(
            "On by default. An ALT or PATCH record is a near-copy of a primary "
            "region, so keeping it rejects good probes for matching a second "
            "copy of their own gene. No effect on a reference without them.")
        form.addRow("Reference records:", self.exclude_alt)

        self.ref_strand = QComboBox()
        self.ref_strand.addItem("Both strands (default)", REFERENCE_STRAND_BOTH)
        self.ref_strand.addItem("Plus strand only", REFERENCE_STRAND_PLUS)
        self.ref_strand.setCurrentIndex(
            self.ref_strand.findData(sr.reference_strand))
        self.ref_strand.setToolTip(
            "Which reference strands the off-target screen searches. Both is "
            "the safe default; plus-only misses matches on the reverse strand.")
        form.addRow("Reference strands:", self.ref_strand)

        self.specificity_note = QLabel()
        self.specificity_note.setWordWrap(True)
        self.specificity_note.setStyleSheet(f"color: {muted_hex()};")
        form.addRow("", self.specificity_note)

        self.screen_engine.currentIndexChanged.connect(
            lambda _i: self._refresh_specificity_availability())
        self._refresh_specificity_availability()
        return box

    def _build_spacing_group(self):
        """The spacing ladder auto-pick and the review colours use.

        Recorded in run_config.json, which is where the review dialog and the
        auto-export read a run's own spacing from, so these apply to both.
        """
        r = RunConfig(fasta_path="", output_dir="")
        box = QGroupBox("Probe spacing")
        form = QFormLayout(box)
        self.min_gap = QSpinBox(); self.min_gap.setRange(1, 5000)
        self.min_gap.setSuffix(" bp"); self.min_gap.setValue(r.min_probe_gap_bp)
        self.relaxed_gap = QSpinBox(); self.relaxed_gap.setRange(1, 5000)
        self.relaxed_gap.setSuffix(" bp")
        self.relaxed_gap.setValue(r.relaxed_probe_gap_bp)
        self.floor_gap = QSpinBox(); self.floor_gap.setRange(1, 5000)
        self.floor_gap.setSuffix(" bp"); self.floor_gap.setValue(r.floor_probe_gap_bp)
        # floor <= relaxed <= target, kept by the boxes themselves.
        self.relaxed_gap.setMaximum(self.min_gap.value())
        self.floor_gap.setMaximum(self.relaxed_gap.value())
        self.min_gap.valueChanged.connect(self.relaxed_gap.setMaximum)
        self.relaxed_gap.valueChanged.connect(self.floor_gap.setMaximum)
        self.min_gap.setToolTip(
            f"Target distance between chosen probes (default "
            f"{r.min_probe_gap_bp} bp, green). A gene that cannot fill "
            "min(probes per gene, 3) at this distance steps down the ladder.")
        self.relaxed_gap.setToolTip(
            f"Second rung (default {r.relaxed_probe_gap_bp} bp, orange). Setting "
            "it equal to the target removes this rung only: a starved gene "
            "then goes straight to the floor.")
        self.floor_gap.setToolTip(
            f"Last rung (default {r.floor_probe_gap_bp} bp, red). No chosen "
            "pair is closer than this. Set it equal to the second rung to "
            "remove the red rung.")
        form.addRow("Target spacing:", self.min_gap)
        form.addRow("Second rung:", self.relaxed_gap)
        form.addRow("Floor:", self.floor_gap)
        return box

    def _build_input_calibration_group(self):
        """FASTA header parsing and the null-model calibration."""
        r = RunConfig(fasta_path="", output_dir="")
        box = QGroupBox("Input parsing and score calibration")
        form = QFormLayout(box)
        self.split_header = QCheckBox(
            "Use only the text before the first '|' as the record id")
        self.split_header.setChecked(bool(r.split_header_on_pipe))
        self.split_header.setToolTip(
            "For NCBI-style headers such as '>gi|123|ref|NM_000001.1|'. Records "
            "sharing that first field then count as ONE record, so turn this on "
            "only when the first field really is a unique id. Off by default.")
        form.addRow("FASTA headers:", self.split_header)
        self.null_shuffles = QSpinBox(); self.null_shuffles.setRange(0, 1000)
        self.null_shuffles.setValue(r.null_model_shuffles)
        self.null_shuffles.setToolTip(
            "Shuffled decoys per sampled window for the null-model calibration "
            "of the specificity score. 0 = off (default). Each shuffle is an "
            "extra screen query, and it is skipped when the screen does not run.")
        self.null_windows = QSpinBox(); self.null_windows.setRange(1, 100000)
        self.null_windows.setValue(r.null_model_max_windows)
        self.null_windows.setToolTip(
            "Windows sampled at random for the calibration (default "
            f"{r.null_model_max_windows}). Cost is at most windows x shuffles "
            "extra screen queries, whatever the run's size.")
        self.null_seed = QSpinBox(); self.null_seed.setRange(0, 2**31 - 1)
        self.null_seed.setValue(r.null_model_seed)
        self.null_seed.setToolTip("Random seed for the calibration sample.")
        self.null_windows.setEnabled(self.null_shuffles.value() > 0)
        self.null_seed.setEnabled(self.null_shuffles.value() > 0)
        self.null_shuffles.valueChanged.connect(
            lambda v: (self.null_windows.setEnabled(v > 0),
                       self.null_seed.setEnabled(v > 0)))
        self.null_shuffles.setSpecialValueText("off")
        form.addRow("Null-model shuffles:", self.null_shuffles)
        form.addRow("Null-model windows:", self.null_windows)
        form.addRow("Null-model seed:", self.null_seed)
        return box

    def _refresh_specificity_availability(self) -> None:
        """Enable only the specificity settings the chosen engine can honour.

        BLAST cannot measure the junction-duplex gate (RunConfig.validate
        refuses it), and the local-alignment count is meaningful only under
        BLAST. Disabled rather than hidden, so the note can say why.
        """
        if not hasattr(self, "ot_counting"):
            return
        is_blast = self.screen_engine.currentData() == ENGINE_BLAST
        self.junction_duplex.setEnabled(not is_blast)
        self.junction_duplex_bp.setEnabled(
            not is_blast and self.junction_duplex.isChecked())
        local_i = self.ot_counting.findData(OT_COUNT_LOCAL)
        self.ot_counting.model().item(local_i).setEnabled(is_blast)
        if not is_blast and self.ot_counting.currentData() == OT_COUNT_LOCAL:
            self.ot_counting.setCurrentIndex(
                self.ot_counting.findData(OT_COUNT_FULL_HALF))
        self.specificity_note.setText(
            "Under BLAST the junction-duplex filter is switched off, because "
            "BLAST cannot measure it." if is_blast else
            "Local-alignment counting is available only under the BLAST "
            "engine.")

    def _build_ranking_group(self):
        # --- Ranking weights (advanced, opt-in) ---
        box = QGroupBox("Ranking weights (advanced)")
        box.setCheckable(True)
        box.setChecked(False)
        self.ranking_box = box
        form = QFormLayout(box)
        self._ranking_spins = {}
        # (key, label, default) — defaults mirror RankingWeights.
        for key, label, default in [
            ("gc_proximity", "GC proximity", 0.25),
            ("gc_balance", "GC balance", 0.15),
            ("specificity_margin", "Specificity margin", 0.30),
            ("homopolymer_penalty", "Homopolymer penalty", 0.15),
            ("central_placement", "Central placement", 0.15),
            ("isoform_coverage", "Isoform coverage", 0.20),
            ("mismatch_placement", "Mismatch placement",
             DEFAULT_WEIGHTS.mismatch_placement),
        ]:
            sp = QDoubleSpinBox()
            sp.setRange(0.0, 1.0)
            sp.setSingleStep(0.05)
            sp.setDecimals(2)
            sp.setValue(default)
            self._ranking_spins[key] = sp
            form.addRow(f"{label}:", sp)
        # Isoform coverage only differentiates probes in gene-aware mode (every
        # probe is coverage 1.0 otherwise), so tie it to the group-by-gene
        # toggle, mirroring the min-coverage spinner. self.group_by_gene already
        # exists — the BLAST group is built before the ranking group.
        iso_sp = self._ranking_spins["isoform_coverage"]
        iso_sp.setToolTip(
            "Weight for isoform coverage in the composite score. Only matters "
            "when 'Group transcripts by gene' is on (otherwise every probe is "
            "coverage 1.0). Higher = prefer probes that bind more of their "
            "gene's isoforms. Left unchecked, gene-aware runs default to 0.20."
        )
        iso_sp.setEnabled(self.group_by_gene.isChecked())
        self.group_by_gene.toggled.connect(iso_sp.setEnabled)
        note = QLabel(
            "Weights are relative importances — they need not sum to 1; the "
            "score is divided by their total."
        )
        note.setWordWrap(True)
        note.setStyleSheet(f"color: {muted_hex()}; font-size: 11px;")
        form.addRow("", note)
        box.setToolTip(
            "Override the seven composite-score weights. Unchecked = use the "
            "defaults. Mismatch placement applies only to conserved homology "
            "runs, where it prefers probes whose on-target mismatches sit away "
            "from the ligation junction."
        )
        return box

    def _build_blast_group(self):
        # --- off-target screening ---
        blast_box = QGroupBox("Off-target screening")
        blast_outer = QVBoxLayout(blast_box)
        blast_l = QFormLayout()
        blast_l.setFieldGrowthPolicy(QFormLayout.ExpandingFieldsGrow)
        blast_outer.addLayout(blast_l)

        # Engine choice. Exposed because the k-mer default cannot search a
        # reference larger than ~2.1 Gbases (a genome), and its error message
        # tells the user to "use --screen-engine blast" — a CLI flag a desktop
        # user cannot reach. Without this control the only escape hatch from the
        # default engine is unavailable to the people most likely to need it.
        self.screen_engine = QComboBox()
        self.screen_engine.addItem("k-mer — exhaustive, end-to-end (default)",
                                   ENGINE_KMER)
        self.screen_engine.addItem("BLAST — legacy; needed for genome-scale "
                                   "references", ENGINE_BLAST)
        self.screen_engine.setCurrentIndex(
            0 if SCREEN_ENGINE_DEFAULT == ENGINE_KMER else 1
        )
        self.screen_engine.setToolTip(
            "k-mer compares every probe half to every candidate position across "
            "all 25 bases, so mismatch counts are exact and no BLAST database "
            "is needed. Point 'Off-target references' in the Inputs section at "
            "the reference FASTA.\n\n"
            "BLAST is the legacy engine. It reports a local alignment, so its "
            "verdicts rest on reconstructed positions. Where the two engines "
            "disagree, the k-mer engine agrees with exhaustive ground truth. "
            "Choose BLAST only for a reference too large for the k-mer index "
            "(over ~2.1 Gbases, i.e. a genome).\n\n"
            "THE ENGINE ALSO CHANGES THE RANKING, substantially and by design. "
            "The composite score's specificity_margin term carries a large "
            "share of the weight. Under BLAST it comes from trimmed-HSP "
            "identity, which the alignment length floor pins near 100% for "
            "almost every probe, so it often carries no information. Under "
            "k-mer it is a real measurement. The two engines' top-N "
            "selections can share few or no probes.\n\n"
            "So runs made under different engines are NOT rank-comparable, and "
            "the same input can give you a different panel. The engine is "
            "recorded as screen_engine in results.csv and the manifest, which "
            "is what tells two such runs apart."
        )
        blast_l.addRow("Engine", self.screen_engine)

        # The ranking consequence above is the kind of thing nobody discovers
        # in a tooltip, and it is invisible in the output: two runs of the same
        # FASTA just quietly disagree about which probes are best. Say it on
        # the form, and only when the user has actually left the default.
        self.engine_warning = QLabel()
        self.engine_warning.setWordWrap(True)
        self.engine_warning.setTextFormat(Qt.RichText)
        self.engine_warning.setStyleSheet("color:#b00020;")
        self.engine_warning.setVisible(False)
        blast_l.addRow("", self.engine_warning)
        self.screen_engine.currentIndexChanged.connect(
            self._refresh_engine_warning)
        self._refresh_engine_warning()

        self.skip_screen = QCheckBox(
            "Skip off-target screening (filter-only smoke test)")
        # Screening runs by default; tick this only for a quick filter-only check.
        self.skip_screen.setChecked(False)
        blast_l.addRow("", self.skip_screen)
        # The references themselves sit in the Inputs group at the top.

        self.cached_hits = QLineEdit()
        ch = QHBoxLayout(); ch.setContentsMargins(0, 0, 0, 0)
        ch.addWidget(self.cached_hits)
        b4 = QPushButton("Browse…"); b4.clicked.connect(self._browse_cached); ch.addWidget(b4)
        blast_l.addRow("Cached hits JSON:", _wrap(ch))

        # Isoform awareness: accessions that must NOT count as off-target (e.g.
        # other transcript isoforms of the target gene that a valid probe will
        # legitimately hit).
        self.extra_on_target = QPlainTextEdit()
        self.extra_on_target.setPlaceholderText(
            "One accession per line — hits to these are NOT treated as "
            "off-target (e.g. other isoforms of the target gene)."
        )
        self.extra_on_target.setFixedHeight(56)
        blast_l.addRow("On-target isoforms:", self.extra_on_target)

        # Gene-aware / whole-transcriptome isoform mode: group input transcripts
        # by gene (Ensembl "gene:" header field) so sibling isoforms count as
        # on-target and each probe is scored on isoform coverage.
        self.group_by_gene = QCheckBox(
            "Group transcripts by gene (isoform-aware: reward probes that hit "
            "all isoforms)"
        )
        self.group_by_gene.setToolTip(
            "Whole-transcriptome mode. Groups input transcripts by their "
            "Ensembl 'gene:' header field, so a probe hitting a sibling isoform "
            "of the same gene is NOT counted as off-target, and probes are "
            "scored higher when they bind more of their gene's isoforms. Use "
            "the same transcriptome FASTA as both input and off-target "
            "reference. On by default: on input with no 'gene:' field "
            "(viral or custom sequences) there is nothing to group, and the "
            "run proceeds ungrouped and says so."
        )
        from flexprobe.constants import GROUP_BY_GENE_DEFAULT
        self.group_by_gene.setChecked(GROUP_BY_GENE_DEFAULT)
        blast_l.addRow("", self.group_by_gene)

        self.min_isoform_coverage = QDoubleSpinBox()
        self.min_isoform_coverage.setRange(0.0, 1.0)
        self.min_isoform_coverage.setSingleStep(0.05)
        self.min_isoform_coverage.setValue(0.0)
        self.min_isoform_coverage.setToolTip(
            "Optional hard gate (0 = off, soft scoring only). Drops probes that "
            "bind less than this fraction of their gene's isoforms — e.g. 0.80 "
            "hides probes covering under 80% of isoforms. Needs gene grouping."
        )
        # Follows the checkbox from the start, which is ticked by default
        # since G15; a hard-coded False left the gate greyed out while on.
        self.min_isoform_coverage.setEnabled(self.group_by_gene.isChecked())
        self.group_by_gene.toggled.connect(self.min_isoform_coverage.setEnabled)
        blast_l.addRow("Min isoform coverage:", self.min_isoform_coverage)

        # Whole-transcriptome speedup: only BLAST the best N candidates per gene.
        self.screen_top_per_gene = QSpinBox()
        self.screen_top_per_gene.setRange(0, 100000)
        self.screen_top_per_gene.setValue(0)
        self.screen_top_per_gene.setToolTip(
            "Whole-transcriptome speedup (0 = off, screen everything). Off-target "
            "screen only the top N candidates per gene by the screen-free "
            "pre-score (GC, balance, homopolymer, central placement), dropping "
            "the rest before the screen runs. Cuts screening volume ~10-50x on a "
            "whole transcriptome — you only select a few probes per gene anyway. "
            "Try ~100. Applies to both engines."
        )
        blast_l.addRow("Screen top N / gene:", self.screen_top_per_gene)

        # Screen parallelism. The dialog set --workers and never
        # --screen-workers, so a GUI whole-transcriptome run screened on ONE
        # core while the HPC bundle it exported set both to the core count.
        self.screen_workers = QSpinBox()
        self.screen_workers.setRange(1, 64)
        self.screen_workers.setValue(1)
        self.screen_workers.setToolTip(
            "Worker processes for the OFF-TARGET SCREEN itself. 1 = serial "
            "(the right choice for a handful of genes).\n\n"
            "On a whole transcriptome the screen dominates the runtime, so "
            "this is the setting that decides how long the run takes. It "
            "changes who does the arithmetic, never the answer: it is not part "
            "of the checkpoint fingerprint, so raising it resumes an "
            "in-progress run instead of invalidating it.\n\n"
            "Workers share the reference and index in memory, so each extra "
            "worker adds little by itself. But the screen starts from the "
            "tiling stage's memory peak, so on a whole transcriptome make "
            "sure the machine's RAM covers that first: an out-of-memory screen "
            "is killed, not slowed. The pipeline "
            "additionally clamps this to the cores it can actually see.\n\n"
            "k-mer engine only. Under BLAST the screen runs serially whatever "
            "this is set to."
        )
        blast_l.addRow("Screen worker processes:", self.screen_workers)

        # Targeting guidance — the single most common reason every probe fails
        # is including a near-identical sequence (e.g. the target's own close
        # relative) in the off-target set, so every candidate looks non-specific.
        help_btn = QPushButton(
            "What should I screen against? (avoid \"all probes failed\")")
        help_btn.setFlat(True)
        help_btn.setStyleSheet(f"text-align:left; color:{link_hex()}; padding:2px 0;")
        help_btn.clicked.connect(self._show_blast_targeting_help)
        blast_outer.addWidget(help_btn)
        return blast_box

    def _refresh_engine_warning(self) -> None:
        """Warn on the form when the non-default screen engine is selected."""
        is_blast = self.screen_engine.currentData() == ENGINE_BLAST
        if is_blast:
            self.engine_warning.setText(
                "\u26a0 BLAST is the legacy engine. Its verdicts rest on "
                "reconstructed alignment positions, and where the two engines "
                "disagree the k-mer engine is the one that matches exhaustive "
                "ground truth. BLAST cannot measure the junction-duplex "
                "filter, so that filter is <b>switched off</b> for BLAST runs "
                "and the screen is weaker. It also <b>changes the ranking</b>: "
                "top-N picks can share few or no probes with a k-mer run of the "
                "same input, so the two are not rank-comparable. Screen worker processes do not apply. Use it "
                "only for a reference the k-mer index cannot address (over "
                "~2.1 Gbases, i.e. a genome)."
            )
        self.engine_warning.setVisible(is_blast)

    def _show_blast_targeting_help(self) -> None:
        QMessageBox.information(
            self, "Off-target reference guide",
            _BLAST_TARGETING_HELP_HTML,
        )

    # ---- Project + run-mode handlers ----

    def _refresh_project_combo(self) -> None:
        # Block signals around the rebuild so currentIndexChanged doesn't
        # fire mid-population.
        self.project_combo.blockSignals(True)
        try:
            current_id = self.project_combo.currentData() if self.project_combo.count() else ""
            self.project_combo.clear()
            self.project_combo.addItem("(no project)", userData="")
            for proj in self._projects.list():
                self.project_combo.addItem(proj.name, userData=proj.project_id)
            if current_id:
                idx = self.project_combo.findData(current_id)
                if idx >= 0:
                    self.project_combo.setCurrentIndex(idx)
        finally:
            self.project_combo.blockSignals(False)

    def _on_new_project(self) -> None:
        from PySide6.QtWidgets import QInputDialog
        name, ok = QInputDialog.getText(
            self, "New project", "Project name:",
        )
        if not ok or not name.strip():
            return
        try:
            proj = self._projects.create(name.strip())
        except Exception as exc:
            QMessageBox.warning(self, "Project", str(exc))
            return
        self._refresh_project_combo()
        idx = self.project_combo.findData(proj.project_id)
        if idx >= 0:
            self.project_combo.setCurrentIndex(idx)

    def _on_run_mode_changed(self, label: str) -> None:
        auto = label.startswith("Auto")
        self.auto_n.setVisible(auto)
        self.auto_n_label.setVisible(auto)

    def _on_strain_mode_changed(self, _label: str) -> None:
        is_multi = self.strain_mode_combo.currentText().startswith("Homology")
        for w in (self.strategy_label, self.strategy_combo):
            w.setVisible(is_multi)
        is_unique = is_multi and self.strategy_combo.currentText().startswith("Sequence-unique")
        for w in (self.unique_target_label, self.unique_target_combo,
                  self.unique_min_mm_label, self.unique_min_mm_row,
                  self.unique_jd_label, self.unique_jd_row,
                  self.unique_per_half):
            w.setVisible(is_unique)
        is_constitutive = is_multi and self.strategy_combo.currentText().startswith("Constitutive")
        for w in (self.isoform_min_label, self.isoform_min_row):
            w.setVisible(is_constitutive)
        # The degenerate tolerance is a property of the conserved rule only:
        # sequence-unique designs from one record's own sequence and
        # constitutive groups by exact shared sequence.
        is_conserved = is_multi and not is_unique and not is_constitutive
        for w in (self.conserved_degen_label, self.conserved_degen_row):
            w.setVisible(is_conserved)
        # The template only means anything where a consensus is emitted;
        # sequence-unique already picks its sequence via "Target sequence".
        for w in (self.template_label, self.template_combo):
            w.setVisible(is_multi and not is_unique)
        if is_multi:
            hint = (
                "Homology mode: every record in the chosen FASTA is treated as a "
                "related sequence (e.g. strain or isoform) of the same fragment. "
                "Enter the actual individual sequences — the pipeline aligns them, "
                "builds a consensus and adds per-sequence mismatch + ΔTm columns "
                "to results.csv. At a variant position a designed probe takes "
                "the base of the chosen Consensus template (the first record by "
                "default), so the exported order files contain no N's."
            )
            if is_constitutive:
                hint += (
                    "\n\nConstitutive (isoform-collapse): designs probes from "
                    "regions shared by ≥ \"Min isoforms shared\" records (0 = all), "
                    "so one probe binds every covered isoform.\n"
                    "⚠ Records within 5 nt of each other in length are trimmed "
                    "to the shortest, not aligned. Otherwise they are aligned "
                    "with MAFFT (bundled with the app, else found on the "
                    "system), then MUSCLE, and only if neither is available a "
                    "built-in progressive aligner that is not a true multiple "
                    "alignment. The aligner used is recorded in the run's "
                    "manifest. For many isoforms or large exon "
                    "insertions/deletions, check it, or supply records already "
                    "aligned to equal length."
                )
            elif is_unique:
                hint += (
                    "\n\nSequence-unique: designs probes that distinguish the "
                    "chosen target sequence from every other by ≥ the set mismatch "
                    "count."
                )
            if not is_unique:
                hint += (
                    "\n\nConsensus template: at a position where the records "
                    "disagree, the ordered probe has to carry ONE concrete base — "
                    "it takes the template's. A single run therefore only ever "
                    "yields that record's flavour of each window, and the other "
                    "records' equally valid probes are never generated. Choose "
                    "\"Every sequence\" to design once per record and queue them "
                    "all; the union is the complete set."
                )
            self.strain_hint.setText(hint)
            self._on_fasta_path_changed(self.fasta_edit.text())
        else:
            self.strain_hint.setText("")

    def _on_fasta_path_changed(self, path: str) -> None:
        # Populate the "Target sequence" / "Consensus template" combos in
        # homology mode by peeking at the FASTA record IDs. Errors are silent —
        # the picker just stays empty and the user will see a clearer error at
        # run time.
        # Gate on the strategy, not isVisible(): this runs from __init__ (via
        # _on_strain_mode_changed) before the dialog is ever shown, so the
        # visibility form left the combo unpopulated for a draft/re-run that
        # opens straight into sequence-unique mode. An empty combo makes
        # currentData() None, which is the second route to a silent target 0.
        _strat = self.strategy_combo.currentText()
        if not self.strain_mode_combo.currentText().startswith("Homology"):
            return
        ids: list = []
        try:
            from flexprobe.errors import FlexprobeError
            from flexprobe.io_fasta import parse_fasta
            ids = [gid for gid, _desc, _seq in parse_fasta(path.strip())]
        except (OSError, ValueError, FlexprobeError) as exc:
            logger.debug("Could not pre-read strain IDs from %s: %s", path, exc)
            ids = []
        self._populate_template_combo(ids)
        if not _strat.startswith("Sequence-unique"):
            return
        keep = self.unique_target_combo.currentData()
        self.unique_target_combo.blockSignals(True)
        try:
            self.unique_target_combo.clear()
            for i, sid in enumerate(ids):
                self.unique_target_combo.addItem(sid, userData=i)
            # The same fan-out entry as the template picker, and build_runs
            # reads the same ALL_TEMPLATES sentinel from either combo. Until
            # now it was never added here, so a strain-unique fan-out could not
            # be selected at all. APPENDED rather than leading: a strain-unique
            # run has always defaulted to the first record, and an untouched
            # form must not start queueing one run per record.
            if len(ids) > 1:
                self.unique_target_combo.addItem(
                    f"Every sequence ({len(ids)} runs, one per target)",
                    userData=ALL_TEMPLATES,
                )
            self._select_after_fill("unique_target_combo", keep)
        finally:
            self.unique_target_combo.blockSignals(False)

    def _select_after_fill(self, attr: str, keep) -> None:
        """Re-select ``keep`` (or a draft's deferred choice) after a refill.

        By userData, not index: indices shift when the fan-out entry appears
        or the record count changes. A deferred draft choice wins over
        ``keep``, which is only the empty combo's None at that point, and is
        consumed once the combo holds records, whether or not it matched.
        """
        combo = getattr(self, attr)
        pending = self._pending_combo_data.get(attr)
        want = pending if pending is not None else keep
        idx = combo.findData(want) if want is not None else -1
        if idx >= 0:
            combo.setCurrentIndex(idx)
        if combo.count():
            self._pending_combo_data.pop(attr, None)

    def _populate_template_combo(self, ids: list) -> None:
        """Rebuild the consensus-template picker from the FASTA record ids.

        The fan-out entry leads because it is the one that returns the COMPLETE
        pan-strain set; picking a single template is the narrower choice, even
        though it is what the engine did unconditionally before this existed.
        Only offered for 2+ records — with one record there is nothing to fan
        out to and every template is the same sequence.
        """
        keep = self.template_combo.currentData()
        self.template_combo.blockSignals(True)
        try:
            self.template_combo.clear()
            if len(ids) > 1:
                self.template_combo.addItem(
                    f"Every sequence — {len(ids)} runs, one per template",
                    userData=ALL_TEMPLATES,
                )
            for i, sid in enumerate(ids):
                self.template_combo.addItem(sid, userData=i)
            self._select_after_fill("template_combo", keep)
        finally:
            self.template_combo.blockSignals(False)

    def _on_mode_changed(self, _label: str) -> None:
        # Toggle the Flex v1 multi-barcode picker when entering/leaving
        # multiplex, then refresh the preview.
        self._refresh_multi_barcode_visibility()
        self._refresh_handle_preview()

    def _is_flex_v1_multiplex(self) -> bool:
        """True when the current modality+mode is Flex v1 multiplex (the only
        case that uses per-sample Probe Barcodes)."""
        modality = MODALITY_PRESETS[self.modality.currentText()]
        return (
            get_chemistry(modality).key == "flex_v1"
            and self.mode.currentText() == "Multiplex"
        )

    def _on_multi_count_changed(self, n: int) -> None:
        bcs = self._barcodes_for_count(n)
        if bcs:
            self.multi_count_hint.setText(f"→ {bcs[0]}–{bcs[-1]}")
        else:
            self.multi_count_hint.setText("")
        # During construction the handle preview may not exist yet; it is
        # refreshed at the end of __init__ via _on_modality_changed.
        if hasattr(self, "handle_preview"):
            self._refresh_handle_preview()

    @staticmethod
    def _barcodes_for_count(n: int) -> list:
        """Probe Barcodes BC001..BC00N for an N-sample multiplex run."""
        n = max(0, min(int(n), len(FLEX_V1_PROBE_BARCODES)))
        return [f"BC{i:03d}" for i in range(1, n + 1)]

    def _refresh_multi_barcode_visibility(self) -> None:
        show_multi = self._is_flex_v1_multiplex()
        self.multi_bc_label.setVisible(show_multi)
        self.multi_bc_widget.setVisible(show_multi)
        # When the sample-count selector is visible, the single-pick combo
        # above is redundant — hide it to avoid confusion.
        if show_multi:
            self.barcode_label.setVisible(False)
            self.barcode_combo.setVisible(False)

    def _selected_multi_barcodes(self) -> list:
        # Only meaningful for Flex v1 multiplex — gate on the modality+mode
        # condition (robust) rather than widget.isVisible() (false when the
        # dialog isn't shown, e.g. in tests). Returns BC001..BC00N for the
        # chosen sample count.
        if not self._is_flex_v1_multiplex():
            return []
        return self._barcodes_for_count(self.multi_count.value())

    # ---- Modality auto-populate + visibility ----

    def _on_modality_changed(self, label: str) -> None:
        defaults = MODALITY_DEFAULTS.get(label) or {}
        modality = MODALITY_PRESETS[label]
        chemistry = get_chemistry(modality)

        # Filter / tile values
        # Locked to the chemistry: every current chemistry assembles 50 bp
        # probes and drops any other tile, so any other size produced no
        # probes at all. The box stays visible so the value is on screen.
        _tl = int(chemistry.tile_length)
        self.tile_size.setRange(_tl, _tl)
        self.tile_size.setValue(_tl)
        self.tile_size.setEnabled(False)
        self.tile_size.setToolTip(
            f"Fixed at {_tl} bp by the chemistry: two {_tl // 2} nt halves. "
            "Other sizes are dropped when the probes are assembled.")
        self.tile_step.setValue(int(defaults.get("tile_step", 1)))
        self.gc_min.setValue(float(defaults.get("gc_min", 0.44)))
        self.gc_max.setValue(float(defaults.get("gc_max", 0.72)))
        self.enforce_gc.setChecked(bool(defaults.get("enforce_gc_range", False)))
        self.tm_min.setValue(float(defaults.get("tm_min", TM_MIN_DEFAULT)))
        self.max_n.setValue(int(defaults.get("max_n_per_half", 0)))
        self.max_hp.setValue(int(defaults.get("max_homopolymer_run", 4)))
        self._set_tn_mode(TN_MODE_STRICT_A
                          if defaults.get("require_tn_junction", True)
                          else TN_MODE_OFF)
        self.primer3.setChecked(bool(defaults.get("enable_primer3", False)))

        # Mode: meaningful when the chemistry can actually distinguish multiplex
        # from singleplex. That's true either when the RHS handles differ
        # (Flex v2) OR when the chemistry embeds a per-sample Probe Barcode
        # (Flex v1) — there multiplex = several barcodes (BC001..BC00N),
        # singleplex = a single barcode, while the RHS handle stays the same.
        mode_enabled = (
            chemistry.has_probe_barcode
            or (
                chemistry.rhs_3prime_singleplex is not None
                and chemistry.rhs_3prime_singleplex != chemistry.rhs_3prime_multiplex
            )
        )
        self.mode.setEnabled(mode_enabled)
        self.mode.setCurrentText("Multiplex" if defaults.get("multiplex", True) else "Singleplex")

        # Probe Barcode: populate dropdown when chemistry uses one (Flex v1).
        self.barcode_combo.blockSignals(True)
        self.barcode_combo.clear()
        if chemistry.has_probe_barcode:
            for bc_label, bc_seq in FLEX_V1_PROBE_BARCODES.items():
                self.barcode_combo.addItem(f"{bc_label}  ({bc_seq})", userData=bc_label)
            preset = defaults.get("probe_barcode", chemistry.default_probe_barcode)
            idx = self.barcode_combo.findData(preset)
            if idx >= 0:
                self.barcode_combo.setCurrentIndex(idx)
            self.barcode_label.setVisible(True)
            self.barcode_combo.setVisible(True)
        else:
            self.barcode_combo.addItem("(not used for this modality)", userData="")
            self.barcode_label.setVisible(False)
            self.barcode_combo.setVisible(False)
        self.barcode_combo.blockSignals(False)

        # Multiplex sample-count selector visibility depends on (modality, mode).
        self._refresh_multi_barcode_visibility()
        self._refresh_tn_mode_availability(chemistry)

        self._refresh_handle_preview()

    # ---- TN junction mode ----

    def _tn_mode(self) -> str:
        """The selected junction gate: strict_a | soft | off."""
        return self.tn_junction_mode.currentData() or TN_MODE_STRICT_A

    def _set_tn_mode(self, mode: str) -> None:
        idx = self.tn_junction_mode.findData(mode)
        if idx >= 0:
            self.tn_junction_mode.setCurrentIndex(idx)

    def _refresh_tn_mode_availability(self, chemistry) -> None:
        """Enable/disable the soft gate for the current chemistry.

        ``chemistry.resolve_tn_mode`` raises ConfigError for soft on Visium
        rather than quietly substituting a different gate, so a dialog that let
        the user pick it there would build a config that cannot run. Disabled
        rather than removed, so the option stays visible and the note says why.
        """
        allowed = bool(getattr(chemistry, "allows_soft_tn_junction", False))
        idx = self.tn_junction_mode.findData(TN_MODE_SOFT)
        if idx >= 0:
            item = self.tn_junction_mode.model().item(idx)
            if item is not None:
                item.setEnabled(allowed)
        if not allowed and self._tn_mode() == TN_MODE_SOFT:
            # Carried in from a draft or the previous modality. Fall back to the
            # conservative gate and say so, rather than silently re-gating.
            self._set_tn_mode(TN_MODE_STRICT_A)
            self.tn_junction_note.setText(
                "Soft is not available for this chemistry — CG000621 Rev D "
                "states the junction adenine as a requirement for spatial "
                "probes. Reset to Hard."
            )
            self.tn_junction_note.setVisible(True)
            return
        if not allowed:
            self.tn_junction_note.setText(
                "Soft is Flex-only; CG000621 Rev D states the junction adenine "
                "as a requirement for spatial probes."
            )
        self.tn_junction_note.setVisible(not allowed)

    def _refresh_handle_preview(self) -> None:
        label = self.modality.currentText()
        modality = MODALITY_PRESETS[label]
        chem = get_chemistry(modality)
        bc_label = self.barcode_combo.currentData() or ""

        multiplex = self.mode.currentText() == "Multiplex"
        # If the Flex v1 multi-barcode picker is visible AND >1 barcodes are
        # checked, show the generic {BC} placeholder rather than picking one.
        multi_labels = self._selected_multi_barcodes()
        if chem.key == "flex_v1" and len(multi_labels) > 1:
            leave_placeholder = True
            preview_bc = None
        else:
            leave_placeholder = False
            preview_bc = bc_label or None

        lhs_regions = chem.lhs_regions()
        rhs_regions = chem.rhs_regions(
            multiplex=multiplex, probe_barcode=preview_bc,
            leave_placeholder=leave_placeholder,
        )
        self.handle_preview.setText(
            _render_regions_html(lhs_regions, rhs_regions)
        )

    # ---- File pickers ----

    def _browse_fasta(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select FASTA", "",
            "FASTA files (*.fa *.fasta *.fna);;All files (*)",
        )
        if path:
            self.fasta_edit.setText(path)

    def _browse_out(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Select output base directory", self.out_edit.text().strip(),
        )
        if path:
            self.out_edit.setText(path)
            self._remember_output_dir(path)

    def _remember_output_dir(self, path: str) -> None:
        """Make ``path`` the default output folder for the next New query.

        Decided 2026-09-15: a folder chosen here sticks. Before, Browse changed
        this query only and the next New query reverted to the saved default,
        which read to users as the app forgetting their choice. Saved through
        AppSettings, so Preferences shows the same folder.
        """
        path = (path or "").strip()
        if not path or path == self._settings.default_output_dir():
            return
        try:
            self._settings.set_default_output_dir(path)
        except OSError as exc:
            # Never block a run over a preference; the run uses the folder
            # either way.
            logger.warning("Could not save the default output folder: %s", exc)

    def _append_db_line(self, line: str) -> None:
        text = self.db_edit.toPlainText().rstrip()
        self.db_edit.setPlainText((text + "\n" if text else "") + line)
        # Adding a reference only makes sense when BLAST will actually run.
        if line:
            self.skip_screen.setChecked(False)

    def _browse_db(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Pick any file from a BLAST DB", "",
            "BLAST DB files (*.nhr *.nin *.nsq);;All files (*)",
        )
        if path:
            base, _ = os.path.splitext(path)
            self._append_db_line(base)

    def _browse_blast_fasta(self) -> None:
        # Add a raw FASTA verbatim — the pipeline auto-indexes it (cached by
        # SHA-256) right before the run via ensure_db_for_path.
        path, _ = QFileDialog.getOpenFileName(
            self, "Pick a FASTA of off-target sequences", "",
            "FASTA files (*.fa *.fasta *.fna *.ffn *.frn);;All files (*)",
        )
        if path:
            self._append_db_line(path)

    def _browse_cached(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select cached hits JSON", "", "JSON (*.json);;All files (*)",
        )
        if path:
            self.cached_hits.setText(path)

    # ---- Validation + build ----

    def _validate_inputs(self) -> bool:
        """Validate the form; warn and return False on the first problem."""
        if self._is_bulk:
            if not self._bulk_paths:
                QMessageBox.warning(
                    self, "Missing FASTAs", "Bulk mode received an empty file list.",
                )
                return False
        elif not self.fasta_edit.text().strip():
            QMessageBox.warning(self, "Missing FASTA", "Please choose a FASTA file.")
            return False
        if not self.out_edit.text().strip():
            QMessageBox.warning(self, "Missing output", "Please choose an output directory.")
            return False
        if self.gc_min.value() > self.gc_max.value():
            QMessageBox.warning(self, "GC range", "GC min must be <= GC max.")
            return False
        # Then ask the engine itself. Every combination it refuses (a Tm
        # ceiling below the floor, for one) is caught here, with the engine's
        # own message, instead of failing the run after Start.
        try:
            runs = self.build_runs()
            if runs:
                runs[0][1].validate()
        except Exception as exc:  # noqa: BLE001 - shown to the user verbatim
            QMessageBox.warning(self, "Settings not accepted", str(exc))
            return False
        return True

    def _validate_and_accept(self) -> None:
        if self._validate_inputs():
            # A folder typed rather than browsed is remembered when the run is
            # started, not per keystroke, so a half-typed path never sticks.
            self._remember_output_dir(self.out_edit.text())
            self.accept()

    def _on_export_hpc(self) -> None:
        """Build the run config(s) and hand off to the HPC-export dialog.

        Reuses ``build_runs()`` so the exported config is byte-for-byte the one
        a local run would use. On success the new-query dialog closes WITHOUT
        queuing a local run — the work happens on the cluster instead.
        """
        if not self._validate_inputs():
            return
        try:
            runs = self.build_runs()
        except Exception as exc:  # build_runs may raise on malformed input
            QMessageBox.critical(
                self, "Export for HPC", f"Could not build the run config:\n{exc}",
            )
            return
        from .hpc_dialog import HpcExportDialog
        dlg = HpcExportDialog(runs, self._settings, self)
        if dlg.exec() == QDialog.Accepted:
            self.reject()

    # ---- Save / load query drafts ----

    # Draft save/load (to_draft/apply_draft/_scalar_*/_SCALAR_FIELDS)
    # lives in DraftMixin — see _nq_draft.py.

    # build_runs/build (assemble RunConfig from the form) live in
    # BuildMixin — see _nq_build.py.


# ---------------------------------------------------------------------------
# BLAST targeting guidance
# ---------------------------------------------------------------------------

_BLAST_TARGETING_HELP_HTML = (
    "<b>What the off-target screen does here</b><br>"
    "Every candidate probe is searched against the references you list — by the "
    "k-mer engine unless you switched to BLAST. A probe "
    "is rejected when it also matches something it <i>shouldn't</i> (an "
    "off-target). So the reference set is your list of sequences the probe must "
    "<b>not</b> cross-react with — not the thing you're designing against.<br><br>"

    "<b>What to include</b><br>"
    "&bull; The host transcriptome/genome (e.g. human tx) so probes don't bind host RNA.<br>"
    "&bull; Other organisms in the sample you want to exclude (e.g. other viruses).<br><br>"

    "<b>What to leave OUT (this is what causes \"all probes failed\")</b><br>"
    "Do <b>not</b> include sequences that are nearly identical to your actual "
    "target. If you do, every probe matches them as well as the target and gets "
    "thrown out as non-specific — leaving zero probes.<br><br>"

    "<b>Example (OROV)</b><br>"
    "Designing probes for the OROV long (L) fragment: screen against the human "
    "transcriptome and all other virus genomes, <i>but exclude the L fragments "
    "of the other OROV strains</i>. Including those near-identical L sequences "
    "makes both long-fragment regions look non-specific and you end up with no "
    "suitable probes. Screen against them separately (or accept the cross-match) "
    "rather than dumping them into the off-target set.<br><br>"

    "<b>Tip</b><br>"
    "You can list a FASTA directly (\"Add FASTA (auto-index)\") — the k-mer "
    "engine reads it as-is, and the BLAST engine turns it into a database for "
    "you. Tick \"Skip off-target screening\" only for a quick filter-only "
    "smoke test; real designs should screen against host + neighbours."
)


# ---------------------------------------------------------------------------
# Colour-coded handle preview rendering
# ---------------------------------------------------------------------------


def _render_regions_html(lhs_regions, rhs_regions) -> str:
    """Render LHS + RHS region lists as HTML for the preview ``QLabel``.

    Each ``HandleRegion`` is a coloured background swatch with the sequence
    written inside it; the second line is a colour-keyed legend.
    """
    from html import escape

    def _swatch(seq: str, color: str) -> str:
        # Slightly lighter background + dark text for contrast across themes.
        return (
            f'<span style="background:{color}; color:white; '
            f'padding:1px 4px; margin:1px; border-radius:3px;">'
            f'{escape(seq)}</span>'
        )

    def _legend_chip(label: str, color: str) -> str:
        return (
            f'<span style="background:{color}; color:white; '
            f'padding:0px 6px; margin:1px 4px 1px 0; border-radius:3px; '
            f'font-size:10px;">{escape(label)}</span>'
        )

    lhs_seq_html = "5'&nbsp;" + "".join(
        _swatch(r.sequence, r.color) for r in lhs_regions
    ) + "&nbsp;3'"
    rhs_seq_html = "5'&nbsp;" + "".join(
        _swatch(r.sequence, r.color) for r in rhs_regions
    ) + "&nbsp;3'"
    # Dedupe legend chips by label so the same colour doesn't appear twice
    # when a chemistry uses (say) two green regions.
    legend_items = []
    seen_labels = set()
    for r in list(lhs_regions) + list(rhs_regions):
        key = (r.label, r.color)
        if key in seen_labels:
            continue
        seen_labels.add(key)
        legend_items.append(_legend_chip(r.label, r.color))
    legend_html = "".join(legend_items)

    return (
        f'<div style="font-family:monospace; font-size:11px; line-height:1.7;">'
        f'<b>LHS:</b>&nbsp;{lhs_seq_html}<br>'
        f'<b>RHS:</b>&nbsp;{rhs_seq_html}<br><br>'
        f'<div style="font-family:sans-serif;">{legend_html}</div>'
        f'</div>'
    )
