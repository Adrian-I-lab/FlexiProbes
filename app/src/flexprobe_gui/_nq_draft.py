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

"""Draft (save/load form state) for NewQueryDialog.

Split out of new_query_dialog.py to keep serialization separate from UI
construction (audit #19). ``DraftMixin`` operates entirely on the host dialog's
own widgets via ``self`` and holds no state of its own; it is mixed into
NewQueryDialog ahead of QDialog in the MRO."""
from __future__ import annotations

import os

from PySide6.QtWidgets import QMessageBox  # noqa: F401  (used by save/load draft)

from flexprobe.constants import TN_MODE_OFF, TN_MODE_STRICT_A
from flexprobe.screen import ENGINE_BLAST

from flexprobe.models import RunConfig, ScreenResources, SpecificityThresholds

from ._nq_build import (
    _OT_COUNT_DEFAULT, _OT_RULE_DEFAULT, strip_generated_run_name,
)

# Engine defaults for the specificity controls, read live, never restated.
_SP_DEFAULT = SpecificityThresholds()
_SR_DEFAULT = ScreenResources()
_R_DEFAULT = RunConfig(fasta_path="", output_dir="")


from typing import TYPE_CHECKING

if TYPE_CHECKING:  # give the mixin the host widget's type for mypy
    from PySide6.QtWidgets import QDialog as _HostBase
else:
    _HostBase = object


class DraftMixin(_HostBase):
    # Single source of truth for the scalar (spin/check) form fields so
    # to_draft and apply_draft can't drift on names or defaults — the drift that
    # once silently flipped primer3 on (audit #12/#11). Each entry is
    # (draft_key, widget_attr, kind, default). "spin" fields are only restored
    # when present in the draft (so a missing key keeps the modality default);
    # "check" fields always restore, using `default` when absent.
    _SCALAR_FIELDS = [
        ("multi_count", "multi_count", "spin_int", None),
        ("unique_min_mm", "unique_min_mm", "spin_int", None),
        ("unique_per_half", "unique_per_half", "check", True),
        ("isoform_min", "isoform_min", "spin_int", None),
        ("conserved_degen", "conserved_degen", "spin_int", None),
        ("auto_n", "auto_n", "spin_int", None),
        ("tile_size", "tile_size", "spin_int", None),
        ("tile_step", "tile_step", "spin_int", None),
        ("gc_min", "gc_min", "spin_float", None),
        ("gc_max", "gc_max", "spin_float", None),
        ("tm_min", "tm_min", "spin_float", None),
        ("enforce_gc", "enforce_gc", "check", False),
        ("max_n", "max_n", "spin_int", None),
        ("max_hp", "max_hp", "spin_int", None),
        ("workers", "workers", "spin_int", None),
        ("min_isoform_coverage", "min_isoform_coverage", "spin_float", 0.0),
        ("screen_top_per_gene", "screen_top_per_gene", "spin_int", 0),
        ("overwrite", "overwrite", "check", True),
        ("screen_workers", "screen_workers", "spin_int", None),
        # DELIBERATELY False while the engine and the modality presets default
        # True. This default applies only when the key is ABSENT, i.e. to a
        # draft saved before primer3 was a field — and such a draft describes a
        # run made WITHOUT structure scoring. Restoring it as True would change
        # the run it claims to reproduce, which is audit #12 exactly. A fresh
        # dialog gets True from the modality preset, never from here.
        ("primer3", "primer3", "check", False),
        ("skip_screen", "skip_screen", "check", False),
        ("group_by_gene", "group_by_gene", "check", False),
        ("ranking_enabled", "ranking_box", "check", False),
        # Off-target specificity (added 2026-09-22). Absent keys mean a draft
        # made before the controls existed, i.e. one built from the engine
        # defaults, so the check defaults below are the engine's.
        ("ot_min_mm", "ot_min_mm", "spin_int", None),
        ("junction_duplex", "junction_duplex", "check",
         _SP_DEFAULT.off_target_junction_duplex_bp > 0),
        ("junction_duplex_bp", "junction_duplex_bp", "spin_int", None),
        ("exclude_alt", "exclude_alt", "check",
         bool(_SR_DEFAULT.exclude_alt_scaffolds)),
        # The rest of the CLI-only design settings (added 2026-09-22).
        # Spacing in THIS order: each box caps the next one's maximum, so the
        # target is restored before the rungs below it.
        ("min_gap", "min_gap", "spin_int", None),
        ("relaxed_gap", "relaxed_gap", "spin_int", None),
        ("floor_gap", "floor_gap", "spin_int", None),
        ("tm_max", "tm_max", "spin_float", None),
        ("split_header", "split_header", "check",
         bool(_R_DEFAULT.split_header_on_pipe)),
        ("null_shuffles", "null_shuffles", "spin_int", None),
        ("null_windows", "null_windows", "spin_int", None),
        ("null_seed", "null_seed", "spin_int", None),
        ("unique_jd", "unique_jd", "check",
         _R_DEFAULT.strain_unique_junction_duplex_bp > 0),
        ("unique_jd_bp", "unique_jd_bp", "spin_int", None),
    ]

    def _scalar_to_draft(self) -> dict:
        out: dict = {}
        for key, attr, kind, _default in self._SCALAR_FIELDS:
            w = getattr(self, attr)
            if kind == "spin_int":
                out[key] = int(w.value())
            elif kind == "spin_float":
                out[key] = float(w.value())
            else:  # check
                out[key] = w.isChecked()
        return out

    def _scalar_from_draft(self, d: dict) -> None:
        for key, attr, kind, default in self._SCALAR_FIELDS:
            w = getattr(self, attr)
            if kind == "check":
                w.setChecked(bool(d.get(key, default)))
            elif key in d:                       # spins: set only if present
                try:
                    w.setValue(int(d[key]) if kind == "spin_int"
                               else float(d[key]))
                except (TypeError, ValueError):
                    pass

    def to_draft(self) -> dict:
        """Capture the current form state as a flat, JSON-serializable dict.

        The bulk FASTA list is intentionally NOT stored — it's tied to the
        files dropped this session. Everything else (chemistry, filters,
        BLAST references, project) round-trips. The output dir is recorded
        but not restored: apply_draft keeps the current default folder.
        """
        draft = {
            # Ordering-sensitive combos + free text + special fields.
            "modality": self.modality.currentText(),
            "mode": self.mode.currentText(),
            "probe_barcode": self.barcode_combo.currentData() or "",
            "strain_mode": self.strain_mode_combo.currentText(),
            "strategy": self.strategy_combo.currentText(),
            # Stored as the userData (int index or the ALL_TEMPLATES sentinel),
            # not the visible text: record ids differ between FASTAs, so a draft
            # re-applied to another input must still resolve.
            "consensus_template": self.template_combo.currentData(),
            # Same encoding for the sequence-unique target picker, which also
            # offers the ALL_TEMPLATES fan-out. None outside homology mode.
            "unique_target": self.unique_target_combo.currentData(),
            "screen_engine": self.screen_engine.currentData(),
            "ot_rule": self.ot_rule.currentData(),
            "ot_counting": self.ot_counting.currentData(),
            "ref_strand": self.ref_strand.currentData(),
            # strict_a | soft | off. `tn_junction` is still written beside it
            # so a draft saved now still loads in an older build, where it is
            # the only key that exists.
            "tn_junction_mode": self._tn_mode(),
            "tn_junction": self._tn_mode() != TN_MODE_OFF,
            "run_mode": self.run_mode.currentText(),
            "project_id": self.project_combo.currentData() or "",
            # Only what the user typed. A name restored from a pre-fix draft
            # may still carry a generated suffix; it is never written back.
            "run_name": strip_generated_run_name(self.run_name.text()),
            "fasta_path": "" if self._is_bulk else self.fasta_edit.text(),
            "output_dir": self.out_edit.text(),
            "blast_dbs": self.db_edit.toPlainText(),
            "cached_hits": self.cached_hits.text(),
            "extra_on_target": self.extra_on_target.toPlainText(),
            "ranking_weights": {
                k: float(sp.value()) for k, sp in self._ranking_spins.items()
            },
        }
        # Scalar spin/check fields from the single shared spec (no drift).
        draft.update(self._scalar_to_draft())
        return draft

    def apply_draft(self, d: dict) -> None:
        """Restore form state from a ``to_draft`` payload.

        Order matters: modality is applied first because changing it resets
        filters/tile to that modality's defaults — we apply the saved filter
        values *after* so they win.
        """
        def _set_combo_text(combo, text):
            if text:
                idx = combo.findText(text)
                if idx >= 0:
                    combo.setCurrentIndex(idx)

        def _set_combo_data(combo, data):
            idx = combo.findData(data)
            if idx >= 0:
                combo.setCurrentIndex(idx)

        # 1. Modality first (rebuilds barcode combo + resets filter defaults).
        _set_combo_text(self.modality, d.get("modality", ""))
        # 2. Mode (toggles multi-barcode picker).
        _set_combo_text(self.mode, d.get("mode", ""))
        # 3. Barcodes.
        _set_combo_data(self.barcode_combo, d.get("probe_barcode", ""))
        if "multi_count" not in d and d.get("multi_barcodes"):
            # Back-compat with drafts saved before the sample-count selector:
            # restore the number of samples from the old checked-barcode list.
            self.multi_count.setValue(max(2, len(d["multi_barcodes"])))
        # 4. Homology mode (toggles strategy / target rows).
        _set_combo_text(self.strain_mode_combo, d.get("strain_mode", ""))
        _set_combo_text(self.strategy_combo, d.get("strategy", ""))
        # After the strategy, because that is what repopulates the record
        # pickers. Before the FASTA is set a picker is usually still EMPTY, so
        # the choice is deferred and applied when the FASTA fills it (see
        # NewQueryDialog._select_after_fill). Selecting straight into the empty
        # combo, as this used to, found nothing: a Re-run of a template-1 run
        # came back as "Every sequence" and queued a run per record.
        for key, attr in (("consensus_template", "template_combo"),
                          ("unique_target", "unique_target_combo")):
            if d.get(key) is None:
                continue
            combo = getattr(self, attr)
            idx = combo.findData(d[key])
            if idx >= 0:
                combo.setCurrentIndex(idx)
                self._pending_combo_data.pop(attr, None)
            else:
                self._pending_combo_data[attr] = d[key]
        # 4b. Off-target engine. A draft with no key predates the setting, so
        # it describes a run made under BLAST — the same reasoning as `primer3`
        # above. Restoring it as the current default would change the run the
        # draft claims to reproduce.
        _set_combo_data(self.screen_engine,
                        d.get("screen_engine", ENGINE_BLAST))
        # 4b'. Specificity rules. After the engine, because the engine decides
        # which counting modes are available; a missing key is a draft from
        # before these controls, built from the engine defaults.
        _set_combo_data(self.ot_rule, d.get("ot_rule", _OT_RULE_DEFAULT))
        _set_combo_data(self.ot_counting,
                        d.get("ot_counting", _OT_COUNT_DEFAULT))
        _set_combo_data(self.ref_strand,
                        d.get("ref_strand", _SR_DEFAULT.reference_strand))
        # 4c. Junction gate. A draft written before the mode existed carries
        # only the boolean, which is exactly "hard or off" — it can never have
        # meant soft, so mapping it that way reproduces the run it describes.
        # Applied before the scalar fields so nothing can overwrite it, and
        # after the modality, whose preset also sets this combo.
        if "tn_junction_mode" in d:
            self._set_tn_mode(d["tn_junction_mode"] or TN_MODE_STRICT_A)
        else:
            self._set_tn_mode(TN_MODE_STRICT_A if d.get("tn_junction", True)
                              else TN_MODE_OFF)
        # Soft is Flex-only, so a soft draft opened against Visium must be
        # caught here rather than becoming a config that cannot run.
        from flexprobe.chemistry import get_chemistry
        from .modality import MODALITY_PRESETS
        self._refresh_tn_mode_availability(
            get_chemistry(MODALITY_PRESETS[self.modality.currentText()]))

        # 5. Workflow combos.
        _set_combo_text(self.run_mode, d.get("run_mode", ""))
        _set_combo_data(self.project_combo, d.get("project_id", ""))
        # 6. Run identity + IO (skip the run name / FASTA in bulk mode).
        if not self._is_bulk:
            # FASTA first: changing it clears a restored run name (see
            # _on_fasta_changed_for_run_name), so the name is restored after
            # it, against the FASTA it belongs to.
            if d.get("fasta_path"):
                self.fasta_edit.setText(d["fasta_path"])
            self._restore_run_name(d.get("run_name", ""), d.get("fasta_path", ""))
        # The output folder is deliberately NOT restored: Re-run and Load draft
        # keep the current default (decided 2026-09-15). A draft's folder is
        # where an old run went, and restoring it silently overrode a default
        # the user had changed since. to_draft still records it, for the record.
        # 7. All scalar spin/check fields via the single shared spec — filters
        #    are applied AFTER modality (step 1) so they override its defaults.
        #    Spin fields restore only if present (keeping modality defaults);
        #    checks always restore with the spec default (so no drift, e.g.
        #    primer3 stays OFF for a legacy draft — audit #12).
        self._scalar_from_draft(d)
        # 8. Ranking weights (advanced, opt-in).
        for k, val in (d.get("ranking_weights") or {}).items():
            if k in self._ranking_spins:
                self._ranking_spins[k].setValue(float(val))
        # 9. BLAST free-text fields.
        self.db_edit.setPlainText(d.get("blast_dbs", ""))
        self.cached_hits.setText(d.get("cached_hits", ""))
        self.extra_on_target.setPlainText(d.get("extra_on_target", ""))

    # ---- Run name: restored names follow the FASTA, typed names stay ----

    def _restore_run_name(self, saved: str, saved_fasta: str) -> None:
        """Put a draft's run name back without letting it go stale.

        Drafts and records written before this fix hold the fully generated
        name (``<stem>__<strategy tag>__<record id>``), so the suffix is cut
        first. What is left is either a name the user typed or, very often,
        just the old FASTA's stem, which build_runs derives anyway; the stem
        is dropped so the box shows its placeholder and the name re-derives
        from whatever FASTA the run ends up using.
        """
        name = strip_generated_run_name(saved)
        stem = os.path.splitext(os.path.basename((saved_fasta or "").strip()))[0]
        if name and name == stem:
            name = ""
        self.run_name.setText(name)
        self._run_name_restored = bool(name)
        self._run_name_fasta = self.fasta_edit.text().strip()

    def _on_run_name_edited(self, _text: str) -> None:
        # textEdited fires for user edits only, never for setText, so a name
        # typed in this dialog is never cleared by a FASTA change.
        self._run_name_restored = False

    def _on_fasta_changed_for_run_name(self, path: str) -> None:
        """Clear a restored run name once the target FASTA changes.

        A restored name describes the draft's input. Keeping it for another
        FASTA is how a Large-segment run came to be named after the Small one.
        """
        if (getattr(self, "_run_name_restored", False)
                and path.strip() != getattr(self, "_run_name_fasta", "")):
            self.run_name.clear()
            self._run_name_restored = False

    def _on_save_draft(self) -> None:
        from PySide6.QtWidgets import QInputDialog
        suggested = self.run_name.text().strip() or self.modality.currentText()
        name, ok = QInputDialog.getText(
            self, "Save query draft", "Draft name:",
            text=suggested,
        )
        if not ok or not name.strip():
            return
        existing = self._drafts.list()
        if name.strip() in existing:
            ans = QMessageBox.question(
                self, "Overwrite draft",
                f"A draft named '{name.strip()}' already exists. Overwrite it?",
                QMessageBox.Yes | QMessageBox.No,
            )
            if ans != QMessageBox.Yes:
                return
        try:
            self._drafts.save(name.strip(), self.to_draft())
        except Exception as exc:
            QMessageBox.warning(self, "Save draft", str(exc))
            return
        QMessageBox.information(
            self, "Draft saved",
            f"Saved '{name.strip()}'. Reload it any time with “Load draft…”.",
        )

    def _on_load_draft(self) -> None:
        from PySide6.QtWidgets import QInputDialog
        names = self._drafts.list()
        if not names:
            QMessageBox.information(
                self, "Load draft",
                "No saved drafts yet. Configure a query and use “Save draft…”.",
            )
            return
        name, ok = QInputDialog.getItem(
            self, "Load query draft", "Draft:", names, 0, False,
        )
        if not ok or not name:
            return
        payload = self._drafts.load(name)
        if payload is None:
            QMessageBox.warning(self, "Load draft", f"Could not read draft '{name}'.")
            return
        self.apply_draft(payload)

