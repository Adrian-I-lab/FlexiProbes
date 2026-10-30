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

"""Modality presets + per-chemistry GUI defaults.

Selecting a modality in the New-Query dialog auto-populates every field
in ``MODALITY_DEFAULTS[label]``. Values come from:
  • CG000621 Rev D — Visium spatial + Chromium Flex v1
  • CG000839 Rev A — GEM-X Flex v2

The chemistry-specific handle sequences (LHS 5' prefix, RHS 3' suffix,
Probe Barcode table) live in ``flexprobe.chemistry`` and are looked up by
the pipeline at run time.
"""
from __future__ import annotations

from dataclasses import asdict
from typing import Any, Dict

from flexprobe.models import FilterThresholds, ModalityDefinition, RunConfig


MODALITY_PRESETS: Dict[str, ModalityDefinition] = {
    "Flex v2 (10X GEM-X)":   ModalityDefinition(name="flex",   version="v2", spatial=False),
    "Flex v1 (Chromium)":    ModalityDefinition(name="flex",   version="v1", spatial=False),
    "Visium spatial":        ModalityDefinition(name="visium", version="v1", spatial=True),
}

DEFAULT_MODALITY = "Flex v2 (10X GEM-X)"


# --- defaults -------------------------------------------------------------
#
# These used to be three hand-written literals, and all 21 of their engine-owned
# values (7 filter thresholds x 3 modalities) were verbatim copies of the
# engine's own defaults. Copies drift, and twice they did: max_n_per_half sat at
# 0 while FilterThresholds said 1, and the conserved-degenerate default did the
# same (3c557e1). Both were silent — the app designed differently from the CLI
# and nothing in the output said so.
#
# So the values are now READ from the engine dataclasses rather than restated.
# There is nothing left to keep in step by hand: a new FilterThresholds field
# arrives here automatically, with the engine's value.
# RunConfig's two path arguments are required but irrelevant to its defaults;
# any placeholder yields the same values.
_ENGINE_RUN = RunConfig(fasta_path="", output_dir="")

_ENGINE_DEFAULTS: Dict[str, Any] = {
    # Every filter threshold, live from the engine's own dataclass.
    **asdict(FilterThresholds()),
    # The remaining form fields the modality preset drives are not filter
    # thresholds, so they come from RunConfig — same principle, different
    # dataclass.
    "tile_size": _ENGINE_RUN.tile_length,
    "tile_step": _ENGINE_RUN.tile_step,
    "enable_primer3": _ENGINE_RUN.primer3_scoring,
    "multiplex": _ENGINE_RUN.multiplex,
}


def _modality_defaults(*, probe_barcode: str,
                       **chemistry_overrides: Any) -> Dict[str, Any]:
    """Engine defaults, with room for a genuine per-chemistry difference.

    ``chemistry_overrides`` is EMPTY for every modality below except Flex v1's
    singleplex flag, and that emptiness is the whole point of this function:
    it is the audit's finding, stated in code. Nothing in CG000621 or CG000839
    asks Flex v1, Flex v2 or Visium to filter differently from one another, so
    nothing here should override a filter threshold.

    If a future tech note DOES — say Visium needs a wider GC band — add it as a
    keyword here and the divergence is one visible line with a citation next to
    it, rather than a re-duplication of all 21 values that hides which one
    actually differs.

    ``probe_barcode`` is separate because it is a chemistry FACT, not an engine
    default: only Flex v1 embeds a per-sample Probe Barcode, and the engine has
    no opinion about which one a preset should open on.
    """
    return {**_ENGINE_DEFAULTS, "probe_barcode": probe_barcode,
            **chemistry_overrides}


# Keys the dialog reads out of each entry (via ``.get``) are a subset of what
# ``_ENGINE_DEFAULTS`` supplies; the extra FilterThresholds fields ride along
# harmlessly and are there for the next form field that needs one. One of them
# needs care before it is wired up: ``tn_junction_mode`` arrives as the engine's
# "auto", which the dialog's three-way picker cannot represent — the picker is
# driven from ``require_tn_junction`` instead, exactly as resolve_tn_mode would
# resolve "auto" anyway.
MODALITY_DEFAULTS: Dict[str, Dict[str, Any]] = {
    # Source: CG000839 Rev A, p.2. TN junction is "preferred" per CG000839 and
    # "must be" per CG000621 — the engine's require_tn_junction=True is the
    # conservative reading of both, so neither modality overrides it.
    "Flex v2 (10X GEM-X)": _modality_defaults(
        probe_barcode="",           # Flex v2 has no embedded Probe Barcode
    ),
    # Source: CG000621 Rev D, Chromium Flex section (p.11+).
    "Flex v1 (Chromium)": _modality_defaults(
        probe_barcode="BC001",
        # The one real chemistry difference in this table: Flex v1 embeds a
        # per-sample Probe Barcode, so its preset opens singleplex on BC001.
        multiplex=False,
    ),
    # Source: CG000621 Rev D, Visium section (p.2+).
    "Visium spatial": _modality_defaults(
        probe_barcode="",
        # `multiplex` is inherited but ignored — Visium has no singleplex
        # variant, so the dialog disables the picker for it.
    ),
}
