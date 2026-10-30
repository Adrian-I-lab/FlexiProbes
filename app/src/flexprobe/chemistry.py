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

"""Per-modality probe chemistry: handle sequences + LHS/RHS oligo assembly.

For a 50-bp mRNA target window (sense strand, as it appears in the input
FASTA), the LHS and RHS probes are built as follows:

    mRNA target (sense, 5'->3'):   [-- 5' half 25 bp --|-- 3' half 25 bp --]

    target_RHS  =  RC(mRNA 5' half)   # 25 bp — RHS probe binds the mRNA 5' end
    target_LHS  =  RC(mRNA 3' half)   # 25 bp — LHS probe binds the mRNA 3' end

    LHS oligo (5'->3') = lhs_5prime_handle + target_LHS
    RHS oligo (5'->3') = "/5Phos/"        + target_RHS + rhs_3prime_handle

The directional convention (LHS sits on the mRNA 3' end, RHS on the 5' end)
was verified against the EGFP example in CG000839 Tables 4 and 5: the EGFP
target sequence "AGCGTGCAGCTCGCCGACCACTACC" appears immediately downstream
(3') of "GATCCGCCACAACATCGAGGACGGC" in EGFP mRNA, and RC of the former is
target_LHS-1 while RC of the latter is target_RHS-1.

Handle sequences and barcode tables come from:
  • CG000621 Rev D — Visium spatial + Chromium Flex v1
  • CG000839 Rev A — GEM-X Flex v2
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, Optional

from .constants import TN_MODE_AUTO, TN_MODE_OFF, TN_MODE_SOFT, TN_MODE_STRICT_A
from .errors import ConfigError
from .models import ModalityDefinition
from .tiling import reverse_complement


# ---------------------------------------------------------------------------
# Flex v1 Probe Barcodes
# ---------------------------------------------------------------------------

# Flex v1 RHS Probe Barcodes (CG000621 Rev D, Table 8). BC001 = singleplex.
FLEX_V1_PROBE_BARCODES: Dict[str, str] = {
    "BC001": "ACTTTAGG",
    "BC002": "AACGGGAA",
    "BC003": "AGTAGGCT",
    "BC004": "ATGTTGAC",
    "BC005": "ACAGACCT",
    "BC006": "ATCCCAAC",
    "BC007": "AAGTAGAG",
    "BC008": "AGCTGTGA",
    "BC009": "ACAGTCTG",
    "BC010": "AGTGAGTG",
    "BC011": "AGAGGCAA",
    "BC012": "ACTACTCA",
    "BC013": "ATACGTCA",
    "BC014": "ATCATGTG",
    "BC015": "AACGCCGA",
    "BC016": "ATTCGGTT",
}


# Placeholder kept in the assembled RHS oligo when more than one Probe
# Barcode is selected (Flex v1 multiplex). The select_export step expands
# each selected probe into one row per barcode at order time.
BARCODE_PLACEHOLDER = "{BC}"


# ---------------------------------------------------------------------------
# Region metadata for the GUI handle preview.
# ---------------------------------------------------------------------------
#
# Each ChemistryConfig exposes ``lhs_regions()`` / ``rhs_regions()`` which
# return the LHS/RHS oligo broken down into labelled, colour-coded chunks
# matching the schematics in CG000621 / CG000839. The GUI concatenates them
# into a coloured HTML preview.

# A consistent colour palette across the chemistries. Background tints
# are chosen for readability against both light and dark Qt themes; matching
# text colours keep the labels under the swatches legible.
REGION_COLORS = {
    "phosphate":     "#bbbbbb",  # /5Phos/
    "read_2s":       "#1f77b4",  # Read 2S (Flex v1/v2/Visium LHS handle)
    "target_lhs":    "#2ca02c",  # target_LHS — the 25 bp probe target
    "target_rhs":    "#2ca02c",  # target_RHS — sense same colour as target_lhs
    "constant":      "#9467bd",  # Constant Sequence (Flex v1) / partial Constant (Flex v2 mux)
    "nn":            "#cccccc",  # NN spacer in Flex v1 RHS
    "barcode":       "#ff7f0e",  # Probe Barcode (Flex v1)
    "pcs1":          "#d62728",  # partial Capture Sequence 1
    "polya":         "#bcbd22",  # Visium 30 × A tail
}


@dataclass(frozen=True)
class HandleRegion:
    """One coloured chunk of an LHS or RHS oligo for the handle preview."""
    label: str           # short label (e.g. "pRead 2S", "Probe Barcode (BC003)")
    sequence: str        # nt string (or placeholder like "[target_LHS 25 bp]")
    color_key: str       # key into ``REGION_COLORS``

    @property
    def color(self) -> str:
        return REGION_COLORS.get(self.color_key, "#999")


@dataclass(frozen=True)
class ChemistryConfig:
    """Per-modality probe-assembly rules.

    Attributes
    ----------
    key : str
        Short canonical key ("flex_v1", "flex_v2", "visium").
    label : str
        Human-readable name used in CSVs / manifests.
    lhs_5prime : str
        Constant sequence prepended to ``target_LHS`` to form the LHS oligo.
    rhs_3prime_multiplex : str
        Constant sequence appended to ``target_RHS`` for multiplex runs.
        May contain ``{probe_barcode}`` placeholder (Flex v1 only).
    rhs_3prime_singleplex : str | None
        Constant sequence appended for singleplex runs. If None, multiplex
        sequence is used for both.
    half_length : int
        Length of each target half (always 25 for current 10X chemistries).
    requires_tn_junction : bool
        Whether this chemistry's 10x tech note asks for the TN junction —
        whether as a requirement (CG000621: "must be a T") or as a
        recommendation (CG000839: "10x Genomics recommends using the T at the
        specified position"). True for all three current chemistries. Drives
        only the "you have relaxed the junction gate" warning; it does not
        select the gate.
    allows_soft_tn_junction : bool
        Whether ``TN_MODE_SOFT`` (reject only ``C``) may be selected for this
        chemistry. True for Flex v1/v2, whose relaxation CG000839 supports;
        False for Visium, where CG000621 p.3 states "The 25th nucleotide of the
        probe (3' most nucleotide of the LHS probe) must be a T. The opposing
        nucleotide in the target RNA must be an A." with no qualifier anywhere
        in either document. Asking for soft on Visium is a ConfigError rather
        than a silent downgrade.

    The gate itself is NOT chemistry-defaulted: the strict adenine is the
    default everywhere it is engaged (see constants.TN_MODE_*), and the soft
    mode is always an explicit opt-in.
    has_probe_barcode : bool
        Whether the RHS template uses a Probe Barcode (Flex v1 only).
    default_probe_barcode : str
        BC label used when has_probe_barcode and the user does not pick one.
    """

    key: str
    label: str
    lhs_5prime: str
    rhs_3prime_multiplex: str
    rhs_3prime_singleplex: Optional[str] = None
    half_length: int = 25
    requires_tn_junction: bool = True
    allows_soft_tn_junction: bool = False
    has_probe_barcode: bool = False
    default_probe_barcode: str = "BC001"

    @property
    def tile_length(self) -> int:
        return self.half_length * 2

    def _resolve_rhs_template(self, *, multiplex: bool) -> str:
        if multiplex or self.rhs_3prime_singleplex is None:
            return self.rhs_3prime_multiplex
        return self.rhs_3prime_singleplex

    def rhs_template_with_placeholder(self, *, multiplex: bool) -> str:
        """Return the RHS 3' template with the ``{BC}`` placeholder in place
        of the embedded Probe Barcode (Flex v1).

        For chemistries that don't carry a barcode (Flex v2, Visium), this
        is identical to ``_resolve_rhs_template``.
        """
        template = self._resolve_rhs_template(multiplex=multiplex)
        if self.has_probe_barcode:
            template = template.replace("{probe_barcode}", BARCODE_PLACEHOLDER)
        return template

    def split_tile(self, mrna_tile: str) -> tuple[str, str, str, str]:
        """Return (mrna_5half, mrna_3half, target_LHS, target_RHS)."""
        if len(mrna_tile) != self.tile_length:
            raise ConfigError(
                f"Expected {self.tile_length}-bp tile, got {len(mrna_tile)} bp"
            )
        mrna_5half = mrna_tile[:self.half_length].upper()
        mrna_3half = mrna_tile[self.half_length:].upper()
        # LHS probe binds the mRNA 3' end; RHS binds the 5' end. See module
        # docstring for the EGFP verification of this convention.
        target_lhs = reverse_complement(mrna_3half)
        target_rhs = reverse_complement(mrna_5half)
        return mrna_5half, mrna_3half, target_lhs, target_rhs

    def build_lhs_oligo(self, target_lhs: str) -> str:
        return f"{self.lhs_5prime}{target_lhs}"

    def build_rhs_oligo(
        self,
        target_rhs: str,
        *,
        multiplex: bool,
        probe_barcode: Optional[str] = None,
        leave_placeholder: bool = False,
    ) -> str:
        """Build an RHS oligo string.

        When ``leave_placeholder=True`` and the chemistry carries a barcode,
        the embedded barcode position is left as ``{BC}`` so a single probe
        record can be expanded into multiple per-barcode rows at order time.
        """
        template = self._resolve_rhs_template(multiplex=multiplex)
        if self.has_probe_barcode:
            if leave_placeholder:
                template = template.replace("{probe_barcode}", BARCODE_PLACEHOLDER)
            else:
                bc_label = probe_barcode or self.default_probe_barcode
                bc_seq = self._lookup_probe_barcode(bc_label)
                template = template.replace("{probe_barcode}", bc_seq)
        return f"/5Phos/{target_rhs}{template}"

    def lhs_regions(self) -> "list[HandleRegion]":
        """Return the LHS oligo broken into labelled, coloured chunks."""
        return [
            HandleRegion(
                label="pRead 2S (partial Read 2S handle)",
                sequence=self.lhs_5prime,
                color_key="read_2s",
            ),
            HandleRegion(
                label="target_LHS (25 bp RC of mRNA 3' half)",
                sequence=f"[target_LHS · {self.half_length} bp]",
                color_key="target_lhs",
            ),
        ]

    def rhs_regions(
        self,
        *,
        multiplex: bool,
        probe_barcode: Optional[str] = None,
        leave_placeholder: bool = False,
    ) -> "list[HandleRegion]":
        """Return the RHS oligo broken into labelled, coloured chunks.

        Per-chemistry segmentation:
          • Flex v1 RHS = /5Phos/ + target_RHS + Constant Seq + NN + BC + pCS1
          • Flex v2 multiplex RHS = /5Phos/ + target_RHS + partial Constant
          • Flex v2 singleplex RHS = /5Phos/ + target_RHS + pCS1
          • Visium RHS = /5Phos/ + target_RHS + 30 × A poly-A tail
        """
        out: list = [
            HandleRegion(label="/5Phos/ (5' phosphate)",
                         sequence="/5Phos/", color_key="phosphate"),
            HandleRegion(
                label="target_RHS (25 bp RC of mRNA 5' half)",
                sequence=f"[target_RHS · {self.half_length} bp]",
                color_key="target_rhs",
            ),
        ]
        if self.key == "flex_v1":
            # /5Phos/-target_RHS-ACGCGGTTAGCACGTA-NN-{BC}-CGGTCCTAGCAA
            if leave_placeholder:
                bc_label = "BC"
                bc_seq = BARCODE_PLACEHOLDER
            else:
                lab = probe_barcode or self.default_probe_barcode
                bc_label = lab
                bc_seq = self._lookup_probe_barcode(lab) if lab else BARCODE_PLACEHOLDER
            out.extend([
                HandleRegion(label="Constant Sequence (16 bp)",
                             sequence="ACGCGGTTAGCACGTA", color_key="constant"),
                HandleRegion(label="NN (2 bp quality-balancing spacer)",
                             sequence="NN", color_key="nn"),
                HandleRegion(label=f"Probe Barcode ({bc_label})",
                             sequence=bc_seq, color_key="barcode"),
                HandleRegion(label="pCS1 (partial Capture Seq 1, 12 bp)",
                             sequence="CGGTCCTAGCAA", color_key="pcs1"),
            ])
            return out
        if self.key == "flex_v2":
            if multiplex:
                out.append(HandleRegion(
                    label="partial Constant Sequence (13 bp)",
                    sequence="CCCATATAAGAAA", color_key="constant",
                ))
            else:
                out.append(HandleRegion(
                    label="pCS1 (partial Capture Seq 1, 12 bp; 4-Samples Kit)",
                    sequence="CGGTCCTAGCAA", color_key="pcs1",
                ))
            return out
        if self.key == "visium":
            out.append(HandleRegion(
                label="poly(A) tail (30 × A)",
                sequence="A" * 30, color_key="polya",
            ))
            return out
        return out

    def _lookup_probe_barcode(self, label: str) -> str:
        # Allow callers to pass either the BC label ("BC003") or the raw 8-mer.
        if re.fullmatch(r"[ACGTacgt]{8}", label):
            return label.upper()
        table = FLEX_V1_PROBE_BARCODES
        if label not in table:
            raise ConfigError(
                f"Unknown Probe Barcode {label!r} for {self.label}. "
                f"Known: {sorted(table)}"
            )
        return table[label]


# ---------------------------------------------------------------------------
# Modality registry
# ---------------------------------------------------------------------------

# Common partial-Read-2S handle shared by Flex v1, Flex v2, and Visium LHS.
_LHS_READ2S = "CCTTGGCACCCGAGAATTCCA"

CHEMISTRIES: Dict[str, ChemistryConfig] = {
    # Source: CG000621 Rev D — Chromium Single Cell Gene Expression Flex.
    "flex_v1": ChemistryConfig(
        key="flex_v1",
        label="Flex v1 (Chromium)",
        lhs_5prime=_LHS_READ2S,
        # /5Phos/-target_RHS-ACGCGGTTAGCACGTA-NN-{BC}-CGGTCCTAGCAA
        rhs_3prime_multiplex="ACGCGGTTAGCACGTANN{probe_barcode}CGGTCCTAGCAA",
        rhs_3prime_singleplex="ACGCGGTTAGCACGTANN{probe_barcode}CGGTCCTAGCAA",
        # CG000621 p.12 states the junction T/A as a requirement, so the strict
        # adenine stays the default. The soft mode is *offered* here because
        # 10x's own Flex v1 panels (human v1.1.0, mouse v1.1.1) carry that A
        # only 51.3% / 52.2% of the time while avoiding C 99.0% / 99.3% of the
        # time — the "must" is demonstrably not what they designed to.
        requires_tn_junction=True,
        allows_soft_tn_junction=True,
        has_probe_barcode=True,
        default_probe_barcode="BC001",
    ),
    # Source: CG000839 Rev A — GEM-X Flex v2 (October 2025).
    "flex_v2": ChemistryConfig(
        key="flex_v2",
        label="Flex v2 (10X GEM-X)",
        lhs_5prime=_LHS_READ2S,
        # Multiplex: partial Constant Sequence
        rhs_3prime_multiplex="CCCATATAAGAAA",
        # Singleplex 4 Samples Kit: partial Capture Sequence 1 (pCS1)
        rhs_3prime_singleplex="CGGTCCTAGCAA",
        # CG000839 p.2, verbatim: "While preferred, the TN motif is not always
        # strictly required, and other nucleotide motifs can also function
        # effectively." Not a requirement for v2 — but the same document then
        # says "for a simplified and robust design process, 10x Genomics
        # recommends using the T at the specified position", so strict stays
        # the default, the relaxation is opt-in, and relaxing it still warns.
        requires_tn_junction=True,
        allows_soft_tn_junction=True,
        has_probe_barcode=False,
    ),
    # Source: CG000621 Rev D — Visium Spatial Gene Expression.
    "visium": ChemistryConfig(
        key="visium",
        label="Visium spatial",
        lhs_5prime=_LHS_READ2S,
        rhs_3prime_multiplex="A" * 30,            # poly-A tail
        rhs_3prime_singleplex=None,                # no distinction
        # CG000621 Rev D p.3 (Visium section): "The 25th nucleotide of the probe
        # (3' most nucleotide of the LHS probe) must be a T. The opposing
        # nucleotide in the target RNA must be an A." Neither document offers
        # any relaxation for spatial — CG000839's "not always strictly
        # required" sentence is inside the GEM-X Flex v2 note and says nothing
        # about Visium — so soft mode is not available here at all.
        requires_tn_junction=True,
        allows_soft_tn_junction=False,
        has_probe_barcode=False,
    ),
}


def resolve_tn_mode(filters, chemistry: ChemistryConfig) -> str:
    """Resolve the effective TN-junction gate for a run.

    The gate is HARD (``strict_a``, mRNA[25] == 'A') wherever it is engaged;
    ``soft`` (reject only 'C') must always be asked for explicitly, and only
    the Flex chemistries offer it.

    Precedence:
      1. An explicit ``filters.tn_junction_mode`` other than ``"auto"`` — the
         user named a gate, so that is the gate. It wins over
         ``require_tn_junction``, which DEFAULTS to True and is therefore a
         default rather than a statement of intent; treating the two as
         contradictory here would make the soft mode unreachable. The CLI,
         which *can* tell an explicitly-passed flag from a default, rejects
         ``--require-tn-junction --tn-junction-mode soft`` as a genuine
         contradiction before it ever reaches this function.
      2. ``"auto"`` + ``filters.require_tn_junction`` True (the default)
         -> ``strict_a``.
      3. ``"auto"`` + ``require_tn_junction`` False -> ``off``, i.e. no
         junction-identity gate.

    Raises ``ConfigError`` if ``soft`` is requested for a chemistry that does
    not allow it (Visium), rather than silently substituting a different gate.
    """
    mode = getattr(filters, "tn_junction_mode", TN_MODE_AUTO) or TN_MODE_AUTO
    if mode == TN_MODE_AUTO:
        mode = (TN_MODE_STRICT_A
                if getattr(filters, "require_tn_junction", True)
                else TN_MODE_OFF)
    if mode == TN_MODE_SOFT and not chemistry.allows_soft_tn_junction:
        raise ConfigError(
            f"tn_junction_mode='soft' is not available for {chemistry.label}. "
            f"CG000621 Rev D states the junction adenine as a requirement for "
            f"spatial probes ('The opposing nucleotide in the target RNA must "
            f"be an A') and neither tech note relaxes it for Visium. Use "
            f"'strict_a' (the default), or a Flex chemistry."
        )
    return mode


def get_chemistry(modality: ModalityDefinition) -> ChemistryConfig:
    """Return the ChemistryConfig matching a ModalityDefinition."""
    if modality.spatial:
        return CHEMISTRIES["visium"]
    if modality.name == "flex":
        if modality.version == "v1":
            return CHEMISTRIES["flex_v1"]
        if modality.version == "v2":
            return CHEMISTRIES["flex_v2"]
    raise ConfigError(
        f"No chemistry registered for modality "
        f"name={modality.name!r} version={modality.version!r}"
    )


# ---------------------------------------------------------------------------
# Convenience: a per-pair assembled view
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AssembledProbePair:
    """Everything you need to write a result row for one tile."""
    mrna_target: str          # 50 bp, sense (input FASTA strand)
    probe_rc: str             # 50 bp, RC(mrna_target) — what the probe pair contains
    target_lhs: str           # 25 bp
    target_rhs: str           # 25 bp
    lhs_oligo: str            # lhs_5prime + target_lhs
    rhs_oligo: str            # /5Phos/ + target_rhs + rhs_3prime (with barcode subbed in)
    tn_junction_ok: bool      # True if mrna_target[25] == 'A'


def assemble_pair(
    mrna_tile: str,
    chemistry: ChemistryConfig,
    *,
    multiplex: bool,
    probe_barcode: Optional[str] = None,
    leave_barcode_placeholder: bool = False,
) -> AssembledProbePair:
    """Assemble LHS + RHS oligos for one 50-bp mRNA tile.

    When ``leave_barcode_placeholder=True`` and the chemistry carries a
    Probe Barcode (Flex v1), the RHS oligo keeps the literal ``{BC}``
    placeholder so the downstream export step can expand it into one row
    per selected barcode. The other chemistries ignore the flag.
    """
    _m5, _m3, target_lhs, target_rhs = chemistry.split_tile(mrna_tile)
    lhs = chemistry.build_lhs_oligo(target_lhs)
    rhs = chemistry.build_rhs_oligo(
        target_rhs, multiplex=multiplex, probe_barcode=probe_barcode,
        leave_placeholder=leave_barcode_placeholder,
    )
    tn_ok = (
        len(mrna_tile) > chemistry.half_length
        and mrna_tile[chemistry.half_length].upper() == "A"
    )
    return AssembledProbePair(
        mrna_target=mrna_tile.upper(),
        probe_rc=reverse_complement(mrna_tile),
        target_lhs=target_lhs,
        target_rhs=target_rhs,
        lhs_oligo=lhs,
        rhs_oligo=rhs,
        tn_junction_ok=tn_ok,
    )

