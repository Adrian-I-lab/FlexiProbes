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

"""IDT oPool order-sheet writer.

Two rows per probe pair (LHS + RHS) in the same pool. Sequences are the full
oligos with the modality-specific handles already attached by the pipeline;
RHS sequences carry the /5Phos/ prefix.

Columns: Pool Name, Sequence Name, Sequence.
"""
from __future__ import annotations
import csv
import os
import tempfile
from typing import List, Optional, Sequence, Tuple

from ..chemistry import BARCODE_PLACEHOLDER
from ..constants import PROBE_ROLE_CONSENSUS
from ..models import ProbePair
from ..errors import WriterError


def _reject_degenerate(pairs: Sequence[ProbePair]) -> None:
    """Refuse to write a pair built from a degenerate consensus parent.

    The pipeline already skips parents when it assembles pairs, and
    ``select_export`` substitutes a selected parent with its variants. This is
    the last line: it makes that invariant a property of the WRITER, so a future
    change to either upstream path fails loudly here instead of quietly emitting
    a mixed-base oligo onto an order sheet — an order that costs real money and
    yields a probe that binds nothing in particular.

    Checked on ``probe_role`` rather than on the oligo characters, because the
    Flex v1 RHS handle legitimately contains a literal 'NN' spacer
    (``ACGCGGTTAGCACGTANN{BC}…``) and a character test would fire on every
    correct Flex v1 row.
    """
    for pair in pairs:
        if (pair.lhs.probe_role == PROBE_ROLE_CONSENSUS
                or pair.rhs.probe_role == PROBE_ROLE_CONSENSUS):
            raise WriterError(
                f"Degenerate consensus probe {pair.lhs.probe_id!r} reached the "
                f"IDT order writer. Consensus probes are never orderable — "
                f"order their concrete variants instead."
            )


def write_idt_opools(
    pairs: List[ProbePair],
    output_path: str,
    overwrite: bool = False,
    barcodes: Optional[Sequence[Tuple[str, str]]] = None,
) -> str:
    """Write the IDT oPool order sheet (Pool Name, Sequence Name, Sequence).

    Single-barcode (default): two rows per pair (LHS + RHS) in one pool per gene.

    Multi-barcode Flex v1 (``barcodes`` = ``[(label, 8-mer), ...]`` with >1
    entry): the RHS oligo carries a ``{BC}`` placeholder; each selected probe is
    emitted once per barcode (shared LHS + one RHS per barcode), grouped into one
    pool per sample/barcode and named after the barcode (``BC001``, ``BC008`` …)
    — mirroring ``select_export.export_selected_probes`` so the one-shot (CLI)
    path matches the discovery/review path.
    """
    if os.path.exists(output_path) and not overwrite:
        raise WriterError(f"Output exists and overwrite=False: {output_path}")
    _reject_degenerate(pairs)

    multi = barcodes is not None and len(barcodes) > 1

    dir_name = os.path.dirname(output_path) or '.'
    os.makedirs(dir_name, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=dir_name, suffix='.tmp')
    try:
        with os.fdopen(fd, 'w', newline='') as fh:
            writer = csv.writer(fh)
            writer.writerow(['Pool Name', 'Sequence Name', 'Sequence'])
            if multi:
                for bc_label, bc_seq in (barcodes or []):
                    pool = bc_label
                    for pair in pairs:
                        writer.writerow(
                            [pool, f"{pair.lhs.probe_id}_LHS", pair.lhs.sequence])
                    for pair in pairs:
                        rhs = pair.rhs.sequence.replace(BARCODE_PLACEHOLDER, bc_seq)
                        writer.writerow(
                            [pool, f"{pair.lhs.probe_id}_RHS_{bc_label}", rhs])
            else:
                for pair in pairs:
                    # Both halves share the same base probe_id; suffix is added
                    # exactly once here, fixing the previous "_RHS_RHS" duplicate.
                    base_id = pair.lhs.probe_id
                    pool = f"{pair.lhs.gene_id}_pool"
                    writer.writerow([pool, f"{base_id}_LHS", pair.lhs.sequence])
                    writer.writerow([pool, f"{base_id}_RHS", pair.rhs.sequence])
        os.replace(tmp_path, output_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise
    return output_path
