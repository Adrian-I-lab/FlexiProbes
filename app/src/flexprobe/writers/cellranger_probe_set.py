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

"""Cell Ranger / Space Ranger probe-set CSV writer (Flex v1/v2 + spatial).

This module is the single source of truth for the CSV's shape (the '#'
metadata header, the `included` column, and the transcript-accession
gene_id check) — `select_export.py`'s two writers import from here so all
three call sites (pipeline export, per-run selection export, project bulk
export) stay in lock-step. The format here was checked against the eight
public 10x Chromium Flex vendor probe-set CSVs (Human/Mouse/Rat), which are
not redistributed with this repository.
"""
from __future__ import annotations
import csv
import logging
import os
import re
import tempfile
from typing import Iterable, List, Optional
from .. import __version__
from ..constants import DEGENERATE_BASES, PROBE_ROLE_CONSENSUS
from ..models import ProbePair
from ..errors import WriterError

logger = logging.getLogger("flexprobe.writers.cellranger_probe_set")


# ---------------------------------------------------------------------------
# '#' metadata header
# ---------------------------------------------------------------------------
# A Cell Ranger / Space Ranger probe-set CSV is not just a table: every real
# 10x vendor file opens with five '#'-prefixed metadata lines before the
# column header. Confirmed against all eight public Chromium Flex vendor
# probe-set CSVs (Human/Mouse/Rat, probe-set file format versions 1.0/2.0/3.0)
# — every one of them opens with exactly:
#   #probe_set_file_format=<spec version>
#   #panel_name=<panel identity>
#   #panel_type=predesigned|custom
#   #reference_genome=<genome/build the probes were designed against>
#   #reference_version=<annotation build/version>
#
# flexiprobes designs directly against an arbitrary user-supplied FASTA, not
# against a curated 10x whole-genome reference — there is no genome build to
# report. So 'reference_genome'/'reference_version' are populated from what
# is actually true of a flexiprobes run (the input FASTA's identity and the
# flexiprobes version that designed against it) rather than a genome-build
# string this tool has no way to know.
#
# The column set below (gene_id, probe_seq, probe_id, included — no region,
# no gene_name) matches vendor format 1.0 exactly, so that is the format
# version declared.
CR_PROBE_SET_FILE_FORMAT = "1.0"
CR_PANEL_TYPE = "custom"


def cr_header_lines(
    panel_name: str,
    reference_genome: str,
    reference_version: Optional[str] = None,
) -> List[str]:
    """The five '#' metadata lines a Cell Ranger probe-set CSV must open with."""
    return [
        f"#probe_set_file_format={CR_PROBE_SET_FILE_FORMAT}",
        f"#panel_name={panel_name}",
        f"#panel_type={CR_PANEL_TYPE}",
        f"#reference_genome={reference_genome}",
        f"#reference_version={reference_version or f'flexiprobes-{__version__}'}",
    ]


def default_panel_identity(fallback: str) -> str:
    """A ``fallback`` panel identity string, cleaned up for the header line.

    Used when a caller has no more specific identity (e.g. an input FASTA
    basename) to hand in.
    """
    return fallback or "flexiprobes custom probe set"


# ---------------------------------------------------------------------------
# gene_id: warn rather than silently pass through a transcript accession
# ---------------------------------------------------------------------------
# A probe's gene_id column is only as good as the input FASTA record id it
# came from (io_fasta.iter_fasta takes the header's first whitespace token
# verbatim — see its docstring for why it is not translated or truncated).
# That is correct for the common case: a viral/custom target name, or an
# actual Ensembl gene id. It is wrong when a user hands the pipeline a
# transcript-level FASTA (ENST.../NM_.../XM_...) — Cell Ranger's gene_id
# column expects a GENE, and nothing in this pipeline carries a GTF or other
# gene model that could translate a transcript accession to its parent gene
# (confirmed: no exon/intron/gene-model input exists anywhere in
# app/src/flexprobe/). Rather than silently emit a wrong-shaped id or invent
# a translation this pipeline cannot actually perform, this
# warns once per export so the person ordering the panel can see the problem.
_TRANSCRIPT_ACCESSION_RE = re.compile(
    r'^(ENS[A-Z]*T\d+(\.\d+)?|[NXY][MR]_\d+(\.\d+)?)$'
)


def warn_if_transcript_like_gene_ids(gene_ids: Iterable[str]) -> None:
    """Log a warning if any ``gene_ids`` look like transcript accessions."""
    bad = sorted({g for g in gene_ids if g and _TRANSCRIPT_ACCESSION_RE.match(g)})
    if bad:
        logger.warning(
            "Cell Ranger probe-set gene_id column expects a gene identifier, "
            "but %d input record id(s) look like transcript accessions "
            "(e.g. %r) rather than gene ids. flexiprobes has no gene model to "
            "translate a transcript accession to its parent gene, so these "
            "are written through unchanged; Cell Ranger will treat each "
            "transcript as its own 'gene'. Use a gene-level FASTA header, or "
            "map these accessions to their parent gene before ordering.",
            len(bad), bad[0],
        )


def _probe_id_schema(gene_id: str, index: int, version: str = "v1") -> str:
    """Generate probe_id per 10x schema."""
    return f"{gene_id}|{version}|{index:04d}"


def write_cellranger_probe_set(
    pairs: List[ProbePair],
    output_path: str,
    version: str = "v1",
    spatial: bool = False,
    overwrite: bool = False,
    panel_name: Optional[str] = None,
    reference_genome: Optional[str] = None,
) -> str:
    """Write a Cell Ranger / Space Ranger custom probe-set CSV.

    Matches real 10x vendor probe-set format 1.0 (the minimal variant
    observed on disk — see module docstring): five '#' metadata lines, then
    columns ``gene_id, probe_seq, probe_id, included``.

    ``panel_name``/``reference_genome`` should be the caller's best
    identification of this design run (e.g. the input FASTA's basename);
    generic fallbacks are used when not supplied (direct/test use). ``spatial``
    does not change the column set: flexiprobes has no exon/intron model to
    populate a trustworthy ``region`` (spliced/unspliced) value for either
    chemistry, so that column is omitted rather than filled with a value that
    would be wrong for most windows. It is kept as a parameter for API
    stability and potential future spatial-specific metadata.
    """
    if os.path.exists(output_path) and not overwrite:
        raise WriterError(f"Output exists and overwrite=False: {output_path}")

    panel_name = default_panel_identity(panel_name or "")
    reference_genome = reference_genome or "input FASTA (no whole-genome reference used)"

    dir_name = os.path.dirname(output_path) or '.'
    os.makedirs(dir_name, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=dir_name, suffix='.tmp')
    try:
        with os.fdopen(fd, 'w', newline='') as fh:
            for line in cr_header_lines(panel_name, reference_genome):
                fh.write(line + "\n")
            writer = csv.writer(fh)
            writer.writerow(['gene_id', 'probe_seq', 'probe_id', 'included'])
            gene_ids_written: List[str] = []
            for idx, pair in enumerate(pairs):
                pid = _probe_id_schema(pair.lhs.gene_id, idx, version)
                # Cell Ranger's probe_seq is the 50-bp ligated probe sequence
                # (RC of the mRNA target window) — NOT the full LHS oligo
                # with its 5' handle.
                probe_seq = (
                    pair.lhs.probe_rc
                    or (pair.lhs.target_lhs + pair.lhs.target_rhs)
                    or pair.lhs.sequence
                )
                # A Cell Ranger probe set is read by the aligner as literal
                # nucleotides; a mixed-base string there is not a permissive
                # probe, it is an invalid record. Unlike the IDT sheet there is
                # no handle here, so probe_seq must be pure ACGT and both the
                # role and the characters are worth asserting. This is the last
                # line — the pipeline never pairs a parent and select_export
                # substitutes one — so reaching it means an upstream regression.
                if (pair.lhs.probe_role == PROBE_ROLE_CONSENSUS
                        or DEGENERATE_BASES & set(probe_seq.upper())):
                    raise WriterError(
                        f"Degenerate probe {pair.lhs.probe_id!r} reached the "
                        f"Cell Ranger probe-set writer. Consensus probes are "
                        f"never orderable — write their variants instead."
                    )
                gene_ids_written.append(pair.lhs.gene_id)
                writer.writerow([pair.lhs.gene_id, probe_seq, pid, 'TRUE'])
            warn_if_transcript_like_gene_ids(gene_ids_written)
        os.replace(tmp_path, output_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise
    return output_path
