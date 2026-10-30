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

"""Sliding-window tiling and reverse-complement utilities."""
from __future__ import annotations
from typing import List, Tuple

_COMP = str.maketrans('ACGTacgtRYSWKMBDHVNryswkmbdhvn',
                       'TGCAtgcaYRSWMKVHDBNyrswmkvhdbn')


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of *seq*."""
    return seq.translate(_COMP)[::-1]


def tile_sequence(
    sequence: str,
    length: int = 50,
    step: int = 1,
) -> List[Tuple[int, int, str]]:
    """Generate (start, end, subseq) tiles across *sequence*.

    Parameters
    ----------
    sequence : str
        The full-length nucleotide sequence.
    length : int
        Tile (probe) length in bases.
    step : int
        Sliding step size.

    Returns
    -------
    list of (start, end, subsequence) tuples (0-based, half-open).
    """
    tiles: List[Tuple[int, int, str]] = []
    seq_len = len(sequence)
    if seq_len < length:
        return tiles
    for i in range(0, seq_len - length + 1, step):
        tiles.append((i, i + length, sequence[i:i + length]))
    return tiles
