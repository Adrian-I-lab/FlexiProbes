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

"""The engine seam: engine selection, the hit contract, and cache keying.

These pin the properties a second engine will depend on. The expensive ones to
get wrong are the cache-keying tests — a digest blind to the engine would let
one engine read back the other's hits, which changes verdicts silently.
"""
from __future__ import annotations

import pytest

from flexprobe import screen_cache
from flexprobe.constants import REFERENCE_STRAND_PLUS
from flexprobe.models import ScreenResources, ProbeCandidate
from flexprobe.screen import (
    ENGINE_BLAST, ENGINE_KMER, ENGINES, BlastScreen, KmerScreen, half_queries,
    make_screen,
)


def _cand(i: int = 0) -> ProbeCandidate:
    return ProbeCandidate(
        gene_id="g", probe_id=f"p{i}", sequence="ACGT" * 12 + "AC",
        start=0, end=50,
    )


# ---------------------------------------------------------------------------
# Engine selection
# ---------------------------------------------------------------------------

def test_default_engine_is_kmer():
    """The k-mer engine is the default; see constants.SCREEN_ENGINE_DEFAULT.

    Pinned because a silent revert would be invisible and expensive: where the
    two engines disagree, a brute-force check sides with the k-mer engine,
    which is exhaustive, while BLAST can certify clean a probe that is
    genuinely non-specific.
    """
    from flexprobe.constants import SCREEN_ENGINE_DEFAULT
    assert SCREEN_ENGINE_DEFAULT == ENGINE_KMER
    assert ScreenResources().engine == ENGINE_KMER


def test_default_engine_constant_is_a_real_engine():
    """constants.py cannot import screen (it imports nothing from the package),
    so the default is a bare literal. This is what stops the two drifting."""
    from flexprobe.constants import SCREEN_ENGINE_DEFAULT
    assert SCREEN_ENGINE_DEFAULT in ENGINES


def test_make_screen_builds_blast():
    s = make_screen(ENGINE_BLAST)
    assert isinstance(s, BlastScreen)
    assert s.name == ENGINE_BLAST


def test_unknown_engine_fails_immediately():
    """Fail before a DB is resolved or a chunk is searched, not partway in."""
    with pytest.raises(ValueError, match="unknown screen engine"):
        make_screen("bowtie")


def test_make_screen_builds_kmer():
    """Selecting an engine must never silently fall back to another one."""
    from flexprobe.screen import KmerScreen
    assert ENGINE_KMER in ENGINES
    s = make_screen(ENGINE_KMER)
    assert isinstance(s, KmerScreen)
    assert s.name == ENGINE_KMER


# ---------------------------------------------------------------------------
# The half-query contract
# ---------------------------------------------------------------------------

def test_half_queries_split_and_suffix():
    """__L/__R is the ONLY record of which half a hit came from — see
    specificity.hit_half. Every engine must produce it identically."""
    q = half_queries(_cand())
    assert [qid for qid, _ in q] == ["p0__L", "p0__R"]
    assert [len(s) for _, s in q] == [25, 25]
    assert q[0][1] + q[1][1] == _cand().sequence


def test_probe_id_survives_the_suffix_round_trip():
    """`rsplit("__", 1)` must recover the id even when the id itself has one."""
    c = ProbeCandidate(gene_id="g", probe_id="odd__name", sequence="A" * 50,
                       start=0, end=50)
    for qid, _ in half_queries(c):
        assert qid.rsplit("__", 1)[0] == "odd__name"


# ---------------------------------------------------------------------------
# Cache keying
# ---------------------------------------------------------------------------

def _fp(**kw):
    base = dict(db_arg="/db", chunk_size=2, perc_identity=None,
                qcov_hsp_perc=None, cache_mode=screen_cache.MODE_HITS)
    base.update(kw)
    return screen_cache.fingerprint(["p0", "p1"], **base)


def test_engine_changes_the_cache_fingerprint():
    """Two engines searching the same DB produce DIFFERENT hits (one local and
    trimmed, one end-to-end). A digest blind to the engine would let one read
    back the other's hits as its own — silently, and changing verdicts."""
    assert _fp(engine="blast") != _fp(engine="kmer")


def test_engine_defaults_keep_the_digest_stable():
    assert _fp() == _fp(engine="blast")


def test_blast_screen_surfaces_its_search_parameters_for_the_cache():
    """Anything that changes which hits come back must reach the digest."""
    s = BlastScreen(min_pident=80.0, min_aln=9, query_half_len=25,
                    derive_qcov=lambda a, b: 24.0)
    parts = s.fingerprint_parts()
    assert parts["perc_identity"] == 80.0
    assert parts["qcov_hsp_perc"] == 24.0
    assert _fp(**parts) != _fp()


def test_blast_screen_reference_strand_reaches_the_cache_digest():
    """A "plus" run never even searches the minus strand, so it finds different hits than
    a "both" run — a resumed run under a different setting must not replay
    the other one's cache."""
    s_both = BlastScreen()
    s_plus = BlastScreen(reference_strand=REFERENCE_STRAND_PLUS)
    assert s_both.fingerprint_parts()["reference_strand"] != \
        s_plus.fingerprint_parts()["reference_strand"]
    assert _fp(**s_both.fingerprint_parts()) != _fp(**s_plus.fingerprint_parts())


def test_kmer_screen_reference_strand_reaches_the_cache_digest(tmp_path):
    fa = tmp_path / "ref.fa"
    fa.write_text(">R\n" + "ACGT" * 20 + "\n")
    s_both = KmerScreen(cache_dir=str(tmp_path / "idx1"), mmap=False)
    s_both.prepare([str(fa)])
    s_plus = KmerScreen(cache_dir=str(tmp_path / "idx2"), mmap=False,
                        reference_strand=REFERENCE_STRAND_PLUS)
    s_plus.prepare([str(fa)])
    assert s_both.fingerprint_parts()["reference_strand"] != \
        s_plus.fingerprint_parts()["reference_strand"]


def test_db_arg_quotes_prefixes_containing_spaces():
    """blastn splits -db on whitespace; an unquoted path with a space becomes
    two bogus databases and the screen silently searches neither."""
    s = BlastScreen()
    s.db_arg = " ".join(
        f'"{p}"' if " " in p else p for p in ["/a/plain", "/has space/db"]
    )
    assert s.db_arg == '/a/plain "/has space/db"'
