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

"""Three ways the off-target screen could go quiet instead of wrong-loudly.

Each one is a place where the screen could return a verdict it had not earned
and say nothing about it:

(1) a single degenerate base must not remove a probe half from the k-mer
    search, or the PAIR could pass while a perfect off-target sat in the
    reference. Reachable with shipped defaults
    (``FilterThresholds.max_n_per_half == 1``).
(2) the junction-duplex gate exists only in the k-mer engine, so selecting
    ``--screen-engine blast`` with the gate on would disable a hard filter the
    config still reported as on. That combination is refused.
(3) on an Ensembl ``cdna.all`` reference no subject id yields a gene, so
    per-gene off-target scoping would silently become per-transcript, the
    permissive direction. That fallback is announced.
"""
from __future__ import annotations

import logging
import random

import pytest

from flexprobe.errors import ConfigError
from flexprobe.models import (
    ProbeCandidate, RunConfig, ScreenResources,
    SpecificityThresholds,
)
from flexprobe.screen import ENGINE_BLAST, ENGINE_KMER, KmerScreen
from flexprobe.specificity import (
    evaluate_specificity, hit_half, split_junction_records,
)

HALF = 25


def _rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _cand(seq, pid="p0", gene="g"):
    return ProbeCandidate(gene_id=gene, probe_id=pid, sequence=seq,
                          start=0, end=len(seq))


def _screen(tmp_path, fasta_text, **kw):
    fa = tmp_path / "ref.fa"
    fa.write_text(fasta_text)
    s = KmerScreen(cache_dir=str(tmp_path / "idx"), mmap=False, **kw)
    s.prepare([str(fa)])
    return s


def _planted(rng, probe_seq, name="OFFGENE"):
    """A reference whose only interesting content is ``probe_seq`` verbatim."""
    return (f">{name}\n"
            f"{_rand_seq(400, rng)}{probe_seq}{_rand_seq(400, rng)}\n")


# ---------------------------------------------------------------------------
# (1) A degenerate base must not remove a half from the search
# ---------------------------------------------------------------------------

def test_degenerate_half_is_still_searched(tmp_path):
    """An 'N' in the LHS must not delete the LHS from the hit list.

    The probe below matches a non-target record over all 50 bp except at the
    single position the 'N' sits on, where every resolution is available to it.
    Before the fix the LHS produced no hits at all, so the pair's verdict rested
    on the RHS alone.
    """
    rng = random.Random(11)
    true_seq = _rand_seq(50, rng)
    s = _screen(tmp_path, _planted(rng, true_seq))
    probe = true_seq[:5] + "N" + true_seq[6:]
    hits = s.search_chunk([_cand(probe)], 0, str(tmp_path))
    aln = split_junction_records(hits["p0"])[0]
    assert {hit_half(h) for h in aln} == {"L", "R"}
    best_left = min(h["mismatch"] for h in aln if hit_half(h) == "L")
    # One resolution of the 'N' IS the reference base, so the best duplex this
    # half can form there is perfect. Anything else understates the off-target.
    assert best_left == 0


def test_degenerate_half_cannot_silently_pass_the_pair(tmp_path):
    """The consequence: a pair with a perfect off-target must be rejected.

    ``off_target_junction_duplex_bp=0`` isolates the count rule — the junction
    gate would reject this pair too, and the point here is that the 10x count
    rule itself must see the LHS.
    """
    rng = random.Random(23)
    true_seq = _rand_seq(50, rng)
    s = _screen(tmp_path, _planted(rng, true_seq))
    probe = true_seq[:5] + "N" + true_seq[6:]
    cand = _cand(probe)
    hits = s.search_chunk([cand], 0, str(tmp_path))
    evaluate_specificity(cand, hits["p0"], set(),
                         off_target_junction_duplex_bp=0)
    assert cand.specificity_pass is False
    assert cand.specificity_status == "failed_off_target"


def test_too_many_degenerate_bases_fails_closed(tmp_path):
    """Above the expansion cap the screen refuses instead of guessing.

    Five 'N's in one half is 1024 concrete oligos; enumerating them is not the
    same search the rest of the pipeline believes it ran. The engine emits an
    explicit unscreenable record and ``specificity`` rejects the pair — the one
    thing it must never do is return "clean".
    """
    rng = random.Random(31)
    true_seq = _rand_seq(50, rng)
    s = _screen(tmp_path, _planted(rng, true_seq))
    probe = "NNNNN" + true_seq[5:]
    cand = _cand(probe)
    hits = s.search_chunk([cand], 0, str(tmp_path))
    evaluate_specificity(cand, hits["p0"], set(),
                         off_target_junction_duplex_bp=0)
    assert cand.specificity_pass is False
    assert cand.specificity_status == "unscreenable_half"


def test_plain_acgt_half_is_unaffected(tmp_path):
    """The hot path must be byte-identical for a probe with no ambiguity."""
    rng = random.Random(37)
    true_seq = _rand_seq(50, rng)
    s = _screen(tmp_path, _planted(rng, true_seq))
    cand = _cand(true_seq)
    hits = s.search_chunk([cand], 0, str(tmp_path))
    aln = split_junction_records(hits["p0"])[0]
    assert {hit_half(h) for h in aln} == {"L", "R"}
    assert all(h["qseq"] in (true_seq[:HALF], true_seq[HALF:]) for h in aln)


# ---------------------------------------------------------------------------
# (2) The junction gate is k-mer-only, so BLAST + gate must be refused
# ---------------------------------------------------------------------------

def _cfg(engine, jd_bp):
    return RunConfig(
        screen=ScreenResources(engine=engine),
        specificity=SpecificityThresholds(off_target_junction_duplex_bp=jd_bp),
    )


def test_blast_engine_with_the_junction_gate_is_refused():
    """A user who selects BLAST must not silently lose a hard filter."""
    with pytest.raises(ConfigError) as ei:
        _cfg(ENGINE_BLAST, 9).validate()
    msg = str(ei.value)
    assert "junction" in msg.lower()
    # Actionable: both ways out have to be in the message.
    assert ENGINE_KMER in msg
    assert "0" in msg


def test_blast_engine_with_the_gate_off_is_accepted():
    _cfg(ENGINE_BLAST, 0).validate()


def test_kmer_engine_keeps_the_gate():
    _cfg(ENGINE_KMER, 9).validate()


# ---------------------------------------------------------------------------
# (3) Per-gene scoping degrading to per-transcript must be announced
# ---------------------------------------------------------------------------

_ENSEMBL = (
    ">ENST00000000001.1 cdna chromosome:GRCh38:1:100:200:1 "
    "gene:ENSG00000000001.1 gene_biotype:protein_coding\n"
    "{a}\n"
    ">ENST00000000002.1 cdna chromosome:GRCh38:1:300:400:1 "
    "gene:ENSG00000000001.1 gene_biotype:protein_coding\n"
    "{b}\n"
)


def test_warns_when_no_reference_id_yields_a_gene(tmp_path, caplog):
    rng = random.Random(41)
    fasta = _ENSEMBL.format(a=_rand_seq(300, rng), b=_rand_seq(300, rng))
    with caplog.at_level(logging.WARNING, logger="flexprobe.screen"):
        _screen(tmp_path, fasta)
    text = "\n".join(r.getMessage() for r in caplog.records)
    assert "gene" in text.lower()
    assert "transcript" in text.lower()


def test_no_warning_when_ids_already_group_by_gene(tmp_path, caplog):
    rng = random.Random(43)
    fasta = (f">ENSG00000000001|ENST00000000001\n{_rand_seq(300, rng)}\n"
             f">ENSG00000000001|ENST00000000002\n{_rand_seq(300, rng)}\n")
    with caplog.at_level(logging.WARNING, logger="flexprobe.screen"):
        _screen(tmp_path, fasta)
    text = "\n".join(r.getMessage() for r in caplog.records)
    assert "per-transcript" not in text.lower()
