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

"""The parallel screen must be an optimisation and nothing else.

Two properties carry the whole feature and both are pinned here:

* **Same answer.** Splitting a chunk across worker processes must produce the
  verdict every probe would have got serially. ``batch_evaluate`` is a pure
  per-candidate loop, so this holds by construction — but "by construction" is
  what every silent-wrong-answer bug was before it shipped.
* **Same cache.** ``screen_workers`` must not reach either fingerprint. If it
  did, raising it on a chain already 111 chunks into a 726-chunk run would
  discard those chunks and re-tile a whole transcriptome — turning a speedup
  into a 27-hour regression.
"""
from __future__ import annotations

import random

import pytest

from flexprobe import candidate_cache, screen_parallel
from flexprobe.screen_cache import VERDICT_FIELDS
from flexprobe.models import ProbeCandidate, RunConfig
from flexprobe.screen import ENGINE_BLAST, ENGINE_KMER, make_screen
from flexprobe.specificity import batch_evaluate


# ---------------------------------------------------------------------------
# Slicing: a partition, not merely a split
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("n,w", [(10, 4), (3, 8), (5000, 12), (1, 4), (7, 7)])
def test_slices_partition_exactly(n, w):
    """Every index exactly once, in order, no gaps — a dropped slice would lose
    probes silently and a duplicated one would double-count the tally."""
    slices = screen_parallel._slice_bounds(n, w)
    seen = [i for r in slices for i in r]
    assert seen == list(range(n))
    assert all(len(r) > 0 for r in slices)


def test_slices_are_balanced():
    """Sizes differ by at most one, so no worker holds a disproportionate share
    of a chunk's hits (the term that sets peak memory)."""
    sizes = [len(r) for r in screen_parallel._slice_bounds(5000, 12)]
    assert max(sizes) - min(sizes) <= 1
    assert sum(sizes) == 5000


def test_zero_candidates_yields_no_slices():
    assert screen_parallel._slice_bounds(0, 8) == []


# ---------------------------------------------------------------------------
# Cache keying: raising screen_workers must RESUME, never invalidate
# ---------------------------------------------------------------------------

def test_screen_workers_excluded_from_candidate_fingerprint():
    """The candidate cache fingerprints by DENYLIST — every RunConfig field is
    folded in unless named — so a new field is invalidating by default. This is
    the test that would have caught the omission."""
    payload = candidate_cache.fingerprint_payload(RunConfig(fasta_path="x"))
    assert "screen_workers" not in payload
    assert "screen_workers" in candidate_cache.FINGERPRINT_IGNORED


def test_raising_screen_workers_keeps_both_digests():
    a = RunConfig(fasta_path="x", screen_workers=1)
    b = RunConfig(fasta_path="x", screen_workers=12)
    assert candidate_cache.fingerprint(a) == candidate_cache.fingerprint(b)


# ---------------------------------------------------------------------------
# Equivalence: the parallel path answers what the serial path answers
# ---------------------------------------------------------------------------

def _reference(tmp_path, n_genes: int = 6, glen: int = 900):
    """A small multi-gene reference that produces real hits (it is its own DB)."""
    rng = random.Random(20260820)
    seqs = {f"G{i}": "".join(rng.choice("ACGT") for _ in range(glen))
            for i in range(n_genes)}
    # Give two genes a shared 300 bp block so cross-gene off-targets exist and
    # the screen has something to actually reject.
    seqs["G1"] = seqs["G0"][100:400] + seqs["G1"][300:]
    fa = tmp_path / "ref.fa"
    fa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
    return str(fa), seqs


def _candidates(seqs, step: int = 7):
    out = []
    for gene, s in seqs.items():
        for start in range(0, len(s) - 50, step):
            out.append(ProbeCandidate(
                gene_id=gene, probe_id=f"{gene}_{start:05d}",
                sequence=s[start:start + 50], start=start, end=start + 50,
            ))
    return out


def _gene_aware(seqs):
    """Per-gene on-target scoping — the production configuration."""
    return {"gene_of": {g: g for g in seqs},
            "group_members": {g: {g} for g in seqs}}


def _screen_for(fa):
    sc = make_screen(ENGINE_KMER, query_half_len=25)
    sc.prepare([fa])
    return sc


def _verdicts(cands):
    return {c.probe_id: {f: getattr(c, f) for f in VERDICT_FIELDS}
            for c in cands}


def _same(a, b) -> bool:
    """Verdict-set equality, treating NaN as equal to NaN.

    Plain ``==`` is wrong here and fails on EVERY probe: a k-mer hit carries
    ``evalue``/``bitscore`` as NaN by design (they are undefinable for an
    end-to-end search, and NaN rather than 0 because 0 is a real bitscore), and
    ``offtarget_min_qcov`` is NaN wherever no off-target was scored. NaN != NaN,
    so a naive comparison reports a difference where the bits are identical.
    """
    import math

    def eq(x, y):
        if isinstance(x, float) and isinstance(y, float):
            return (math.isnan(x) and math.isnan(y)) or x == y
        return x == y

    if set(a) != set(b):
        return False
    return all(eq(a[p][f], b[p][f]) for p in a for f in VERDICT_FIELDS)


@pytest.mark.skipif(not screen_parallel.available(),
                    reason="requires the fork start method")
def test_parallel_verdicts_match_serial(tmp_path):
    """The load-bearing test: same verdicts, probe for probe, field for field."""
    import copy

    fa, seqs = _reference(tmp_path)
    cands = _candidates(seqs)
    assert len(cands) > 100, "fixture must be big enough to span slices"
    target_ids = set(seqs)
    # Gene-aware, as every whole-transcriptome run is (group_by_gene: true):
    # each probe's on-target set is its OWN gene, so the block G0 and G1 share
    # makes real cross-gene off-targets for the screen to reject.
    spec_kwargs = _gene_aware(seqs)

    # Serial: one search over the whole chunk, then score it.
    sc = _screen_for(fa)
    serial_cands = copy.deepcopy(cands)
    hits = sc.search_chunk(serial_cands, 0, str(tmp_path))
    batch_evaluate(serial_cands, hits, target_ids, subjects=None, **spec_kwargs)
    expected = _verdicts(serial_cands)
    sc.close()

    # Parallel: four slices, each searched AND scored in its own process.
    sc2 = _screen_for(fa)
    par = screen_parallel.ParallelScreen(
        sc2, workers=4, output_dir=str(tmp_path),
        target_ids=target_ids, spec_kwargs=spec_kwargs,
    )
    with par:
        got = par.verdicts_for(0, copy.deepcopy(cands))
    sc2.close()

    assert set(got) == set(expected)
    assert _same(got, expected)


@pytest.mark.skipif(not screen_parallel.available(),
                    reason="requires the fork start method")
def test_parallel_actually_rejects_something(tmp_path):
    """Guards the test above from passing vacuously: if the fixture produced no
    rejections, 'identical verdicts' would prove nothing."""
    import copy

    fa, seqs = _reference(tmp_path)
    cands = _candidates(seqs)
    sc = _screen_for(fa)
    par = screen_parallel.ParallelScreen(
        sc, workers=3, output_dir=str(tmp_path),
        target_ids=set(seqs), spec_kwargs=_gene_aware(seqs),
    )
    with par:
        got = par.verdicts_for(0, copy.deepcopy(cands))
    sc.close()
    statuses = {v["specificity_status"] for v in got.values()}
    assert len(statuses) > 1, f"fixture produced only {statuses}"


# ---------------------------------------------------------------------------
# The factory refuses loudly rather than running serially in silence
# ---------------------------------------------------------------------------

def _factory(cfg, **kw):
    from flexprobe.pipeline import _make_parallel_factory
    kw.setdefault("needs_subjects", False)
    kw.setdefault("output_dir", "/tmp")
    return _make_parallel_factory(cfg, set(), {}, **kw)


def test_serial_by_default():
    assert _factory(RunConfig(fasta_path="x")) is None


def test_refused_for_blast_engine(caplog):
    cfg = RunConfig(fasta_path="x", screen_workers=8)
    cfg.screen.engine = ENGINE_BLAST
    assert _factory(cfg) is None
    assert "screen_workers" in caplog.text


def test_refused_when_subjects_are_needed(caplog):
    cfg = RunConfig(fasta_path="x", screen_workers=8)
    cfg.screen.engine = ENGINE_KMER
    assert _factory(cfg, needs_subjects=True) is None
    assert "screen_workers" in caplog.text


def test_enabled_for_kmer():
    cfg = RunConfig(fasta_path="x", screen_workers=4)
    cfg.screen.engine = ENGINE_KMER
    assert _factory(cfg) is not None


def test_screen_workers_must_be_positive():
    from flexprobe.errors import ConfigError

    with pytest.raises(ConfigError, match="screen_workers"):
        RunConfig(fasta_path="x", screen_workers=0).validate()
