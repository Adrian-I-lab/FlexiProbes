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

"""Correctness of the k-mer search engine.

Two tiers, deliberately:

* the fast tests below pin the engine's CONTRACTS — record isolation, N
  handling, persistence, the seeding plans' arithmetic. They run every time.
* :func:`test_seeded_search_equals_brute_force` is the ORACLE, marked ``slow``
  and deselected by default. Pigeonhole seeding is a *method*; its answer is
  worth nothing unless it equals the *definition*, which is
  ``kmer_index.brute``. Same discipline, and the same reasoning, as
  ``test_oracle_exhaustive.py`` — run it whenever the search path changes.
"""
from __future__ import annotations

import numpy as np
import pytest

from flexprobe import kmer_index as K


def _rand(n: int, rng: np.random.Generator) -> np.ndarray:
    return rng.integers(0, 4, n, dtype=np.uint8)


def _mutate(q: np.ndarray, n_mm: int, rng: np.random.Generator) -> np.ndarray:
    q = q.copy()
    for p in rng.choice(len(q), n_mm, replace=False):
        q[p] = (q[p] + 1 + rng.integers(0, 3)) % 4
    return q


# ---------------------------------------------------------------------------
# The oracle
# ---------------------------------------------------------------------------

@pytest.mark.slow
@pytest.mark.parametrize("plan", list(K.iter_seed_plans()), ids=lambda p: p[0])
def test_seeded_search_equals_brute_force(plan):
    """Every seeding plan must return exactly what exhaustive comparison returns.

    Queries are PLANTED in the reference before being mutated: a random 25-mer
    almost never hits a random 400 kb reference, so testing with random queries
    would let a completely broken engine pass by agreeing that there is nothing
    to find.
    """
    _name, qlen, max_mm, offs, seed_mm = plan
    rng = np.random.default_rng(20260810)
    for _ in range(4):
        ref = _rand(400_000, rng)
        idx = K.KmerIndex(ref, 12)
        for _ in range(30):
            s = int(rng.integers(0, len(ref) - qlen))
            q = _mutate(ref[s:s + qlen], int(rng.integers(0, max_mm + 1)), rng)
            got_pos, got_mm = K.search(ref, idx, q, max_mm, offs, seed_mm)
            want_pos, want_mm = K.brute(ref, q, max_mm)
            order_g = np.argsort(got_pos)
            order_w = np.argsort(want_pos)
            assert np.array_equal(got_pos[order_g], want_pos[order_w])
            assert np.array_equal(got_mm[order_g], want_mm[order_w])


# ---------------------------------------------------------------------------
# Contracts (fast)
# ---------------------------------------------------------------------------

def test_pad_prevents_a_window_spanning_two_records(tmp_path):
    """A hit must never be assembled from two different records.

    Without the PAD separator, the tail of record A concatenated with the head
    of record B would present a 50 nt window that exists in no sequence at all —
    an off-target the screen would report and no molecule could produce.
    """
    fa = tmp_path / "two.fa"
    fa.write_text(">a\n" + "A" * 60 + "\n>b\n" + "C" * 60 + "\n")
    ref = K.load_reference(str(fa))
    idx = K.KmerIndex(ref.seq, 12)
    # 25 A's followed by 25 C's exists only across the record boundary.
    q = np.concatenate([np.zeros(25, dtype=np.uint8),
                        np.ones(25, dtype=np.uint8)])
    pos, _ = K.search(ref.seq, idx, q, 4, (0, 12, 25, 37), 2)
    assert len(pos) == 0


def test_query_longer_than_pad_is_refused():
    with pytest.raises(ValueError, match="record separator"):
        K.assert_pad_covers(K.PAD + 1)


# ---------------------------------------------------------------------------
# Record containment — the PAD is a cost, not a boundary
# ---------------------------------------------------------------------------

def _two_record_reference(tmp_path):
    """recA (300 A) then recB, whose opening 21 nt are distinctive.

    Sized so a 25 nt query of ``TTTT`` + ``recB[:21]`` matches at exactly 4
    mismatches when it starts 4 bases inside recA's trailing pad — the budget
    the 25 nt half search actually runs at.
    """
    rec_a = "A" * 300
    rec_b = ("GTTACCGATCAGTTACCGATCAGTTACCGATCAGGGCCTTAACCGGTTAA"
             + "C" * 250)
    fa = tmp_path / "two_records.fa"
    fa.write_text(f">recA_geneA\n{rec_a}\n>recB_geneB\n{rec_b}\n")
    ref = K.load_reference(str(fa), exclude_non_primary=False)
    return ref, rec_a, rec_b


def test_a_window_starting_in_a_pad_is_not_a_hit(tmp_path):
    """REGRESSION. The sentinels cost mismatches; they do not disqualify.

    A window can therefore START inside record A's trailing pad, spend its
    whole budget on the four sentinels it covers, and still pass — while 21 of
    its 25 bases are record B's. ``record_of`` maps its start back to record A,
    because the pad sits on A's side of ``starts[1]``, so the site is reported
    against the WRONG GENE. Measured at 0.0197% of verdicts over 10,175 probes,
    and one-directional: it can excuse a real off-target (when the preceding
    record is the probe's own gene) or invent one.

    Such a window is a chimaera of sentinel and neighbour, so it is DISCARDED,
    not re-attributed — re-attributing it to record B would still be claiming a
    site four of whose bases are not sequence. The genuine in-record hit at
    3 mismatches must survive untouched.
    """
    ref, _rec_a, rec_b = _two_record_reference(tmp_path)
    idx = K.KmerIndex(ref.seq, 12)
    start_b = int(ref.starts[1])
    q = K.encode_seq("TTTT" + rec_b[:21])

    pos, mm = K.search(ref.seq, idx, q, 4, (0, 13), 2, bounds=ref.bounds)
    reported = sorted(zip(pos.tolist(), mm.tolist()))
    assert reported == [(start_b + 7, 3)], reported
    # And it is not merely absent, it was absent for the right reason: the
    # unbounded search still proposes it, so the filter is what removes it.
    unbounded = sorted(K.search(ref.seq, idx, q, 4, (0, 13), 2)[0].tolist())
    assert start_b - 4 in unbounded
    assert int(ref.record_of(np.array([start_b - 4]))[0]) == 0


def test_a_hit_at_the_very_first_base_of_a_record_is_kept(tmp_path):
    """Off-by-one guard: the containment test must not eat real hits.

    ``starts[i]`` is a legal window start. If ``contains_window`` used ``<``
    where it uses ``<=``, or tested ``pos > starts[rec]``, every record's
    opening window would vanish and the screen would quietly stop seeing
    off-targets at transcript 5' ends.
    """
    ref, rec_a, rec_b = _two_record_reference(tmp_path)
    idx = K.KmerIndex(ref.seq, 12)
    for rec_i, seq in enumerate((rec_a, rec_b)):
        start = int(ref.starts[rec_i])
        q = K.encode_seq(seq[:25])
        pos, mm = K.search(ref.seq, idx, q, 4, (0, 13), 2, bounds=ref.bounds)
        assert start in pos.tolist(), f"record {rec_i} lost its opening window"
        assert mm[pos.tolist().index(start)] == 0


def test_the_last_record_ends_a_pad_before_the_array_does(tmp_path):
    """The final record's arithmetic differs and is easy to get wrong.

    ``ends[i] = starts[i+1] - pad`` has no ``starts[i+1]`` for the last record;
    its end is ``len(seq) - pad``, because ``load_reference`` appends a pad
    block after EVERY record, the final one included. Using ``len(seq)`` there
    would leave the last record's own trailing pad searchable — exactly the
    defect, surviving in the one place it is hardest to notice.
    """
    ref, _rec_a, rec_b = _two_record_reference(tmp_path)
    assert int(ref.bounds.ends[-1]) == int(ref.starts[-1]) + len(rec_b)
    assert int(ref.bounds.ends[-1]) == len(ref.seq) - K.PAD

    idx = K.KmerIndex(ref.seq, 12)
    last_legal = len(ref.seq) - K.PAD - 25       # last window wholly in recB
    q = ref.seq[last_legal:last_legal + 25].copy()
    pos, _ = K.search(ref.seq, idx, q, 4, (0, 13), 2, bounds=ref.bounds)
    assert last_legal in pos.tolist()
    assert max(pos.tolist()) <= last_legal, "a window overhung the last record"


def test_record_bounds_agree_with_a_per_base_validity_mask(tmp_path):
    """Independent path: the O(1)-per-record arithmetic vs an O(1)-per-base map.

    ``RecordBounds`` is index arithmetic, and index arithmetic is exactly what
    was wrong to begin with, so it is checked against a mask built a different
    way — mark every real base True, then ask whether a window is all-True.
    The mask is the definition; the arithmetic is the method.

    The assertion below is only worth anything if the mask can say both things,
    so that is asserted first. A check that returns one verdict everywhere is
    measuring nothing.
    """
    ref, rec_a, rec_b = _two_record_reference(tmp_path)
    real = np.zeros(len(ref.seq), dtype=bool)
    for start, length in zip(ref.starts.tolist(), (len(rec_a), len(rec_b))):
        real[start:start + length] = True
    assert real.any() and not real.all()         # positive control on the mask

    L = 25
    pos = np.arange(len(ref.seq) - L + 1, dtype=np.int64)
    want = np.array([bool(real[p:p + L].all()) for p in pos.tolist()])
    assert want.any() and not want.all()         # ... and on the answer
    got = ref.bounds.contains_window(pos, L)
    assert np.array_equal(got, want)


@pytest.mark.slow
def test_search_equals_brute_force_over_a_padded_reference(tmp_path):
    """The oracle, re-run where the padding actually exists.

    ``test_seeded_search_equals_brute_force`` searches a BARE array: no
    records, no separators, nothing to be misattributed. It is the right test
    for the pigeonhole guarantee and it cannot see this defect at all — both
    sides of it agreed on the fictional sites, because brute force over a
    padded array counts them too.

    So the definition is restated for a real reference — every offset whose
    window lies inside one record — and the method must equal THAT. Neither
    side was changed to make them agree: ``brute`` gained the same ``bounds``
    argument, and passing it to only one side fails.
    """
    rng = np.random.default_rng(20260910)
    seqs = ["".join("ACGT"[c] for c in rng.integers(0, 4, int(n)))
            for n in rng.integers(200, 900, 12)]
    fa = tmp_path / "many.fa"
    fa.write_text("".join(f">rec{i}|gene{i}\n{s}\n"
                          for i, s in enumerate(seqs)))
    ref = K.load_reference(str(fa), exclude_non_primary=False)
    idx = K.KmerIndex(ref.seq, 12)

    queries = []
    for _ in range(60):                          # planted in a record, mutated
        r = int(rng.integers(0, len(seqs)))
        s = int(rng.integers(0, max(1, len(seqs[r]) - 25)))
        p = int(ref.starts[r]) + s
        queries.append(_mutate(ref.seq[p:p + 25].copy(),
                               int(rng.integers(0, 5)), rng))
    # Planted ACROSS a boundary: four filler bases then a record's opening 21.
    # Its best window starts four bases inside the PREVIOUS record's pad and
    # costs exactly the budget, which is the case the whole fix is about. A
    # purely random query almost never lands one, so the positive control at
    # the end of this test would be vacuous without these.
    for r in range(1, len(seqs)):
        queries.append(K.encode_seq("TTTT" + seqs[r][:21]))

    disagreements_without_bounds = 0
    for q in queries:
        got_pos, got_mm = K.search(ref.seq, idx, q, 4, (0, 13), 2,
                                   bounds=ref.bounds)
        want_pos, want_mm = K.brute(ref.seq, q, 4, bounds=ref.bounds)
        g, w = np.argsort(got_pos), np.argsort(want_pos)
        assert np.array_equal(got_pos[g], want_pos[w])
        assert np.array_equal(got_mm[g], want_mm[w])
        loose, _ = K.brute(ref.seq, q, 4)
        disagreements_without_bounds += len(loose) - len(want_pos)
    # Positive control: the unbounded definition really does admit sites the
    # bounded one refuses, so the agreement above was not vacuous.
    assert disagreements_without_bounds > 0


def test_reference_n_counts_as_a_mismatch(tmp_path):
    """An N must make the screen stricter, never invisible.

    INVALID (4) can equal no query code, so a window covering an N accrues a
    mismatch there. The alternative — treating it as a wildcard match — would
    let an off-target hide behind an unsequenced base.
    """
    fa = tmp_path / "n.fa"
    fa.write_text(">a\n" + "ACGT" * 10 + "N" + "ACGT" * 10 + "\n")
    ref = K.load_reference(str(fa))
    assert (ref.seq == K.INVALID).sum() >= 1
    idx = K.KmerIndex(ref.seq, 12)
    exact = ref.seq[0:25].copy()
    assert (exact != K.INVALID).all()
    pos, mm = K.search(ref.seq, idx, exact, 0, (0, 13), 2)
    assert len(pos) >= 1 and mm.min() == 0


def test_seed_containing_n_proposes_nothing_but_others_still_do(tmp_path):
    """An N in one seed must not abort the search — the other seed carries it."""
    rng = np.random.default_rng(7)
    ref = _rand(50_000, rng)
    idx = K.KmerIndex(ref, 12)
    s = 1234
    q = ref[s:s + 25].copy()
    q[2] = K.INVALID                      # inside the first seed only
    pos, _ = K.search(ref, idx, q, 4, (0, 13), 2)
    assert s in set(pos.tolist())


def test_search_rejects_a_seed_that_does_not_fit():
    rng = np.random.default_rng(1)
    ref = _rand(20_000, rng)
    idx = K.KmerIndex(ref, 12)
    with pytest.raises(ValueError, match="does not fit"):
        K.search(ref, idx, ref[0:25].copy(), 4, (0, 20), 2)


@pytest.mark.parametrize("max_mm,expected", [(0, 1), (1, 37), (2, 631)])
def test_neighbourhood_sizes_are_exact(max_mm, expected):
    """1 + 12*3 + C(12,2)*9 for a 12-mer. Wrong sizes mean lost or wasted hits."""
    seed = np.array([0, 1, 2, 3] * 3, dtype=np.uint8)
    codes = K.neighbourhood(seed, max_mm)
    assert len(codes) == expected
    assert len(np.unique(codes)) == expected


def test_neighbourhood_refuses_unsupported_budget():
    seed = np.array([0, 1, 2, 3] * 3, dtype=np.uint8)
    with pytest.raises(NotImplementedError):
        K.neighbourhood(seed, 3)


def test_seed_plans_satisfy_the_pigeonhole_bound():
    """The guarantee is arithmetic; assert it rather than trusting the comment.

    If some seed budget were ever set below ``max_mm // n_seeds`` the search
    would still run and would still look fine — it would just silently stop
    being exhaustive, which is the one failure mode this engine exists to avoid.
    """
    for name, qlen, max_mm, offs, seed_mm in K.iter_seed_plans():
        assert max_mm // len(offs) <= seed_mm, name
        for off in offs:                       # seeds must fit the query
            assert 0 <= off and off + 12 <= qlen, name
        ranges = sorted((o, o + 12) for o in offs)
        for (_, end), (start, _) in zip(ranges, ranges[1:]):
            assert start >= end, f"{name}: seeds overlap"


def test_plan_for_known_and_unknown():
    assert K.plan_for(25, 4) == ((0, 13), 2)
    with pytest.raises(KeyError, match="no seeding plan"):
        K.plan_for(30, 4)


def test_revcomp_codes_roundtrips_and_keeps_invalid():
    codes = np.array([0, 1, 2, 3, K.INVALID], dtype=np.uint8)
    rc = K.revcomp_codes(codes)
    assert rc.tolist() == [K.INVALID, 0, 1, 2, 3]
    assert K.revcomp_codes(rc).tolist() == codes.tolist()


def test_lookup_many_matches_repeated_lookup():
    rng = np.random.default_rng(3)
    ref = _rand(30_000, rng)
    idx = K.KmerIndex(ref, 8)
    codes = np.array([code_i for code_i in rng.integers(0, 4 ** 8, 40)],
                     dtype=np.int64)
    want = np.concatenate([idx.lookup(int(c)) for c in codes]
                          + [np.empty(0, dtype=np.int32)])
    got = idx.lookup_many(codes)
    assert sorted(got.tolist()) == sorted(want.tolist())


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------

def test_index_roundtrips_through_disk(tmp_path):
    rng = np.random.default_rng(11)
    ref = _rand(40_000, rng)
    idx = K.KmerIndex(ref, 12)
    prefix = str(tmp_path / "idx")
    idx.save(prefix, source_sha="deadbeef")
    back = K.KmerIndex.load(prefix)
    assert back.k == idx.k
    assert np.array_equal(np.asarray(back.pos), idx.pos)
    assert np.array_equal(back.offsets, idx.offsets)
    q = _mutate(ref[500:525], 3, rng)
    a = K.search(ref, idx, q, 4, (0, 13), 2)
    b = K.search(ref, back, q, 4, (0, 13), 2)
    assert np.array_equal(a[0], b[0]) and np.array_equal(a[1], b[1])


def test_build_or_load_caches_then_rebuilds_when_the_fasta_changes(tmp_path):
    """The sidecar must invalidate on content, like blast.build_db's .sha256.

    An index silently reused after its reference changed would screen against a
    sequence that is no longer there — wrong in the unsafe direction, and
    invisible.
    """
    fa = tmp_path / "r.fa"
    fa.write_text(">a\n" + "ACGTACGTACGTAAGGCCTTACGTACGT" * 4 + "\n")
    idx1, _ = K.KmerIndex.build_or_load(str(fa), 8, cache_dir=str(tmp_path))
    meta = tmp_path / "r.k8.meta.json"
    assert meta.exists()
    stamp = meta.stat().st_mtime_ns
    idx2, _ = K.KmerIndex.build_or_load(str(fa), 8, cache_dir=str(tmp_path))
    assert meta.stat().st_mtime_ns == stamp          # cache hit, not rebuilt
    assert np.array_equal(np.asarray(idx2.pos), np.asarray(idx1.pos))

    fa.write_text(">a\n" + "TTTTGGGGCCCCAAAATTTTGGGGCCCC" * 4 + "\n")
    idx3, _ = K.KmerIndex.build_or_load(str(fa), 8, cache_dir=str(tmp_path))
    assert meta.stat().st_mtime_ns != stamp          # rebuilt
    assert not np.array_equal(np.asarray(idx3.pos), np.asarray(idx1.pos))


def test_load_refuses_a_foreign_format(tmp_path):
    rng = np.random.default_rng(5)
    idx = K.KmerIndex(_rand(5_000, rng), 8)
    prefix = str(tmp_path / "idx")
    idx.save(prefix)
    import json
    p = tmp_path / "idx.meta.json"
    meta = json.loads(p.read_text())
    meta["format"] = K.INDEX_FORMAT + 1
    p.write_text(json.dumps(meta))
    with pytest.raises(ValueError, match="rebuild"):
        K.KmerIndex.load(prefix)


# ---------------------------------------------------------------------------
# Reference parsing
# ---------------------------------------------------------------------------

def test_reference_maps_positions_back_to_records(tmp_path):
    fa = tmp_path / "m.fa"
    fa.write_text(">first|x\nACGTACGTAC\n>second|y\nTTTTGGGGCC\n")
    ref = K.load_reference(str(fa))
    assert ref.names == ["first|x", "second|y"]
    assert ref.name_of(0) == "first|x"
    assert ref.name_of(int(ref.starts[1])) == "second|y"
    assert ref.record_of(np.array([0, int(ref.starts[1])])).tolist() == [0, 1]


def test_empty_fasta_is_an_error(tmp_path):
    fa = tmp_path / "e.fa"
    fa.write_text("")
    with pytest.raises(ValueError, match="no FASTA records"):
        K.load_reference(str(fa))


# ---------------------------------------------------------------------------
# Reference-size limits (Phase 5)
# ---------------------------------------------------------------------------

def test_oversized_reference_is_refused_not_silently_wrapped():
    """Above the int32 range, positions would wrap to NEGATIVE and the search
    would return silent garbage. That is a correctness limit, not a memory
    preference, so it must raise rather than warn."""
    with pytest.raises(ValueError, match="32-bit position index"):
        K.check_indexable(K.MAX_INDEXABLE_BASES + 1)


def test_refusal_names_the_genome_confusion_and_the_escape_hatch():
    """The limit falls exactly where the intended scope does, so the message has
    to say so — the user's mistake is almost always 'genome, not transcriptome'."""
    with pytest.raises(ValueError) as e:
        K.check_indexable(3_100_000_000, name="GRCh38.fa")
    msg = str(e.value)
    assert "GRCh38.fa" in msg
    assert "GENOME" in msg
    assert "--screen-engine blast" in msg


def test_mammalian_genomes_exceed_the_limit_and_transcriptomes_do_not():
    """Pins the coincidence the error message relies on."""
    for genome in (3_100_000_000, 2_700_000_000):        # human, mouse
        with pytest.raises(ValueError):
            K.check_indexable(genome)
    for tx in (372_000_000, 150_000_000):                # human cDNA, smaller
        K.check_indexable(tx)


def test_large_but_indexable_reference_warns_with_a_footprint(caplog):
    import logging
    with caplog.at_level(logging.WARNING, logger="flexprobe.kmer_index"):
        K.check_indexable(K.LARGE_REFERENCE_WARN_BASES + 1, name="big.fa")
    assert caplog.records, "a multi-Gbase index should announce its footprint"


def test_index_construction_enforces_the_limit(monkeypatch):
    """The guard must sit on the constructor, not only on the helper — that is
    the path every caller actually takes."""
    monkeypatch.setattr(K, "MAX_INDEXABLE_BASES", 100)
    rng = np.random.default_rng(2)
    with pytest.raises(ValueError, match="32-bit position index"):
        K.KmerIndex(_rand(500, rng), 8)


# ---------------------------------------------------------------------------
# Where the index is cached
# ---------------------------------------------------------------------------

def test_index_is_cached_beside_the_reference_by_default(tmp_path):
    """An index is a property of the REFERENCE, not of the run that needed it.

    Scoping it to an output directory silently defeats the SHA-256 cache key:
    every run rebuilds from scratch and writes its own copy. Measured on GRCh38
    cdna.all that is ~5 min per run against 1.4 s for a cache hit.
    """
    refdir = tmp_path / "refs"
    refdir.mkdir()
    fa = refdir / "r.fa"
    fa.write_text(">a\n" + "ACGTACGTACGTAAGGCCTTACGTACGT" * 4 + "\n")
    K.KmerIndex.build_or_load(str(fa), 8)
    assert (refdir / "r.k8.meta.json").exists(), "index not beside the FASTA"


def test_a_second_run_reuses_the_index_rather_than_rebuilding(tmp_path):
    """The property the location exists to deliver."""
    refdir = tmp_path / "refs"
    refdir.mkdir()
    fa = refdir / "r.fa"
    fa.write_text(">a\n" + "ACGTACGTACGTAAGGCCTTACGTACGT" * 4 + "\n")
    K.KmerIndex.build_or_load(str(fa), 8)
    stamp = (refdir / "r.k8.meta.json").stat().st_mtime_ns
    # A different "run" — different output dir, same reference.
    K.KmerIndex.build_or_load(str(fa), 8,
                              fallback_dir=str(tmp_path / "run2"))
    assert (refdir / "r.k8.meta.json").stat().st_mtime_ns == stamp


def test_read_only_reference_dir_falls_back_without_failing(tmp_path):
    """A read-only reference mount must still run — just without sharing."""
    import os
    import stat
    refdir = tmp_path / "refs"
    refdir.mkdir()
    fa = refdir / "r.fa"
    fa.write_text(">a\n" + "ACGTACGTACGTAAGGCCTTACGTACGT" * 4 + "\n")
    fallback = tmp_path / "out" / "_kmer_index"
    os.chmod(refdir, stat.S_IRUSR | stat.S_IXUSR)
    try:
        K.KmerIndex.build_or_load(str(fa), 8, fallback_dir=str(fallback))
        assert (fallback / "r.k8.meta.json").exists()
    finally:
        os.chmod(refdir, stat.S_IRWXU)


def test_explicit_cache_dir_still_wins(tmp_path):
    refdir = tmp_path / "refs"; refdir.mkdir()
    fa = refdir / "r.fa"
    fa.write_text(">a\n" + "ACGTACGTACGTAAGGCCTTACGTACGT" * 4 + "\n")
    explicit = tmp_path / "explicit"
    K.KmerIndex.build_or_load(str(fa), 8, cache_dir=str(explicit))
    assert (explicit / "r.k8.meta.json").exists()
    assert not (refdir / "r.k8.meta.json").exists()
