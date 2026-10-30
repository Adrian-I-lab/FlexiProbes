"""The per-gene candidate count must be ONE pass over the pairs.

The original comprehension rescanned every pair once per gene — O(genes x
pairs). It produced the right answer, so no correctness test could have caught
it; what it cost was 12.0 h on rat (46,856 genes x 1.73M pairs) and 23.5 h on
mouse (87,587 x 1.81M), invisible in `elapsed_seconds` because the stage timer
snapshots before this runs. So the guard here counts attribute reads rather
than seconds: quadratic and linear differ by a factor of `len(gene_seqs)`
exactly, with no timing flake.
"""
from flexprobe.pipeline import _candidates_per_gene


class _Lhs:
    """Counts every read of `gene_id` — the access the quadratic form repeats."""

    def __init__(self, gene_id, counter):
        self._gene_id = gene_id
        self._counter = counter

    @property
    def gene_id(self):
        self._counter[0] += 1
        return self._gene_id


class _Pair:
    def __init__(self, lhs):
        self.lhs = lhs


def _pairs(gene_ids, counter):
    return [_Pair(_Lhs(g, counter)) for g in gene_ids]


GENE_SEQS = {"g1": "A", "g2": "C", "g3": "G", "empty": "T"}
PAIR_GENES = ["g1", "g2", "g1", "g3", "g1", "g2"]


def test_counts_match_a_brute_force_reference():
    counter = [0]
    pairs = _pairs(PAIR_GENES, counter)
    got = _candidates_per_gene(pairs, GENE_SEQS)
    assert got == {"g1": 3, "g2": 2, "g3": 1, "empty": 0}


def test_every_gene_is_a_key_even_with_no_pairs():
    # Consumers index by gene id; a missing key is a KeyError, not a zero.
    counter = [0]
    got = _candidates_per_gene(_pairs(PAIR_GENES, counter), GENE_SEQS)
    assert set(got) == set(GENE_SEQS)
    assert got["empty"] == 0


def test_key_order_follows_gene_seqs():
    # The dict is serialised into design_summary.json; a reordering would show
    # up as a spurious diff between two runs of the same input.
    counter = [0]
    got = _candidates_per_gene(_pairs(PAIR_GENES, counter), GENE_SEQS)
    assert list(got) == list(GENE_SEQS)


def test_reads_each_pair_exactly_once():
    # THE regression. Quadratic reads len(pairs) * len(gene_seqs) = 24 here.
    counter = [0]
    pairs = _pairs(PAIR_GENES, counter)
    _candidates_per_gene(pairs, GENE_SEQS)
    assert counter[0] == len(PAIR_GENES), (
        f"read gene_id {counter[0]} times for {len(PAIR_GENES)} pairs and "
        f"{len(GENE_SEQS)} genes — the screen is rescanning pairs per gene"
    )


def test_empty_inputs():
    assert _candidates_per_gene([], {}) == {}
    assert _candidates_per_gene([], GENE_SEQS) == {g: 0 for g in GENE_SEQS}
