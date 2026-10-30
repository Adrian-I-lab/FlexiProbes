"""The composite-score defaults are pinned, and nothing may drift them.

Without a test that names these values, they could change with only the golden
digests noticing, and a golden re-record would hide it. The provenance of the
numbers is in ``RankingWeights``'s docstring.
"""
from __future__ import annotations

import pytest

from flexprobe.ranking import DEFAULT_WEIGHTS, RankingWeights

G13 = {
    "gc_proximity": 0.35,
    "gc_balance": 0.34,
    "specificity_margin": 0.28,
    "homopolymer_penalty": 0.03,
    "central_placement": 0.00,
    "isoform_coverage": 0.0,
    "mismatch_placement": 0.10,
}


@pytest.mark.parametrize("name,value", sorted(G13.items()))
def test_default_weight(name, value):
    assert getattr(DEFAULT_WEIGHTS, name) == pytest.approx(value, abs=1e-12)
    assert getattr(RankingWeights(), name) == pytest.approx(value, abs=1e-12)


def test_five_main_weights_sum_to_one():
    """The fit is normalised over five components; keeping their sum at 1.00
    keeps mismatch_placement's 0.10 on the same scale as the other weights."""
    five = ("gc_proximity", "gc_balance", "specificity_margin",
            "homopolymer_penalty", "central_placement")
    assert sum(getattr(DEFAULT_WEIGHTS, k) for k in five) == pytest.approx(1.0)
