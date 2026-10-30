"""MAFFT versus the built-in progressive aligner: do they yield the same probes?

``consensus._pairwise_align`` aligns every record against the FIRST record
rather than performing a true multiple alignment, and is documented as a
fallback only. It nevertheless runs whenever MAFFT and MUSCLE are absent, which
can happen without notice on a Finder-launched .app.

On two-strain segmented-virus inputs and on synthetic multi-record cases, the
two aligners give the same conserved-window sets, with at most small
differences in discriminatory windows. This module builds sequences of that
character deterministically and pins the agreement here. The day an input
class breaks it, that is a test failure rather than a surprise inside a result.

Skipped when no MAFFT is installed. The comparison is meaningless without a
second aligner to compare against, and a skipped test says so honestly where a
vacuous pass would not.
"""
import random

import pytest

import flexprobe.consensus as C


def _panel(seed: int = 20260909):
    """Two diverged strains plus a third carrying an internal deletion.

    Mirrors the OROV panel's structure: ~2% substitutions between strains, one
    record truncated at both termini (the partial-CDS case, which produces
    terminal gap columns), and a third record with a 30 nt internal deletion so
    the alignment cannot be a simple column-for-column comparison.
    """
    rng = random.Random(seed)
    a = "".join(rng.choice("ACGT") for _ in range(3000))
    b = list(a)
    for i in rng.sample(range(len(b)), k=int(len(b) * 0.02)):
        b[i] = rng.choice([x for x in "ACGT" if x != b[i]])
    b = "".join(b)[60:-60]              # terminal truncation, as for a partial CDS
    c = a[:1500] + a[1530:]             # 30 nt internal deletion
    return [("strain_a", a), ("strain_b", b), ("strain_c", c)]


def _windows(aligned, ids):
    """Conserved and discriminatory window sets for one alignment."""
    cons, var = C._consensus_from_aligned(aligned)
    aln = C.StrainAlignment(
        fragment_name="parity", strain_ids=list(ids),
        aligned_seqs=list(aligned), consensus=cons, variant_positions=var,
    )
    conserved = {w.start for w in C.degenerate_conserved_windows(
        aln, tile_length=50, tile_step=1, max_degenerate_per_half=1)}
    unique = {s for s, _ in C.strain_unique_start_positions(
        aln, target_strain_idx=0, tile_length=50, tile_step=1)}
    gaps = sum(1 for i in range(len(cons))
               if any(s[i] == "-" for s in aligned))
    return conserved, unique, gaps, len(var)


@pytest.mark.skipif(C._find_msa_tool("mafft") is None,
                    reason="MAFFT not installed; nothing to compare against")
def test_builtin_aligner_matches_mafft_on_conserved_windows():
    records = _panel()
    ids = [sid for sid, _ in records]

    msa = C._external_msa(records)
    assert msa is not None, "MAFFT resolved but produced no alignment"
    mafft_aligned, tool = msa
    assert tool.startswith("MAFFT"), tool

    builtin_aligned = C._pairwise_align(records)

    m_cons, m_uniq, m_gaps, m_var = _windows(mafft_aligned, ids)
    b_cons, b_uniq, b_gaps, b_var = _windows(builtin_aligned, ids)

    # Positive control: the comparison must be non-trivial. Without gaps,
    # variant positions and a non-empty window set, two aligners would agree
    # vacuously and this test would pass while measuring nothing.
    assert m_gaps > 0, "no gap columns: the truncation/indel did not survive"
    assert m_var > 0, "no variant positions: the strains did not diverge"
    assert len(m_cons) > 0, "no conserved windows to compare"

    assert m_cons == b_cons, (
        f"conserved windows diverge: MAFFT-only {len(m_cons - b_cons)}, "
        f"built-in-only {len(b_cons - m_cons)}")


@pytest.mark.skipif(C._find_msa_tool("mafft") is None,
                    reason="MAFFT not installed; nothing to compare against")
def test_builtin_aligner_agrees_on_alignment_geometry():
    """Same alignment length and same gap columns, not merely the same probes.

    Equal window sets could in principle arise from two different alignments
    that happen to admit the same tiles. Checking the geometry as well makes
    the agreement a statement about the alignment itself.
    """
    records = _panel()
    msa = C._external_msa(records)
    assert msa is not None
    mafft_aligned, _ = msa
    builtin_aligned = C._pairwise_align(records)

    assert len({len(s) for s in mafft_aligned}) == 1
    assert len(mafft_aligned[0]) == len(builtin_aligned[0])
    for m, b in zip(mafft_aligned, builtin_aligned):
        assert m.upper() == b.upper()
