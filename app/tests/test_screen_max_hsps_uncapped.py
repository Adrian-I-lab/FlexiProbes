"""The BLAST screen must not cap HSPs per subject.

``-max_hsps 1`` keeps the best-SCORING alignment per subject, and best-scoring
is not most-disqualifying. A 5-mismatch alignment scoring 20 masks a 4-mismatch
alignment scoring 17 in the same record, so the screen never sees the hit that
should have rejected the probe.

With the vendored blastn 2.17.0 and real probe halves against a CDS
reference, a per-subject cap hides an alignment with FEWER mismatches than the
one kept in a few percent of (query, subject) pairs, occasionally an exact
25 nt match behind a 2-mismatch hit.

``run_blastn`` defaults ``max_hsps`` to ``None`` and omits the flag entirely,
which is correct. This guard pins the CALL SITE, because a value passed here
overrides ``run_blastn``'s own defaults (see the evalue/word_size note in
``BlastScreen.search_chunk``).
"""
from types import SimpleNamespace

from flexprobe import screen as S


class _Cand(SimpleNamespace):
    pass


def _capture_kwargs(monkeypatch, tmp_path):
    """Drive one ``search_chunk`` and return the kwargs it gave ``run_blastn``."""
    seen = {}

    def fake_run_blastn(*args, **kwargs):
        seen.update(kwargs)
        tsv = tmp_path / "hits.tsv"
        tsv.write_text("")
        return str(tsv)

    monkeypatch.setattr(S.blast, "run_blastn", fake_run_blastn)
    monkeypatch.setattr(S.blast, "parse_hits_tsv", lambda _p: [])

    sc = S.BlastScreen(min_pident=80.0, min_aln=20, num_threads=1)
    sc._tools = {"blastn": "blastn"}
    sc.db_arg = "db"
    cand = _Cand(probe_id="p1", sequence="ACGT" * 12 + "AA")
    sc.search_chunk([cand], 0, str(tmp_path))
    return seen


def test_screen_does_not_cap_hsps_per_subject(tmp_path, monkeypatch):
    kwargs = _capture_kwargs(monkeypatch, tmp_path)
    assert kwargs.get("max_hsps") is None, (
        "BlastScreen passed max_hsps=%r. Capping HSPs per subject hides "
        "lower-scoring but more-disqualifying alignments." % (kwargs.get("max_hsps"),)
    )


def test_screen_still_pushes_the_filters_it_should(tmp_path, monkeypatch):
    """Negative control: the same capture DOES see the flags that belong there.

    Without this, a capture that silently saw nothing would make the test above
    pass for the wrong reason.
    """
    kwargs = _capture_kwargs(monkeypatch, tmp_path)
    assert kwargs.get("perc_identity") == 80.0
    assert "max_target_seqs" in kwargs
    assert kwargs.get("num_threads") == 1
