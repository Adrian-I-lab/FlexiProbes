"""Integration: the pre-screen per-gene cap reduces what reaches BLAST."""
from flexprobe.models import (
    ScreenResources, ExportSettings, ModalityDefinition, RunConfig,
)


# A ~200 bp, GC-balanced, low-homopolymer sequence so step-1 tiling yields many
# candidates that clear the filters (one gene).
_SEQ = ("ACGTACGTGCATGCATACGTACGTGCATGCATTGCATGCAACGTTGCAAC"
        "GTACGTGCATGCATACGTGTCAGCATGCATACGTACGTGCATTGCATGCA"
        "ACGTACGTGCATGCATACGTACGTGCATGCATTGCAACGTTGCATGCATA"
        "CGTACGTGCATGCATACGTGTCAGCATGCATACGTACGTGCATTGCATGC")


def _make_cfg(tmp_path, top_n):
    fa = tmp_path / "g.fasta"
    fa.write_text(f">GENE1 cdna\n{_SEQ}\n")
    return RunConfig(
        fasta_path=str(fa), output_dir=str(tmp_path / "out"),
        modality=ModalityDefinition(name="flex", version="v2"),
        # A real (nonexistent-prefix) DB would fail; we monkeypatch BLAST away,
        # so just need db_paths non-empty + skip_screen False to take the path.
        screen=ScreenResources(skip_screen=False, db_paths=["/fake/db"]),
        export=ExportSettings(overwrite=True), discovery_only=True,
        screen_top_per_gene=top_n,
    )


def _patch_blast(monkeypatch):
    import flexprobe.pipeline as pl
    seen = {}

    def fake_blast(passing, *_a, **_k):
        seen["n"] = len(passing)
        return {c.probe_id: [] for c in passing}

    monkeypatch.setattr(pl, "_run_screen", fake_blast)
    return seen


def test_cap_limits_candidates_reaching_blast(tmp_path, monkeypatch):
    from flexprobe.pipeline import run_pipeline
    seen = _patch_blast(monkeypatch)
    run_pipeline(_make_cfg(tmp_path, top_n=5))
    assert seen["n"] == 5, "one gene, cap 5 -> exactly 5 screened"


def test_no_cap_screens_all(tmp_path, monkeypatch):
    from flexprobe.pipeline import run_pipeline
    seen = _patch_blast(monkeypatch)
    run_pipeline(_make_cfg(tmp_path, top_n=0))   # 0 = screen everything
    assert seen["n"] > 5, "without a cap, far more than 5 candidates are screened"


def test_cap_ignored_when_blast_skipped(tmp_path, monkeypatch):
    # With skip_screen, the cap must NOT drop candidates (nothing to screen).
    from flexprobe.pipeline import run_pipeline
    cfg = _make_cfg(tmp_path, top_n=5)
    cfg.screen.skip_screen = True
    result = run_pipeline(cfg)
    # discovery results.csv keeps every candidate (cap is BLAST-only)
    assert result["stats"]["filter_passed"] > 5
