"""The manifest must be able to reproduce the run.

The failure this guards against: if ``pipeline`` built ``manifest.params`` from
a hand-curated allowlist, anything not explicitly listed would be silently
dropped. The settings checked here are the ones such a manifest loses: every
``filters.*`` threshold (the whole filter funnel), all four ``specificity.*``
thresholds (the off-target gate itself), ``ranking_weights``,
``group_by_gene``, ``min_isoform_coverage``,
``extra_on_target_ids``, ``workers``, ``screen_cache_mode``, and
``blast.db_paths`` (i.e. WHICH DATABASE was screened against). ``inputs`` must
carry a hash, not a bare path, and ``db_sha256s`` must be filled, not ``{}``.

Net effect: a reader asking "what GC range and off-target identity produced
this panel?" can answer it from the outputs.
"""
from flexprobe import blast as fp_blast
from flexprobe.screen import ENGINE_BLAST
import dataclasses
import hashlib
import json
import os

import pytest

from flexprobe.config_io import load_config
from flexprobe.models import (
    ScreenResources, ExportSettings, FilterThresholds, ModalityDefinition,
    RunConfig, SpecificityThresholds,
)
from flexprobe.pipeline import run_pipeline

from conftest import SYNTH_200, write_fasta


def _leaf_paths(obj, prefix=""):
    for f in dataclasses.fields(obj):
        val = getattr(obj, f.name)
        path = f"{prefix}{f.name}"
        if dataclasses.is_dataclass(val):
            yield from _leaf_paths(val, prefix=f"{path}.")
        else:
            yield path, val


ALL_LEAVES = [p for p, _v in _leaf_paths(RunConfig())]


def _dig(params, path):
    """Follow a dotted path into the manifest's params dict."""
    cur = params
    for part in path.split("."):
        if not isinstance(cur, dict) or part not in cur:
            raise KeyError(path)
        cur = cur[part]
    return cur


def _cfg(tmp_path, fa, **over):
    """A config with every field set to a NON-default value where possible, so
    a manifest that merely echoed the defaults would not pass."""
    kw = dict(
        fasta_path=fa,
        output_dir=str(tmp_path / "out"),
        modality=ModalityDefinition(name="flex", version="v2", spatial=False),
        filters=FilterThresholds(
            gc_min=0.31, gc_max=0.81, max_n_per_half=2,
            max_homopolymer_run=6, reject_junction_base_n=False,
            require_tn_junction=False,
        ),
        specificity=SpecificityThresholds(
            off_target_pident_min=77.0, off_target_aln_min=19,
            on_target_pident_min=91.0, on_target_aln_min=22,
        ),
        screen=ScreenResources(skip_screen=True),
        export=ExportSettings(overwrite=True),
        tile_length=50, tile_step=10,
        multiplex=False, probe_barcodes=["BC003"],
        seed=1234,
        workers=1, screen_cache_mode="hits",
        extra_on_target_ids=["NM_999999.1"],
        min_isoform_coverage=0.25,
        group_by_gene=False,        # non-default: True is the default
        screen_top_per_gene=0,
        ranking_weights={"gc_proximity": 0.42},
        split_header_on_pipe=False,
        resume_screen=False,
        primer3_scoring=False,
    )
    kw.update(over)
    return RunConfig(**kw)


@pytest.fixture
def run_manifest(tmp_path):
    fa = write_fasta(tmp_path / "g.fa", [("G1", SYNTH_200), ("G2", SYNTH_200[::-1])])
    cfg = _cfg(tmp_path, fa)
    run_pipeline(cfg)
    with open(os.path.join(cfg.output_dir, "manifest.json")) as fh:
        return cfg, json.load(fh)


# ---------------------------------------------------------------------------
# (a) every RunConfig field is recorded
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", ALL_LEAVES)
def test_every_run_config_field_appears_in_the_manifest(run_manifest, path):
    """Derived from RunConfig, not hand-listed: a field added later gets its own
    case automatically and must be recorded rather than silently dropped."""
    cfg, manifest = run_manifest
    try:
        recorded = _dig(manifest["params"], path)
    except KeyError:
        pytest.fail(
            f"manifest.params has no {path!r} — the run that produced these "
            f"outputs cannot be reproduced from them."
        )
    # Absolute paths are recorded relative to the output directory.
    from flexprobe.config_io import portable_config_dict
    expected = _dig(portable_config_dict(cfg, cfg.output_dir), path)
    assert recorded == expected, f"{path}: manifest says {recorded!r}, config was {expected!r}"


def test_the_specific_fields_the_allowlist_dropped(run_manifest):
    """Named explicitly so the reported defects can't quietly come back."""
    _cfg_, manifest = run_manifest
    p = manifest["params"]
    assert p["filters"]["gc_min"] == 0.31
    assert p["filters"]["gc_max"] == 0.81
    assert p["filters"]["max_n_per_half"] == 2
    assert p["filters"]["max_homopolymer_run"] == 6
    assert p["filters"]["reject_junction_base_n"] is False
    # Subset, not exact equality: SpecificityThresholds grows over time (the
    # corrected 10x off-target rule added four fields), and the point of
    # deriving params from the whole config is that new fields are recorded
    # without anyone updating a list. Asserting equality here would turn that
    # success into a failure. The generic per-field test above still checks
    # every field, including any added later.
    assert p["specificity"].items() >= {
        "off_target_pident_min": 77.0, "off_target_aln_min": 19,
        "on_target_pident_min": 91.0, "on_target_aln_min": 22,
    }.items()
    # The off-target rule switches decide which probes survive, so a Methods
    # section needs them on the record.
    for _k in ("off_target_full_half", "off_target_min_mismatch",
               "off_target_require_both_halves", "off_target_same_subject"):
        assert _k in p["specificity"], f"specificity.{_k} missing from manifest"
    assert p["ranking_weights"] == {"gc_proximity": 0.42}
    assert p["group_by_gene"] is False
    assert p["min_isoform_coverage"] == 0.25
    assert p["extra_on_target_ids"] == ["NM_999999.1"]
    assert p["workers"] == 1
    assert p["screen_cache_mode"] == "hits"
    assert "db_paths" in p["screen"]


def test_legacy_param_keys_survive(run_manifest):
    """Readers of older manifests (GUI review dialog, HPC import, coverage
    selection) must keep working."""
    _cfg_, manifest = run_manifest
    p = manifest["params"]
    for key in ("modality_label", "probe_barcode", "probe_barcode_seq",
                "probe_barcodes", "probe_barcodes_applied",
                "gene_lengths", "gene_isoforms",
                "discovery_only", "multi_strain", "consensus_strategy",
                "strain_ids", "variant_positions", "variant_position_list",
                "require_tn_junction", "seed", "primer3_scoring",
                "screen_top_per_gene", "blast_saw_prefiltered_subset",
                "blast_db_sizes_bytes", "tile_length", "tile_step"):
        assert key in p, f"legacy manifest key {key!r} disappeared"
    assert p["modality_label"] == "flex_v2"
    assert p["gene_lengths"].keys() == {"G1", "G2"}
    # 'probe_barcodes' is the config request; '_applied' is what chemistry used
    # (Flex v2 has no probe barcode, so it drops the request entirely).
    assert p["probe_barcodes"] == ["BC003"]
    assert p["probe_barcodes_applied"] == []


def test_applied_probe_barcodes_helper_reads_old_and_new_manifests():
    from flexprobe.select_export import applied_probe_barcodes

    # legacy: only the applied list, under the plain key
    assert applied_probe_barcodes({"probe_barcodes": ["BC001", "BC003"]}) \
        == ["BC001", "BC003"]
    # new: the explicit applied list wins over the config request
    assert applied_probe_barcodes(
        {"probe_barcodes": ["BC003"], "probe_barcodes_applied": []}
    ) == []
    assert applied_probe_barcodes({}) == []


def test_modality_label_helper_reads_old_and_new_manifests():
    from flexprobe.select_export import modality_label

    assert modality_label({"modality": "flex_v1"}) == "flex_v1"          # legacy
    assert modality_label(
        {"modality": {"name": "flex", "version": "v2"}, "modality_label": "flex_v2"}
    ) == "flex_v2"
    # a structured modality with no explicit label still yields one
    assert modality_label(
        {"modality": {"name": "visium", "version": "v1", "spatial": True}}
    ) == "visium_v1"
    assert modality_label({}) == ""


# ---------------------------------------------------------------------------
# (c) the input FASTA is identified by content, not just by path
# ---------------------------------------------------------------------------

def test_input_fasta_is_hashed(run_manifest):
    cfg, manifest = run_manifest
    inputs = manifest["inputs"]
    # Path still recorded, relative to the output directory, and it
    # resolves to the input; the file name is kept beside the hash.
    assert os.path.normpath(os.path.join(cfg.output_dir, inputs["fasta"])) \
        == os.path.abspath(cfg.fasta_path)
    assert inputs["fasta_name"] == os.path.basename(cfg.fasta_path)
    with open(cfg.fasta_path, "rb") as fh:
        expected = hashlib.sha256(fh.read()).hexdigest()
    assert inputs["fasta_sha256"] == expected
    assert int(inputs["fasta_bytes"]) == os.path.getsize(cfg.fasta_path)


# ---------------------------------------------------------------------------
# (d) which build of the engine ran
# ---------------------------------------------------------------------------

def _assert_code_version_identifies_a_build(cv):
    """Either a git checkout (commit + dirty flag) or an installed version.

    Guarding this with ``if cv:`` would let an empty stamp pass. A run from an
    installed package in a container has no checkout to describe, so it must
    still record the installed version rather than ``code_version: {}``.
    """
    from flexprobe import __version__

    assert cv, "code_version is empty: this run cannot be traced to a build"
    assert cv.get("source") in ("git", "installed"), cv
    if cv["source"] == "git":
        assert len(cv["commit"]) == 40
        assert cv["dirty"] in ("true", "false")
    else:
        assert cv["installed_version"] == __version__


def test_code_version_stamp(run_manifest):
    _cfg_, manifest = run_manifest
    _assert_code_version_identifies_a_build(manifest["code_version"])


def test_code_version_outside_a_checkout_is_the_installed_version(monkeypatch):
    """No git (an installed wheel, the container, the packaged app) must still
    name the build, not return ``{}``."""
    import subprocess

    import flexprobe.pipeline as pl

    def no_git(*_a, **_k):
        raise FileNotFoundError("git")

    monkeypatch.setattr(subprocess, "run", no_git)
    cv = pl._code_version()
    _assert_code_version_identifies_a_build(cv)
    assert cv["source"] == "installed"


def test_code_version_ignores_an_unrelated_enclosing_repo(monkeypatch, tmp_path):
    """A package installed into a venv that happens to sit inside SOME git repo
    must not claim that repo's commit. The checkout has to track this package."""
    import subprocess

    import flexprobe.pipeline as pl

    real_run = subprocess.run

    def other_repo(cmd, *a, **k):
        if cmd[:1] == ["git"] and "ls-files" in cmd:
            return subprocess.CompletedProcess(cmd, 1, "", "")
        return real_run(cmd, *a, **k)

    monkeypatch.setattr(subprocess, "run", other_repo)
    assert pl._code_version()["source"] == "installed"


# ---------------------------------------------------------------------------
# (e) the config rides along with the outputs
# ---------------------------------------------------------------------------

def test_run_config_json_round_trips(run_manifest):
    """The CLI writes only design_summary.json and the HPC export tars only
    output/, so without this an imported HPC archive carried no configuration
    at all."""
    cfg, _manifest = run_manifest
    path = os.path.join(cfg.output_dir, "run_config.json")
    assert os.path.exists(path)
    assert load_config(path) == cfg          # re-executable, not just readable


def test_run_config_json_is_in_the_hpc_results_tarball(tmp_path):
    """hpc_export tars `-C output .`, so anything the pipeline writes into the
    output dir rides home. Guard that run_config.json lands there."""
    fa = write_fasta(tmp_path / "g.fa", [("G1", SYNTH_200)])
    cfg = _cfg(tmp_path, fa)
    run_pipeline(cfg)
    assert "run_config.json" in os.listdir(cfg.output_dir)


# ---------------------------------------------------------------------------
# (b) db_sha256s: which database was actually screened against
# ---------------------------------------------------------------------------

def _fake_db(tmp_path, name, payload):
    prefix = str(tmp_path / name)
    for ext, extra in ((".nin", b"i"), (".nsq", b"s"), (".nhr", b"h")):
        with open(prefix + ext, "wb") as fh:
            fh.write(payload + extra)
    return prefix


def test_db_sha256s_populated_for_every_db_used(tmp_path, monkeypatch):
    """The field existed and was serialized from day one but was never assigned
    anywhere in the codebase — every shipped manifest had ``db_sha256s: {}``."""
    import flexprobe.pipeline as pl

    monkeypatch.setattr(fp_blast, "discover_blast",
                        lambda: {"blastn": "b", "makeblastdb": "m"})
    monkeypatch.setattr(fp_blast, "check_version", lambda *a, **k: "2.13.0")
    monkeypatch.setattr(pl, "is_db_prefix", lambda _p: True)

    def fake_run_blastn(_b, query_fa, _db, **_k):
        tsv = query_fa + ".tsv"
        open(tsv, "w").close()
        return tsv

    monkeypatch.setattr(fp_blast, "run_blastn", fake_run_blastn)
    monkeypatch.setattr(fp_blast, "parse_hits_tsv", lambda _t: [])

    db_a = _fake_db(tmp_path, "dbA", b"AAAA")
    db_b = _fake_db(tmp_path, "dbB", b"BBBB")
    fa = write_fasta(tmp_path / "g.fa", [("G1", SYNTH_200)])
    cfg = _cfg(tmp_path, fa,
               screen=ScreenResources(db_paths=[db_a, db_b],
                                    engine=ENGINE_BLAST),
               # k-mer-only gate; a BLAST config must switch it off itself.
               specificity=SpecificityThresholds(
                   off_target_junction_duplex_bp=0))
    run_pipeline(cfg)

    with open(os.path.join(cfg.output_dir, "manifest.json")) as fh:
        digests = json.load(fh)["db_sha256s"]

    # Keyed by the path relative to the output directory.
    digests = {os.path.normpath(os.path.join(cfg.output_dir, k)): v
               for k, v in digests.items()}
    assert set(digests) == {db_a, db_b}
    assert all(len(v) == 64 for v in digests.values())
    assert digests[db_a] != digests[db_b]        # different DBs, different hash

    # Deterministic: rehashing the same DB reproduces the digest, so a
    # published claim is checkable years later.
    assert pl._blast_db_sha256s([db_a])[db_a] == digests[db_a]

    # ...and it is a CONTENT hash — editing an index file changes it, unlike
    # blast_db_sizes_bytes, which any two same-size databases share.
    with open(db_a + ".nsq", "wb") as fh:
        fh.write(b"XXXXs")
    assert pl._blast_db_sha256s([db_a])[db_a] != digests[db_a]


def test_skip_screen_run_claims_no_db_hashes(run_manifest):
    """A --skip-screen run screened against nothing; hashing a DB it never
    opened would be the same fabricated claim as a placeholder version string."""
    _cfg_, manifest = run_manifest
    assert manifest["db_sha256s"] == {}


def test_db_hashing_never_raises_on_a_missing_db(tmp_path):
    import flexprobe.pipeline as pl

    assert pl._blast_db_sha256s([str(tmp_path / "nope")]) == {}
    assert pl._blast_db_sha256s([]) == {}


# ---------------------------------------------------------------------------
# the filter report suffered the same allowlist bug
# ---------------------------------------------------------------------------

def test_filter_report_records_every_threshold(run_manifest):
    cfg, _manifest = run_manifest
    with open(os.path.join(cfg.output_dir, "filter_report.json")) as fh:
        report = json.load(fh)
    assert report["thresholds"] == dataclasses.asdict(cfg.filters)
    # the one the hand-listed version omitted
    assert "reject_junction_base_n" in report["thresholds"]
