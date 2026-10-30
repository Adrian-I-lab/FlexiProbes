"""Round-trip tests for RunConfig JSON serialization (HPC-export backbone)."""
import json

from flexprobe.config_io import (
    config_from_dict, config_to_dict, load_config, save_config,
)
from flexprobe.models import (
    ScreenResources, ExportSettings, FilterThresholds, ModalityDefinition,
    RunConfig, SpecificityThresholds,
)


def _rich_config() -> RunConfig:
    return RunConfig(
        fasta_path="in.fa",
        output_dir="out",
        modality=ModalityDefinition(name="flex", version="v1", spatial=False),
        filters=FilterThresholds(gc_min=0.5, gc_max=0.7, require_tn_junction=True),
        screen=ScreenResources(db_paths=["/scratch/human", "/scratch/virus"]),
        export=ExportSettings(write_idt_opools=False),
        specificity=SpecificityThresholds(off_target_pident_min=0.85),
        group_by_gene=True,
        discovery_only=True,
        auto_select_per_gene=3,
        min_isoform_coverage=0.8,
        ranking_weights={"gc_balance": 2.0, "isoform_coverage": 1.5},
        extra_on_target_ids=["ENST1", "ENST2"],
        probe_barcodes=["BC001", "BC002"],
    )


def test_round_trip_preserves_everything():
    cfg = _rich_config()
    assert config_from_dict(config_to_dict(cfg)) == cfg


def test_round_trip_default_config():
    cfg = RunConfig(fasta_path="x.fa")
    assert config_from_dict(config_to_dict(cfg)) == cfg


def test_ranking_weights_none_round_trips():
    cfg = RunConfig(fasta_path="x.fa", ranking_weights=None)
    out = config_from_dict(config_to_dict(cfg))
    assert out.ranking_weights is None


def test_dict_is_json_serializable():
    # The whole point is to write it to disk for the cluster.
    json.dumps(config_to_dict(_rich_config()))


def test_from_dict_ignores_unknown_keys():
    d = config_to_dict(RunConfig(fasta_path="x.fa"))
    d["a_future_field"] = 123
    d["filters"]["another_future_field"] = "ignored"
    cfg = config_from_dict(d)
    assert cfg.fasta_path == "x.fa"


def test_from_dict_tolerates_missing_blocks():
    cfg = config_from_dict({"fasta_path": "x.fa"})
    # Missing nested blocks fall back to dataclass defaults.
    assert cfg.filters.gc_min == FilterThresholds().gc_min
    assert cfg.screen.db_paths == []


def test_save_and_load_file(tmp_path):
    cfg = _rich_config()
    path = tmp_path / "config.json"
    save_config(cfg, str(path))
    assert path.exists()
    assert load_config(str(path)) == cfg


def test_legacy_blast_block_still_loads(tmp_path):
    """0.6.0 renamed the config block "blast" -> "screen".

    The key moved; the meaning did not. A config written by 0.5.x must still
    load, because a saved GUI config is not a run artefact the user expects to
    have to recreate.
    """
    import json

    from flexprobe.config_io import load_config

    p = tmp_path / "legacy.json"
    p.write_text(json.dumps({
        "fasta_path": "x", "output_dir": "y",
        "modality": {"name": "flex", "version": "v2"},
        "blast": {"engine": "kmer", "db_paths": ["/db"], "query_chunk": 123},
    }))
    cfg = load_config(str(p))
    assert cfg.screen.engine == "kmer"
    assert cfg.screen.db_paths == ["/db"]
    assert cfg.screen.query_chunk == 123


def test_new_screen_block_wins_over_a_stale_blast_block(tmp_path):
    """Both present means a half-migrated file; the current key is the truth."""
    import json

    from flexprobe.config_io import load_config

    p = tmp_path / "both.json"
    p.write_text(json.dumps({
        "fasta_path": "x", "output_dir": "y",
        "modality": {"name": "flex", "version": "v2"},
        "blast": {"engine": "blast", "query_chunk": 1},
        "screen": {"engine": "kmer", "query_chunk": 999},
    }))
    cfg = load_config(str(p))
    assert cfg.screen.engine == "kmer"
    assert cfg.screen.query_chunk == 999
