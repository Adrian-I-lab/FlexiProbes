"""Config single-source-of-truth defaults + RunConfig.validate()."""
import pytest

from flexprobe import constants
from flexprobe.errors import ConfigError
from flexprobe.models import FilterThresholds, RunConfig


def test_filter_defaults_track_constants():
    # Bare defaults must match the documented 10x spec in constants.py (the
    # single source of truth), not drift from the CLI/GUI defaults.
    thr = FilterThresholds()
    assert thr.gc_min == constants.GC_MIN_DEFAULT == 0.44
    assert thr.gc_max == constants.GC_MAX_DEFAULT == 0.72
    assert thr.max_n_per_half == 1
    assert thr.max_homopolymer_run == constants.MAX_HOMOPOLYMER_RUN_DEFAULT


def test_default_runconfig_is_valid():
    RunConfig().validate()  # must not raise


@pytest.mark.parametrize("mutate", [
    lambda c: setattr(c.filters, "gc_min", 0.8),                 # gc_min > gc_max
    lambda c: setattr(c.filters, "gc_max", 1.5),                 # > 1.0
    lambda c: setattr(c.filters, "max_n_per_half", -1),
    lambda c: setattr(c, "tile_length", 51),                     # odd
    lambda c: setattr(c, "tile_length", 0),
    lambda c: setattr(c, "tile_step", 0),
    lambda c: setattr(c, "consensus_strategy", "bogus"),
    lambda c: setattr(c, "auto_select_per_gene", -3),
])
def test_validate_rejects_bad_config(mutate):
    cfg = RunConfig()
    mutate(cfg)
    with pytest.raises(ConfigError):
        cfg.validate()


def test_validate_message_lists_problems():
    cfg = RunConfig()
    cfg.filters.gc_min = 0.9
    cfg.tile_length = 7
    with pytest.raises(ConfigError) as exc:
        cfg.validate()
    msg = str(exc.value)
    assert "GC range" in msg and "tile_length" in msg
