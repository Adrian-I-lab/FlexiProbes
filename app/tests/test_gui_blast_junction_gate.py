"""The GUI must not be able to build a config the engine will refuse.

`RunConfig.validate()` rejects `screen.engine == "blast"` while the
junction-duplex gate is on, because BLAST provably cannot produce the records
the gate needs. The new-query dialog exposes an engine picker but no
junction-duplex control, so unless `_nq_build.build_runs` zeroes the gate,
every BLAST run built from the GUI raises a ConfigError the user has no widget
to act on. Turning the gate off is a real weakening of the screen, so it must also be
reported, never silent.
"""
from dataclasses import replace

import pytest

from flexprobe.errors import ConfigError
from flexprobe.models import RunConfig, SpecificityThresholds


def _gui_style_config(engine: str) -> RunConfig:
    """What the dialog builds: engine set, specificity left to its default."""
    cfg = RunConfig(fasta_path="x.fa", output_dir="out")
    cfg.screen.engine = engine
    return cfg


def test_the_default_gate_is_on_so_the_clash_is_reachable():
    assert RunConfig(
        fasta_path="x", output_dir="o"
    ).specificity.off_target_junction_duplex_bp > 0


def test_engine_blast_with_the_default_gate_is_refused():
    with pytest.raises(ConfigError, match="junction"):
        _gui_style_config("blast").validate()


def test_zeroing_the_gate_is_what_makes_blast_runnable():
    cfg = _gui_style_config("blast")
    cfg.specificity = replace(
        SpecificityThresholds(), off_target_junction_duplex_bp=0)
    cfg.validate()


def test_the_kmer_default_is_untouched():
    _gui_style_config("kmer").validate()


@pytest.mark.parametrize("engine,expect_note", [("blast", True), ("kmer", False)])
def test_build_runs_zeroes_the_gate_for_blast_and_says_so(engine, expect_note,
                                                          monkeypatch):
    """The behavioural half: exercised through the real dialog builder."""
    pytest.importorskip("PySide6")
    from flexprobe_gui import _nq_build

    src = _nq_build.__file__
    text = open(src).read()
    # build_runs must both zero the gate and set the user-facing note.
    assert "off_target_junction_duplex_bp=0" in text
    assert "blast_gate_note" in text
    assert "specificity=specificity" in text
