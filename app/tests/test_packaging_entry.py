"""Guard the frozen-app entry point.

The packaged .app spawns ProcessPoolExecutor workers via multiprocessing's
"spawn" method, which re-launches the frozen binary. The entry point MUST call
multiprocessing.freeze_support() before launching the GUI, or every worker opens
a new app window (the "--workers 8 spawns 8 apps" bug). This test fails loudly
if that call is ever dropped.
"""
import os


_PKG = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "packaging",
)
ENTRY = os.path.join(_PKG, "flexprobe_main.py")
SPEC = os.path.join(_PKG, "flexprobe.spec")


def test_entry_calls_freeze_support_before_launch():
    with open(ENTRY) as fh:
        src = fh.read()
    assert "multiprocessing.freeze_support()" in src, (
        "frozen entry point must call multiprocessing.freeze_support() or "
        "ProcessPoolExecutor workers will each open a new GUI window"
    )
    # It must run before the GUI is launched: freeze_support() appears before
    # the _launch() call site (rindex skips the `def _launch()` line above it).
    fs = src.index("multiprocessing.freeze_support()")
    launch_call = src.rindex("_launch()")
    assert fs < launch_call, (
        "freeze_support() must be called before launching the GUI"
    )


def test_spec_bundles_primer3_compiled_libs():
    # primer3-py ships two C extensions (p3helpers + thermoanalysis); collecting
    # only submodules left p3helpers.so out, so `import primer3` failed in the
    # frozen app and every Tm/hairpin came back n/a. The spec must collect
    # primer3's dynamic libs (and the thermo-config data for hairpin scoring).
    with open(SPEC) as fh:
        src = fh.read()
    assert 'collect_dynamic_libs("primer3")' in src, (
        "spec must collect primer3's compiled extensions or Tm/hairpin break "
        "in the packaged app"
    )
    assert "primer3_binaries" in src and "binaries=blast_binaries + primer3_binaries" in src
    assert 'collect_data_files(' in src and "primer3_config" in src
