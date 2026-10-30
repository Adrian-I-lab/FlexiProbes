"""Byte-identical output gate for the ``run_pipeline`` break-up.

Two layers, deliberately:

**Layer 1 — determinism (always runs, everywhere).** Each scenario is executed
TWICE in the same interpreter and the two digests must match. This needs no
committed baseline, so it cannot be invalidated by a different primer3 build,
and it catches the refactor failure mode that unit tests miss almost by
definition: state that leaks between stages, or an extraction that makes output
depend on dict/set iteration order.

**Layer 2 — the golden baseline (REQUIRED; must be recorded for this
environment).** Digests committed in ``golden/baseline.json`` are compared
against this run. This is the layer that actually gates the refactor: it says
"the code you just moved produces the same bytes it did before you moved it".

A MISSING BASELINE IS A FAILURE, NOT A SKIP. If it skipped, a CI leg whose
Python version has no recorded baseline would report green having verified
nothing, and sail past the only byte-level gate the project has. A gate that
silently disables itself on an unrecognised environment is not a gate. The
failure names the environment and the command that records it.

If you genuinely need to run without a baseline — bisecting on a throwaway
interpreter, or a CI leg whose baseline has not been generated YET and is
tracked as an outstanding action — set ``FLEXPROBE_GOLDEN_OPTIONAL=1``. That
downgrades the missing baseline to a skip. It is deliberately an explicit,
greppable opt-out: someone has to write it down, in a file that gets
reviewed, rather than have the suite quietly decide for them.

WHY THE BASELINE IS KEYED BY ENVIRONMENT. Tm comes from ``primer3.calc_tm``
(``_utils.py:112``), so the numbers in ``results.csv`` are a function of the
installed primer3 build, and float formatting is a function of the Python minor
version. A digest is therefore only comparable within the environment that
produced it. Rather than scrub the science columns out — which would gut the
gate — the baseline stores **one block of digests per environment**
(``py3.11+primer3-2.3.0``, ``py3.13+primer3-2.3.0``, …). A run compares against
its own block if one exists, and skips with a loud message if not. Recording a
new environment never disturbs the ones already there, so the same file can gate
the developer's venv, a conda interpreter and each CI leg at once.

REGENERATING / RECORDING AN ENVIRONMENT::

    cd app && FLEXPROBE_UPDATE_GOLDEN=1 .venv/bin/python -m pytest \\
        tests/test_golden_output.py

That MERGES this interpreter's block into the file. Recording a *new*
environment is routine — and, since the gate now fails without one, required
for every interpreter the CI matrix runs. Overwriting an *existing* one means you are asserting
that an output change is intended and reviewed — during the ``run_pipeline``
break-up the correct number of those is zero, because a digest change there
means the refactor was not a refactor.
"""
import json
import os
import pathlib
import sys

import pytest

from golden_scenarios import SCENARIOS, run_scenario

BASELINE_PATH = pathlib.Path(__file__).parent / "golden" / "baseline.json"

_UPDATE = os.environ.get("FLEXPROBE_UPDATE_GOLDEN") == "1"

# Explicit, greppable escape hatch for "no baseline recorded for this
# interpreter yet". Absent it, a missing baseline FAILS — see the module
# docstring for why the old silent skip was the actual defect.
_OPTIONAL = os.environ.get("FLEXPROBE_GOLDEN_OPTIONAL") == "1"


def _no_baseline(message: str) -> None:
    """Fail on a missing baseline; skip only under the explicit opt-out."""
    if _OPTIONAL:
        pytest.skip(
            "GOLDEN GATE DISABLED for this run "
            "(FLEXPROBE_GOLDEN_OPTIONAL=1). Nothing was byte-compared. "
            + message
        )
    pytest.fail(message)


def _env_key() -> str:
    """Identify what the digests are a function of, beyond flexprobe's code.

    Used as the baseline dict key, so it must be a stable, readable string:
    the file is reviewed by eye in a PR.
    """
    try:
        import primer3
        p3 = getattr(primer3, "__version__", "unknown")
    except Exception:                                    # pragma: no cover
        p3 = "absent"
    return f"py{sys.version_info.major}.{sys.version_info.minor}+primer3-{p3}"


def _load_baseline() -> dict:
    if not BASELINE_PATH.exists():
        return {}
    with open(BASELINE_PATH) as fh:
        return json.load(fh)


@pytest.fixture(scope="module")
def baseline() -> dict:
    return _load_baseline()


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_scenario_is_deterministic(name, tmp_path):
    """Layer 1: the same scenario run twice gives byte-identical outputs.

    Environment-independent, so this is the check that keeps running on every
    CI leg even when the golden baseline is skipped.
    """
    first = run_scenario(name, str(tmp_path / "a"))
    second = run_scenario(name, str(tmp_path / "b"))

    assert first, f"scenario {name!r} produced no output files at all"
    assert first == second, (
        f"{name}: two runs of identical code produced different output.\n"
        + _describe(first, second, "run A", "run B")
    )


@pytest.mark.parametrize("name", sorted(SCENARIOS))
def test_scenario_matches_golden_baseline(name, tmp_path, baseline):
    """Layer 2: outputs match the committed baseline byte for byte."""
    if _UPDATE:
        pytest.skip("regenerating baseline (FLEXPROBE_UPDATE_GOLDEN=1)")

    if not baseline:
        pytest.fail(
            f"{BASELINE_PATH} is missing. Generate it with:\n"
            "  cd app && FLEXPROBE_UPDATE_GOLDEN=1 python -m pytest "
            "tests/test_golden_output.py"
        )

    env = _env_key()
    block = baseline.get("environments", {}).get(env)
    if block is None:
        _no_baseline(
            f"no golden block recorded for {env!r} (recorded: "
            f"{sorted(baseline.get('environments', {}))}).\n"
            "Tm comes from primer3 and float formatting from the "
            "interpreter, so digests are not comparable across environments "
            "— this one has to record its own. On a machine running THIS "
            "interpreter:\n"
            "  cd app && FLEXPROBE_UPDATE_GOLDEN=1 python -m pytest "
            "tests/test_golden_output.py\n"
            "then commit tests/golden/baseline.json. Recording a new "
            "environment MERGES; it never disturbs an existing block.\n"
            "To run without the gate anyway set FLEXPROBE_GOLDEN_OPTIONAL=1 "
            "— that leaves only the determinism layer, which cannot detect "
            "an output CHANGE, only non-determinism."
        )

    expected = block.get(name)
    if expected is None:
        pytest.fail(
            f"scenario {name!r} has no entry under {env!r} in "
            f"{BASELINE_PATH.name}; re-record that environment."
        )

    actual = run_scenario(name, str(tmp_path))
    assert actual == expected, (
        f"{name}: output no longer matches the golden baseline.\n"
        + _describe(expected, actual, "baseline", "now")
        + "\n\nIf this change is intended and reviewed, regenerate with "
          "FLEXPROBE_UPDATE_GOLDEN=1. A pure refactor must not move a "
          "byte."
    )


def _describe(left: dict, right: dict, lname: str, rname: str) -> str:
    """Name the files that differ, so a failure points at a stage."""
    lines = []
    for key in sorted(set(left) | set(right)):
        lv, rv = left.get(key), right.get(key)
        if lv == rv:
            continue
        if lv is None:
            lines.append(f"  + {key}  (only in {rname})")
        elif rv is None:
            lines.append(f"  - {key}  (only in {lname})")
        else:
            lines.append(f"  ~ {key}  {lname}={lv[:12]} {rname}={rv[:12]}")
    return "\n".join(lines) or "  (no per-file difference — dict compare bug?)"


def test_baseline_covers_every_scenario(baseline):
    """A scenario added without regenerating the baseline is a silent hole."""
    if _UPDATE or not baseline:
        pytest.skip("baseline being (re)generated")
    block = baseline.get("environments", {}).get(_env_key())
    if block is None:
        _no_baseline(
            f"no golden block recorded for {_env_key()!r}; record it with "
            "FLEXPROBE_UPDATE_GOLDEN=1 python -m pytest "
            "tests/test_golden_output.py"
        )
    missing = sorted(set(SCENARIOS) - set(block))
    assert not missing, (
        f"scenarios with no golden entry under {_env_key()!r}: {missing}. "
        "Re-record that environment."
    )


@pytest.mark.skipif(not _UPDATE, reason="set FLEXPROBE_UPDATE_GOLDEN=1")
def test_zz_regenerate_baseline(tmp_path):
    """Not a test: the baseline writer, parked behind an env var.

    Named ``zz`` so alphabetical collection runs it after the scenarios above
    have had their chance to report. Lives here rather than in a script so it
    shares the suite's ``FLEXPROBE_HOME`` isolation fixture — a standalone
    script wrote into the developer's real home once already (see conftest).

    MERGES rather than overwrites: recording a second interpreter must not
    silently delete the block that is currently gating someone else's work.
    """
    payload = _load_baseline() or {}
    payload["_comment"] = (
        "Byte-identical output baseline, keyed "
        "by environment because Tm comes from primer3 and float formatting "
        "from the interpreter. Record an environment with "
        "FLEXPROBE_UPDATE_GOLDEN=1 python -m pytest tests/test_golden_output.py "
        "-- adding a new key is routine; CHANGING an existing block asserts "
        "that an output change is intended and reviewed."
    )
    env = _env_key()
    existing = payload.setdefault("environments", {}).get(env)
    fresh = {
        name: run_scenario(name, str(tmp_path / name))
        for name in sorted(SCENARIOS)
    }
    payload["environments"][env] = fresh

    BASELINE_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(BASELINE_PATH, "w") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")

    verb = "added" if existing is None else "REPLACED"
    print(f"\n{verb} block {env!r} in {BASELINE_PATH} "
          f"({len(fresh)} scenarios; file now holds "
          f"{len(payload['environments'])} environment(s))")
    if existing is not None:
        # Only scenarios that were ALREADY recorded can have "changed" — a
        # newly added scenario has no prior digest to move away from, and
        # flagging it would cry wolf on the routine case.
        changed = sorted(k for k in set(existing) & set(fresh)
                         if existing[k] != fresh[k])
        added = sorted(set(fresh) - set(existing))
        if added:
            print(f"  new scenarios recorded: {added}")
        if changed:
            print(f"  WARNING: digests CHANGED for existing scenarios: "
                  f"{changed}\n"
                  "  That is an output change, not a re-record. Make sure it "
                  "was intended and reviewed.")
