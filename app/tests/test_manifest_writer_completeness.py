"""Every RunManifest field must actually reach manifest.json.

WHY THIS EXISTS. ``writers.manifest.write_manifest`` builds an EXPLICIT dict of
keys rather than calling ``asdict(manifest)``. A field added to the dataclass and
not added to that dict is silently dropped, and nothing notices: a new field
can be populated at every manifest site while the golden tests stay green,
precisely BECAUSE the output has not changed.

``_StageTimer._write`` follows the same fixed-key-set pattern, and the same
risk applies there.

This test compares the dataclass's field names against the written JSON's keys,
so the next field to be added cannot go missing quietly.
"""
from __future__ import annotations

import dataclasses
import json

from flexprobe.models import RunManifest
from flexprobe.writers.manifest import write_manifest


def test_every_manifest_field_is_written(tmp_path):
    """No RunManifest field may be missing from the written JSON."""
    manifest = RunManifest()
    out = write_manifest(manifest, str(tmp_path / "manifest.json"))
    written = json.loads(open(out).read())

    declared = {f.name for f in dataclasses.fields(RunManifest)}
    present = set(written)

    dropped = sorted(declared - present)
    assert not dropped, (
        f"RunManifest declares {sorted(declared)} but manifest.json contains "
        f"{sorted(present)}.\nDROPPED SILENTLY: {dropped}\n"
        f"write_manifest() uses an explicit allow-list dict, not "
        f"asdict(manifest), so a new field must be added there too. A green "
        f"golden suite does NOT catch this: the output simply does not change."
    )


def test_no_unexpected_extra_keys():
    """The written JSON may not invent keys the dataclass does not declare.

    The other direction matters too. A key written but not declared cannot be
    set by a caller, so it is either dead or computed in the writer, and both
    belong on the dataclass where they are discoverable.
    """
    import tempfile, os
    with tempfile.TemporaryDirectory() as d:
        out = write_manifest(RunManifest(), os.path.join(d, "m.json"))
        written = set(json.loads(open(out).read()))
    declared = {f.name for f in dataclasses.fields(RunManifest)}
    extra = sorted(written - declared)
    assert not extra, (
        f"manifest.json contains keys RunManifest does not declare: {extra}. "
        f"Add them to the dataclass so they are discoverable, or stop writing "
        f"them."
    )


def test_consensus_field_specifically_round_trips():
    """The field whose absence motivated this file.

    Guards the value, not just the key: a writer that emitted ``{}`` regardless
    of what was passed would satisfy the key check above and still lose the
    aligner.
    """
    import tempfile, os
    payload = {"aligner": "MAFFT v7.526 (bundled)", "n_strains": 3}
    with tempfile.TemporaryDirectory() as d:
        out = write_manifest(RunManifest(consensus=payload),
                             os.path.join(d, "m.json"))
        written = json.loads(open(out).read())
    assert written["consensus"] == payload, (
        f"consensus did not round-trip: wrote {payload!r}, read back "
        f"{written.get('consensus')!r}"
    )
