"""The engine is frozen, and an engine edit must change the manifest too.

``flexprobe/ENGINE_FREEZE.json`` records the SHA-256 of every file in
``engine_identity.ENGINE_FILES`` and the combined ``engine_sha256`` that every
cache stamp carries. An engine file that changes without the manifest changing
in the same commit fails here, naming the file. Unfreezing is deliberate: it
means regenerating the manifest, which a reviewer sees in the diff::

    cd app && python -c "from flexprobe.engine_identity import write_freeze_manifest; write_freeze_manifest()"

GUI, docs, writers and other non-engine files are outside the manifest, so
changing them never touches it.
"""
import hashlib
import json

from flexprobe import engine_identity
from flexprobe.engine_identity import ENGINE_FILES, FREEZE_MANIFEST, PACKAGE_ROOT


def _manifest():
    with open(FREEZE_MANIFEST, encoding="utf-8") as fh:
        return json.load(fh)


def test_the_manifest_lists_exactly_the_engine_files():
    assert set(_manifest()["files"]) == set(ENGINE_FILES)


def test_no_engine_file_changed_without_the_manifest():
    recorded = _manifest()["files"]
    checked, changed = 0, []
    for rel in sorted(ENGINE_FILES):
        data = (PACKAGE_ROOT / rel).read_bytes()
        checked += 1
        if hashlib.sha256(data).hexdigest() != recorded[rel]:
            changed.append(rel)
    # A comparison that examined nothing would pass. Refuse that.
    assert checked == len(ENGINE_FILES) and checked > 20, checked
    assert not changed, (
        f"engine file(s) changed but ENGINE_FREEZE.json did not: {changed}. "
        "The engine is frozen. If this change is a deliberate unfreeze, "
        "regenerate the manifest in the same commit (see this file's "
        "docstring) and bump the engine version.")


def test_the_combined_hash_is_the_one_every_cache_stamp_carries():
    engine_identity.engine_sha256.cache_clear()
    assert _manifest()["engine_sha256"] == engine_identity.engine_sha256()


def test_the_manifest_names_an_engine_version():
    version = _manifest()["engine_version"]
    parts = version.split(".")
    assert len(parts) == 3 and all(p.isdigit() for p in parts), version
