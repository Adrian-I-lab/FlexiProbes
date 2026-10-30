"""Tests for the pure-Python GUI stores (no Qt needed).

Covers the persistence layers behind the home-window features:
  • AppSettings        — default output folder preference
  • QueryDraftStore    — save/resume a half-configured query
  • ProjectHistory.update — rename / re-describe a project
"""
import pytest

from flexprobe_gui.settings import AppSettings
from flexprobe_gui.query_drafts import QueryDraftStore
from flexprobe_gui.project_history import ProjectHistory
from flexprobe_gui.run_history import RunHistory, RunRecord


# --- AppSettings -----------------------------------------------------------

def test_settings_default_output_falls_back_then_persists(tmp_path):
    path = str(tmp_path / "settings.json")
    s = AppSettings(path)
    # Falls back to the built-in default when unset.
    assert s.default_output_dir().endswith("flexprobe_out")
    s.set_default_output_dir("/tmp/my_outputs")
    assert s.default_output_dir() == "/tmp/my_outputs"
    # Survives a reload from disk.
    assert AppSettings(path).default_output_dir() == "/tmp/my_outputs"


def test_settings_blank_value_uses_fallback(tmp_path):
    s = AppSettings(str(tmp_path / "settings.json"))
    s.set_default_output_dir("   ")
    assert s.default_output_dir().endswith("flexprobe_out")


def test_settings_a_stale_copy_does_not_overwrite_another_copys_key(tmp_path):
    """Two open copies of the app each hold an AppSettings.

    If each rewrote the whole file from its start-up snapshot, the second
    save would silently revert the first copy's change.
    """
    path = str(tmp_path / "settings.json")
    first = AppSettings(path)
    second = AppSettings(path)             # snapshot taken before any write
    first.set_default_output_dir("/data/first_choice")
    second.set("some_other_key", 7)
    fresh = AppSettings(path)
    assert fresh.default_output_dir() == "/data/first_choice"
    assert fresh.get("some_other_key") == 7


def test_settings_reads_see_another_copys_write(tmp_path):
    path = str(tmp_path / "settings.json")
    reader = AppSettings(path)
    AppSettings(path).set_default_output_dir("/data/changed_elsewhere")
    assert reader.default_output_dir() == "/data/changed_elsewhere"


def test_settings_an_unreadable_file_is_not_treated_as_empty(tmp_path):
    """A read failure must not become 'no settings' and be written back."""
    path = tmp_path / "settings.json"
    s = AppSettings(str(path))
    s.set_default_output_dir("/data/kept")
    path.write_text("{ not json")
    s.set("other", 1)                      # falls back to the last good cache
    fresh = AppSettings(str(path))
    assert fresh.default_output_dir() == "/data/kept"
    assert fresh.get("other") == 1
    assert not [p for p in tmp_path.iterdir() if p.name.endswith(".tmp")]


# --- QueryDraftStore -------------------------------------------------------

def test_draft_save_list_load_roundtrip(tmp_path):
    store = QueryDraftStore(str(tmp_path))
    assert store.list() == []
    store.save("OROV L", {"modality": "Flex v1", "gc_min": 0.3})
    assert store.list() == ["OROV L"]
    loaded = store.load("OROV L")
    assert loaded["modality"] == "Flex v1"
    assert loaded["gc_min"] == 0.3
    assert loaded["name"] == "OROV L"


def test_draft_overwrite_and_delete(tmp_path):
    store = QueryDraftStore(str(tmp_path))
    store.save("draft", {"v": 1})
    store.save("draft", {"v": 2})              # same name overwrites
    assert store.list() == ["draft"]
    assert store.load("draft")["v"] == 2
    assert store.delete("draft") is True
    assert store.delete("draft") is False
    assert store.load("draft") is None


def test_draft_blank_name_rejected(tmp_path):
    store = QueryDraftStore(str(tmp_path))
    with pytest.raises(ValueError):
        store.save("   ", {"v": 1})


def test_draft_unsafe_name_does_not_collide_silently(tmp_path):
    store = QueryDraftStore(str(tmp_path))
    store.save("a/b", {"v": 1})
    # Listing reflects the user-visible name, not the sanitized filename.
    assert store.list() == ["a/b"]
    assert store.load("a/b")["v"] == 1


# --- ProjectHistory.update -------------------------------------------------

def test_project_update_rename_and_describe(tmp_path):
    ph = ProjectHistory(str(tmp_path / "projects.json"))
    a = ph.create("Alpha")
    ph.update(a.project_id, name="Alpha2", description="hello")
    got = ph.get(a.project_id)
    assert got.name == "Alpha2"
    assert got.description == "hello"
    # Persists across reload.
    assert ProjectHistory(str(tmp_path / "projects.json")).get(a.project_id).name == "Alpha2"


def test_project_update_rejects_blank_and_collisions(tmp_path):
    ph = ProjectHistory(str(tmp_path / "projects.json"))
    a = ph.create("Alpha")
    ph.create("Beta")
    with pytest.raises(ValueError):
        ph.update(a.project_id, name="   ")
    with pytest.raises(ValueError):
        ph.update(a.project_id, name="Beta")
    # Renaming to its own (differently-cased) name is fine.
    ph.update(a.project_id, name="ALPHA")
    assert ph.get(a.project_id).name == "ALPHA"


def test_project_update_missing_returns_none(tmp_path):
    ph = ProjectHistory(str(tmp_path / "projects.json"))
    assert ph.update("does-not-exist", name="x") is None


# --- RunHistory.delete -----------------------------------------------------

def test_run_delete_terminal_persists(tmp_path):
    path = str(tmp_path / "runs.json")
    rh = RunHistory(path)
    a = RunRecord(name="a", status="failed", project_id="p1")
    b = RunRecord(name="b", status="complete", project_id="p1")
    rh.append(a)
    rh.append(b)
    assert {r.name for r in rh.list_completed()} == {"a", "b"}
    assert rh.delete(a.run_id) is True
    assert rh.delete("missing") is False
    # Deletion of a terminal run is flushed to disk.
    assert [r.name for r in RunHistory(path).list_completed()] == ["b"]


def test_run_delete_active_not_persisted(tmp_path):
    path = str(tmp_path / "runs.json")
    rh = RunHistory(path)
    running = RunRecord(name="live", status="running")
    rh.append(running)
    assert rh.delete(running.run_id) is True
    assert rh.list_active() == []


def test_run_delete_many_persists_terminal_only(tmp_path):
    path = str(tmp_path / "runs.json")
    rh = RunHistory(path)
    a = RunRecord(name="a", status="complete")
    b = RunRecord(name="b", status="failed")
    c = RunRecord(name="c", status="cancelled")
    live = RunRecord(name="live", status="running")
    for r in (a, b, c, live):
        rh.append(r)
    # Bulk delete two terminal runs + the active one in a single call.
    assert rh.delete_many([a.run_id, c.run_id, live.run_id]) == 3
    assert rh.delete_many([]) == 0
    assert rh.delete_many(["missing"]) == 0
    # Terminal deletions flushed to disk; the surviving terminal run is "b".
    assert [r.name for r in RunHistory(path).list_completed()] == ["b"]
    assert rh.list_active() == []


def test_run_set_project_many_reassigns_and_persists(tmp_path):
    path = str(tmp_path / "runs.json")
    rh = RunHistory(path)
    a = RunRecord(name="a", status="complete", project_id="old")
    b = RunRecord(name="b", status="complete", project_id="")
    rh.append(a)
    rh.append(b)
    assert rh.set_project_many([a.run_id, b.run_id], "new") == 2
    # Re-running with the same target is a no-op (no records change).
    assert rh.set_project_many([a.run_id, b.run_id], "new") == 0
    reloaded = {r.name: r.project_id for r in RunHistory(path).list_completed()}
    assert reloaded == {"a": "new", "b": "new"}
    # Empty target unfiles the runs.
    assert rh.set_project_many([a.run_id], "") == 1
    assert RunHistory(path).get(a.run_id).project_id == ""


def test_project_delete_many(tmp_path):
    ph = ProjectHistory(str(tmp_path / "projects.json"))
    a = ph.create("Alpha")
    b = ph.create("Beta")
    ph.create("Gamma")
    assert ph.delete_many([a.project_id, b.project_id]) == 2
    assert ph.delete_many([]) == 0
    assert ph.delete_many(["missing"]) == 0
    remaining = [p.name for p in ProjectHistory(str(tmp_path / "projects.json")).list()]
    assert remaining == ["Gamma"]


def test_run_record_query_draft_roundtrips(tmp_path):
    path = str(tmp_path / "runs.json")
    rh = RunHistory(path)
    rec = RunRecord(
        name="r", status="complete",
        query_draft={"modality": "Flex v1", "blast_dbs": "/refs/human_tx"},
    )
    rh.append(rec)
    loaded = RunHistory(path).get(rec.run_id)
    assert loaded.query_draft["modality"] == "Flex v1"
    assert loaded.query_draft["blast_dbs"] == "/refs/human_tx"


def test_review_table_surfaces_the_evidence_behind_a_verdict():
    """Pass/fail alone can hide defects, so the evidence columns must be present.

    A specificity bug can produce a results table that looks normal: the
    verdict is there, what produced it is not.
    """
    from flexprobe_gui.review_dialog import _BASE_COLS
    keys = {c.source for c in _BASE_COLS}
    for k in ("offtarget_recount_basis", "offtarget_min_qcov",
              "offtarget_worst_mismatches", "offtarget_worst_junction_mm",
              "probe_role"):
        assert k in keys, f"review table no longer shows {k}"


def test_rerun_draft_preserves_the_pan_strain_tolerance():
    """A re-run must not silently revert to exact-identity conserved windows.

    On the OROV L segment that reset is the difference between 44 orderable
    probes and none at all, and it would happen without a visible change.
    """
    from flexprobe_gui.home_window import HomeWindow
    from flexprobe_gui.run_history import RunRecord

    class _Stub:                       # _record_to_draft only needs settings
        class settings:
            @staticmethod
            def default_output_dir(): return "/tmp/out"

    rec = RunRecord(
        name="r", modality_label="Flex v1", fasta_path="in.fa",
        output_dir="out", project_id=None,
        config_summary={"conserved_max_degenerate_per_half": 1,
                        "max_n_per_half": 1},
    )
    draft = HomeWindow._record_to_draft(_Stub(), rec)
    assert draft["conserved_degen"] == 1
