"""The container recipes pin what a rebuild could otherwise change.

A rebuilt image could change BLAST silently if a recipe started from the
moving tag ``python:3.11-slim`` and installed apt's ``ncbi-blast+``, whose
version is whatever that Debian release carries on build day. So both recipes
pin the base by digest and ship the vendored BLAST 2.17.0 the desktop app uses.

These tests read the recipe text, so they run anywhere. They cannot prove a
build works; only a build can. The digest is the multi-platform index digest,
so the platform is pinned separately.
"""
import hashlib
import importlib.util
import json
import re
from pathlib import Path

APP = Path(__file__).resolve().parents[1]
CONTAINER = APP / "packaging" / "container"
VENDOR_BLAST = APP / "packaging" / "vendor" / "blast" / "linux-x64" / "bin"
RECIPES = {"Dockerfile": CONTAINER / "Dockerfile",
           "flexiprobe.def": CONTAINER / "flexiprobe.def"}

BASE_RE = re.compile(r"^(?:FROM|From:)\s+(?:--platform=(\S+)\s+)?(\S+)", re.M)
PINNED = re.compile(r"^python:3\.11-slim@sha256:[0-9a-f]{64}$")


def _base(text):
    found = BASE_RE.findall(text)
    assert len(found) == 1, found
    return found[0][1]


def _platform(text):
    return BASE_RE.findall(text)[0][0]


def _code(text):
    """Recipe lines that do something: comments dropped."""
    return "\n".join(l for l in text.splitlines() if not l.lstrip().startswith("#"))


def check(text):
    """Every way a recipe can fail to pin. Empty list means pinned."""
    out = []
    base = _base(text)
    if not PINNED.match(base):
        out.append(f"base image not pinned by digest: {base}")
    code = _code(text)
    if "ncbi-blast+" in code:
        out.append("installs apt's ncbi-blast+")
    if "packaging/vendor/blast/linux-x64/bin" not in code:
        out.append("does not copy the vendored BLAST")
    for tool in ("blastn", "makeblastdb", "blastdbcmd"):
        digest = hashlib.sha256((VENDOR_BLAST / tool).read_bytes()).hexdigest()
        if f"{digest} /opt/blast/bin/{tool}" not in code:
            out.append(f"does not check the sha256 of the vendored {tool}")
    if "grep -q 'blastn: 2.17.0+'" not in code:
        out.append("does not assert BLAST 2.17.0")
    if "write_build_id.py" not in code:
        out.append("does not write a build ID")
    # Without blastdbcmd the BLAST engine silently falls back to conservative
    # counts, so an image lacking it must not build.
    if "blastdbcmd -version | grep -q 'blastdbcmd: 2.17.0+'" not in code:
        out.append("does not require blastdbcmd 2.17.0")
    return out


def test_both_recipes_pin_base_and_blast():
    bad = {name: check(p.read_text()) for name, p in RECIPES.items()}
    assert not any(bad.values()), bad


def test_the_dockerfile_pins_the_platform_of_the_index_digest():
    """The digest names a multi-platform index; amd64 is what BLAST is built for."""
    text = RECIPES["Dockerfile"].read_text()
    assert _platform(text) == "linux/amd64"
    assert "--platform linux/amd64" in text


def test_a_placeholder_digest_is_not_pinned():
    assert check("FROM python:3.11-slim@sha256:REPLACE_WITH_DIGEST\n")[0].startswith(
        "base image not pinned")


def test_both_recipes_use_the_same_base_digest():
    bases = {name: _base(p.read_text()) for name, p in RECIPES.items()}
    assert len(set(bases.values())) == 1, bases
    # The build ID must name the base the recipe actually starts from.
    for name, p in RECIPES.items():
        assert f"write_build_id.py {name} \\\n        {bases[name]}" in p.read_text()


def test_the_check_fails_an_unpinned_recipe():
    old = ("FROM python:3.11-slim\n"
           "RUN apt-get install -y --no-install-recommends ncbi-blast+\n")
    problems = check(old)
    assert any("not pinned" in p for p in problems)
    assert any("ncbi-blast+" in p for p in problems)
    assert any("sha256 of the vendored blastn" in p for p in problems)
    assert any("build ID" in p for p in problems)


def test_the_vendored_blast_is_2_17_0():
    info = (VENDOR_BLAST.parent / "ncbi_package_info").read_text()
    assert "PACKAGE_VERSION: 2.17.0" in info


def _build_id_module():
    spec = importlib.util.spec_from_file_location(
        "write_build_id", CONTAINER / "write_build_id.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_build_id_content_id_ignores_recipe_and_time(tmp_path, monkeypatch):
    mod = _build_id_module()
    for d in ("src", "blast", "mafft"):
        (tmp_path / d).mkdir()
        (tmp_path / d / "f").write_text(d)
    real_tree = mod.tree_sha256
    roots = {"/opt/flexiprobe-src/src": tmp_path / "src",
             "/opt/blast/bin": tmp_path / "blast", "/opt/mafft": tmp_path / "mafft"}
    monkeypatch.setattr(mod, "tree_sha256", lambda r: real_tree(str(roots[r])))
    monkeypatch.setattr(mod.subprocess, "run", lambda *a, **k: type(
        "R", (), {"stdout": "flexiprobes 0.7.0\n"})())
    # EVERY output of main() goes to tmp: it writes /opt/flexiprobe/BUILD_ID
    # and a .pth into the running interpreter's site-packages otherwise.
    monkeypatch.setattr(mod, "OUT_PLAIN", str(tmp_path / "BUILD_ID"))
    monkeypatch.setattr(mod, "site_dir", lambda: str(tmp_path / "site"))
    (tmp_path / "site").mkdir()
    ids = []
    for recipe in ("Dockerfile", "flexiprobe.def"):
        monkeypatch.setattr(mod, "OUT", str(tmp_path / f"{recipe}.json"))
        mod.main(recipe, "python:3.11-slim@sha256:" + "0" * 64)
        ids.append(json.loads((tmp_path / f"{recipe}.json").read_text()))
    assert ids[0]["content_id"] == ids[1]["content_id"]
    assert ids[0]["recipe"] != ids[1]["recipe"]
    assert (tmp_path / "BUILD_ID").read_text().strip() == ids[1]["content_id"]
    assert ids[1]["content_id"] in (tmp_path / "site" / mod.PTH_NAME).read_text()

    # ...and it moves when an input does.
    (tmp_path / "blast" / "f").write_text("other bytes")
    monkeypatch.setattr(mod, "OUT", str(tmp_path / "changed.json"))
    mod.main("Dockerfile", "python:3.11-slim@sha256:" + "0" * 64)
    changed = json.loads((tmp_path / "changed.json").read_text())
    assert changed["content_id"] != ids[0]["content_id"]


def test_the_build_id_reaches_every_python_process(tmp_path):
    """Positive control for the .pth: a fresh interpreter sees the variable."""
    import subprocess
    import sys

    mod = _build_id_module()
    mod.write_pth("abc123", str(tmp_path))
    code = ("import site, os; site.addsitedir(%r); "
            "print(os.environ.get(%r))" % (str(tmp_path), mod.ENV_NAME))
    env = {k: v for k, v in __import__("os").environ.items() if k != mod.ENV_NAME}
    out = subprocess.run([sys.executable, "-c", code], capture_output=True,
                         text=True, env=env, check=True).stdout.strip()
    assert out == "abc123"
    # An operator's explicit value wins over the baked one.
    env[mod.ENV_NAME] = "operator"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True,
                         text=True, env=env, check=True).stdout.strip()
    assert out == "operator"


def test_both_recipes_set_the_build_id_variable_name_the_stamps_read():
    mod = _build_id_module()
    assert mod.ENV_NAME == "FLEXIPROBES_CONTAINER_BUILD_ID"
    assert "FLEXIPROBES_CONTAINER_BUILD_ID" in RECIPES["flexiprobe.def"].read_text()
