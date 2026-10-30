"""Every version declaration must name the same version.

``pyproject.toml``, both Python packages and ``CITATION.cff`` each carry a
version string. If they disagree, a fresh citation or the GUI About box names a
different version from the engine that ships.

``CITATION.cff`` is not Python, so it is checked with a plain regex, which is
as strong a guarantee as the hand-written string it guards deserves.
"""
from __future__ import annotations

import pathlib
import re

from flexprobe import __version__ as engine_version
from flexprobe_gui import __version__ as gui_version

APP_DIR = pathlib.Path(__file__).resolve().parents[1]
REPO_ROOT = APP_DIR.parent
PYPROJECT = APP_DIR / "pyproject.toml"
CITATION_CFF = REPO_ROOT / "CITATION.cff"


def _pyproject_version() -> str:
    """Pull ``version = "X.Y.Z"`` out of ``[project]`` without a TOML parser.

    The project's floor is Python 3.11, so ``tomllib`` would also work. A
    single-purpose regex is enough for one field.
    """
    text = PYPROJECT.read_text(encoding="utf-8")
    match = re.search(r'(?m)^version\s*=\s*"([^"]+)"', text)
    assert match, "pyproject.toml has no top-level version = \"...\" line"
    return match.group(1)



def test_pyproject_engine_and_gui_versions_agree():
    """The three Python-visible version declarations must be identical."""
    pyproject_version = _pyproject_version()
    assert engine_version == pyproject_version, (
        f"flexprobe.__version__ ({engine_version!r}) != "
        f"pyproject.toml version ({pyproject_version!r})"
    )
    assert gui_version == pyproject_version, (
        f"flexprobe_gui.__version__ ({gui_version!r}) != "
        f"pyproject.toml version ({pyproject_version!r})"
    )


def test_citation_cff_version_matches_engine():
    """CITATION.cff must name a version a reader can actually check out."""
    text = CITATION_CFF.read_text(encoding="utf-8")
    match = re.search(r'(?m)^version:\s*"([^"]+)"', text)
    assert match, "CITATION.cff has no version: \"...\" line"
    cff_version = match.group(1)
    assert cff_version == engine_version, (
        f"CITATION.cff version ({cff_version!r}) != "
        f"flexprobe.__version__ ({engine_version!r})"
    )
