"""Guards for the redistribution-licensing invariants.

These are cheap file-level checks, but each one corresponds to a defect that
actually shipped: a GPLv3 §6 written offer naming the wrong version, and a
binary that redistributed primer3-py (GPL-2.0-or-later) and Qt without their
license texts because the build relied on PyInstaller carrying each wheel's
``*.dist-info`` (it does not — it collects modules, not metadata).
"""
from __future__ import annotations

import pathlib
import re

import pytest

from flexprobe import __version__

APP_DIR = pathlib.Path(__file__).resolve().parents[1]
TPL_DIR = APP_DIR / "THIRD-PARTY-LICENSES"
SOURCE_OFFER = APP_DIR / "SOURCE-OFFER.txt"
SPEC = APP_DIR / "packaging" / "flexprobe.spec"

# Every component the desktop binary redistributes as compiled code, with the
# vendored text that must accompany it.
BUNDLED_LICENCE_TEXTS = {
    "primer3-py": "primer3-py-GPL-2.0.txt",
    "biopython": "biopython-BIOPYTHON-LICENSE.txt",
    "numpy": "numpy-BSD-3-Clause.txt",
    "qt/pyside6 (LGPL option)": "qt-pyside6-LGPL-3.0.txt",
    "ncbi blast+": "PUBLIC-DOMAIN-NCBI-BLAST.txt",
}


@pytest.mark.parametrize(
    "component,filename", sorted(BUNDLED_LICENCE_TEXTS.items())
)
def test_bundled_component_licence_text_is_vendored(component, filename):
    """Each bundled component's licence text exists and is non-trivial."""
    path = TPL_DIR / filename
    assert path.is_file(), (
        f"{component}: {filename} missing from THIRD-PARTY-LICENSES/. "
        f"Run packaging/vendor_licenses.py."
    )
    assert len(path.read_text(encoding="utf-8").strip()) > 500, (
        f"{component}: {filename} is suspiciously short for a licence text."
    )


def test_primer3_text_is_the_gpl2():
    """primer3-py is the copyleft the whole project inherits — check the text."""
    text = (TPL_DIR / "primer3-py-GPL-2.0.txt").read_text(encoding="utf-8")
    assert "GNU GENERAL PUBLIC LICENSE" in text
    assert "Version 2, June 1991" in text


def test_qt_lgpl_text_is_the_lgpl3():
    text = (TPL_DIR / "qt-pyside6-LGPL-3.0.txt").read_text(encoding="utf-8")
    assert "GNU LESSER GENERAL PUBLIC LICENSE" in text
    assert "Version 3, 29 June 2007" in text


def test_source_offer_is_a_template_not_a_hardcoded_version():
    """GPLv3 §6 offer must name the version of the object code it covers.

    It is rendered from ``flexprobe.__version__`` by the PyInstaller spec, so
    the tree copy must carry the placeholder and must NOT carry a stale
    hand-typed version.
    """
    text = SOURCE_OFFER.read_text(encoding="utf-8")
    assert "@FLEXIPROBES_VERSION@" in text, (
        "SOURCE-OFFER.txt lost its version placeholder; the spec will refuse "
        "to build."
    )
    stale = re.findall(r"\bversion \d+\.\d+\.\d+\b", text)
    assert not stale, f"hand-typed version(s) in SOURCE-OFFER.txt: {stale}"


def test_source_offer_renders_to_the_current_version():
    """Reproduce the spec's substitution and check the result."""
    text = SOURCE_OFFER.read_text(encoding="utf-8")
    rendered = text.replace("@FLEXIPROBES_VERSION@", __version__)
    assert f"version {__version__}" in rendered
    assert "@FLEXIPROBES_VERSION@" not in rendered


def test_source_offer_makes_no_release_tag_promise():
    """The offer must not depend on a route that does not resolve.

    A promise that "a release tag matches each binary version" cannot be kept
    unless every binary version has a tag.
    """
    text = SOURCE_OFFER.read_text(encoding="utf-8").lower()
    assert "a release tag matches each binary version" not in text


def test_spec_ships_the_licence_texts():
    """The spec must copy THIRD-PARTY-LICENSES/ into the bundle."""
    text = SPEC.read_text(encoding="utf-8")
    assert 'datas.append((_tpl, "licenses/THIRD-PARTY-LICENSES"))' in text
    assert '"@FLEXIPROBES_VERSION@" not in _offer_text' in text, (
        "spec no longer fails closed when the offer template is broken"
    )


def test_notice_and_third_party_readme_agree_on_the_qt_option():
    """Task 6: NOTICE, THIRD-PARTY-LICENSES/README and the spec told three
    different stories about which Qt licence option is in force."""
    notice = (APP_DIR / "NOTICE").read_text(encoding="utf-8")
    readme = (TPL_DIR / "README.md").read_text(encoding="utf-8")
    spec = SPEC.read_text(encoding="utf-8")
    for name, text in (("NOTICE", notice), ("TPL README", readme), ("spec", spec)):
        assert "GPL-3.0 option" in text, (
            f"{name} does not state that Qt is used under the GPL-3.0 option"
        )


def test_idt_is_attributed():
    """IDT trademarks are used throughout (oPools exports) — must be noticed."""
    notice = (APP_DIR / "NOTICE").read_text(encoding="utf-8")
    assert "Integrated DNA Technologies" in notice
    assert "oPools" in notice


def test_notice_dependency_table_matches_constraints():
    """Every pinned distribution appears in the NOTICE component table."""
    notice = (APP_DIR / "NOTICE").read_text(encoding="utf-8").lower()
    constraints = (APP_DIR / "constraints.txt").read_text(encoding="utf-8")
    pinned = [
        line.split("==")[0].strip().lower()
        for line in constraints.splitlines()
        if "==" in line and not line.strip().startswith("#")
    ]
    assert pinned, "constraints.txt parsed to nothing"
    missing = [p for p in pinned if p.replace("_", "-") not in notice]
    assert not missing, f"pinned but absent from NOTICE table: {missing}"


def test_gui_about_box_surfaces_the_trademark_disclaimer():
    """A GUI-only user never opens NOTICE, so the About box must carry it.

    The nominative-fair-use position for the 10x and IDT marks, and the "no
    patent rights granted" position, both depend on the disclaimer reaching the
    people actually using the tool.
    """
    pytest.importorskip("PySide6")
    from unittest.mock import patch

    from PySide6.QtWidgets import QApplication, QMessageBox

    from flexprobe_gui import tools_dialog

    QApplication.instance() or QApplication([])
    with patch.object(QMessageBox, "about") as mock_about:
        tools_dialog.show_about(None)
    text = mock_about.call_args[0][2]

    assert __version__ in text
    for required in (
        "10x Genomics",
        "Integrated DNA Technologies",
        "oPools",
        "No patent rights are granted",
        "NOT a medical device",
        "GPL-3.0-or-later",
    ):
        assert required in text, f"About box does not mention {required!r}"
