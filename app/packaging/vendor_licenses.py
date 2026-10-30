#!/usr/bin/env python3
"""Vendor third-party licence texts into ``app/THIRD-PARTY-LICENSES/``.

Why this exists
---------------
The desktop binaries redistribute compiled third-party code. GPLv2 §1 and
GPLv3 §4 require that a recipient of a conveyed copy receives the licence text
with it; BSD-3-Clause and MIT require the copyright notice and permission
notice be reproduced in binary redistributions. Relying on PyInstaller to carry
each wheel's ``*.dist-info`` into the bundle does NOT satisfy that: PyInstaller
collects imported *modules*, not distribution metadata, so in practice only
packages whose own hooks happen to request the metadata survive (in the
0.4.0 build tree, only ``numpy-2.4.6.dist-info/``).

So we vendor the texts as ordinary files in the source tree. The PyInstaller
spec already ships the whole ``THIRD-PARTY-LICENSES/`` directory into the
bundle under ``licenses/``, which makes the obligation robust and auditable.

Usage
-----
Run from the repo root (or anywhere; paths resolve off this file)::

    app/.venv-build/bin/python app/packaging/vendor_licenses.py
    app/.venv-build/bin/python app/packaging/vendor_licenses.py --check

``--check`` re-derives every file and fails if the tree is out of date, so CI
can catch a dependency bump that silently changes a licence.

Texts are taken from the *installed distributions* in the build virtualenv, so
they match the pinned versions in ``app/constraints.txt`` exactly. Two texts
cannot come from a wheel and are handled explicitly:

* **PySide6 / shiboken6** ship no licence file in their PyPI wheels at all.
  Qt is tri-licensed (LGPL-3.0 / GPL-2.0 / GPL-3.0); flexiprobes uses it under
  the **GPL-3.0** option, whose text is ``app/LICENSE``. The LGPL-3.0 text is
  vendored as well because it is the option a downstream re-user is most
  likely to want to exercise, and it is not otherwise present in the tree.
  Fetched once from https://www.gnu.org/licenses/lgpl-3.0.txt and committed;
  ``--check`` verifies its SHA-256 rather than re-fetching.
* **NCBI BLAST+** is a US Government work in the public domain; its notice was
  already vendored by hand as ``PUBLIC-DOMAIN-NCBI-BLAST.txt``.

Why ``--check`` normalises ``Files:`` lines
-------------------------------------------
numpy's ``LICENSE.txt`` is not identical across wheels, and the difference is
the *paths of the bundled shared libraries*, not the licence terms. The same
file reads, on a Linux wheel and a macOS wheel respectively::

    Files: numpy.libs/libscipy_openblas*.so
    Files: numpy/.dylibs/libscipy_openblas*.so

    Files: numpy.libs/libgfortran*.so
    Files: numpy/.dylibs/libgfortran*, numpy/.dylibs/libgcc*

Four lines in a 935-line file, all of them Debian-copyright-style file
manifests naming which bundled binary each following sub-licence covers. Every
word of every licence is byte-identical.

That made ``--check`` unsatisfiable. The vendored copy matches whichever OS
last ran the generator, so a text written on macOS failed CI's ubuntu runner --
and regenerating it on Linux would only have inverted the failure onto every
macOS developer. The ``lockfile`` job was red on every run for that reason
alone, which is worse than useless: a check that can never pass is a check
nobody reads.

So the comparison, and *only* the comparison, blurs ``Files:`` lines. What gets
written to disk stays a verbatim copy -- these are third-party licence
documents and this script does not edit them. A wheel that genuinely changes
its terms, or that bundles a new library under a new licence, still trips the
check, because that shows up as changed licence *text* rather than as a changed
path. numpy is the only vendored text with any ``Files:`` line at all, so
nothing else is affected.
"""
from __future__ import annotations

import argparse
import hashlib
import re
import sys
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent.parent
OUT_DIR = APP_DIR / "THIRD-PARTY-LICENSES"

# (distribution name, path within the installed dist-info, output filename)
# The dist-info path differs between wheels: newer metadata puts licence files
# under ``licenses/``, older ones sit at the dist-info root, and the spelling
# alternates between LICENSE / LICENCE / *.rst.
WHEEL_TEXTS = [
    ("primer3_py", "licenses/LICENSE", "primer3-py-GPL-2.0.txt"),
    ("biopython", "licenses/LICENSE.rst", "biopython-BIOPYTHON-LICENSE.txt"),
    ("numpy", "licenses/LICENSE.txt", "numpy-BSD-3-Clause.txt"),
    ("openpyxl", "LICENCE.rst", "openpyxl-MIT.txt"),
    ("et_xmlfile", "LICENCE.rst", "et-xmlfile-MIT.txt"),
]

# Hand-fetched texts that no wheel provides. name -> (filename, sha256).
FETCHED_TEXTS = {
    "LGPL-3.0": (
        "qt-pyside6-LGPL-3.0.txt",
        "e3a994d82e644b03a792a930f574002658412f62407f5fee083f2555c5f23118",
    ),
}


def _dist_info_dir(dist_name: str) -> Path:
    """Locate the installed .dist-info directory for ``dist_name``."""
    from importlib.metadata import distribution

    dist = distribution(dist_name)
    # ``_path`` is the .dist-info dir for a normal filesystem install.
    path = getattr(dist, "_path", None)
    if path is None:
        raise RuntimeError(f"cannot locate dist-info for {dist_name}")
    return Path(path)


# A Debian-copyright-style file manifest, e.g.
#     Files: numpy.libs/libscipy_openblas*.so
# The globs name platform-specific bundled binaries and differ between a Linux
# wheel and a macOS one; the licence stanza that follows does not. See
# "Why ``--check`` normalises ``Files:`` lines" in the module docstring.
_FILES_MANIFEST = re.compile(r"^Files:.*$", re.MULTILINE)
_FILES_PLACEHOLDER = "Files: <bundled binaries - path varies by platform>"


def _normalise_for_compare(text: str) -> str:
    """``text`` with file manifests blurred, for equality checks only.

    Never used on the way to disk. The vendored files stay verbatim.
    """
    return _FILES_MANIFEST.sub(_FILES_PLACEHOLDER, text)


def _header(dist_name: str, version: str, rel_path: str) -> str:
    return (
        f"Licence text for {dist_name} {version}, redistributed with flexiprobes.\n"
        f"Copied verbatim (below the line) from the installed distribution at\n"
        f"    {dist_name}-{version}.dist-info/{rel_path}\n"
        f"by app/packaging/vendor_licenses.py. Do not edit by hand.\n"
        f"{'=' * 78}\n"
    )


def build() -> dict[str, str]:
    """Return {output filename: full file content} for every wheel-sourced text."""
    out: dict[str, str] = {}
    for dist_name, rel_path, out_name in WHEEL_TEXTS:
        info = _dist_info_dir(dist_name)
        src = info / rel_path
        if not src.exists():
            raise FileNotFoundError(
                f"{dist_name}: expected licence text at {src} — the wheel layout "
                f"changed; update WHEEL_TEXTS in {Path(__file__).name}."
            )
        # ``numpy-2.4.6.dist-info`` -> ``2.4.6``. Strip the suffix BEFORE
        # splitting: ``split("-")[1]`` on the raw name yields "2.4.6.dist",
        # because ".dist-info" contains the separator, and the trailing
        # ``removesuffix`` then has nothing left to match, giving a header
        # that reads "numpy 2.4.6.dist".
        version = info.name.removesuffix(".dist-info").split("-", 1)[1]
        out[out_name] = _header(dist_name, version, rel_path) + src.read_text(
            encoding="utf-8"
        )
    return out


def check_fetched() -> list[str]:
    """Verify the hand-fetched texts are present and unmodified."""
    problems = []
    for label, (fname, want) in FETCHED_TEXTS.items():
        path = OUT_DIR / fname
        if not path.exists():
            problems.append(f"missing {fname} ({label})")
            continue
        got = hashlib.sha256(path.read_bytes()).hexdigest()
        if got != want:
            problems.append(f"{fname}: sha256 {got} != expected {want}")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--check",
        action="store_true",
        help="verify the vendored texts are current; do not write.",
    )
    args = ap.parse_args()

    texts = build()
    problems = check_fetched()

    if args.check:
        for name, content in texts.items():
            path = OUT_DIR / name
            if not path.exists():
                problems.append(f"missing {name}")
            elif _normalise_for_compare(
                path.read_text(encoding="utf-8")
            ) != _normalise_for_compare(content):
                problems.append(f"stale {name}")
        if problems:
            print("Vendored licence texts are out of date:", file=sys.stderr)
            for p in problems:
                print(f"  - {p}", file=sys.stderr)
            print(
                f"\nRe-run: {sys.executable} {Path(__file__)}",
                file=sys.stderr,
            )
            return 1
        print(f"OK — {len(texts) + len(FETCHED_TEXTS)} vendored licence texts current.")
        return 0

    OUT_DIR.mkdir(exist_ok=True)
    for name, content in texts.items():
        (OUT_DIR / name).write_text(content, encoding="utf-8")
        print(f"wrote {name}")
    if problems:
        print("\nWARNING — hand-fetched texts need attention:", file=sys.stderr)
        for p in problems:
            print(f"  - {p}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
