#!/usr/bin/env python3
"""Assert what PyInstaller ACTUALLY collected into a built bundle.

Why this exists
---------------
``packaging/build.py`` already refuses a bundle containing a foreign Qt5,
because that failure is silent: the build reports success and the app dies at
launch. Every other way the bundle can come out wrong is silent in exactly the
same way, and several of them have already happened here:

* ``primer3`` ships TWO C extensions. An earlier build collected
  ``thermoanalysis`` but not ``p3helpers``, so ``import primer3`` failed inside
  the frozen app and **every probe's Tm and hairpin came back n/a** — a
  scientifically wrong output from a bundle that built cleanly (see the comment
  in ``flexprobe.spec``).
* The vendored NCBI BLAST+ binaries are added by the spec only ``if
  os.path.exists(...)``; a missing one prints ``[spec] WARNING`` and the build
  continues. A warning in thousands of lines of PyInstaller output is not a
  gate, and the result is a bundle whose BLAST fallback engine simply is not
  there.
* The licence texts under ``licenses/`` are what makes redistributing the
  binary lawful at all (GPLv2 §1 / GPLv3 §4 for primer3-py and Qt, the
  BSD/MIT notice requirements for NumPy and Biopython). A build that loses them
  is a licence violation that looks like a green build.
* PyInstaller collects whatever is importable in the build interpreter. A dirty
  environment can silently sweep in packages such as ``PIL``,
  ``charset_normalizer`` or a compiled ``__mypyc`` extension, none of them
  dependencies and none with licence text. ``vendor_licenses.py --check`` does
  NOT catch this, because it checks constraints.txt and not what PyInstaller
  actually collected. This script is the missing half.

So: read-only, no build side effects, run it on the artefact in CI right after
``build.py`` and before anything is uploaded or released.

Usage
-----
    python packaging/verify_bundle.py ../dist/flexiprobes.app     # macOS
    python packaging/verify_bundle.py ../dist/flexiprobes         # Linux/Windows

Layout-agnostic on purpose. PyInstaller has moved things between releases
(``_internal/`` appeared in 6.x) and a macOS ``.app`` splits datas into
``Contents/Resources`` and binaries into ``Contents/Frameworks``. Every check
below therefore matches on a path SUFFIX or a filename anywhere in the tree,
so it keeps working across PyInstaller versions and all three platforms rather
than encoding one layout that will quietly stop matching.
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
from pathlib import Path

APP_DIR = Path(__file__).resolve().parent.parent
THIRD_PARTY = APP_DIR / "THIRD-PARTY-LICENSES"

# Packages that must NEVER appear: none is a dependency, and each one showed up
# in a real artefact built from a dirty interpreter.
FOREIGN_DIRS = {"PIL", "charset_normalizer", "mypy", "mypyc", "tkinter"}
FOREIGN_FILE_RE = re.compile(r"(?i)(^__mypyc|^_?mypyc)")
# The anaconda-Qt5 contamination build.py already guards at build time; re-check
# it here so the assertion also covers an artefact someone hands you.
QT5_RE = re.compile(r"(?i)^(lib)?qt5")


def _source_version() -> str:
    init = APP_DIR / "src" / "flexprobe" / "__init__.py"
    for line in init.read_text(encoding="utf-8").splitlines():
        if line.startswith("__version__"):
            return line.split("=", 1)[1].strip().strip('"').strip("'")
    raise SystemExit(f"__version__ not found in {init}")


class Report:
    def __init__(self) -> None:
        self.failures: list[str] = []

    def ok(self, what: str) -> None:
        print(f"  OK    {what}")

    def fail(self, what: str, detail: str) -> None:
        print(f"  FAIL  {what}\n        {detail}")
        self.failures.append(what)

    def require(self, cond: bool, what: str, detail: str) -> bool:
        if cond:
            self.ok(what)
        else:
            self.fail(what, detail)
        return bool(cond)


def collect(root: Path) -> tuple[list[str], list[str]]:
    """Return (relative posix paths of files, of directories).

    ``followlinks=False`` matters: a macOS ``.app`` symlinks
    ``Contents/Resources/blast`` at ``Contents/Frameworks/blast`` and following
    those both double-counts and can loop.
    """
    files: list[str] = []
    dirs: list[str] = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        rel = Path(dirpath).relative_to(root)
        for d in dirnames:
            dirs.append((rel / d).as_posix())
        for f in filenames:
            files.append((rel / f).as_posix())
    return files, dirs


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("bundle", type=Path, help="the built .app / onedir folder")
    ap.add_argument(
        "--expect-version",
        default=None,
        help="version the bundled GPLv3 §6 source offer must name "
             "(default: flexprobe.__version__ in this source tree)",
    )
    args = ap.parse_args(argv)

    root: Path = args.bundle
    if not root.is_dir():
        print(f"FAIL  bundle not found (or not a directory): {root}")
        return 1
    want_version = args.expect_version or _source_version()

    files, dirs = collect(root)
    names = {Path(p).name for p in files}
    print(f"Verifying {root}  ({len(files)} files, {len(dirs)} directories)")
    print(f"Expecting the bundled source offer to name version {want_version}")

    r = Report()

    # --- Legal: the texts that make redistribution lawful --------------------
    for legal in ("LICENSE", "NOTICE", "SOURCE-OFFER.txt"):
        r.require(
            any(p.endswith(f"licenses/{legal}") for p in files),
            f"licenses/{legal} is in the bundle",
            "GPLv3 §4/§6 require the licence and the written source offer to "
            "travel with the object code; the spec adds these from app/.",
        )

    expected_texts = sorted(p.name for p in THIRD_PARTY.iterdir() if p.is_file())
    if not expected_texts:
        r.fail(
            "app/THIRD-PARTY-LICENSES/ has files to check",
            f"{THIRD_PARTY} is empty or missing in the source tree.",
        )
    missing_texts = [
        t for t in expected_texts
        if not any(p.endswith(f"licenses/THIRD-PARTY-LICENSES/{t}") for p in files)
    ]
    r.require(
        not missing_texts,
        f"all {len(expected_texts)} vendored third-party licence texts shipped",
        f"missing from the bundle: {missing_texts}",
    )

    # The GPLv3 §6 written offer must name the version of THIS object code.
    # The spec substitutes @FLEXIPROBES_VERSION@ at build time; a bundle that
    # still contains the placeholder, or names another release, is the exact
    # defect that shipped in 0.4.0 (an offer declaring 0.3.6).
    offers = [p for p in files if p.endswith("licenses/SOURCE-OFFER.txt")]
    if offers:
        text = (root / offers[0]).read_text(encoding="utf-8", errors="replace")
        r.require(
            "@FLEXIPROBES_VERSION@" not in text,
            "bundled source offer has its version substituted",
            "the @FLEXIPROBES_VERSION@ placeholder is still literal — the "
            "GPLv3 §6 offer does not name any version at all.",
        )
        r.require(
            want_version in text,
            f"bundled source offer names {want_version}",
            f"the offer does not mention {want_version}; it must cover the "
            f"Corresponding Source for THIS object code.",
        )

    # --- Vendored NCBI BLAST+ (the fallback off-target engine) ---------------
    # flexprobe.blast._bundled_blast_dir() looks for these at "blast/bin".
    for tool in ("blastn", "makeblastdb"):
        r.require(
            any(
                re.search(rf"(^|/)blast/bin/{tool}(\.exe)?$", p) for p in files
            ),
            f"vendored BLAST+ {tool} shipped under blast/bin/",
            "the spec only prints '[spec] WARNING' when a vendored BLAST+ "
            "binary is missing and builds anyway, so the packaged app would "
            "ship with no BLAST fallback at all.",
        )

    # --- Vendored MAFFT (the multi-strain consensus aligner) -----------------
    # flexprobe.consensus._bundled_mafft() probes for a driver named
    # "mafft.bat" at the top of a relocatable folder, and ONLY that: the
    # docstring records that requiring "mafftdir" as well would have silently
    # failed on Windows, whose portable package puts its helpers under "usr/"
    # inside a partial MSYS tree instead. So the driver check below matches
    # exactly what the runtime matches, and the helper check accepts either
    # layout rather than encoding the POSIX one.
    #
    # Without MAFFT the pipeline falls back to a progressive merge that
    # consensus.py itself documents as NOT a true multiple alignment, and the
    # spec only gates on the driver existing in the vendor tree, so a bundle
    # built from an incomplete vendor folder would ship and align differently
    # from every other build of the same input.
    mafft_driver = [p for p in files if re.search(r"(^|/)mafft/mafft\.bat$", p)]
    r.require(
        bool(mafft_driver),
        "vendored MAFFT driver shipped as mafft/mafft.bat",
        "consensus._bundled_mafft() probes this exact path and returns None "
        "without it, so multi-strain designs would silently downgrade to the "
        "built-in progressive merge.",
    )
    r.require(
        any(re.search(r"(^|/)mafft/(mafftdir|usr)(/|$)", p) for p in files + dirs),
        "vendored MAFFT helper tree shipped (mafftdir/ or usr/)",
        "the driver locates its helpers relative to itself; a driver with no "
        "helper tree beside it exits non-zero on the first real alignment.",
    )

    # Functional check, not just a content check. The driver is shipped as
    # DATA, so PyInstaller strips its executable bit and the runtime invokes it
    # through an explicit `bash` — this mirrors that exactly, including popping
    # MAFFT_BINARIES so an inherited value cannot point the driver at a
    # DIFFERENT MAFFT install. Version goes to STDERR and the exit status is
    # not dependable across versions (see packaging/container/Dockerfile), so
    # the streams are merged and the returncode is deliberately ignored.
    if mafft_driver and os.name == "posix":
        driver = root / mafft_driver[0]
        env = dict(os.environ)
        env.pop("MAFFT_BINARIES", None)
        try:
            out = subprocess.run(
                ["bash", str(driver), "--version"],
                capture_output=True, text=True, timeout=60, env=env,
            )
            blob = (out.stdout or "") + (out.stderr or "")
        except (OSError, subprocess.SubprocessError) as exc:
            blob = ""
            print(f"  NOTE  could not execute the bundled MAFFT driver: {exc}")
        found = re.search(r"v\d+\.\d+", blob)
        r.require(
            bool(found),
            f"bundled MAFFT driver runs and reports a version "
            f"({found.group(0) if found else 'none'})",
            "the driver is present but did not report a version when run the "
            "way the app runs it. Its helper tree is incomplete or was built "
            "for another platform, and every multi-strain design in this "
            "bundle would fall back to the progressive merge.",
        )

    # --- primer3: BOTH C extensions, plus the thermodynamic parameters -------
    for ext in ("p3helpers", "thermoanalysis"):
        r.require(
            any(n.startswith(ext + ".") for n in names),
            f"primer3 C extension {ext} shipped",
            "collect_submodules() only takes the pure-Python modules. Missing "
            "p3helpers is what made `import primer3` fail inside the frozen "
            "app and every Tm/hairpin come back n/a.",
        )
    r.require(
        any("/primer3_config/" in p or p.endswith("/primer3_config") for p in files + dirs),
        "primer3_config/ thermodynamic parameters shipped",
        "primer3 cannot compute Tm without its parameter files.",
    )

    # --- Qt / PySide6: the product is a GUI ----------------------------------
    for lib in ("QtCore", "QtWidgets"):
        r.require(
            any(re.search(rf"(?i)(^|[^a-z]){lib}", n) for n in names),
            f"Qt {lib} shipped",
            "PySide6's Qt libraries are missing; the GUI cannot start.",
        )

    # --- our own data file ---------------------------------------------------
    r.require(
        any(p.endswith("flexprobe_gui/assets/icon.png") for p in files),
        "flexprobe_gui/assets/icon.png shipped",
        "the spec adds this so QIcon can find the window icon at runtime; if "
        "it is absent, the datas list did not make it into the bundle.",
    )

    # --- Negative assertions -------------------------------------------------
    qt5 = sorted(n for n in names if QT5_RE.match(n))
    r.require(
        not qt5,
        "no foreign Qt5 libraries",
        f"{len(qt5)} Qt5 librar(y/ies) collected alongside PySide6's Qt6 "
        f"({qt5[:4]}). This bundle WILL abort at launch. You are building "
        f"with the wrong interpreter — almost certainly anaconda.",
    )

    strays = sorted(
        {Path(p).parts[i] for p in files for i in range(len(Path(p).parts) - 1)}
        & FOREIGN_DIRS
    )
    stray_files = sorted(n for n in names if FOREIGN_FILE_RE.match(n))
    r.require(
        not strays and not stray_files,
        "no foreign packages swept in from the build interpreter",
        f"found {strays + stray_files}. None is a dependency and none carries "
        f"licence text. Build from a clean, constraints-pinned venv (see "
        f"app/README.md, Package a desktop binary).",
    )

    print()
    if r.failures:
        print(f"FAILED {len(r.failures)} check(s): {r.failures}")
        return 1
    print("All bundle content checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
