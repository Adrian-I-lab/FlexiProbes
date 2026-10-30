# Third-party licenses

flexiprobes itself is licensed **GPL-3.0-or-later** (see `../LICENSE`). It depends
on, and its desktop binaries bundle, the third-party components below. Each keeps
its own license; this folder holds the **full verbatim texts**.

The texts are vendored into this folder as ordinary files by
`../packaging/vendor_licenses.py`, which copies them out of the installed
distributions so they match the versions pinned in `../constraints.txt`. The
PyInstaller spec ships this whole folder inside the binary under
`licenses/THIRD-PARTY-LICENSES/`. Re-run the script after any dependency bump;
`vendor_licenses.py --check` fails if the tree has drifted.

## Bundled in the desktop binaries

These are redistributed as compiled code inside the `.app` / onedir bundle, so
their license texts must travel with the binary.

| Component | Version | License | Full text |
|---|---|---|---|
| **flexiprobes** (this app) | the `version` in `../pyproject.toml` | GPL-3.0-or-later | `../LICENSE` |
| **NCBI BLAST+** (`blastn`, `makeblastdb`, `blastdbcmd`) | 2.17.0 | Public domain (US Government work) | `PUBLIC-DOMAIN-NCBI-BLAST.txt` (also `../packaging/vendor/blast/LICENSE.ncbi-blast`) |
| **PySide6 / Qt** (incl. `shiboken6`) | 6.11.1 | Tri-licensed LGPL-3.0 / GPL-2.0 / GPL-3.0 — **used here under the GPL-3.0 option** | `../LICENSE` (GPL-3.0, the option in force); `qt-pyside6-LGPL-3.0.txt` (LGPL-3.0, for downstream re-users who prefer that option) |
| **primer3-py** (incl. vendored libprimer3) | 2.3.1 | GPL-2.0-or-later | `primer3-py-GPL-2.0.txt` |
| **Biopython** (optional alignment) | 1.87 | Biopython License Agreement (BSD-style) | `biopython-BIOPYTHON-LICENSE.txt` |
| **NumPy** | 2.4.6 | BSD-3-Clause | `numpy-BSD-3-Clause.txt` |
| **MAFFT** (multiple alignment for homology runs) | 7.526 | BSD (the optional `extensions` directory carries its own terms, see the file) | `mafft-BSD.txt` (also `../packaging/vendor/mafft/LICENSE.mafft`) |
| **MSYS2 tools inside Windows MAFFT** (Windows build only) | see `msys2-README.txt` | GPL-3.0, LGPL-3.0, LGPL-2.1, GPL-2.0, GCC Runtime Library Exception, BSD, MIT, zlib (per component) | `msys2-README.txt` and the `msys2-*.txt` texts it lists |

## Runtime dependencies, not bundled

Declared in `../pyproject.toml` and pinned in `../constraints.txt`, but **not**
present in the desktop bundle — pip installs them into the user's own
environment when flexiprobes is installed from source. Their texts are vendored
anyway so a source distribution is self-describing.

| Component | Version | License | Full text |
|---|---|---|---|
| **openpyxl** (optional `.xlsx` export) | 3.1.5 | MIT | `openpyxl-MIT.txt` |
| **et-xmlfile** (openpyxl transitive) | 2.0.0 | MIT | `et-xmlfile-MIT.txt` |

Why these are not in the bundle: `openpyxl` is imported lazily inside
`flexprobe/writers/results_table.py` and degrades gracefully when absent, so
PyInstaller's static analysis does not see it — the frozen app writes CSV but
not XLSX. `et-xmlfile` is openpyxl's own dependency.

## Notes

- **Qt is used under its GPL-3.0 option.** The combined work is
  GPL-3.0-or-later, and the Qt Company offers GPL-3.0 as one of the three
  options for Qt and PySide6. Under that option the relinking /
  replaceable-Qt requirement that LGPL-3.0 §4 would impose **does not apply**;
  the obligation that does apply is GPLv3's — convey the Corresponding Source,
  which `../SOURCE-OFFER.txt` and the public repository satisfy.

  The binaries are nevertheless built as PyInstaller **onedir**, not onefile.
  That is a deliberate belt-and-braces choice rather than an LGPL obligation:
  it leaves the Qt shared libraries as separate, replaceable files on disk,
  which keeps the LGPL-3.0 option genuinely exercisable by anyone downstream
  who would rather take Qt that way, and it keeps the bundled `licenses/`
  folder visible to the end user instead of sealed inside a self-extracting
  blob. (It also avoids onefile's per-launch extraction overhead and macOS
  code-signing friction.)
- **NCBI BLAST+ is public domain**, so bundling and redistributing it carries no
  license obligation. Its usage-reporting network call is disabled by the app
  (`BLAST_USAGE_REPORT=false`); see `../NOTICE`.
- **primer3-py is GPL-2.0-or-later**, which is why flexiprobes as a whole is
  GPL. Taken under the "or later" grant it is usable under GPL-3.0, which is how
  it combines with Qt-under-GPL-3.0 and with this project's own license. See the
  "WHY GPL: COPYLEFT CHAIN" section of `../NOTICE`.
- Trademark and patent notices for 10x Genomics and Integrated DNA Technologies
  are in `../NOTICE`; neither company is affiliated with, endorses, or supports
  this tool.
