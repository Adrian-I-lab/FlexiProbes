<!--
Thanks for contributing. Keep this short. The checklist matters more than
the prose.
-->

## What and why

<!-- What does this change, and what problem does it solve? Link the issue. -->

Closes #

## Does this change probe design output?

- [ ] **No** — refactor, docs, tests, packaging, or a fix with no effect on
      which probes are produced from an unchanged input.
- [ ] **Yes** — an unchanged input will now produce a different probe set.

If **yes**, all of these are required:

- [ ] An opt-out flag restores the previous behaviour exactly
- [ ] A release-notes paragraph in this pull request describing the change,
      its rationale, and the opt-out
- [ ] The rationale cites its source (10x tech note + revision + page, a
      published panel, or an explicit statement that it is a heuristic)
- [ ] Defaults in `app/src/flexprobe/constants.py` carry the citation inline

## Checks

- [ ] `QT_QPA_PLATFORM=offscreen python -m pytest -q` passes (run from `app/`)
- [ ] `python -m mypy` clean
- [ ] `python -m pyflakes src/flexprobe src/flexprobe_gui` clean
- [ ] Added a test that would have failed before this change

## Housekeeping

- [ ] New files under `app/src/` have the GPL header
      (`python packaging/add_license_headers.py`)
- [ ] If dependencies changed: `python packaging/vendor_licenses.py` re-run so
      the vendored licence texts match `app/constraints.txt`
- [ ] User-visible behaviour is documented in `README.md` / `app/README.md`

## Licensing

- [ ] I agree my contribution is licensed **GPL-3.0-or-later**, and it does not
      include code from an incompatibly licensed project
