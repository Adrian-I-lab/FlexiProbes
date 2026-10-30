# flexiprobes — probe design for 10x Genomics Flex v1 / Flex v2 / Visium.
# Copyright (C) 2026 Adrian Ilich and contributors
#
# This file is part of flexiprobes.
#
# flexiprobes is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# flexiprobes is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Run manifest writer (JSON provenance record)."""
from __future__ import annotations
import json
import os
import tempfile
from ..models import RunManifest
from ..errors import WriterError
import flexprobe


def write_manifest(
    manifest: RunManifest,
    output_path: str,
    overwrite: bool = False,
) -> str:
    """Write JSON manifest with provenance information.

    Uses atomic temp-then-rename with overwrite guard.
    """
    if os.path.exists(output_path) and not overwrite:
        raise WriterError(f"Output exists and overwrite=False: {output_path}")

    dir_name = os.path.dirname(output_path) or '.'
    os.makedirs(dir_name, exist_ok=True)

    data = {
        'app_version': manifest.app_version or flexprobe.__version__,
        'code_version': manifest.code_version,
        'timestamp': manifest.timestamp,
        'inputs': manifest.inputs,
        'params': manifest.params,
        'tool_versions': manifest.tool_versions,
        'elapsed_seconds': manifest.elapsed_seconds,
        'db_sha256s': manifest.db_sha256s,
        # Beside db_sha256s deliberately: which reference, and what was in it.
        'reference_composition': manifest.reference_composition,
        # Which aligner built the consensus, on a multi-strain run. ALWAYS
        # emitted, as `{}` on a single-species run, rather than omitted when
        # empty: an absent key cannot be told apart from an engine too old to
        # record it, and that ambiguity is the exact failure this field exists
        # to end. `params['multi_strain']` says which case an empty dict is.
        'consensus': manifest.consensus,
        # Which aligner built the consensus, on a multi-strain run. ALWAYS
        # emitted, as `{}` on a single-species run, rather than omitted when
        # empty: an absent key cannot be told apart from an engine too old to
        # record it, and that ambiguity is the exact failure this field exists
        # to end. `params['multi_strain']` says which case an empty dict is.
        'seed': manifest.seed,
        # Keys the loaded config did not set, filled from this build's
        # defaults. ALWAYS emitted, for the same reason as `consensus`: null
        # means the config was built in memory (no file, nothing to lack), []
        # means a file that set every key. An absent key would only say "an
        # engine too old to record it".
        'config_filled_from_defaults': manifest.config_filled_from_defaults,
    }

    # ⚠️ THIS DICT IS AN EXPLICIT ALLOW-LIST, NOT asdict(manifest). A field
    # added to RunManifest and NOT added here is silently dropped, and the test
    # suite stays green because every golden digest still matches. A field can
    # therefore be populated in memory and absent from manifest.json. The same
    # fixed-key pattern applies to ``_StageTimer._write``. If you add a field
    # above, add it here, and CHECK THE WRITTEN FILE rather than trusting a
    # green suite.

    fd, tmp_path = tempfile.mkstemp(dir=dir_name, suffix='.tmp')
    try:
        with os.fdopen(fd, 'w') as fh:
            json.dump(data, fh, indent=2, sort_keys=True)
            fh.write('\n')
        os.replace(tmp_path, output_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise
    return output_path


def record_stage_elapsed(
    manifest_path: str,
    stage: str,
    seconds: float,
) -> bool:
    """Fold a POST-pipeline stage's wall clock into an existing manifest.

    ``run_pipeline`` writes ``elapsed_seconds`` when it finishes, but the HPC
    bundle then runs ``select-coverage`` as a *separate* invocation. On the July
    2026 mouse run that step took 3 h 05 m of a 44 h 14 m job (~7%), which the
    manifest could not see — so a Methods table built from the manifest alone
    understated the run. This lets the later stage add its own share.

    Both ``stage`` and ``total`` accumulate, so ``total`` remains >= the sum of
    the individual stages and stays comparable to the scheduler's walltime.

    Returns False (never raises) if the manifest is missing or unreadable: the
    selection step is explicitly non-fatal, so a provenance update must not be
    the thing that fails it.
    """
    try:
        with open(manifest_path) as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return False

    elapsed = data.setdefault('elapsed_seconds', {})
    if not isinstance(elapsed, dict):       # hand-edited / foreign manifest
        return False
    elapsed[stage] = round(elapsed.get(stage, 0.0) + seconds, 3)
    elapsed['total'] = round(elapsed.get('total', 0.0) + seconds, 3)

    dir_name = os.path.dirname(manifest_path) or '.'
    fd, tmp_path = tempfile.mkstemp(dir=dir_name, suffix='.tmp')
    try:
        with os.fdopen(fd, 'w') as fh:
            json.dump(data, fh, indent=2, sort_keys=True)
            fh.write('\n')
        os.replace(tmp_path, manifest_path)
    except Exception:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        return False
    return True
