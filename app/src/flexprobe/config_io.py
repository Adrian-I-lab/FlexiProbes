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

"""Serialize a :class:`~flexprobe.models.RunConfig` to/from JSON.

This is the backbone of the HPC-export feature: a run configured in the GUI is
serialized verbatim, carried to a cluster, and reconstructed by the headless CLI
(``flexiprobes design --config config.json``). Round-tripping the whole dataclass
guarantees the cluster run is identical to what the GUI would have launched,
instead of re-deriving every field through CLI flags (which drift out of sync).

``RunConfig`` is a tree of plain dataclasses, so ``dataclasses.asdict`` handles
serialization. The inverse needs to rebuild the nested dataclasses by hand;
:func:`config_from_dict` does so defensively — unknown keys are ignored and
missing keys fall back to dataclass defaults — so a config written by an older or
newer build still loads (mirroring ``RunRecord.from_json``).
"""
from __future__ import annotations

import json
import os
from dataclasses import asdict, fields
from typing import Any, Dict, List, Optional, Type, TypeVar

import logging

from .constants import SCREEN_ENGINE_DEFAULT
from .models import (
    ScreenResources,
    ExportSettings,
    FilterThresholds,
    ModalityDefinition,
    RunConfig,
    SpecificityThresholds,
)
from .screen import ENGINE_BLAST

logger = logging.getLogger("flexprobe.config_io")

_T = TypeVar("_T")


def config_to_dict(config: RunConfig) -> Dict[str, Any]:
    """Return a JSON-serializable dict for ``config`` (nested dataclasses flattened)."""
    return asdict(config)


def _build(cls: Type[_T], data: Any) -> _T:
    """Construct dataclass ``cls`` from ``data``, ignoring unknown keys.

    Returns a default-constructed instance when ``data`` is not a mapping (e.g. a
    missing nested block), so a partial config still yields a valid object.
    """
    if not isinstance(data, dict):
        return cls()
    known = {f.name for f in fields(cls)}  # type: ignore[arg-type]
    return cls(**{k: v for k, v in data.items() if k in known})


def config_from_dict(data: Dict[str, Any]) -> RunConfig:
    """Reconstruct a :class:`RunConfig` from :func:`config_to_dict` output.

    Tolerant of missing/extra keys for forward/backward compatibility.
    """
    if not isinstance(data, dict):
        raise ValueError("config payload must be a JSON object")

    nested = {
        "modality": (ModalityDefinition, data.get("modality")),
        "filters": (FilterThresholds, data.get("filters")),
        # Older configs call this block "blast"; it is "screen" because the
        # default engine is the k-mer screen. A legacy config is still read:
        # the key moved, the meaning did not.
        "screen": (ScreenResources, data.get("screen", data.get("blast"))),
        "export": (ExportSettings, data.get("export")),
        "specificity": (SpecificityThresholds, data.get("specificity")),
    }

    known = {f.name for f in fields(RunConfig)}
    scalar = {
        k: v for k, v in data.items()
        if k in known and k not in nested
    }
    cfg = RunConfig(
        **scalar,
        **{name: _build(cls, raw) for name, (cls, raw) in nested.items()},
    )
    _pin_legacy_screen_engine(data, cfg)

    filled = sorted(
        [f.name for f in fields(RunConfig)
         if f.name not in nested and f.name not in data]
        + [f"{name}.{f.name}" for name, (cls, raw) in nested.items()
           for f in fields(cls)
           if not (isinstance(raw, dict) and f.name in raw)]
    )
    # Stored as an instance attribute, NOT a dataclass field: it is a fact
    # about where this config came from, not a setting, so asdict (and so
    # save_config and run_config.json) never sees it and equality ignores it.
    cfg.__dict__[_FILLED_ATTR] = filled
    if filled:
        logger.warning(
            "This config does not set %d key(s), so they take TODAY's "
            "defaults, which may not be what the original run used: %s. "
            "They are recorded in the run manifest as "
            "config_filled_from_defaults.", len(filled), ", ".join(filled))
    return cfg


_FILLED_ATTR = "_config_filled_from_defaults"


def filled_from_defaults(cfg: RunConfig) -> Optional[List[str]]:
    """Dotted keys :func:`config_from_dict` filled from defaults.

    ``[]`` means the payload set every key. ``None`` means ``cfg`` was not
    built from a payload at all (constructed in memory by the GUI or the CLI
    flags), so there is nothing it could have lacked.
    """
    value = cfg.__dict__.get(_FILLED_ATTR)
    return list(value) if value is not None else None


def _pin_legacy_screen_engine(data: Dict[str, Any], cfg: RunConfig) -> None:
    """A config written before ``blast.engine`` existed must keep BLAST.

    A config file is a reproducibility artefact: it records what a run DID. The
    field arrived together with the k-mer engine as the default, so a payload
    without it was authored under BLAST, and letting it inherit the current
    default would silently re-run it on a different engine, producing different
    probes in a different order with nothing in the file to say so. Such a
    config typically points ``db_paths`` at a FASTA, which the k-mer engine
    accepts happily, so there would be no error to notice.

    So the absent key is pinned to ``blast`` and the substitution is logged with
    the flag that opts in. Configs written by this version always carry the
    field (it is a dataclass member), so this only ever fires for legacy
    payloads.
    """
    blast_raw = data.get("screen", data.get("blast"))
    if isinstance(blast_raw, dict) and "engine" in blast_raw:
        return                       # explicit — honour whatever it says
    if cfg.screen.skip_screen or not cfg.screen.db_paths:
        return                       # no screen runs; the engine is moot
    cfg.screen.engine = ENGINE_BLAST
    logger.warning(
        "This config predates the off-target screen engine setting, so it was "
        "written for BLAST; pinning engine='%s' to reproduce it faithfully "
        "rather than silently adopting the current default ('%s'), which would "
        "change both which probes pass and their ranking. Re-export the config, "
        "or pass --screen-engine %s, to run it on the new engine deliberately.",
        ENGINE_BLAST, SCREEN_ENGINE_DEFAULT, SCREEN_ENGINE_DEFAULT,
    )


#: The RunConfig fields that hold filesystem paths.
PATH_FIELDS = ("fasta_path", "output_dir", "screen.db_paths",
               "screen.cached_hits_path")

#: Key a portable record carries: the dotted fields (list items as
#: ``field[i]``) its writer turned from absolute into relative-to-this-file.
#: Only those are resolved on load, so a hand-written relative path keeps its
#: usual current-directory meaning. Not a RunConfig field, so the loader's
#: unknown-key rule ignores it and asdict never writes it.
RELATIVE_MARKER = "_paths_relative_to_this_file"


def portable_path(path: Any, base: str) -> Optional[str]:
    """``path`` relative to directory ``base``, or None to leave it as is.

    Only absolute paths are rewritten. A path on another Windows drive has no
    relative form and stays whole. Relative to the record's own directory
    rather than the repository, because
    an engine record has to reload to the same files wherever the run tree is
    moved, and an installed package has no repository.
    """
    if not isinstance(path, str) or not path or not os.path.isabs(path):
        return None
    try:
        return os.path.relpath(path, base)
    except ValueError:
        return None


def _portable_fields(d: Dict[str, Any], base: str) -> List[str]:
    """Rewrite PATH_FIELDS in ``d`` in place; return the dotted names changed."""
    converted: List[str] = []
    for dotted in PATH_FIELDS:
        *parents, leaf = dotted.split(".")
        holder: Any = d
        for part in parents:
            holder = holder.get(part) if isinstance(holder, dict) else None
        if not isinstance(holder, dict) or leaf not in holder:
            continue
        value = holder[leaf]
        if isinstance(value, list):
            for i, item in enumerate(value):
                rel = portable_path(item, base)
                if rel is not None:
                    value[i] = rel
                    converted.append(f"{dotted}[{i}]")
        else:
            rel = portable_path(value, base)
            if rel is not None:
                holder[leaf] = rel
                converted.append(dotted)
    return converted


def portable_config_dict(config: RunConfig, base: str) -> Dict[str, Any]:
    """``config_to_dict`` with absolute paths made relative to ``base``.

    Used for records written into a run's output directory (run_config.json
    and the manifest's ``params``), which are committed and must not publish a
    home directory or a cluster layout.
    """
    d = config_to_dict(config)
    converted = _portable_fields(d, base)
    if converted:
        d[RELATIVE_MARKER] = sorted(converted)
    return d


def _resolve_portable_fields(d: Dict[str, Any], base: str) -> None:
    """Inverse of :func:`_portable_fields`, for the fields the marker names."""
    for entry in d.pop(RELATIVE_MARKER, None) or []:
        dotted, _, index = entry.partition("[")
        *parents, leaf = dotted.split(".")
        holder: Any = d
        for part in parents:
            holder = holder.get(part) if isinstance(holder, dict) else None
        if not isinstance(holder, dict) or leaf not in holder:
            continue
        if index:
            i = int(index.rstrip("]"))
            items = holder[leaf]
            if isinstance(items, list) and i < len(items) and isinstance(items[i], str):
                items[i] = os.path.normpath(os.path.join(base, items[i]))
        elif isinstance(holder[leaf], str):
            holder[leaf] = os.path.normpath(os.path.join(base, holder[leaf]))


def save_config(config: RunConfig, path: str, *, portable: bool = False,
                extra: Optional[Dict[str, Any]] = None) -> str:
    """Write ``config`` to ``path`` as indented JSON. Returns ``path``.

    ``extra`` adds top-level keys that are not RunConfig fields, such as the
    run's ``effective_ranking_weights``. :func:`load_config` ignores
    unknown keys, so they are a record, never an input.

    ``portable=True`` writes absolute paths relative to the file's directory
    and records which (see RELATIVE_MARKER); the pipeline uses it for
    run_config.json. The default writes paths as given: the HPC bundle's
    config.json holds CLUSTER paths, written on another machine, which must
    not be rewritten against a local directory.
    """
    data = (portable_config_dict(config, os.path.dirname(os.path.abspath(path)))
            if portable else config_to_dict(config))
    if extra:
        data.update(extra)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, sort_keys=True)
    return path


def load_config(path: str) -> RunConfig:
    """Load a :class:`RunConfig` from a JSON file written by :func:`save_config`.

    Fields a portable writer made relative are resolved against the file's
    directory, so a run tree can be moved and its config still names its
    inputs.
    """
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    if isinstance(data, dict) and RELATIVE_MARKER in data:
        _resolve_portable_fields(data, os.path.dirname(os.path.abspath(path)))
    return config_from_dict(data)
