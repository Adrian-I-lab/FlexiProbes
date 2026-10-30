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

"""Structured logging setup for flexprobe."""
import logging
import sys
import json
from typing import Any, Optional


class StructuredFormatter(logging.Formatter):
    """Emit log records as JSON lines for machine consumption."""

    def format(self, record: logging.LogRecord) -> str:
        log_entry = {
            "ts": self.formatTime(record, self.datefmt),
            "level": record.levelname,
            "logger": record.name,
            "msg": record.getMessage(),
        }
        if record.exc_info and record.exc_info[1]:
            log_entry["exception"] = self.formatException(record.exc_info)
        extra = getattr(record, "extra_data", None)
        if extra:
            log_entry["data"] = extra
        return json.dumps(log_entry)


def setup_logging(
    level: str = "INFO",
    structured: bool = False,
    stream: Optional[Any] = None,
) -> logging.Logger:
    """Configure and return the root flexprobe logger.

    Parameters
    ----------
    level : str
        Logging level name (DEBUG, INFO, WARNING, ERROR, CRITICAL).
    structured : bool
        If True, emit JSON-lines instead of human-readable messages.
    stream : file-like, optional
        Output stream. Defaults to sys.stderr.
    """
    logger = logging.getLogger("flexprobe")
    logger.setLevel(getattr(logging, level.upper(), logging.INFO))
    # Remove existing handlers to avoid duplicates on re-init
    logger.handlers.clear()

    handler = logging.StreamHandler(stream or sys.stderr)
    if structured:
        handler.setFormatter(StructuredFormatter())
    else:
        handler.setFormatter(
            logging.Formatter(
                "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                datefmt="%Y-%m-%dT%H:%M:%S",
            )
        )
    logger.addHandler(handler)
    return logger
