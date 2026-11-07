#!/usr/bin/env python3
"""Deterministic Matplotlib SVG export for published repository artifacts."""

from __future__ import annotations

import sys as _export_sys
from pathlib import Path as _ExportPath
_export_sys.path.insert(0, str(_ExportPath(__file__).resolve().parents[1] / "plots"))
from vector_export import save_figure


from pathlib import Path
from typing import Any


SVG_HASH_SALT = "ocean-imu-published-svg-v1"
SVG_METADATA = {
    "Date": None,
    "Creator": "ocean-imu deterministic Matplotlib export",
}


def configure_svg(matplotlib_module: Any) -> None:
    """Fix Matplotlib's SVG id salt before any SVG figures are rendered."""

    matplotlib_module.rcParams["svg.hashsalt"] = SVG_HASH_SALT


def save_svg(figure: Any, path: str | Path, **kwargs: Any) -> None:
    """Save one SVG without volatile creation metadata."""

    metadata = dict(kwargs.pop("metadata", {}) or {})
    metadata.update(SVG_METADATA)
    save_figure(figure, path, format="svg", metadata=metadata, **kwargs)
