#!/usr/bin/env python3
"""Deterministic Matplotlib SVG export for published repository artifacts."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
from vector_export import save_figure


SVG_HASH_SALT = "ocean-imu-published-svg-v1"
SVG_METADATA = {
    "Date": None,
    "Creator": "ocean-imu deterministic Matplotlib export",
}


def configure_svg(matplotlib_module: Any) -> None:
    """Fix Matplotlib's SVG id salt before any SVG figures are rendered, and
    decimate long lines of every saved figure (figure_decimation.py)."""

    matplotlib_module.rcParams["svg.hashsalt"] = SVG_HASH_SALT
    spec = importlib.util.spec_from_file_location("figure_decimation", Path(__file__).with_name("figure_decimation.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.install(matplotlib_module)


def save_svg(figure: Any, path: str | Path, **kwargs: Any) -> None:
    """Save one SVG without volatile creation metadata."""

    metadata = dict(kwargs.pop("metadata", {}) or {})
    metadata.update(SVG_METADATA)
    save_figure(figure, path, format="svg", metadata=metadata, **kwargs)
