#!/usr/bin/env python3
"""Shared sampling/downsampling settings for plot generation."""

import os
from typing import Optional

BASE_SAMPLE_RATE_HZ = 200.0
DEFAULT_PLOT_RATE_HZ = BASE_SAMPLE_RATE_HZ / 40


def get_plot_rate_hz() -> float:
    """Return target plot sample rate from env, defaulting to 200 Hz."""
    raw = os.getenv("PLOT_SAMPLE_RATE_HZ", "").strip()
    if not raw:
        return DEFAULT_PLOT_RATE_HZ
    try:
        value = float(raw)
    except ValueError:
        return DEFAULT_PLOT_RATE_HZ
    if value <= 0.0:
        return DEFAULT_PLOT_RATE_HZ
    return min(value, BASE_SAMPLE_RATE_HZ)


def get_decimation_step(base_rate_hz: float = BASE_SAMPLE_RATE_HZ, target_rate_hz: Optional[float] = None) -> int:
    """Compute integer keep-every-N-samples step for plotting."""
    rate = get_plot_rate_hz() if target_rate_hz is None else target_rate_hz
    if rate <= 0.0 or rate >= base_rate_hz:
        return 1
    step = int(round(base_rate_hz / rate))
    return max(1, step)
