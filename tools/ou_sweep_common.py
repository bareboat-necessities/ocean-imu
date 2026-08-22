"""Shared record inventory and simulator-output parsing for OU sweep utilities.

This module contains no regularizer-law assumptions.  It exists so current
research tools can share the eight reference records and common RMS parsers
without depending on a retired OU-III tuning study.
"""

from __future__ import annotations

import re
from typing import Any


# The eight scored records, in the order used by the OU publication tables.
RECORDS = (
    ("JONSWAP", 0.27, "wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv"),
    ("JONSWAP", 1.50, "wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv"),
    ("JONSWAP", 4.00, "wave_data_jonswap_H4.000_L112.766_A30.00_P30.00.csv"),
    ("JONSWAP", 8.50, "wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv"),
    ("PM-Stokes", 0.27, "wave_data_pmstokes_H0.270_L14.047_A30.00_P60.00.csv"),
    ("PM-Stokes", 1.50, "wave_data_pmstokes_H1.500_L50.710_A-30.00_P120.00.csv"),
    ("PM-Stokes", 4.00, "wave_data_pmstokes_H4.000_L112.766_A30.00_P30.00.csv"),
    ("PM-Stokes", 8.50, "wave_data_pmstokes_H8.500_L202.839_A-30.00_P72.00.csv"),
)

FLOAT = r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
PATTERNS = {
    "z_rms_m": re.compile(rf"XYZ RMS \(m\): X={FLOAT} Y={FLOAT} Z=({FLOAT})"),
    "z_pct_hs": re.compile(rf"XYZ RMS \(%Hs\): X={FLOAT}% Y={FLOAT}% Z=({FLOAT})%"),
    "rms_3d_m": re.compile(rf"3D RMS \(m\):\s*({FLOAT})"),
    "roll_deg": re.compile(rf"Angles RMS \(deg\): Roll=({FLOAT})"),
    "pitch_deg": re.compile(rf"Angles RMS \(deg\): Roll={FLOAT} Pitch=({FLOAT})"),
    "yaw_deg": re.compile(rf"Angles RMS \(deg\): Roll={FLOAT} Pitch={FLOAT} Yaw=({FLOAT})"),
    "rs_applied": re.compile(rf"RS_applied=({FLOAT})"),
    "tau_applied": re.compile(rf"tau_applied=({FLOAT})"),
    "sigma_applied": re.compile(rf"sigma_applied=({FLOAT})"),
}


def summarize(rows: list[dict[str, Any]]) -> dict[str, float]:
    """Return common vertical and 3-D RMS summary metrics for a sweep cell."""
    z = [float(row["z_pct_hs"]) for row in rows]
    d3 = [float(row["rms_3d_m"]) for row in rows]
    return {
        "mean_z_pct_hs": sum(z) / len(z),
        "max_z_pct_hs": max(z),
        "mean_rms_3d_m": sum(d3) / len(d3),
        "max_rms_3d_m": max(d3),
    }
