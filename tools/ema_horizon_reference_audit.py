#!/usr/bin/env python3
"""Summarize dynamically computed EMA horizons on the reference simulations.

The three filter families print the final tuning operating point for each of the
same eight versioned wave records. With the deployed coefficients:

  T_sea = tau_target                    (all three currently use c_tau = 1)
  tau_common = 0.40 * T_sea
  tau_freq   = 0.10 * T_sea = 0.05 * T_z
  tau_var    = 2 * T_z = 4 * T_sea
  tau_reg    = c_reg * tau_target
  tau_period_estimator = 8 * T_z

where c_reg is 1.5 for OU-III and 3.0 for OU-II/TFG.

This tool intentionally derives the *unclamped requested* horizons from the
reference logs. It is an evidence tool for choosing safety bounds, not part of
the runtime implementation.
"""

from __future__ import annotations

import argparse
import math
import pathlib
import re
from dataclasses import dataclass

PROCESS_RE = re.compile(r"Processing\s+(.+?\.csv)")
TAU_RE = re.compile(r"tau_target=([+\-0-9.eE]+)")
WAVE_RE = re.compile(r"wave_period_s=([+\-0-9.eE]+)")

REG_MULT = {
    "ou3": 1.5,
    "ou2": 3.0,
    "tfg": 3.0,
}


@dataclass
class Row:
    family: str
    record: str
    tau_target: float
    wave_period: float

    @property
    def common(self) -> float:
        return 0.40 * self.tau_target

    @property
    def freq(self) -> float:
        return 0.05 * self.wave_period

    @property
    def variance(self) -> float:
        return 4.0 * self.tau_target

    @property
    def regularizer(self) -> float:
        return REG_MULT[self.family] * self.tau_target

    @property
    def period_estimator(self) -> float:
        return 8.0 * self.wave_period


def parse_log(family: str, path: pathlib.Path) -> list[Row]:
    current_record: str | None = None
    current_tau: float | None = None
    out: list[Row] = []

    for line in path.read_text(errors="replace").splitlines():
        m = PROCESS_RE.search(line)
        if m:
            current_record = pathlib.Path(m.group(1)).name
            current_tau = None
            continue

        m = TAU_RE.search(line)
        if m:
            try:
                current_tau = float(m.group(1))
            except ValueError:
                current_tau = None
            continue

        m = WAVE_RE.search(line)
        if m and current_record is not None and current_tau is not None:
            try:
                wave = float(m.group(1))
            except ValueError:
                continue
            if math.isfinite(current_tau) and current_tau > 0.0 and math.isfinite(wave) and wave > 0.0:
                out.append(Row(family, current_record, current_tau, wave))
                current_tau = None

    if len(out) != 8:
        raise SystemExit(f"{family}: expected 8 reference records, parsed {len(out)} from {path}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--log",
        action="append",
        required=True,
        metavar="FAMILY=PATH",
        help="family is one of ou3, ou2, tfg",
    )
    args = ap.parse_args()

    rows: list[Row] = []
    for spec in args.log:
        family, sep, raw_path = spec.partition("=")
        if not sep or family not in REG_MULT:
            raise SystemExit(f"bad --log {spec!r}; expected ou3=..., ou2=..., or tfg=...")
        rows.extend(parse_log(family, pathlib.Path(raw_path)))

    print(
        "family,record,tau_target_s,wave_period_s,freq_ema_s,common_ema_s,"
        "variance_ema_s,regularizer_ema_s,period_estimator_ema_s"
    )
    for r in rows:
        print(
            f"{r.family},{r.record},{r.tau_target:.6f},{r.wave_period:.6f},"
            f"{r.freq:.6f},{r.common:.6f},{r.variance:.6f},{r.regularizer:.6f},"
            f"{r.period_estimator:.6f}"
        )

    channels = {
        "frequency": [r.freq for r in rows],
        "common_tau_sigma": [r.common for r in rows],
        "variance": [r.variance for r in rows],
        "regularizer": [r.regularizer for r in rows],
        "wave_period_estimator": [r.period_estimator for r in rows],
    }
    print("\nGLOBAL_HORIZON_RANGES")
    all_values: list[float] = []
    for name, values in channels.items():
        all_values.extend(values)
        print(f"{name}: min={min(values):.6f}s max={max(values):.6f}s")
    print(f"all_dynamic: min={min(all_values):.6f}s max={max(all_values):.6f}s")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
