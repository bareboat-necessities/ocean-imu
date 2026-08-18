#!/usr/bin/env python3
"""Two-dimensional (c_sigma, C_R) calibration of the OU-III integral regularizer.

The experiment this drives replaces the acceleration-amplitude multiplier in
the integral-regularizer base schedule with the accelerometer noise scale,

    old:  r_S,base = c_R sigma_aw tau^3
    new:  r_S,base = C_R sqrt(R_a) tau^3

following Eq. (adapt-rs-base) of the OU-III adaptation derivation.  The
pseudo-update cadence renormalization r_S,filter = r_S,base sqrt(T_S,0/T_S)
with T_S ~ tau is unchanged, so the applied law stays r_S,filter ~ tau^(5/2) in
both cases.  The OU-prior mapping sigma_aw = c_sigma sigma_a,B is untouched;
sigma_aw simply no longer sets pseudo-measurement strength.

Both coefficients have analytical reference points the derivation supplies, and
this sweep exists to measure the complete-MEKF optimum against them:

  C_R      = sqrt(2h/T_S,0) / kappa^3  ~= 17.11  for h = 5 ms, T_S,0 = 15 ms,
             kappa = 0.3627.  C_R is a pole placement, not a gain: it sets the
             normalized regularizer corner kappa = omega_R tau directly, and at
             the analytical value the cadence-normalized cubic base is
             identical to the StrongRiccati law.
  c_sigma  = F_OU^(-1/2) ~= 1.80 from band-energy matching, with
             F_OU = (2/pi)[atan(4 pi) - atan(pi/2)] ~= 0.3104 for the applied
             [0.5, 4] f_z band at c_tau = 1.  The deployed empirical value is
             0.9, so the grid brackets both.

The scalar reduction behind those references omits attitude/gravity leakage,
residual bias, three-axis covariance coupling and the cadence clamps, so the
sweep is the arbiter and the references are what the outcome is reported
against.  c_tau is held at 1 throughout.

This driver runs the deterministic eight-record protocol (default seeds,
last-900 s RMS window, one realization per record) over a grid of
(c_sigma, C_R) and reports the primary endpoint of the OU-III paper, the
vertical displacement RMS in percent of H_s.  The deterministic protocol is the
screening instrument only; the winning pair is confirmed against main with the
paired multi-seed harness (tools/ou_validation.py).

Typical use:

    python3 tools/ou_rs_amplitude_retune_sweep.py --stage coarse --jobs 4
    python3 tools/ou_rs_amplitude_retune_sweep.py --stage fine   --jobs 4
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BINARY = REPO_ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim"
DEFAULT_RECORD_DIR = REPO_ROOT / "plots" / "kalman_ou_iii"

# The eight scored records, in the order the paper's simulation table uses.
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

# Analytical references the derivation supplies; see the module docstring.
C_R_ANALYTICAL = 17.11
C_SIGMA_BAND_MATCHED = 1.80
C_SIGMA_DEPLOYED = 0.90

# Coarse grid, half-octave in both axes.  C_R spans a factor of sixteen around
# the analytical 17.11.  c_sigma is laid out so that both the deployed 0.90 and
# the band-matching 1.80 are grid points rather than interpolated, and the grid
# runs a full octave past each of them.
COARSE_SIGMA = (0.32, 0.45, 0.64, 0.90, 1.27, 1.80, 2.55, 3.60)
COARSE_R = (4.28, 6.05, 8.56, 12.10, 17.11, 24.20, 34.22, 48.40, 68.44)

# Fine grid: quarter-octave refinement around the coarse optimum, which put
# C_R at the analytical 17.11 in every c_sigma row.  The C_R axis brackets it
# on both sides; the c_sigma axis is refined below the coarse floor, because
# the coarse mean endpoint was still falling at its smallest value, while
# keeping the deployed 0.90 and the band-matching 1.80 as anchors so the
# refinement can still be read against both references.
FINE_SIGMA = (0.22, 0.27, 0.32, 0.38, 0.45, 0.64, 0.90, 1.80)
FINE_R = (12.1, 14.4, 17.11, 20.3, 24.2)


def run_record(binary: Path, record: Path, env_extra: dict[str, str]) -> dict[str, float]:
    env = dict(os.environ)
    env["W3D_WRITE_TIMESERIES"] = "0"
    # Score every record, including any that trips a historical gate: the
    # sweep reports what the filter did, not whether it passed.
    env["W3D_COLLECT_ALL_GATES"] = "1"
    env.update(env_extra)
    completed = subprocess.run(
        [str(binary), "--input", str(record)],
        cwd=binary.parent, env=env, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    out = completed.stdout
    row: dict[str, float] = {}
    for name, pattern in PATTERNS.items():
        match = pattern.search(out)
        if match is None:
            tail = "\n".join(out.splitlines()[-20:])
            raise RuntimeError(f"{record.name}: no {name} in simulator output\n{tail}")
        row[name] = float(match.group(1))
    return row


def run_config(binary: Path, record_dir: Path, c_sigma: float, c_R: float,
               jobs: int) -> list[dict[str, Any]]:
    env_extra = {
        "OU_SIGMA_COEFF": repr(c_sigma),
        "OU_III_R_S_COEFF": repr(c_R),
    }
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [
            pool.submit(run_record, binary, record_dir / name, env_extra)
            for _, _, name in RECORDS
        ]
        results = [f.result() for f in futures]
    rows = []
    for (family, hs, name), metrics in zip(RECORDS, results):
        rows.append({"c_sigma": c_sigma, "c_R": c_R, "family": family,
                     "hs_m": hs, "record": name, **metrics})
    return rows


def summarize(rows: list[dict[str, Any]]) -> dict[str, float]:
    z = [r["z_pct_hs"] for r in rows]
    d3 = [r["rms_3d_m"] for r in rows]
    return {
        "mean_z_pct_hs": sum(z) / len(z),
        "max_z_pct_hs": max(z),
        "mean_rms_3d_m": sum(d3) / len(d3),
        "max_rms_3d_m": max(d3),
    }


def parse_grid(text: str | None, default: tuple[float, ...]) -> tuple[float, ...]:
    if not text:
        return default
    return tuple(float(v) for v in text.replace(",", " ").split())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("coarse", "fine"), default="coarse")
    ap.add_argument("--sigma-grid", help="explicit c_sigma values (overrides --stage)")
    ap.add_argument("--r-grid", help="explicit C_R values (overrides --stage)")
    ap.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    ap.add_argument("--record-dir", type=Path, default=DEFAULT_RECORD_DIR)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--output-dir", type=Path,
                    default=REPO_ROOT / "reports" / "results" / "ou_rs_amplitude_retune")
    args = ap.parse_args(argv)

    if args.stage == "coarse":
        sigmas = parse_grid(args.sigma_grid, COARSE_SIGMA)
        radii = parse_grid(args.r_grid, COARSE_R)
    else:
        sigmas = parse_grid(args.sigma_grid, FINE_SIGMA)
        radii = parse_grid(args.r_grid, FINE_R)

    binary = args.binary.resolve()
    if not binary.exists():
        raise SystemExit(f"simulator not built: {binary}")
    record_dir = args.record_dir.resolve()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.output_dir / f"{args.stage}_raw.csv"
    grid_path = args.output_dir / f"{args.stage}_grid.csv"

    all_rows: list[dict[str, Any]] = []
    grid_rows: list[dict[str, Any]] = []

    total = len(sigmas) * len(radii)
    done = 0
    for c_sigma in sigmas:
        for c_R in radii:
            rows = run_config(binary, record_dir, c_sigma, c_R, args.jobs)
            all_rows.extend(rows)
            summary = {"c_sigma": c_sigma, "c_R": c_R, **summarize(rows)}
            grid_rows.append(summary)
            done += 1
            print(f"[{done:3d}/{total}] c_sigma={c_sigma:<6g} C_R={c_R:<7g} "
                  f"mean_z={summary['mean_z_pct_hs']:.4f}%Hs "
                  f"max_z={summary['max_z_pct_hs']:.4f}%Hs "
                  f"mean_3d={summary['mean_rms_3d_m']:.5f} m", flush=True)

    with raw_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)
    with grid_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(grid_rows[0].keys()))
        writer.writeheader()
        writer.writerows(grid_rows)

    grid_rows.sort(key=lambda r: r["mean_z_pct_hs"])
    print(f"\nwrote {raw_path}\nwrote {grid_path}\n")
    print("best 10 by mean vertical RMS (%Hs) over the eight records:")
    print(f"{'c_sigma':>8} {'C_R':>8} {'mean_z':>9} {'max_z':>9} {'mean_3d':>9}")
    for r in grid_rows[:10]:
        print(f"{r['c_sigma']:8g} {r['c_R']:8g} {r['mean_z_pct_hs']:9.4f} "
              f"{r['max_z_pct_hs']:9.4f} {r['mean_rms_3d_m']:9.5f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
