#!/usr/bin/env python3
"""Calibrate the one free scale of the OU-II physical-MSE pseudo law.

tools/ou2_dual_mse_coefficients.py supplies C_P and C_P/C_V analytically, by
solving the derivation's own stationarity conditions on the eight reference
spectra.  Those two numbers come out of a *reduced* model: a scalar
double-integrator with one residual-acceleration intensity, no attitude or
gravity leakage, no accelerometer bias, no three-axis covariance coupling and
no cadence clamps.  The note that derives them says as much, and asks for a
full-state physical-H_2 calculation before the coefficients are called optimal.

This driver is the empirical stand-in for that calculation.  It sweeps the
single overall scale C_P while holding everything the derivation determines --
the 4/5 amplitude power, the 12/5 and 7/5 period powers, and the spectral
channel ratio C_P/C_V -- fixed.  The analytical value is a grid point, so the
outcome is reported against it rather than merely near it.

The deterministic eight-record protocol is the screening instrument here, as it
is for the other OU sweep utilities: default seeds, the trailing 900 s of each
replay, one realization per record, W3D_COLLECT_ALL_GATES=1 so a breach scores
the remaining records instead of exiting at the first.  The verdict is taken on
the paired multi-seed harness (tools/ou2_pseudo_law_compare.py).

Both displacement endpoints are printed, because they do not move together
under this parameter: the pseudo channels regularize all three axes, but the
vertical axis carries the gravity-leakage term the reduced model drops.

Typical use:

    python3 tools/ou2_pseudo_mse_scale_sweep.py --jobs 4
    python3 tools/ou2_pseudo_mse_scale_sweep.py --cp-grid 0.14,0.174,0.21
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ou_sweep_common import PATTERNS, RECORDS, summarize  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BINARY = REPO_ROOT / "tests" / "kalman_ou_ii" / "kalman_ou_ii-sim"
DEFAULT_RECORD_DIR = REPO_ROOT / "plots" / "kalman_ou_ii"

# Analytical prediction of the reduced dual-channel derivation.
C_P_ANALYTICAL = 0.1116
# The scale at which the law's r_p meets the empirical schedule's at the
# H_s = 1.5 m JONSWAP anchor; included so the sweep brackets both references.
# The two sit 3 % apart once q_eff is referred to the acceleration noise floor
# rather than to the accelerometer bench spec -- which is what the first run of
# this sweep established, by putting the optimum a factor of 1.52 above the
# analytical value while q_eff was still the bench figure.  See
# docs/ou-ii-dual-mse-adaptation.md.
C_P_EMPIRICAL_ANCHOR = 0.1146
C_P_GRID = (0.070, 0.085, 0.098, C_P_ANALYTICAL, C_P_EMPIRICAL_ANCHOR,
            0.130, 0.150, 0.175)


def run_record(binary: Path, record: Path,
               env_extra: dict[str, str]) -> dict[str, float]:
    env = dict(os.environ)
    env["W3D_WRITE_TIMESERIES"] = "0"
    env["W3D_COLLECT_ALL_GATES"] = "1"
    env.update(env_extra)
    completed = subprocess.run(
        [str(binary), "--input", str(record)],
        cwd=binary.parent, env=env, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
    out = completed.stdout
    row: dict[str, float] = {}
    for name, pattern in PATTERNS.items():
        # rs_applied is an OU-III channel; OU-II reports p0_S_applied instead.
        if name == "rs_applied":
            continue
        match = pattern.search(out)
        if match is None:
            tail = "\n".join(out.splitlines()[-20:])
            raise RuntimeError(
                f"{record.name}: no {name} in simulator output\n{tail}")
        row[name] = float(match.group(1))
    return row


def parse_grid(text: str | None, default: tuple[float, ...]) -> tuple[float, ...]:
    if not text:
        return default
    return tuple(float(v) for v in text.replace(",", " ").split())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cp-grid")
    ap.add_argument("--ratio", type=float, default=None,
                    help="override the channel ratio C_P/C_V (default: the "
                         "compiled-in spectral value)")
    ap.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    ap.add_argument("--record-dir", type=Path, default=DEFAULT_RECORD_DIR)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--output-dir", type=Path,
                    default=REPO_ROOT / "reports" / "results" / "ou2_pseudo_mse_scale")
    args = ap.parse_args(argv)

    grid = parse_grid(args.cp_grid, C_P_GRID)
    binary = args.binary.resolve()
    if not binary.exists():
        raise SystemExit(f"simulator not built: {binary}")
    record_dir = args.record_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict[str, Any]] = []
    grid_rows: list[dict[str, Any]] = []
    # The empirical law is run as arm zero so every C_P is read against the
    # schedule it replaces on identical records.
    arms: list[tuple[str, dict[str, str]]] = [
        ("empirical", {"OU_II_PSEUDO_LAW": "0"})]
    for c_p in grid:
        env = {"OU_II_PSEUDO_LAW": "1", "OU_II_PSEUDO_MSE_COEFF": repr(c_p)}
        if args.ratio is not None:
            env["OU_II_PSEUDO_MSE_RATIO"] = repr(args.ratio)
        arms.append((f"C_P={c_p:g}", env))

    for label, env_extra in arms:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(run_record, binary, record_dir / name, env_extra)
                       for _, _, name in RECORDS]
            results = [f.result() for f in futures]
        rows = [{"arm": label, "family": fam, "hs_m": hs, "record": name, **m}
                for (fam, hs, name), m in zip(RECORDS, results)]
        all_rows.extend(rows)
        summary = {"arm": label, **summarize(rows)}
        summary["mean_roll_deg"] = sum(r["roll_deg"] for r in rows) / len(rows)
        summary["mean_pitch_deg"] = sum(r["pitch_deg"] for r in rows) / len(rows)
        grid_rows.append(summary)
        print(f"{label:<14} mean_z={summary['mean_z_pct_hs']:.4f} "
              f"max_z={summary['max_z_pct_hs']:.4f} "
              f"mean_3d={summary['mean_rms_3d_m']:.5f} "
              f"roll={summary['mean_roll_deg']:.4f} "
              f"pitch={summary['mean_pitch_deg']:.4f}", flush=True)

    for name, rows in (("raw.csv", all_rows), ("grid.csv", grid_rows)):
        with (args.output_dir / name).open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    print(f"\nwrote {args.output_dir}/raw.csv and grid.csv")

    base = next(r for r in grid_rows if r["arm"] == "empirical")
    print("\nagainst the empirical law (negative is better):")
    print(f"{'arm':<14}{'d mean_z':>10}{'d max_z':>10}{'d mean_3d':>11}"
          f"{'d roll':>9}{'d pitch':>9}")
    for r in grid_rows:
        if r["arm"] == "empirical":
            continue
        print(f"{r['arm']:<14}"
              f"{r['mean_z_pct_hs'] - base['mean_z_pct_hs']:+10.4f}"
              f"{r['max_z_pct_hs'] - base['max_z_pct_hs']:+10.4f}"
              f"{r['mean_rms_3d_m'] - base['mean_rms_3d_m']:+11.5f}"
              f"{r['mean_roll_deg'] - base['mean_roll_deg']:+9.4f}"
              f"{r['mean_pitch_deg'] - base['mean_pitch_deg']:+9.4f}")
    print(f"\n  analytical C_P = {C_P_ANALYTICAL}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
