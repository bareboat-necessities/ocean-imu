#!/usr/bin/env python3
"""Sensitivity sweep of the deployed OU-II PhysicalMSE coefficient.

This is the analytical successor to the retired empirical-vs-MSE comparison.
It keeps every exponent and the spectral channel ratio fixed at the deployed
PhysicalMSE derivation and perturbs only the common coefficient C_P.  The
analytical C_P=0.1116 arm is the reference; no legacy adaptation law is run.
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
C_P_ANALYTICAL = 0.1116
C_P_GRID = (0.070, 0.085, 0.098, C_P_ANALYTICAL, 0.130, 0.150, 0.175)


def run_record(binary: Path, record: Path, env_extra: dict[str, str]) -> dict[str, float]:
    env = dict(os.environ)
    env["W3D_WRITE_TIMESERIES"] = "0"
    env["W3D_COLLECT_ALL_GATES"] = "1"
    env.update(env_extra)
    completed = subprocess.run(
        [str(binary), "--input", str(record)],
        cwd=binary.parent,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    out = completed.stdout
    row: dict[str, float] = {}
    for name, pattern in PATTERNS.items():
        if name == "rs_applied":
            continue
        match = pattern.search(out)
        if match is None:
            tail = "\n".join(out.splitlines()[-20:])
            raise RuntimeError(f"{record.name}: no {name} in simulator output\n{tail}")
        row[name] = float(match.group(1))
    return row


def parse_grid(text: str | None) -> tuple[float, ...]:
    if not text:
        return C_P_GRID
    values = tuple(float(v) for v in text.replace(",", " ").split())
    if C_P_ANALYTICAL not in values:
        raise ValueError(f"grid must include analytical reference {C_P_ANALYTICAL}")
    return values


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--cp-grid")
    parser.add_argument(
        "--ratio",
        type=float,
        default=None,
        help="override analytical C_P/C_V for a dedicated ratio sensitivity check",
    )
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--record-dir", type=Path, default=DEFAULT_RECORD_DIR)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "reports" / "results" / "ou2_physical_mse_scale",
    )
    args = parser.parse_args(argv)

    grid = parse_grid(args.cp_grid)
    binary = args.binary.resolve()
    if not binary.exists():
        raise SystemExit(f"simulator not built: {binary}")
    record_dir = args.record_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict[str, Any]] = []
    grid_rows: list[dict[str, Any]] = []
    for c_p in grid:
        env = {"OU_II_PSEUDO_MSE_COEFF": repr(c_p)}
        if args.ratio is not None:
            env["OU_II_PSEUDO_MSE_RATIO"] = repr(args.ratio)
        label = f"C_P={c_p:g}"
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [
                pool.submit(run_record, binary, record_dir / name, env)
                for _, _, name in RECORDS
            ]
            results = [future.result() for future in futures]
        rows = [
            {"arm": label, "family": fam, "hs_m": hs, "record": name, **metrics}
            for (fam, hs, name), metrics in zip(RECORDS, results)
        ]
        all_rows.extend(rows)
        summary = {"arm": label, **summarize(rows)}
        summary["mean_roll_deg"] = sum(r["roll_deg"] for r in rows) / len(rows)
        summary["mean_pitch_deg"] = sum(r["pitch_deg"] for r in rows) / len(rows)
        grid_rows.append(summary)
        print(
            f"{label:<14} mean_z={summary['mean_z_pct_hs']:.4f} "
            f"max_z={summary['max_z_pct_hs']:.4f} "
            f"mean_3d={summary['mean_rms_3d_m']:.5f} "
            f"roll={summary['mean_roll_deg']:.4f} "
            f"pitch={summary['mean_pitch_deg']:.4f}",
            flush=True,
        )

    for name, rows in (("raw.csv", all_rows), ("grid.csv", grid_rows)):
        with (args.output_dir / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)

    reference = next(
        row for row in grid_rows if row["arm"] == f"C_P={C_P_ANALYTICAL:g}"
    )
    print("\nagainst analytical C_P=0.1116 (negative is better):")
    print(
        f"{'arm':<14}{'d mean_z':>10}{'d max_z':>10}{'d mean_3d':>11}"
        f"{'d roll':>9}{'d pitch':>9}"
    )
    for row in grid_rows:
        print(
            f"{row['arm']:<14}"
            f"{row['mean_z_pct_hs'] - reference['mean_z_pct_hs']:+10.4f}"
            f"{row['max_z_pct_hs'] - reference['max_z_pct_hs']:+10.4f}"
            f"{row['mean_rms_3d_m'] - reference['mean_rms_3d_m']:+11.5f}"
            f"{row['mean_roll_deg'] - reference['mean_roll_deg']:+9.4f}"
            f"{row['mean_pitch_deg'] - reference['mean_pitch_deg']:+9.4f}"
        )
    print(f"\nwrote {args.output_dir}/raw.csv and grid.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
