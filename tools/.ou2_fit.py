#!/usr/bin/env python3
"""Scratch one-factor-at-a-time refit of the OU-II wave-band operating point.

Scores the four stationary JONSWAP records on the 900 s window and reports the
seed-mean normalized vertical error (primary) and mean 3D RMS (secondary).
"""
import argparse
import itertools
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "tests" / "kalman_ou_ii" / "kalman_ou_ii-sim"
DATA = ROOT / "plots" / "kalman_ou_ii"

RECORDS = [
    "wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv",
    "wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv",
    "wave_data_jonswap_H4.000_L112.766_A30.00_P30.00.csv",
    "wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv",
]

SEEDS = [(0, 0), (101, 1009), (211, 1103)]


def run(record, overrides, imu_seed, init_seed):
    env = os.environ.copy()
    env.update(
        {
            "W3D_VALIDATION_WINDOW_SEC": "900",
            "W3D_COLLECT_ALL_GATES": "1",
            "W3D_WRITE_TIMESERIES": "0",
            "W3D_IMU_SEED": str(imu_seed),
            "W3D_INIT_SEED": str(init_seed),
        }
    )
    env.update({k: f"{v:.9g}" for k, v in overrides.items()})
    out = subprocess.run(
        [str(BIN), "--input", str(DATA / record)],
        cwd=str(BIN.parent),
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    ).stdout
    line = next((l for l in out.splitlines() if l.startswith("VALIDATION_METRICS")), None)
    if line is None:
        raise RuntimeError(f"no metrics for {record} {overrides}\n{out[-2000:]}")
    fields = dict(
        kv.split("=", 1) for kv in line.split() if "=" in kv
    )
    return {
        "z": float(fields["disp_z_pct_hs"]),
        "d3": float(fields["disp_3d_rms_m"]),
        "roll": float(fields["roll_rms_deg"]),
        "tau": float(fields["tau_applied_s"]),
        "sigma": float(fields["sigma_applied_mps2"]),
    }


def score(overrides, pool):
    jobs = [
        pool.submit(run, rec, overrides, imu, init)
        for rec in RECORDS
        for imu, init in SEEDS
    ]
    rows = [j.result() for j in jobs]
    n = len(rows)
    return (
        sum(r["z"] for r in rows) / n,
        sum(r["d3"] for r in rows) / n,
        sum(r["roll"] for r in rows) / n,
        [round(r["tau"], 2) for r in rows[::len(SEEDS)]],
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--param", required=True)
    ap.add_argument("--values", required=True)
    ap.add_argument("--base", default="")
    args = ap.parse_args()

    base = {}
    for item in filter(None, args.base.split(",")):
        k, v = item.split("=")
        base[k] = float(v)

    with ThreadPoolExecutor(max_workers=4) as pool:
        for value in [float(v) for v in args.values.split(",")]:
            overrides = dict(base)
            overrides[args.param] = value
            z, d3, roll, taus = score(overrides, pool)
            print(f"{args.param}={value:<8g} meanZ%Hs={z:7.4f}  mean3D={d3:7.4f}  roll={roll:6.4f}  tau={taus}", flush=True)


if __name__ == "__main__":
    main()
