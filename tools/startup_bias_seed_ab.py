#!/usr/bin/env python3
"""A/B of the startup vertical accelerometer-bias seed in OU-II and OU-III.

The seed (``SF_STARTUP_BIAS_SEED=1``, off by default) runs a PII observer on
the startup Mahony proxy's vertical acceleration while the filter is not yet
Live and, at go-live, sets the world-down component of the MEKF's
accelerometer-bias mean to the observer's implicit bias ``-(b + ks * S)``.
See ``src/tuner/StartupVerticalBiasSeed.h``.

Every run is a full 20-minute replay of a pinned 28 ft sailboat RAO record
(all 20: Gerstner, cnoidal, Fenton, JONSWAP, PM-Stokes at four sea states),
with the deployed configuration otherwise, over five noise/initialization
seeds, with the nominal diesel vibration off and on.  Scored per run:

* startup window (35-180 s): displacement, roll, pitch and yaw error RMS;
* the gated trailing 900 s window: the same channels;
* the peak heave error in the first 300 s;
* the seed value against the true world-down bias at go-live.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "reports" / "results" / "startup_bias_seed"

FAMILIES = {
    "OU-II": ROOT / "tests/kalman_ou_ii/kalman_ou_ii-sim",
    "OU-III": ROOT / "tests/kalman_ou_iii/kalman_ou_iii-sim",
}
ARMS = {"deployed": {}, "seed": {"SF_STARTUP_BIAS_SEED": "1"}}
ENGINES = {"off": {}, "on": {"W3D_ENGINE_RPM": "2400", "W3D_ENGINE_LEVEL_MPS2": "0.6",
                             "W3D_ENGINE_BANDWIDTH_HZ": "80"}}
SEEDS = ("default", "11", "23", "202", "3033")
WAVES = ("gerstner", "cnoidal", "fenton", "jonswap", "pmstokes")
HEIGHTS = ("0.270", "1.500", "4.000", "8.500")
STARTUP_S = (35.0, 180.0)
GATED_S = 900.0
PEAK_S = 300.0

LIVE_RE = re.compile(r"STARTUP live_s=([-+0-9.eEnaif]+)")
SEED_RE = re.compile(r"STARTUP vertical_bias_seed_mps2=([-+0-9.eEnaif]+)")

COLS = ["time", "disp_ref_x", "disp_ref_y", "disp_ref_z", "disp_est_x", "disp_est_y",
        "disp_est_z", "roll_ref", "pitch_ref", "yaw_ref", "roll_est", "pitch_est",
        "yaw_est", "acc_bias_x", "acc_bias_y", "acc_bias_z"]


def find_records(data_dir: Path) -> list[Path]:
    out = []
    for wave in WAVES:
        for h in HEIGHTS:
            hits = sorted(data_dir.rglob(f"wave_data_{wave}_H{h}_*.csv"))
            if len(hits) != 1:
                raise FileNotFoundError(f"{wave} H{h}: found {len(hits)} records")
            out.append(hits[0])
    return out


def angle_err(est: np.ndarray, ref: np.ndarray) -> np.ndarray:
    return (est - ref + 180.0) % 360.0 - 180.0


def rms(x: np.ndarray) -> float:
    return float(np.sqrt(np.nanmean(x * x)))


def run_one(task, work: Path) -> dict:
    family, arm, engine, seed, record = task
    run_dir = work / family / arm / engine / seed / record.stem
    run_dir.mkdir(parents=True, exist_ok=True)
    for old in run_dir.iterdir():
        old.unlink()
    (run_dir / record.name).symlink_to(record.resolve())
    env = dict(os.environ, W3D_COLLECT_ALL_GATES="1", W3D_ALL_WAVE_TYPES="1",
               **ARMS[arm], **ENGINES[engine])
    if seed != "default":
        env["W3D_INIT_SEED"] = seed
    done = subprocess.run([str(FAMILIES[family])], cwd=run_dir, env=env,
                          capture_output=True, text=True)
    outs = list(run_dir.glob("w3d_*.csv"))
    if done.returncode not in (0, 1) or len(outs) != 1:
        raise RuntimeError(f"{task}: exit {done.returncode}\n{done.stderr[-2000:]}")
    frame = pd.read_csv(outs[0], usecols=COLS)
    for f in run_dir.iterdir():
        f.unlink()

    live = LIVE_RE.search(done.stderr)
    seed_val = SEED_RE.search(done.stderr)
    t = frame["time"].to_numpy()
    startup = (t >= STARTUP_S[0]) & (t <= STARTUP_S[1])
    gated = t >= t[-1] - GATED_S
    head = t <= PEAK_S
    row = {"family": family, "arm": arm, "engine": engine, "seed": seed,
           "record": record.stem.replace("wave_data_", "").split("_L")[0],
           "live_s": float(live.group(1)) if live else float("nan"),
           "seed_mps2": float(seed_val.group(1)) if seed_val else float("nan")}
    if np.isfinite(row["live_s"]):
        k = int(np.argmin(np.abs(t - row["live_s"])))
        row["true_bz_at_live_mps2"] = float(frame["acc_bias_z"].iloc[k])
    for axis in "xyz":
        err = (frame[f"disp_est_{axis}"] - frame[f"disp_ref_{axis}"]).to_numpy()
        row[f"{axis}_startup_rms_m"] = rms(err[startup])
        row[f"{axis}_gated_rms_m"] = rms(err[gated])
        row[f"{axis}_ref_gated_rms_m"] = rms(frame[f"disp_ref_{axis}"].to_numpy()[gated])
        if axis == "z":
            row["z_peak_300_m"] = float(np.nanmax(np.abs(err[head])))
            row["z_ref_startup_rms_m"] = rms(frame["disp_ref_z"].to_numpy()[startup])
    for ang in ("roll", "pitch", "yaw"):
        err = angle_err(frame[f"{ang}_est"].to_numpy(), frame[f"{ang}_ref"].to_numpy())
        row[f"{ang}_startup_rms_deg"] = rms(err[startup])
        row[f"{ang}_gated_rms_deg"] = rms(err[gated])
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, required=True)
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    records = find_records(args.data_dir)
    tasks = [(fam, arm, eng, seed, rec) for fam in FAMILIES for arm in ARMS
             for eng in ENGINES for seed in SEEDS for rec in records]
    with ThreadPoolExecutor(args.jobs) as pool:
        rows = list(pool.map(lambda task: run_one(task, args.work_dir), tasks))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with (args.output_dir / "startup_bias_seed_runs.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {len(rows)} runs to {args.output_dir}")


if __name__ == "__main__":
    main()
