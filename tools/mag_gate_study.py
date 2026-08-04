#!/usr/bin/env python3
"""One-factor-at-a-time study of the magnetometer start gating, scored on yaw.

Yaw error is dominated by the per-seed hard-iron draw, which no gating setting
can touch, so absolute yaw is not comparable across seeds.  Everything here is
therefore reported as a change relative to the shipped configuration *within*
each seed, and aggregated over seeds only after that normalization.

What it found, OU-III over five seeds and the eight scored records:

  - No setting of any gate improves yaw usefully.  The largest reductions are
    about 3%, which on a mean of 1.65 deg is 0.05 deg, against the roughly
    2.1 deg the hard iron puts there.  Opening the gravity gate, shortening the
    hold, or halving the gate's accelerometer low-pass all land in that band.

  - Every one of those settings costs roll, heavily.  Opening the gravity gate
    fully trades 3.0% of yaw for 127% more roll error; a hold of zero trades
    1.1% for 96%.  The gate is not merely admitting fewer samples, it is
    protecting tilt.

The mechanism is that the magnetometer update constrains all three attitude
degrees of freedom against the learned reference, not heading alone.  The
reference's dip is only as good as the tilt estimate at the moment it was
learned, so a reference learned through a poor-tilt window drags roll and pitch
for the rest of the run.  The gravity gate exists to stop that, and the numbers
above are what happens when it is relaxed.

The obvious remedy -- let the accelerometer own tilt and the magnetometer own
heading -- was prototyped by projecting the attitude block of the mag Kalman
gain onto the world vertical.  That is not sufficient: the residual the
attitude can no longer absorb is taken up by the gyro bias instead, whose error
went from 35% to 166% of the true bias, and six of the eight records failed
their gates.  Doing it properly means a one-degree-of-freedom heading
measurement with its own Jacobian and innovation covariance, so that the whole
state update stays consistent, rather than a projection applied after the gain.
That is a change to the measurement model and is not attempted here.

Usage:

    python3 tools/mag_gate_study.py --sim ou3
"""
import argparse
import concurrent.futures
import os
import re
import statistics
import subprocess
import sys
from pathlib import Path

REPO = Path("/home/user/ocean-imu")
DATA = REPO / "plots" / "kalman_ou_ii"
SIMS = {
    "ou2": REPO / "tests" / "kalman_ou_ii" / "kalman_ou_ii-sim",
    "ou3": REPO / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim",
}

PAT_ANG = re.compile(
    r"Angles RMS \(deg\): Roll=([-\d.eE+]+) Pitch=([-\d.eE+]+) Yaw=([-\d.eE+]+)"
)
PAT_Z = re.compile(r"XYZ RMS \(%Hs\): X=[-\d.eE+]+% Y=[-\d.eE+]+% Z=([-\d.eE+]+)%")
PAT_3D = re.compile(r"3D % of max \|disp_ref\|_3D = ([-\d.eEnan+]+)%")

SEEDS = ["", "11", "202", "3033", "40404"]


def records():
    return sorted(
        p for p in DATA.glob("wave_data_*.csv")
        if "jonswap" in p.name or "pmstokes" in p.name
    )


def run_one(sim, record, env_extra):
    env = dict(os.environ)
    env["W3D_WRITE_TIMESERIES"] = "0"
    env["W3D_COLLECT_ALL_GATES"] = "1"
    env.update(env_extra)
    out = subprocess.run(
        [str(sim), "--input", str(record)],
        cwd=sim.parent, env=env, capture_output=True, text=True,
    ).stdout
    m = PAT_ANG.search(out)
    if not m:
        return None
    roll, pitch, yaw = (float(v) for v in m.groups())
    z = float(PAT_Z.search(out).group(1))
    d3 = float(PAT_3D.search(out).group(1))
    return {"roll": roll, "pitch": pitch, "yaw": yaw, "z": z, "d3": d3}


def score(sim, env_extra, seed, pool):
    env = dict(env_extra)
    if seed:
        env["W3D_SEED"] = seed
    rows = [r for r in pool.map(lambda p: run_one(sim, p, env), records()) if r]
    if not rows:
        return None
    return {
        k: statistics.fmean(r[k] for r in rows) for k in ("roll", "pitch", "yaw", "z", "d3")
    } | {"yaw_max": max(r["yaw"] for r in rows)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim", choices=list(SIMS), default="ou3")
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    args = ap.parse_args()
    sim = SIMS[args.sim]

    # Shipped configuration, then one knob moved at a time.
    factors = [
        ("baseline", {}),
        ("grav_max_sin=0.02", {"SF_MAG_GRAV_ALIGN_MAX_SIN": "0.02"}),
        ("grav_max_sin=0.04", {"SF_MAG_GRAV_ALIGN_MAX_SIN": "0.04"}),
        ("grav_max_sin=0.15", {"SF_MAG_GRAV_ALIGN_MAX_SIN": "0.15"}),
        ("grav_max_sin=0.30", {"SF_MAG_GRAV_ALIGN_MAX_SIN": "0.30"}),
        ("grav_max_sin=1.0 (no gate)", {"SF_MAG_GRAV_ALIGN_MAX_SIN": "1.0"}),
        ("fallback=0 (gate bypassed)", {"SF_MAG_TILT_FALLBACK_SEC": "0"}),
        ("fallback=5", {"SF_MAG_TILT_FALLBACK_SEC": "5"}),
        ("fallback=60", {"SF_MAG_TILT_FALLBACK_SEC": "60"}),
        ("hold=0", {"SF_MAG_GRAV_ALIGN_HOLD_SEC": "0"}),
        ("hold=1", {"SF_MAG_GRAV_ALIGN_HOLD_SEC": "1"}),
        ("hold=5", {"SF_MAG_GRAV_ALIGN_HOLD_SEC": "5"}),
        ("lpf_tau=0.25", {"SF_MAG_GRAV_ALIGN_LPF_TAU": "0.25"}),
        ("lpf_tau=0.5", {"SF_MAG_GRAV_ALIGN_LPF_TAU": "0.5"}),
        ("lpf_tau=2", {"SF_MAG_GRAV_ALIGN_LPF_TAU": "2"}),
        ("lpf_tau=4", {"SF_MAG_GRAV_ALIGN_LPF_TAU": "4"}),
        ("window=5", {"SF_MAG_MIN_WINDOW_SEC": "5"}),
        ("window=10", {"SF_MAG_MIN_WINDOW_SEC": "10"}),
        ("window=30", {"SF_MAG_MIN_WINDOW_SEC": "30"}),
        ("window=60", {"SF_MAG_MIN_WINDOW_SEC": "60"}),
        ("window=120", {"SF_MAG_MIN_WINDOW_SEC": "120"}),
        ("gyro_veto=10", {"SF_MAG_EXTREME_GYRO_DPS": "10"}),
        ("gyro_veto=100", {"SF_MAG_EXTREME_GYRO_DPS": "100"}),
        ("mag_delay=3", {"SF_MAG_DELAY_SEC": "3"}),
        ("mag_delay=15", {"SF_MAG_DELAY_SEC": "15"}),
    ]

    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        base = {s: score(sim, {}, s, pool) for s in SEEDS}

        print(f"### {args.sim}: mean yaw RMS (deg) over 8 records, per seed")
        hdr = f"{'config':30s}" + "".join(f"{('seed ' + (s or 'dflt')):>11s}" for s in SEEDS)
        print(hdr + f"{'mean d%':>9s}{'roll d%':>9s}{'3D d%':>8s}")
        for name, env in factors:
            deltas, rolls, d3s, cells = [], [], [], []
            for s in SEEDS:
                cur = score(sim, env, s, pool) if env else base[s]
                b = base[s]
                cells.append(f"{cur['yaw']:11.3f}")
                deltas.append(100.0 * (cur["yaw"] - b["yaw"]) / b["yaw"])
                rolls.append(100.0 * (cur["roll"] - b["roll"]) / b["roll"])
                d3s.append(100.0 * (cur["d3"] - b["d3"]) / b["d3"])
            mark = ""
            md = statistics.fmean(deltas)
            if md < -2.0:
                mark = "  <== better"
            elif md > 2.0:
                mark = "  (worse)"
            print(f"{name:30s}" + "".join(cells)
                  + f"{md:+9.1f}{statistics.fmean(rolls):+9.1f}"
                  + f"{statistics.fmean(d3s):+8.1f}{mark}")


if __name__ == "__main__":
    sys.exit(main())
