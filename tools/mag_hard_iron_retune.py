#!/usr/bin/env python3
"""Re-tune the continuous magnetic hard-iron correction against the filter as
it currently stands.

The correction's defaults (`mag_hi_*`) were fitted when it was introduced, on
a filter that has since changed underneath it: the r_S adaptation law, the
S_factor = 1 anisotropy constant, and the MEKF sensor-variance sweep all moved
the attitude solution the estimator's answer is applied to.  This tool
re-measures the knobs one factor at a time, and then over the shortlist of
joint settings the one-factor pass suggests.

Scoring follows `tools/mag_gate_study.py`: absolute yaw is dominated by the
per-draw hard-iron offset, which no setting of these knobs can touch, so every
number is a change relative to the shipped configuration *within* a draw, and
draws are pooled only after that normalization.

Two seed axes matter here and they are not interchangeable:

  --axis init   redraws the magnetometer calibration (offset and distortion)
                via W3D_INIT_SEED.  This is the axis the correction is
                actually tuned against -- it changes the thing being
                estimated.
  --axis imu    redraws the sensor noise via W3D_IMU_SEED, the correction's
                answer held fixed.  Use it to check a difference is not one
                noise realization.

What it found, OU-III over five calibration draws and the eight scored
records (see docs/continuous-mag-hard-iron.md for the full write-up):

  - The absolute floor `model_ridge` dominated the shrinkage and should not
    have.  At 4e-3 it was ten times the calmest record's excitation, so it was
    a fixed ridge rather than a floor, shrinking hardest in exactly the calm
    seas where the standing yaw error is largest.  Cutting it to 5e-4 is worth
    10% of pooled yaw and 10% of the worst record's, and it replicates on the
    IMU-noise axis (-11.8%).  That is the change that shipped.

  - `model_ridge_relative` looked worth cutting too -- the joint grid put the
    pooled optimum at 0.3 rather than 0.5 -- and it did not replicate.  On the
    noise axis that point is worth -1.4% against the floor cut's -11.8%: it
    was fitting the five calibration draws.  Two axes are not a luxury here.

  - `apply_fraction` is the same lever as the ridge seen from outside, and
    points the same way; `memory_sec`, `slew_tau_sec` and the two admission
    gates are flat or safety-critical and are left alone.

Usage:

    python3 tools/mag_hard_iron_retune.py --sim ou3 --axis init --seeds 5
    python3 tools/mag_hard_iron_retune.py --sim ou3 --axis imu --only-configs \
        --config ridge=5e-4
"""
import argparse
import concurrent.futures
import json
import os
import re
import statistics
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
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
# What the applied correction leaves of the injected offset.  The simulators
# report the wrapper's applied hard iron in this slot, so this is the estimator
# scored against truth directly rather than through the heading it moves.
PAT_MAGB = re.compile(
    r"Bias error RMS \(mag, uT\): X=[-\d.eE+]+ Y=[-\d.eE+]+ Z=[-\d.eE+]+ "
    r"\|3D\|=([-\d.eE+]+)")

# Draws of the magnetometer calibration.  "" is the shipped default draw; the
# rest are the same four the robustness study uses, so the numbers here can be
# read next to it.
INIT_SEEDS = ["", "11", "23", "202", "3033"]
IMU_SEEDS = ["", "11", "202", "3033", "40404"]

METRICS = ("roll", "pitch", "yaw", "z", "d3", "magb")


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
    return {
        "roll": roll,
        "pitch": pitch,
        "yaw": yaw,
        "z": float(PAT_Z.search(out).group(1)),
        "d3": float(PAT_3D.search(out).group(1)),
        "magb": float(PAT_MAGB.search(out).group(1)) if PAT_MAGB.search(out) else float("nan"),
        "record": record.name,
    }


def score(sim, env_extra, seed, seed_var, pool):
    env = dict(env_extra)
    if seed:
        env[seed_var] = seed
    rows = [r for r in pool.map(lambda p: run_one(sim, p, env), records()) if r]
    if not rows:
        return None
    out = {k: statistics.fmean(r[k] for r in rows) for k in METRICS}
    out["yaw_max"] = max(r["yaw"] for r in rows)
    out["rows"] = rows
    return out


ONE_FACTOR = [
    ("baseline (shipped)", {}),
    ("off (SF_MAG_CONT_HI=0)", {"SF_MAG_CONT_HI": "0"}),

    # The ridge that scales with the excitation.  This is the knob the
    # estimator's own documentation calls the one that matters.
    ("ridge_rel=0.0", {"SF_MAG_HI_RIDGE_REL": "0.0"}),
    ("ridge_rel=0.15", {"SF_MAG_HI_RIDGE_REL": "0.15"}),
    ("ridge_rel=0.25", {"SF_MAG_HI_RIDGE_REL": "0.25"}),
    ("ridge_rel=0.35", {"SF_MAG_HI_RIDGE_REL": "0.35"}),
    ("ridge_rel=0.5 (shipped)", {"SF_MAG_HI_RIDGE_REL": "0.5"}),
    ("ridge_rel=0.7", {"SF_MAG_HI_RIDGE_REL": "0.7"}),
    ("ridge_rel=1.0", {"SF_MAG_HI_RIDGE_REL": "1.0"}),
    ("ridge_rel=1.5", {"SF_MAG_HI_RIDGE_REL": "1.5"}),

    # The absolute floor, which only binds on a window with little excitation.
    ("ridge=0", {"SF_MAG_HI_RIDGE": "0"}),
    ("ridge=1e-3", {"SF_MAG_HI_RIDGE": "1e-3"}),
    ("ridge=8e-3", {"SF_MAG_HI_RIDGE": "8e-3"}),
    ("ridge=2e-2", {"SF_MAG_HI_RIDGE": "2e-2"}),

    ("memory=150", {"SF_MAG_HI_MEMORY_SEC": "150"}),
    ("memory=300", {"SF_MAG_HI_MEMORY_SEC": "300"}),
    ("memory=1200", {"SF_MAG_HI_MEMORY_SEC": "1200"}),
    ("memory=cumulative", {"SF_MAG_HI_MEMORY_SEC": "-1"}),

    ("fraction=0.5", {"SF_MAG_HI_FRACTION": "0.5"}),
    ("fraction=0.75", {"SF_MAG_HI_FRACTION": "0.75"}),
    ("fraction=1.25", {"SF_MAG_HI_FRACTION": "1.25"}),
    ("fraction=1.5", {"SF_MAG_HI_FRACTION": "1.5"}),

    ("slew_tau=10", {"SF_MAG_HI_SLEW_TAU": "10"}),
    ("slew_tau=20", {"SF_MAG_HI_SLEW_TAU": "20"}),
    ("slew_tau=90", {"SF_MAG_HI_SLEW_TAU": "90"}),
    ("slew_tau=180", {"SF_MAG_HI_SLEW_TAU": "180"}),

    ("min_info=0.5", {"SF_MAG_HI_MIN_INFO": "0.5"}),
    ("min_info=8", {"SF_MAG_HI_MIN_INFO": "8"}),
    ("min_weight=200", {"SF_MAG_HI_MIN_WEIGHT": "200"}),
    ("min_weight=2000", {"SF_MAG_HI_MIN_WEIGHT": "2000"}),
    ("max_resid=2", {"SF_MAG_HI_MAX_RESID": "2"}),
    ("max_resid=6", {"SF_MAG_HI_MAX_RESID": "6"}),
]


def parse_joint(spec):
    """`ridge_rel=0.3,fraction=0.8` -> (label, env)."""
    keys = {
        "ridge_rel": "SF_MAG_HI_RIDGE_REL",
        "ridge": "SF_MAG_HI_RIDGE",
        "memory": "SF_MAG_HI_MEMORY_SEC",
        "fraction": "SF_MAG_HI_FRACTION",
        "slew_tau": "SF_MAG_HI_SLEW_TAU",
        "min_info": "SF_MAG_HI_MIN_INFO",
        "min_weight": "SF_MAG_HI_MIN_WEIGHT",
        "max_resid": "SF_MAG_HI_MAX_RESID",
    }
    env = {}
    for part in spec.split(","):
        k, _, v = part.partition("=")
        k = k.strip()
        if k not in keys:
            raise SystemExit(f"unknown knob {k!r}; known: {', '.join(keys)}")
        env[keys[k]] = v.strip()
    return spec, env


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim", choices=list(SIMS), default="ou3")
    ap.add_argument("--axis", choices=("init", "imu"), default="init")
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--jobs", type=int, default=os.cpu_count())
    ap.add_argument(
        "--config", action="append", default=[],
        help="extra joint setting, e.g. ridge_rel=0.3,fraction=0.8; repeatable")
    ap.add_argument(
        "--only-configs", action="store_true",
        help="score the baseline, the ablation and --config only")
    ap.add_argument("--json", type=Path, help="write the per-record rows here")
    args = ap.parse_args()

    sim = SIMS[args.sim]
    if not sim.exists():
        raise SystemExit(f"{sim} not built; run make -C {sim.parent}")

    seed_var = "W3D_INIT_SEED" if args.axis == "init" else "W3D_IMU_SEED"
    seeds = (INIT_SEEDS if args.axis == "init" else IMU_SEEDS)[: args.seeds]

    factors = list(ONE_FACTOR)
    if args.only_configs:
        factors = factors[:2]
    factors += [parse_joint(c) for c in args.config]

    dump = {}
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        base = {s: score(sim, {}, s, seed_var, pool) for s in seeds}

        print(f"### {args.sim}: continuous hard-iron retune, axis={args.axis}, "
              f"{len(seeds)} draws x {len(records())} records")
        print(f"### yaw RMS (deg), mean over records, per draw; "
              f"d% is the within-draw change, pooled")
        hdr = f"{'config':28s}" + "".join(
            f"{(seed_var.split('_')[1].lower() + ' ' + (s or 'dflt')):>11s}" for s in seeds)
        print(hdr + f"{'yaw d%':>8s}{'ywmx d%':>8s}{'roll d%':>8s}"
                    f"{'pitch d%':>9s}{'3D d%':>7s}{'worse':>7s}")
        for name, env in factors:
            cells, dy, dymx, dr, dp, d3, mb = [], [], [], [], [], [], []
            worse = 0
            for s in seeds:
                cur = score(sim, env, s, seed_var, pool) if env else base[s]
                b = base[s]
                cells.append(f"{cur['yaw']:11.3f}")
                dy.append(100.0 * (cur["yaw"] - b["yaw"]) / b["yaw"])
                dymx.append(100.0 * (cur["yaw_max"] - b["yaw_max"]) / b["yaw_max"])
                dr.append(100.0 * (cur["roll"] - b["roll"]) / b["roll"])
                dp.append(100.0 * (cur["pitch"] - b["pitch"]) / b["pitch"])
                d3.append(100.0 * (cur["d3"] - b["d3"]) / b["d3"])
                mb.append(cur["magb"])
                worse += sum(
                    1 for c, bb in zip(cur["rows"], b["rows"]) if c["yaw"] > bb["yaw"])
                dump.setdefault(name, {})[s or "dflt"] = [
                    {k: r[k] for k in (*METRICS, "record")} for r in cur["rows"]]
            md = statistics.fmean(dy)
            mark = "  <== better" if md < -2.0 else ("  (worse)" if md > 2.0 else "")
            print(f"{name:28s}" + "".join(cells)
                  + f"{md:+8.1f}{statistics.fmean(dymx):+8.1f}"
                  + f"{statistics.fmean(dr):+8.1f}{statistics.fmean(dp):+9.1f}"
                  + f"{statistics.fmean(d3):+7.1f}{statistics.fmean(mb):7.3f}"
                  + f"{worse:4d}/{len(seeds) * len(records()):<3d}" + mark)

    if args.json:
        args.json.write_text(json.dumps(dump, indent=1))
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    sys.exit(main())
