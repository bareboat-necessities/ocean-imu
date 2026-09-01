#!/usr/bin/env python3
"""Cross-family check that the front-end vibration guard behaves the same way.

``tools/ou3_engine_noise_mitigation.py`` is the full mitigation study: it sweeps
engine conditions, both guard stages, release behaviour and the corner choice,
and it is the generator of the committed evidence in
``reports/results/engine_noise_mitigation/``.  It runs OU-III, because that is
where the guard was designed.

This is the much smaller question that follows from arming the guard everywhere:
having established on OU-III what the guard does, does it do the same thing in
the other two families at the same operating point?  It replays the eight
stationary records through each family with the engine off and at the nominal
cruise condition, guard off and guard+R, and pools 3-D displacement error the
same way the degradation study does.

Two of the four cells per family are already committed evidence -- the
engine-off and unguarded-cruise columns are the degradation study's own numbers
-- so the run checks itself: if those do not come back equal to
``reports/results/engine_noise_degradation/engine_noise_summary.csv``, the two
studies are not on the same protocol and the comparison is void.

Nothing is written to disk.  The table it prints is the one at the end of the
mitigation section of ``docs/engine-noise-degradation.md``.

Usage:

    make -C tests/kalman_ou_ii kalman_ou_ii-sim
    make -C tests/kalman_ou_iii kalman_ou_iii-sim
    make -C tests/kalman_tfg kalman_tfg-sim
    python3 tools/engine_noise_guard_families.py --data-dir <wave-csv-dir>
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_WINDOW_SEC = 900.0
DEGRADATION_SUMMARY = (
    ROOT / "reports" / "results" / "engine_noise_degradation" / "engine_noise_summary.csv"
)

# The nominal cruise condition of the degradation study: the one condition it
# reports as a table rather than a sweep point.
CRUISE_RPM = "2400"

# Tolerance on the self-check against the committed degradation summary.  The
# replays are deterministic and the two studies pool the same way, so what is
# left is float summation order; a genuine protocol mismatch moves these cells
# by percent, not by parts per million.
SELF_CHECK_REL_TOL = 1e-4


@dataclass(frozen=True)
class Family:
    label: str
    simulator: Path
    env_prefix: str


FAMILIES = (
    Family("OU-II", ROOT / "tests" / "kalman_ou_ii" / "kalman_ou_ii-sim", "OU_II"),
    Family("OU-III", ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim", "OU_III"),
    Family("TFG", ROOT / "tests" / "kalman_tfg" / "kalman_tfg-sim", "TFG"),
)

# The same eight stationary records the degradation and mitigation studies use.
RECORDS = tuple(
    f"wave_data_{family}_H{tail}.csv"
    for family in ("jonswap", "pmstokes")
    for tail in (
        "0.270_L14.047_A30.00_P60.00",
        "1.500_L50.710_A-30.00_P120.00",
        "4.000_L112.766_A30.00_P30.00",
        "8.500_L202.839_A-30.00_P72.00",
    )
)

CONDITIONS = {"engine off": {}, "2400 rpm": {"W3D_ENGINE_RPM": CRUISE_RPM}}

RMS_3D_RE = re.compile(r"^3D RMS \(m\): ([0-9.eE+-]+)", re.MULTILINE)
YAW_RMS_RE = re.compile(r"^Angles RMS \(deg\): \S+ \S+ Yaw=([0-9.eE+-]+)", re.MULTILINE)
GUARD_RE = re.compile(r"^ACC_GUARD (.*)$", re.MULTILINE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True,
                        help="directory containing the versioned wave CSVs")
    parser.add_argument("--window-sec", type=float, default=DEFAULT_WINDOW_SEC,
                        help="trailing scoring window in seconds (default: 900)")
    parser.add_argument("--jobs", type=int, default=4,
                        help="parallel simulator processes (default: 4)")
    parser.add_argument("--no-self-check", action="store_true",
                        help="skip the comparison against the committed degradation summary")
    return parser.parse_args()


def find_record(data_dir: Path, name: str) -> Path:
    matches = list(data_dir.rglob(name))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"expected exactly one {name} under {data_dir}, found {len(matches)}"
        )
    return matches[0]


def arm_env(family: Family, arm: str) -> dict[str, str]:
    """Explicit in both arms, so neither inherits whatever the filter defaults to."""
    if arm == "off":
        return {
            f"{family.env_prefix}_ACC_GUARD_HZ": "0",
            f"{family.env_prefix}_ACC_GUARD_RACC_GAIN": "0",
        }
    # guard+R is the shipped configuration, so it is left to the filter
    # defaults rather than restated here: restating them would let this study
    # keep passing after the deployed point moved.
    return {}


def run_one(family: Family, input_path: Path, condition: str, arm: str,
            window_sec: float) -> dict[str, Any]:
    env = os.environ.copy()
    env.update({
        "W3D_WRITE_TIMESERIES": "0",
        "W3D_VALIDATION_WINDOW_SEC": f"{window_sec:.9g}",
        # The deployed gates were fitted to a vibration-free input; collect the
        # metrics rather than stopping on the first one this study crosses.
        "W3D_COLLECT_ALL_GATES": "1",
    })
    env.update(CONDITIONS[condition])
    env.update(arm_env(family, arm))

    completed = subprocess.run(
        [str(family.simulator), "--input", str(input_path.resolve())],
        cwd=family.simulator.parent, env=env, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False,
    )
    if completed.returncode not in (0, 1):
        tail = "\n".join(completed.stdout.splitlines()[-30:])
        raise RuntimeError(
            f"{family.label} {input_path.name} {condition} guard={arm}: "
            f"simulator exit {completed.returncode}\n{tail}"
        )
    if condition != "engine off" and "ENGINE_VIBRATION" not in completed.stdout:
        raise RuntimeError(f"{family.label}: engine model did not engage")
    if arm == "off" and "ACC_GUARD" in completed.stdout:
        raise RuntimeError(f"{family.label}: guard was not switched off")
    if arm != "off" and "ACC_GUARD" not in completed.stdout:
        raise RuntimeError(f"{family.label}: guard is not armed by default")

    disp = RMS_3D_RE.search(completed.stdout)
    yaw = YAW_RMS_RE.search(completed.stdout)
    if not disp or not yaw:
        raise RuntimeError(f"{family.label} {input_path.name}: no scoring summary")
    guard = GUARD_RE.search(completed.stdout)
    return {
        "family": family.label,
        "record": input_path.name,
        "condition": condition,
        "arm": arm,
        "disp_3d_rms_m": float(disp.group(1)),
        "yaw_rms_deg": float(yaw.group(1)),
        "guard": guard.group(1) if guard else "",
    }


def pooled(rows: list[dict[str, Any]], family: str, condition: str, arm: str,
           key: str) -> float:
    values = [
        row[key] for row in rows
        if row["family"] == family and row["condition"] == condition and row["arm"] == arm
    ]
    if not values:
        raise KeyError(f"no rows for {family} {condition} {arm}")
    return math.sqrt(sum(v * v for v in values) / len(values))


def committed_degradation() -> dict[tuple[str, str], float]:
    """The engine-off and unguarded-cruise cells this study must reproduce."""
    wanted = {"baseline": "engine off", "cruise": "2400 rpm"}
    out: dict[tuple[str, str], float] = {}
    with DEGRADATION_SUMMARY.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            if row["arm"] == "baseline":
                out[(row["family"], wanted["baseline"])] = float(row["disp_3d_rms_m"])
            elif row["arm"] == "speed" and row["rpm"] == f"{float(CRUISE_RPM):.1f}":
                out[(row["family"], wanted["cruise"])] = float(row["disp_3d_rms_m"])
    return out


def main() -> int:
    args = parse_args()
    if not (math.isfinite(args.window_sec) and args.window_sec > 0.0):
        raise SystemExit("--window-sec must be positive and finite")

    missing = [f.label for f in FAMILIES if not f.simulator.exists()]
    if missing:
        raise SystemExit(
            "build the simulators first; missing: " + ", ".join(missing)
        )

    inputs = {name: find_record(args.data_dir, name) for name in RECORDS}
    jobs = [
        (family, inputs[name], condition, arm)
        for family in FAMILIES
        for name in RECORDS
        for condition in CONDITIONS
        for arm in ("off", "guard+R")
    ]
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        rows = list(pool.map(
            lambda job: run_one(*job, args.window_sec), jobs))

    if not args.no_self_check:
        committed = committed_degradation()
        for (family, condition), expected in committed.items():
            got = pooled(rows, family, condition, "off", "disp_3d_rms_m")
            if abs(got - expected) > SELF_CHECK_REL_TOL * max(abs(expected), 1.0):
                raise SystemExit(
                    f"{family} {condition} unguarded: {got:.9g} m against the "
                    f"committed {expected:.9g} m -- this study and the "
                    f"degradation study are not on the same protocol"
                )

    print(f"Pooled over {len(RECORDS)} stationary records, "
          f"{args.window_sec:.0f} s scoring window\n")
    header = (f"{'Family':8s} {'engine off':>11s} {'2400 off':>11s} "
              f"{'2400 guard+R':>13s} {'residual':>9s} {'yaw off':>9s} {'yaw +R':>8s}")
    print(header)
    print("-" * len(header))
    for family in FAMILIES:
        base = pooled(rows, family.label, "engine off", "off", "disp_3d_rms_m")
        cruise_off = pooled(rows, family.label, "2400 rpm", "off", "disp_3d_rms_m")
        cruise_on = pooled(rows, family.label, "2400 rpm", "guard+R", "disp_3d_rms_m")
        yaw_off = pooled(rows, family.label, "2400 rpm", "off", "yaw_rms_deg")
        yaw_on = pooled(rows, family.label, "2400 rpm", "guard+R", "yaw_rms_deg")
        print(f"{family.label:8s} {base:11.3f} {cruise_off:11.3f} "
              f"{cruise_on:13.3f} {cruise_on / base:8.3f}x {yaw_off:9.2f} {yaw_on:8.2f}")

    print("\nEngine-off replays, guarded against unguarded:")
    for family in FAMILIES:
        guarded = pooled(rows, family.label, "engine off", "guard+R", "disp_3d_rms_m")
        base = pooled(rows, family.label, "engine off", "off", "disp_3d_rms_m")
        same = "identical" if guarded == base else f"DIFFER by {guarded - base:.3e} m"
        print(f"  {family.label:8s} {same}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
