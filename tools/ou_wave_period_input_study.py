#!/usr/bin/env python3
"""Compare the possible inputs to the OU-II / OU-III wave-period estimator.

The adaptation tuner takes its operating point from the zero-crossing wave
period, and that period is estimated by ``WavePeriodEstimator`` from a vertical
acceleration.  Which vertical acceleration it gets is a design choice with a
stability consequence:

  ``leveled``            the heading-frame up component, i.e. the accelerometer
                         rotated by the filter's own attitude solution.  It
                         tracks the sea state well, but attitude is a filter
                         state, so the tuner is inside a loop.  This is what
                         the filter shipped before ``complementary``.

  ``body_z``             the raw body-Z proxy ``-(acc.z + g)`` that the
                         frequency tracker already runs on.  It never touches
                         the attitude solution, so the loop is open, but a
                         tilting platform leaks gravity into it below the wave
                         band, where double integration weights the spectrum by
                         1/omega^4.

  ``complementary``      (shipped) levelled with a private Mahony observer
                         running on the raw gyro and accelerometer.  Also a
                         pure function of the measurements, so the loop is
                         equally open, but levelled, so the leakage is removed
                         rather than accepted.

This tool replays the eight reference records under every input, for both
filters, and reports the RMS errors side by side.  The protocol is the
deterministic one of ``tools/ou_sim_table.py``: one realization per record,
default seeds, the final 900 s of a 20-minute replay.  It is not the ten-seed
ensemble study and must not be quoted interchangeably with it.

Usage:
  tools/ou_wave_period_input_study.py
  tools/ou_wave_period_input_study.py --family OU_III --csv out.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# Trailing window the simulators score their own quality gates over.
WINDOW_SEC = 900.0

FAMILIES = {
    "OU_II": {"subdir": "kalman_ou_ii", "binary": "kalman_ou_ii-sim"},
    "OU_III": {"subdir": "kalman_ou_iii", "binary": "kalman_ou_iii-sim"},
}

# The eight records, in the order the historical table reports them.
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

INPUTS = ("leveled", "body_z", "complementary")

# The baseline every ratio is taken against.  This is deliberately not the
# shipped input, which is "complementary": the question the tool exists to
# answer is what the shipped input changed relative to attitude levelling, so
# the old behaviour stays the denominator.
BASELINE = "leveled"

# Fields pulled out of the VALIDATION_METRICS line.
FIELDS = (
    "disp_z_rms_m",
    "disp_3d_rms_m",
    "disp_z_pct_hs",
    "disp_x_rms_m",
    "disp_y_rms_m",
    "roll_rms_deg",
    "pitch_rms_deg",
    "yaw_rms_deg",
    "wave_period_s",
    "period_s",
    "tau_applied_s",
)


def run_one(family: str, record: str, wave_period_input: str) -> dict[str, float]:
    """Replay one record under one input choice and parse its metrics line."""
    spec = FAMILIES[family]
    workdir = (ROOT / "tests" / spec["subdir"]).resolve()
    binary = workdir / spec["binary"]
    if not binary.exists():
        raise SystemExit(f"missing {binary}; run `make -C {workdir} build` first")
    path = workdir / record
    if not path.exists():
        raise SystemExit(f"missing record {path}; run `make fetch-sim-data` first")

    env = dict(os.environ)
    env["W3D_WRITE_TIMESERIES"] = "0"
    # Every record has to be scored, including any that trips a historical gate;
    # this reports what the filter did, not whether it passed.
    env["W3D_COLLECT_ALL_GATES"] = "1"
    env["W3D_VALIDATION_WINDOW_SEC"] = str(WINDOW_SEC)
    env["W3D_WAVE_PERIOD_INPUT"] = wave_period_input

    completed = subprocess.run(
        [str(binary), "--input", str(path)],
        cwd=workdir,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    line = next(
        (l for l in completed.stdout.splitlines() if l.startswith("VALIDATION_METRICS ")),
        None,
    )
    if line is None:
        sys.stderr.write(completed.stdout[-2000:] + completed.stderr[-2000:])
        raise SystemExit(f"no VALIDATION_METRICS for {family} {record} {wave_period_input}")

    parsed = dict(
        token.split("=", 1) for token in line.split() if "=" in token
    )
    out: dict[str, float] = {}
    for field in FIELDS:
        try:
            out[field] = float(parsed[field])
        except (KeyError, ValueError):
            out[field] = float("nan")
    return out


def pct(new: float, old: float) -> str:
    if not (old > 0.0):
        return "   n/a"
    return f"{100.0 * (new / old - 1.0):+6.1f}%"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", choices=sorted(FAMILIES), action="append")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--csv", type=Path, help="write the raw rows here")
    args = ap.parse_args()

    families = args.family or ["OU_II", "OU_III"]

    jobs = [
        (family, family_label, hs, record, source)
        for family in families
        for family_label, hs, record in RECORDS
        for source in INPUTS
    ]

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {
            (family, record, source): pool.submit(run_one, family, record, source)
            for family, _, _, record, source in jobs
        }
        results = {key: future.result() for key, future in futures.items()}

    others = [s for s in INPUTS if s != BASELINE]

    rows = []
    for family in families:
        for metric, label in (("disp_z_rms_m", "vertical"), ("disp_3d_rms_m", "3D")):
            print()
            print(f"=== {family}: last {int(WINDOW_SEC)} s {label} RMS vs {BASELINE} ===")
            head = f"{'record':<12}{'Hs':>6}{BASELINE:>12}"
            for source in others:
                head += f"{source:>14}{'delta':>8}"
            print(head)
            for family_label, hs, record in RECORDS:
                base = results[(family, record, BASELINE)][metric]
                line = f"{family_label:<12}{hs:>6.2f}{base:>12.4f}"
                for source in others:
                    value = results[(family, record, source)][metric]
                    line += f"{value:>14.4f}{pct(value, base):>8}"
                print(line)
            # Geometric mean of the per-record ratio, a scale-free pooling of
            # relative change across sea states that differ by 30x in Hs.
            tail = f"{'geometric mean':<18}{'':>12}"
            for source in others:
                logs = [
                    math.log(results[(family, r, source)][metric]
                             / results[(family, r, BASELINE)][metric])
                    for _, _, r in RECORDS
                    if results[(family, r, BASELINE)][metric] > 0.0
                ]
                ratio = math.exp(sum(logs) / len(logs)) if logs else float("nan")
                tail += f"{ratio:>13.3f}x{'':>8}"
            print(tail)

        print()
        print(f"=== {family}: reported zero-crossing period T_z [s] ===")
        head = f"{'record':<12}{'Hs':>6}{BASELINE:>12}"
        for source in others:
            head += f"{source:>14}"
        print(head)
        for family_label, hs, record in RECORDS:
            line = (f"{family_label:<12}{hs:>6.2f}"
                    f"{results[(family, record, BASELINE)]['wave_period_s']:>12.2f}")
            for source in others:
                line += f"{results[(family, record, source)]['wave_period_s']:>14.2f}"
            print(line)

        for family_label, hs, record in RECORDS:
            for source in INPUTS:
                row = dict(results[(family, record, source)])
                row.update(
                    family=family, record=record, sea_state=family_label,
                    hs_m=hs, wave_period_input=source,
                )
                rows.append(row)

    if args.csv:
        with args.csv.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "family", "sea_state", "hs_m", "record", "wave_period_input", *FIELDS
                ],
            )
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nwrote {args.csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
