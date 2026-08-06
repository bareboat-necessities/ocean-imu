#!/usr/bin/env python3
"""Compare the vertical-acceleration inputs the OU-II / OU-III filters can use.

Two places in the filter consume a vertical acceleration, and each has more than
one signal available.  This tool replays the eight reference records under every
choice, for both filters, and reports the RMS errors and the direction metrics
side by side.

``--axis wave_period`` (default) varies what drives ``WavePeriodEstimator``,
which sets the OU operating point:

  ``leveled``            the heading-frame up component, i.e. the accelerometer
                         rotated by the filter's own attitude solution.  It
                         tracks the sea state well, but attitude is a filter
                         state, so the tuner is inside a loop.  This is what
                         the filter shipped before ``complementary``.

  ``body_z``             the raw body-Z proxy ``-(acc.z + g)`` that the
                         frequency tracker runs on.  It never touches the
                         attitude solution, so the loop is open, but a tilting
                         platform leaks gravity into it below the wave band,
                         where double integration weights the spectrum by
                         1/omega^4.

  ``complementary``      (shipped) levelled with a private Mahony observer
                         running on the raw gyro and accelerometer.  Also a
                         pure function of the measurements, so the loop is
                         equally open, but levelled, so the leakage is removed
                         rather than accepted.

``--axis freq_tracker`` varies what drives the frequency tracker, and with it
the stillness detector that shares its input.  The tracker frequency is the
direction demodulator's carrier, so this axis moves the direction estimates as
well as the displacement RMS:

  ``body_z``             (shipped) the raw proxy.
  ``complementary``      the same levelled signal as above.

Both choices on this axis are measurement-only, so neither closes a loop.
Levelling is not obviously right here the way it is for the period estimator:
the tracker output is not integrated, so it does not suffer the 1/omega^4
weighting that makes sub-band gravity leakage fatal downstream.

The default protocol is the deterministic one of ``tools/ou_sim_table.py``: one
realization per record, default seeds, the final 900 s of a 20-minute replay.
It is not the ten-seed ensemble study and must not be quoted interchangeably
with it.  ``--seeds N`` replicates each record over N seeds and reports a
paired confidence interval alongside the pooled ratio, which is necessary
whenever the effect being measured is smaller than the record-to-record
scatter -- as it is on the ``freq_tracker`` axis.

Usage:
  tools/ou_wave_period_input_study.py
  tools/ou_wave_period_input_study.py --axis freq_tracker --seeds 6
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

# axis -> (env var, inputs, baseline).  The baseline is the denominator of every
# ratio; on the wave_period axis it is deliberately *not* the shipped value,
# because the question that axis answers is what the shipped value changed
# relative to attitude levelling.
AXES = {
    "wave_period": (
        "W3D_WAVE_PERIOD_INPUT",
        ("leveled", "body_z", "complementary"),
        "leveled",
    ),
    "freq_tracker": (
        "W3D_FREQ_TRACKER_INPUT",
        ("body_z", "complementary"),
        "body_z",
    ),
}

# Fields pulled out of the VALIDATION_METRICS line.
FIELDS = (
    "disp_z_rms_m",
    "disp_3d_rms_m",
    "dir_travel_rmse_deg",
    "dir_travel_correct_pct",
    "dir_axis_rmse_deg",
    "dir_axis_circ_std_deg",
    "dir_sense_dominant_pct",
    "disp_x_rms_m",
    "disp_y_rms_m",
    "disp_z_pct_hs",
    "roll_rms_deg",
    "pitch_rms_deg",
    "yaw_rms_deg",
    "dir_travel_error_deg",
    "dir_travel_unresolved_pct",
    "wave_period_s",
    "period_s",
    "tau_applied_s",
)

# Metrics where a ratio is meaningful (errors, lower is better).
RATIO_METRICS = (
    ("disp_z_rms_m", "vertical RMS [m]"),
    ("disp_3d_rms_m", "3D RMS [m]"),
    ("dir_travel_rmse_deg", "travel-sense RMSE [deg]"),
    ("dir_axis_rmse_deg", "axis RMSE [deg]"),
)

# Metrics already bounded, where a ratio would mislead.
PLAIN_METRICS = (
    ("dir_travel_correct_pct", "travel-sense correct [%]"),
    ("dir_axis_circ_std_deg", "axis circular std [deg]"),
)


def run_one(family: str, record: str, env_var: str, value: str,
            seed: int | None = None) -> dict[str, float]:
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
    env[env_var] = value
    if seed is not None:
        env["W3D_SEED"] = str(seed)

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
        raise SystemExit(
            f"no VALIDATION_METRICS for {family} {record} {value} seed={seed}")

    parsed = dict(token.split("=", 1) for token in line.split() if "=" in token)
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


def geometric_mean_ratio(results, family, metric, source, baseline) -> float:
    """Scale-free pooling of relative change across sea states 30x apart in Hs."""
    logs = [
        math.log(results[(family, r, source)][metric]
                 / results[(family, r, baseline)][metric])
        for _, _, r in RECORDS
        if results[(family, r, baseline)][metric] > 0.0
        and results[(family, r, source)][metric] > 0.0
    ]
    return math.exp(sum(logs) / len(logs)) if logs else float("nan")


def paired_log_ratios(raw, family, metric, source, baseline, seeds):
    """Per-(record, seed) log ratios, paired on the realization.

    Pairing matters: the same record replayed under two inputs shares its wave
    realization and its IMU noise, so the difference between them is far better
    resolved than either against the spread across seeds.
    """
    out = []
    for _, _, record in RECORDS:
        for seed in seeds:
            num = raw[(family, record, source, seed)][metric]
            den = raw[(family, record, baseline, seed)][metric]
            if num > 0.0 and den > 0.0:
                out.append(math.log(num / den))
    return out


def summarize_ratio(logs):
    """Geometric mean of a paired ratio with a normal-approximation 95% CI."""
    n = len(logs)
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    mean = sum(logs) / n
    if n < 2:
        return math.exp(mean), float("nan"), float("nan")
    var = sum((x - mean) ** 2 for x in logs) / (n - 1)
    half = 1.96 * math.sqrt(var / n)
    return math.exp(mean), math.exp(mean - half), math.exp(mean + half)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--axis", choices=sorted(AXES), default="wave_period")
    ap.add_argument("--family", choices=sorted(FAMILIES), action="append")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument(
        "--seeds", type=int, default=1,
        help="replay each record under this many seeds and pool the paired "
             "per-seed ratios.  1 (default) keeps the deterministic protocol; "
             "more is needed whenever the effect is smaller than the "
             "record-to-record scatter.")
    ap.add_argument("--csv", type=Path, help="write the raw rows here")
    args = ap.parse_args()

    families = args.family or ["OU_II", "OU_III"]
    env_var, inputs, baseline = AXES[args.axis]
    others = [s for s in inputs if s != baseline]

    # seed=None reproduces the simulator's default seeds exactly, so --seeds 1
    # stays bit-identical to the deterministic protocol.
    seeds = [None] if args.seeds <= 1 else list(range(1, args.seeds + 1))

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {
            (family, record, source, seed):
                pool.submit(run_one, family, record, env_var, source, seed)
            for family in families
            for _, _, record in RECORDS
            for source in inputs
            for seed in seeds
        }
        raw = {key: future.result() for key, future in futures.items()}

    # Per-record display value is the mean over seeds; the pooled ratio below is
    # computed from the paired per-seed ratios instead, which is tighter.
    results = {
        (family, record, source): {
            field: sum(raw[(family, record, source, s)][field] for s in seeds) / len(seeds)
            for field in FIELDS
        }
        for family in families
        for _, _, record in RECORDS
        for source in inputs
    }

    rows = []
    for family in families:
        for metric, label in RATIO_METRICS:
            print()
            print(f"=== {family}: {label}, last {int(WINDOW_SEC)} s, vs {baseline} ===")
            head = f"{'record':<12}{'Hs':>6}{baseline:>14}"
            for source in others:
                head += f"{source:>15}{'delta':>8}"
            print(head)
            for family_label, hs, record in RECORDS:
                base = results[(family, record, baseline)][metric]
                line = f"{family_label:<12}{hs:>6.2f}{base:>14.4f}"
                for source in others:
                    value = results[(family, record, source)][metric]
                    line += f"{value:>15.4f}{pct(value, base):>8}"
                print(line)
            tail = f"{'geometric mean':<18}{'':>14}"
            for source in others:
                ratio = geometric_mean_ratio(results, family, metric, source, baseline)
                tail += f"{ratio:>14.3f}x{'':>8}"
            print(tail)
            if len(seeds) > 1:
                ci = f"{'  paired 95% CI':<18}{'':>14}"
                for source in others:
                    _, lo, hi = summarize_ratio(
                        paired_log_ratios(raw, family, metric, source, baseline, seeds))
                    ci += f"{'[' + f'{lo:.3f}, {hi:.3f}' + ']':>23}"
                print(ci)

        for metric, label in PLAIN_METRICS:
            print()
            print(f"=== {family}: {label} ===")
            head = f"{'record':<12}{'Hs':>6}{baseline:>14}"
            for source in others:
                head += f"{source:>15}"
            print(head)
            for family_label, hs, record in RECORDS:
                line = (f"{family_label:<12}{hs:>6.2f}"
                        f"{results[(family, record, baseline)][metric]:>14.3f}")
                for source in others:
                    line += f"{results[(family, record, source)][metric]:>15.3f}"
                print(line)

        for family_label, hs, record in RECORDS:
            for source in inputs:
                row = dict(results[(family, record, source)])
                row.update(
                    family=family, record=record, sea_state=family_label,
                    hs_m=hs, axis=args.axis, input=source, seeds=len(seeds),
                )
                rows.append(row)

    if args.csv:
        with args.csv.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["family", "sea_state", "hs_m", "record",
                            "axis", "input", "seeds", *FIELDS],
            )
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nwrote {args.csv}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
