#!/usr/bin/env python3
"""Amplitude-exponent ablation of the OU-III integral-regularizer schedule.

Section "Amplitude-exponent ablation of the regularizer schedule" of the
OU-III paper asks a single question: does the drift-band regularizer need its
factor of sigma_aw?  The reduced Riccati model, read with the accelerometer
sensor floor, says no --- the deployed filter is deep in the strong
acceleration-observation branch (zeta >> 1), where the pole-preserving
schedule is r_S ~ tau^(5/2) with no leading-order amplitude dependence.  Read
with a sea-state-proportional model-mismatch term instead, it says yes, and
predicts the deployed r_S ~ sigma_aw tau^(5/2).

The filter exposes exactly that one degree of freedom,

    r_S = sqrt(2 r_a) tau^3 / (sqrt(T_S) kappa^3) * (sigma_aw/sigma_ref)^p,

with kappa calibrated so sigma_ref is the sigma_aw of the nominal Hs = 1.5 m
sea.  The common tau^3/sqrt(T_S) factor cancels between family members, so
p = 1 reproduces the deployed schedule exactly, p = 0 is the strong-branch
asymptote, and every member agrees at sigma_ref for every tau.  The comparison
is therefore of the amplitude exponent alone and cannot be confounded by an
overall change of regularizer gain.

Stationary scenarios still come from tools/ou_validation.py.  The ablation's
non-stationary instrument is deliberately local to this driver so changing it
does not rewrite the committed general-validation evidence protocol: a 120 s
C2 transition starts at 400 s, holds the high sea, and a second 120 s C2
transition starts at 800 s and returns to the *same phase-randomized low-sea
realization*.  The 120 s interval is the same one used by the previous one-way
instrument; only its placement and the added return leg change.

Typical use:

    python3 tools/ou_rs_law_ablation.py --exponents 0,0.5,1,1.25 --jobs 4

Add --posterior to include the full transition law
Eq. (riccati-pole-target-rs) as an extra column; it is expected to be
numerically indistinguishable from p = 0 wherever zeta >> 1.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import statistics as st
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Sequence

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
VALIDATION = REPO_ROOT / "tools" / "ou_validation.py"
OU_III_DATA = REPO_ROOT / "tests" / "kalman_ou_iii"

# Metrics reported by the ablation table.
METRICS = ("disp_z_pct_hs", "disp_3d_rms_m")

# Law selector understood by tests/kalman_ou_iii/kalman_ou_iii-sim.cpp.
LAW_CUBIC = "0"
LAW_STRONG = "1"
LAW_POSTERIOR = "2"

# The former transition occupied 540--660 s: 120 s.  Keep that interval
# unchanged, move the low->high onset to 400 s, then make the same 120 s
# transition back to the original low realization from 800 s.
TRANSITION_DURATION_SEC = 120.0
DEFAULT_TRANSITION_START_SEC = 400.0
DEFAULT_TRANSITION_RETURN_START_SEC = 800.0
ROUNDTRIP_DURATION_SEC = 1200.0
ROUNDTRIP_WINDOW_SEC = 900.0
ROUNDTRIP_SCENARIO = "nonstationary_H1_5_to_H4_0_to_H1_5_Tp5_7_to_11_4"
# Name the ablation already gives its generated round-trip surrogate; the
# diagnostic figure reuses it so the plotted realization is the scored one.
ROUNDTRIP_WAVE_NAME = "wave_data_jonswap_H4.000_L202.839_A30.00_P120.00.csv"

# Endpoint heights of the round trip.  The high sea is the H_s = 8.5 m source
# record scaled to 4.0 m, the same endpoint the one-way instrument uses.
DIAGNOSTIC_LOW_HEIGHT_M = 1.5
DIAGNOSTIC_HIGH_HEIGHT_M = 4.0
DIAGNOSTIC_SOURCE_HEIGHT_M = 8.5


def _validation_module():
    tools = str(REPO_ROOT / "tools")
    if tools not in sys.path:
        sys.path.insert(0, tools)
    import ou_validation as ov
    return ov


def transition_bounds(
    start_sec: float = DEFAULT_TRANSITION_START_SEC,
    return_start_sec: float = DEFAULT_TRANSITION_RETURN_START_SEC,
) -> tuple[float, float, float, float]:
    """Return low->high and high->low bounds with the fixed 120 s duration."""
    end_sec = start_sec + TRANSITION_DURATION_SEC
    return_end_sec = return_start_sec + TRANSITION_DURATION_SEC
    if not (0.0 <= start_sec < end_sec <= return_start_sec < return_end_sec <= ROUNDTRIP_DURATION_SEC):
        raise ValueError(
            "round-trip transition requires 0 <= rise < high hold < fall <= record"
        )
    return start_sec, end_sec, return_start_sec, return_end_sec


def roundtrip_profile(
    ov,
    times: np.ndarray,
    start_sec: float = DEFAULT_TRANSITION_START_SEC,
    return_start_sec: float = DEFAULT_TRANSITION_RETURN_START_SEC,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """High-sea weight and its first two derivatives for one round trip."""
    start, end, return_start, return_end = transition_bounds(
        start_sec, return_start_sec
    )
    up, up_rate, up_accel = ov.smoothstep_profile(times, start, end)
    down, down_rate, down_accel = ov.smoothstep_profile(
        times, return_start, return_end
    )
    return up - down, up_rate - down_rate, up_accel - down_accel


def make_roundtrip_wave(
    ov,
    columns: Sequence[str],
    low_data: np.ndarray,
    high_data: np.ndarray,
    seed: int,
    high_scale: float,
    start_sec: float = DEFAULT_TRANSITION_START_SEC,
    return_start_sec: float = DEFAULT_TRANSITION_RETURN_START_SEC,
) -> np.ndarray:
    """Build a kinematically closed low->high->same-low realization.

    The same randomized ``low`` array is used on both sides of the high-sea
    excursion.  The quintic weight's exact derivative cross terms are retained
    in v and a, matching ou_validation.make_nonstationary_wave rather than
    independently crossfading p/v/a and breaking their derivative relations.
    """
    count = min(low_data.shape[0], high_data.shape[0])
    low = ov.phase_randomize_wave(columns, low_data[:count], seed * 2 + 1)
    high = ov.phase_randomize_wave(columns, high_data[:count], seed * 2 + 2)
    high = ov.scale_wave_motion(columns, high, high_scale)
    times = low[:, columns.index("time")]
    weight, weight_rate, weight_acceleration = roundtrip_profile(
        ov, times, start_sec, return_start_sec
    )

    def indices(names: Sequence[str]) -> list[int]:
        return [columns.index(name) for name in names]

    displacement = indices(("disp_x", "disp_y", "disp_z"))
    velocity = indices(("vel_x", "vel_y", "vel_z"))
    acceleration = indices(("acc_x", "acc_y", "acc_z"))
    attitude = indices(("roll_deg", "pitch_deg", "yaw_deg"))

    result = low.copy()
    w = weight[:, None]
    dw = weight_rate[:, None]
    ddw = weight_acceleration[:, None]
    displacement_delta = high[:, displacement] - low[:, displacement]
    velocity_delta = high[:, velocity] - low[:, velocity]

    result[:, displacement] = (
        (1.0 - w) * low[:, displacement] + w * high[:, displacement]
    )
    result[:, velocity] = (
        (1.0 - w) * low[:, velocity]
        + w * high[:, velocity]
        + dw * displacement_delta
    )
    result[:, acceleration] = (
        (1.0 - w) * low[:, acceleration]
        + w * high[:, acceleration]
        + 2.0 * dw * velocity_delta
        + ddw * displacement_delta
    )
    result[:, attitude] = (
        (1.0 - w) * low[:, attitude] + w * high[:, attitude]
    )
    return ov.rebuild_body_imu(columns, result)


def roundtrip_segments(
    start_sec: float = DEFAULT_TRANSITION_START_SEC,
    return_start_sec: float = DEFAULT_TRANSITION_RETURN_START_SEC,
) -> tuple[tuple[str, float, float], ...]:
    """Expose both transitions, their run-on, and both settled sea states."""
    start, end, return_start, return_end = transition_bounds(
        start_sec, return_start_sec
    )
    window_start = ROUNDTRIP_DURATION_SEC - ROUNDTRIP_WINDOW_SEC
    high_recover_end = end + TRANSITION_DURATION_SEC
    low_recover_end = return_end + TRANSITION_DURATION_SEC
    return (
        ("low_start", window_start, start),
        ("rise", start, end),
        ("high_recover", end, high_recover_end),
        ("high", high_recover_end, return_start),
        ("fall", return_start, return_end),
        ("low_recover", return_end, low_recover_end),
        ("low_return", low_recover_end, ROUNDTRIP_DURATION_SEC),
    )


def run_one(label: str, env_extra: dict[str, str], out_dir: Path,
            jobs: int, extra_args: list[str]) -> Path:
    """Run only the stationary validation scenarios for one law setting."""
    out_dir.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(env_extra)
    # The ablation supplies its own round-trip non-stationary scenario below.
    # Keeping the general validation driver's one-way transition out of this
    # subprocess prevents the old and new instruments from being pooled.
    cmd = [
        sys.executable, str(VALIDATION),
        "--mode", "full",
        "--families", "OU_III",
        "--adaptation-modes", "Adaptive",
        "--no-plots",
        "--skip-nonstationary",
        "--jobs", str(jobs),
        "--output-dir", str(out_dir),
        *extra_args,
    ]
    print(f"[{label}] {' '.join(f'{k}={v}' for k, v in sorted(env_extra.items()))}",
          flush=True)
    log = out_dir.with_suffix(".log")
    with log.open("w") as fh:
        proc = subprocess.run(cmd, env=env, stdout=fh, stderr=subprocess.STDOUT,
                              cwd=REPO_ROOT, check=False)
    raw = out_dir / "ou_validation_raw.csv"
    if not raw.exists():
        raise SystemExit(
            f"[{label}] produced no raw rows (exit {proc.returncode}); see {log}")
    # A nonzero exit only means the publication table could not be built,
    # which needs paired FixedNominal/FixedOracle modes this study does not
    # run.  The per-run rows this study consumes are already written.
    return raw


def load(raw: Path) -> dict[tuple[str, int], dict[str, str]]:
    with raw.open() as fh:
        return {(r["scenario"], int(r["repetition"])): r
                for r in csv.DictReader(fh)}


def stationary_rows(rows: dict[tuple[str, int], dict[str, str]]) -> dict:
    """Discard an old cached one-way transition when --reuse reads legacy raw data."""
    return {
        key: row for key, row in rows.items()
        if str(row.get("scenario_kind", "stationary")) == "stationary"
    }


def seed_triplets(rows: dict[tuple[str, int], dict[str, str]]) -> list[tuple[int, int, int, int]]:
    """Recover exactly the paired seeds selected by the validation invocation."""
    by_repetition: dict[int, tuple[int, int, int, int]] = {}
    for (_, repetition), row in sorted(rows.items()):
        by_repetition.setdefault(
            repetition,
            (
                repetition,
                int(row["wave_phase_seed"]),
                int(row["imu_noise_seed"]),
                int(row["initialization_seed"]),
            ),
        )
    if not by_repetition:
        raise SystemExit("stationary validation produced no seed triplets")
    return list(by_repetition.values())


def _write_rows(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run_roundtrip(
    label: str,
    env_extra: dict[str, str],
    stationary: dict[tuple[str, int], dict[str, str]],
    out_dir: Path,
    reuse: bool,
    start_sec: float,
    return_start_sec: float,
) -> dict[tuple[str, int], dict[str, str]]:
    """Run the ablation-only low->high->same-low scenario on matched seeds."""
    raw = out_dir / "roundtrip_raw.csv"
    if reuse and raw.exists():
        return load(raw)

    ov = _validation_module()
    low_path = ov.find_default_input(OU_III_DATA, "1.500", "50.710")
    high_path = ov.find_default_input(OU_III_DATA, "8.500", "202.839")
    columns, low_data = ov.read_wave_csv(low_path, ROUNDTRIP_DURATION_SEC)
    _, high_data = ov.read_wave_csv(high_path, ROUNDTRIP_DURATION_SEC)
    high_scale = 4.0 / 8.5
    segments = roundtrip_segments(start_sec, return_start_sec)

    old_env = {key: os.environ.get(key) for key in env_extra}
    os.environ.update(env_extra)
    rows: list[dict[str, Any]] = []
    try:
        with tempfile.TemporaryDirectory(prefix="ou_rs_roundtrip_") as tmp:
            root = Path(tmp)
            for repetition, wave_seed, imu_seed, init_seed in seed_triplets(stationary):
                generated = make_roundtrip_wave(
                    ov, columns, low_data, high_data, wave_seed, high_scale,
                    start_sec, return_start_sec,
                )
                sub = root / f"wave_{wave_seed}"
                path = sub / ROUNDTRIP_WAVE_NAME
                ov.write_wave_csv(path, columns, generated)
                metrics, gate_pass, returncode = ov.run_simulator(
                    "OU_III",
                    path,
                    ROUNDTRIP_WINDOW_SEC,
                    imu_seed=imu_seed,
                    initialization_seed=init_seed,
                    tuning_mode="adaptive",
                    aw_cov_sync="periodic",
                    segments=segments,
                )
                rows.append(
                    {
                        "scenario": ROUNDTRIP_SCENARIO,
                        "scenario_kind": "nonstationary",
                        "spectrum": "JONSWAP",
                        "family": "OU_III",
                        "mode": "Adaptive",
                        "repetition": repetition,
                        "wave_phase_seed": wave_seed,
                        "imu_noise_seed": imu_seed,
                        "initialization_seed": init_seed,
                        "score_window_sec": ROUNDTRIP_WINDOW_SEC,
                        "quality_gate_pass": int(gate_pass),
                        "simulator_returncode": returncode,
                        **metrics,
                    }
                )
    finally:
        for key, value in old_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value

    _write_rows(raw, rows)
    print(f"[{label}] round-trip transition rows -> {raw}", flush=True)
    return load(raw)


def write_roundtrip_diagnostic(
    svg_path: Path,
    csv_path: Path,
    wave_seed: int,
    imu_seed: int,
    init_seed: int,
    start_sec: float = DEFAULT_TRANSITION_START_SEC,
    return_start_sec: float = DEFAULT_TRANSITION_RETURN_START_SEC,
    decimation: int = 20,
) -> dict[str, Any]:
    """Plot one round-trip realization the way the one-way figure is plotted.

    The scored tables cannot show whether the tuner follows the sea in both
    directions or only up, so the ablation gets the same four-panel view the
    general-validation transition figure uses: the sea state the record
    actually carries, the wave and its estimate, the vertical error, and the
    applied tuning against the two frozen operating points.  Only the protocol
    differs -- two crossfades and a return to the same low realization.
    """

    ov = _validation_module()
    low_path = ov.find_default_input(OU_III_DATA, "1.500", "50.710")
    high_path = ov.find_default_input(OU_III_DATA, "8.500", "202.839")
    columns, low_data = ov.read_wave_csv(low_path, ROUNDTRIP_DURATION_SEC)
    _, high_data = ov.read_wave_csv(high_path, ROUNDTRIP_DURATION_SEC)
    high_scale = DIAGNOSTIC_HIGH_HEIGHT_M / DIAGNOSTIC_SOURCE_HEIGHT_M

    _, _, _, return_end_sec = transition_bounds(start_sec, return_start_sec)
    window_start_sec = ROUNDTRIP_DURATION_SEC - ROUNDTRIP_WINDOW_SEC

    with tempfile.TemporaryDirectory(prefix="ou_rs_diagnostic_") as tmp:
        root = Path(tmp)

        # The two frozen points are the reference lines of the bottom panel:
        # each is calibrated once from its own noise-free record, exactly as
        # the general validation calibrates its fixed-tuning modes.
        low_reference = root / "calibration" / "low" / low_path.name
        ov.write_wave_csv(low_reference, columns, low_data)
        high_reference = root / "calibration" / "high" / ROUNDTRIP_WAVE_NAME
        ov.write_wave_csv(
            high_reference,
            columns,
            ov.scale_wave_motion(columns, high_data, high_scale),
        )
        low_point = ov.calibrate_tuning_point(
            "OU_III", low_reference, ROUNDTRIP_WINDOW_SEC
        )
        high_point = ov.calibrate_tuning_point(
            "OU_III", high_reference, ROUNDTRIP_WINDOW_SEC
        )

        generated = make_roundtrip_wave(
            ov, columns, low_data, high_data, wave_seed, high_scale,
            start_sec, return_start_sec,
        )
        surrogate = root / f"wave_{wave_seed}" / ROUNDTRIP_WAVE_NAME
        ov.write_wave_csv(surrogate, columns, generated)
        metrics, _, _ = ov.run_simulator(
            "OU_III",
            surrogate,
            ROUNDTRIP_WINDOW_SEC,
            imu_seed=imu_seed,
            initialization_seed=init_seed,
            tuning_mode="adaptive",
            aw_cov_sync="periodic",
            segments=roundtrip_segments(start_sec, return_start_sec),
            write_timeseries=True,
        )
        series = ov.read_diagnostic_timeseries("OU_III", surrogate)

    time = series["time"]
    weight, _, _ = roundtrip_profile(ov, time, start_sec, return_start_sec)
    reference = series["disp_ref_z"]
    estimate = series["disp_est_z"]
    error = estimate - reference

    rolling_hs = ov.rolling_significant_height(reference)
    mixture_hs = np.asarray(
        [
            ov.mixture_significant_height_m(
                DIAGNOSTIC_LOW_HEIGHT_M, DIAGNOSTIC_HIGH_HEIGHT_M, value
            )
            for value in weight
        ]
    )
    linear_hs = (
        (1.0 - weight) * DIAGNOSTIC_LOW_HEIGHT_M
        + weight * DIAGNOSTIC_HIGH_HEIGHT_M
    )

    step = max(1, int(decimation))
    _write_rows(
        csv_path,
        [
            {
                "time_s": float(time[index]),
                "blend_weight": float(weight[index]),
                "rolling_hs_m": float(rolling_hs[index]),
                "mixture_hs_m": float(mixture_hs[index]),
                "disp_ref_z_m": float(reference[index]),
                "disp_est_z_m": float(estimate[index]),
                "disp_err_z_m": float(error[index]),
                "tau_applied_s": float(series["tau_applied"][index]),
                "sigma_aw_applied_mps2": float(series["sigma_a_applied"][index]),
                "r_s_applied_ms": float(series["R_p0_applied"][index]),
                "freq_tracker_hz": float(series["freq_tracker_hz"][index]),
            }
            for index in range(0, time.size, step)
        ],
    )

    sliced = slice(None, None, step)
    svg_path.parent.mkdir(parents=True, exist_ok=True)
    ov.plot_transition_diagnostic(
        svg_path,
        time[sliced],
        reference[sliced],
        estimate[sliced],
        rolling_hs[sliced],
        mixture_hs[sliced],
        linear_hs[sliced],
        series["tau_applied"][sliced],
        series["sigma_a_applied"][sliced],
        series["R_p0_applied"][sliced],
        (
            (r"low-sea fixed point ($H_s=1.5$ m)", low_point, ":"),
            (r"high-sea fixed point ($H_s=4.0$ m)", high_point, "--"),
        ),
        (
            (start_sec, start_sec + TRANSITION_DURATION_SEC),
            (return_start_sec, return_end_sec),
        ),
        window_start_sec,
    )
    print(f"round-trip transition diagnostic -> {svg_path}", flush=True)

    return {
        "family": "OU_III",
        "wave_phase_seed": wave_seed,
        "imu_noise_seed": imu_seed,
        "initialization_seed": init_seed,
        "decimation": step,
        "disp_z_pct_hs": float(metrics["disp_z_pct_hs"]),
        "low_tau_s": low_point.tau_s,
        "low_sigma_aw_mps2": low_point.sigma_a_mps2,
        "low_r_s_ms": low_point.RS_ms,
        "high_tau_s": high_point.tau_s,
        "high_sigma_aw_mps2": high_point.sigma_a_mps2,
        "high_r_s_ms": high_point.RS_ms,
    }


def paired_delta(runs: dict[str, dict], keys: list, label: str,
                 base: str, metric: str) -> tuple[float, float, float]:
    """Return (mean of label, mean paired delta vs base, 95% half-width)."""
    vals = [float(runs[label][k][metric]) for k in keys]
    diffs = [float(runs[label][k][metric]) - float(runs[base][k][metric])
             for k in keys]
    mean = st.mean(vals)
    dmean = st.mean(diffs)
    half = (1.96 * st.stdev(diffs) / math.sqrt(len(diffs))
            if len(diffs) > 1 and st.stdev(diffs) > 0 else 0.0)
    return mean, dmean, half


def print_roundtrip_segments(runs: dict[str, dict], order: Sequence[str]) -> None:
    """Show rise/fall separately so directional adaptation lag cannot cancel."""
    names = [name for name, _, _ in roundtrip_segments()]
    print("roundtrip transition disp_z_pct_refrms")
    print(f"  {'schedule':12} " + " ".join(f"{name:>13}" for name in names))
    for label in order:
        rows = [
            row for (scenario, _), row in runs[label].items()
            if scenario == ROUNDTRIP_SCENARIO
        ]
        cells = []
        for name in names:
            metric = f"seg_{name}_disp_z_pct_refrms"
            values = [float(row[metric]) for row in rows if metric in row]
            cells.append(st.mean(values) if values else math.nan)
        print(f"  {label:12} " + " ".join(f"{value:13.3f}" for value in cells))
    print()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exponents", default="0,0.5,1,1.25",
                    help="comma-separated amplitude exponents p (default 0,0.5,1,1.25)")
    ap.add_argument("--baseline", default="1",
                    help="exponent treated as the deployed baseline (default 1)")
    ap.add_argument("--posterior", action="store_true",
                    help="also run the full posterior transition law")
    ap.add_argument("--transition-start-sec", type=float,
                    default=DEFAULT_TRANSITION_START_SEC,
                    help="low->high transition start (default 400 s)")
    ap.add_argument("--transition-return-start-sec", type=float,
                    default=DEFAULT_TRANSITION_RETURN_START_SEC,
                    help="high->same-low transition start (default 800 s)")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--output-dir", type=Path,
                    default=REPO_ROOT / "reports" / "results" / "ou_rs_law")
    ap.add_argument("--reuse", action="store_true",
                    help="reuse existing stationary and round-trip raw CSVs")
    ap.add_argument("--diagnostic-seeds", default="11,101,1009",
                    help="wave,imu,init seeds of the plotted round-trip "
                         "realization (default 11,101,1009, the first triplet "
                         "of the full validation ensemble)")
    ap.add_argument("--diagnostic-only", action="store_true",
                    help="write only the round-trip transition figure and exit")
    ap.add_argument("rest", nargs=argparse.REMAINDER,
                    help="extra args forwarded to ou_validation.py")
    args = ap.parse_args(argv)

    # Validate the requested placement before any expensive simulator work.
    transition_bounds(args.transition_start_sec, args.transition_return_start_sec)

    diagnostic_seeds = [int(x) for x in args.diagnostic_seeds.split(",") if x.strip()]
    if len(diagnostic_seeds) != 3:
        raise SystemExit("--diagnostic-seeds takes wave,imu,init")

    def diagnostic() -> None:
        write_roundtrip_diagnostic(
            args.output_dir / "ou_rs_roundtrip_transition.svg",
            args.output_dir / "ou_rs_roundtrip_transition.csv",
            *diagnostic_seeds,
            start_sec=args.transition_start_sec,
            return_start_sec=args.transition_return_start_sec,
        )

    if args.diagnostic_only:
        diagnostic()
        return 0

    extra = [a for a in args.rest if a != "--"]
    exponents = [e.strip() for e in args.exponents.split(",") if e.strip()]
    if args.baseline not in exponents:
        raise SystemExit(f"baseline p={args.baseline} must be among --exponents")

    runs: dict[str, dict] = {}
    order: list[str] = []
    for p in exponents:
        label = f"p={p}"
        order.append(label)
        out = args.output_dir / f"p_{p}"
        raw = out / "ou_validation_raw.csv"
        env_extra = {"OU_III_RS_LAW": LAW_STRONG, "OU_III_RS_SIGMA_EXP": p}
        # p == 1 is reproduced exactly by the Riccati path, so the whole family
        # is run through the same code path for a like-for-like comparison.
        if not (args.reuse and raw.exists()):
            raw = run_one(label, env_extra, out, args.jobs, extra)
        stationary = stationary_rows(load(raw))
        roundtrip = run_roundtrip(
            label, env_extra, stationary, out, args.reuse,
            args.transition_start_sec, args.transition_return_start_sec,
        )
        runs[label] = {**stationary, **roundtrip}

    if args.posterior:
        label = "posterior"
        order.append(label)
        out = args.output_dir / "posterior"
        raw = out / "ou_validation_raw.csv"
        env_extra = {"OU_III_RS_LAW": LAW_POSTERIOR}
        if not (args.reuse and raw.exists()):
            raw = run_one(label, env_extra, out, args.jobs, extra)
        stationary = stationary_rows(load(raw))
        roundtrip = run_roundtrip(
            label, env_extra, stationary, out, args.reuse,
            args.transition_start_sec, args.transition_return_start_sec,
        )
        runs[label] = {**stationary, **roundtrip}

    base = f"p={args.baseline}"
    keys = sorted(set.intersection(*(set(r) for r in runs.values())))
    if not keys:
        raise SystemExit("no (scenario, repetition) units are common to all runs")

    print(f"\nPaired amplitude-exponent ablation, n={len(keys)} units, "
          f"baseline {base}\n")
    width = max(len(x) for x in order) + 2
    for metric in METRICS:
        print(f"{metric}")
        print(f"  {'schedule':{width}} {'mean':>9} {'delta':>9} {'95% CI':>18}")
        for label in order:
            mean, dmean, half = paired_delta(runs, keys, label, base, metric)
            if label == base:
                print(f"  {label:{width}} {mean:9.4f} {'---':>9} {'':>18}")
            else:
                print(f"  {label:{width}} {mean:9.4f} {dmean:+9.4f}"
                      f"   [{dmean - half:+7.4f},{dmean + half:+7.4f}]")
        print()

    print("per-scenario disp_z_pct_hs")
    scenarios = sorted({k[0] for k in keys})
    print(f"  {'scenario':54} " + " ".join(f"{x:>12}" for x in order))
    for scenario in scenarios:
        sub = [k for k in keys if k[0] == scenario]
        cells = " ".join(
            f"{st.mean([float(runs[x][k]['disp_z_pct_hs']) for k in sub]):12.3f}"
            for x in order)
        print(f"  {scenario[:54]:54} {cells}")
    print()

    print_roundtrip_segments(runs, order)
    diagnostic()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
