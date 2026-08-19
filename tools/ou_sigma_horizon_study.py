#!/usr/bin/env python3
"""Measure what the sigma_a averaging horizon is worth on the wave-band
operating point.

Both OU families take the whole operating point -- the sigma band's corners,
the sigma_a averaging horizon, and tau -- from the wave-band zero-crossing
period of ``WavePeriodEstimator``, adopting the estimator's own value as soon
as it exists and standing a fixed wave-band prior in before that.  The
acceleration-band frequency tracker has no path into the adaptation, so the
``source`` axis this harness used to sweep (``W3D_TUNER_FREQ_SOURCE``, with its
``tracker_fallback`` and ``wave_band_gated`` arms) is gone with it.

``--axis k_periods`` varies the sigma_a averaging horizon
(``OU_SIGMA_VAR_K_PERIODS``), which ``SeaStateAutoTuner`` expresses in periods
of whatever tuning frequency it is fed: tau_var = K * T_eff, and the variance is
a two-stage EWMA at that horizon, so the memory is about 2*K*T_eff.  K has not
moved since the operating point did, and the move rescaled it: fed the
acceleration band, T_eff was ~2 s whatever the sea was doing, so K = 2 meant a
~4 s horizon; fed the wave band it is 2*T_z, i.e. 5 s on the shortest reference
sea and 17 s on the longest.

``--axis prior`` varies the constant that stands in before the period estimator
has a value (``OU_TUNE_FREQ_PRIOR_HZ``), which is what makes the wave-band
source exogenous.  It is a sensitivity check, not a tuning surface.

``--axis horizon_max`` varies the absolute ceiling on the horizon
(``OU_SIGMA_VAR_HORIZON_MAX_S``), which is what would bind first if K grew.

Scenarios
---------
The eight released stationary records, plus -- unless ``--skip-transition`` --
the controlled non-stationary record of ``tools/ou_validation.py``: a C2 quintic
crossfade from H_s = 1.5 m, T_p = 5.7 s to H_s = 4.0 m, T_p = 11.4 s between 540
and 660 s.  That record is scored four more times, over its start, blend,
recover and end intervals, so a horizon that is merely slow can be told apart
from one that is wrong: a longer sigma_a memory should cost in ``blend``, where
it lags the moving sea, and again in ``recover``, where it is still carrying the
start sea into the new one.

Protocol
--------
``--seeds N`` replays every scenario under N IMU/initialization seeds and pools
the *paired* per-seed log ratios, which is what resolves effects smaller than
the record-to-record scatter.  ``--seeds 1`` keeps the deterministic
single-realization protocol of ``tools/ou_sim_table.py`` and is bit-identical to
it.  Every ratio is against the axis baseline on the identical realization.

Usage:
  tools/ou_sigma_horizon_study.py --seeds 6
  tools/ou_sigma_horizon_study.py --axis k_periods --seeds 6 --csv k.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))

import ou_validation as ouv  # noqa: E402  (path set above)

# Trailing window the simulators score their own quality gates over.
WINDOW_SEC = 900.0
DURATION_SEC = 1200.0
# Crossfade of the controlled transition record, as tools/ou_validation.py
# builds it: 120 s centred on the 1200 s replay.
TRANSITION_START_SEC = 540.0
TRANSITION_END_SEC = 660.0
# The endpoint sea is scored twice, split one crossfade length past the blend.
# "recover" is where a schedule whose averaging horizon is too long is still
# carrying the start sea; "end" is the settled endpoint sea.
TRANSITION_RECOVER_END_SEC = TRANSITION_END_SEC + (
    TRANSITION_END_SEC - TRANSITION_START_SEC
)
TRANSITION_END_HEIGHT_M = 4.0

FAMILIES = {
    "OU_II": {"subdir": "kalman_ou_ii", "binary": "kalman_ou_ii-sim"},
    "OU_III": {"subdir": "kalman_ou_iii", "binary": "kalman_ou_iii-sim"},
}

# The eight released records, in the order the historical table reports them.
STATIONARY = (
    ("JONSWAP", 0.27, "wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv"),
    ("JONSWAP", 1.50, "wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv"),
    ("JONSWAP", 4.00, "wave_data_jonswap_H4.000_L112.766_A30.00_P30.00.csv"),
    ("JONSWAP", 8.50, "wave_data_jonswap_H8.500_L202.839_A-30.00_P72.00.csv"),
    ("PM-Stokes", 0.27, "wave_data_pmstokes_H0.270_L14.047_A30.00_P60.00.csv"),
    ("PM-Stokes", 1.50, "wave_data_pmstokes_H1.500_L50.710_A-30.00_P120.00.csv"),
    ("PM-Stokes", 4.00, "wave_data_pmstokes_H4.000_L112.766_A30.00_P30.00.csv"),
    ("PM-Stokes", 8.50, "wave_data_pmstokes_H8.500_L202.839_A-30.00_P72.00.csv"),
)

TRANSITION_LABEL = "transition"
TRANSITION_KEY = "nonstationary_H1_5_to_H4_0"

# Segments of the transition record, keyed to the same 900 s window.
SEGMENTS = (
    ("start", DURATION_SEC - WINDOW_SEC, TRANSITION_START_SEC),
    ("blend", TRANSITION_START_SEC, TRANSITION_END_SEC),
    ("recover", TRANSITION_END_SEC, TRANSITION_RECOVER_END_SEC),
    ("end", TRANSITION_RECOVER_END_SEC, DURATION_SEC),
)

# axis -> (env var, values, baseline).  The baseline is the denominator of every
# ratio.
AXES = {
    "k_periods": (
        "OU_SIGMA_VAR_K_PERIODS",
        ("0.5", "1", "2", "3", "4", "6", "8"),
        "2",
    ),
    "prior": (
        "OU_TUNE_FREQ_PRIOR_HZ",
        ("0.1", "0.14", "0.2", "0.28", "0.4"),
        "0.2",
    ),
    "horizon_max": (
        "OU_SIGMA_VAR_HORIZON_MAX_S",
        ("10", "20", "30", "60", "120"),
        "60",
    ),
}

SEGMENT_FIELDS = tuple(
    f"seg_{name}_{metric}"
    for name, _, _ in SEGMENTS
    for metric in ("disp_z_rms_m", "disp_3d_rms_m", "disp_z_pct_refrms")
)

FIELDS = (
    "disp_z_rms_m",
    "disp_3d_rms_m",
    "disp_x_rms_m",
    "disp_y_rms_m",
    "disp_z_pct_hs",
    "disp_z_pct_refrms",
    "roll_rms_deg",
    "pitch_rms_deg",
    "yaw_rms_deg",
    "dir_travel_rmse_deg",
    "dir_axis_rmse_deg",
    "dir_travel_correct_pct",
    "wave_period_s",
    "period_s",
    "tau_applied_s",
    "sigma_applied_mps2",
    "accel_variance_m2ps4",
    *SEGMENT_FIELDS,
)

# Metrics where a ratio is meaningful (errors, lower is better).  The attitude
# channels are here because two of the simulators' quality-gate sentinels are
# roll and pitch, and they are fitted with about half a percent of margin: a
# change can be invisible in displacement and still move a gate.
RATIO_METRICS = (
    ("disp_z_rms_m", "vertical RMS [m]"),
    ("disp_3d_rms_m", "3D RMS [m]"),
    ("roll_rms_deg", "roll RMS [deg]"),
    ("pitch_rms_deg", "pitch RMS [deg]"),
    ("yaw_rms_deg", "yaw RMS [deg]"),
)

TRANSITION_RATIO_METRICS = (
    ("seg_start_disp_z_rms_m", "transition start-sea vertical RMS [m]"),
    ("seg_blend_disp_z_rms_m", "transition blend vertical RMS [m]"),
    ("seg_recover_disp_z_rms_m", "transition run-on vertical RMS [m]"),
    ("seg_end_disp_z_rms_m", "transition settled-sea vertical RMS [m]"),
)


def transition_record(wave_seed: int, data_dir: Path, tmpdir: Path) -> Path:
    """Synthesize the crossfade record for one wave-phase seed."""
    start_path = ouv.find_default_input(data_dir, "1.500", "50.710")
    end_path = ouv.find_default_input(data_dir, "8.500", "202.839")
    columns, start_data = ouv.read_wave_csv(start_path, DURATION_SEC)
    _, end_data = ouv.read_wave_csv(end_path, DURATION_SEC)
    generated = ouv.make_nonstationary_wave(
        columns,
        start_data,
        end_data,
        wave_seed,
        TRANSITION_END_HEIGHT_M / 8.5,
        TRANSITION_START_SEC,
        TRANSITION_END_SEC,
    )
    # The simulator reads H_s, wavelength and azimuth out of the filename, so
    # the surrogate has to carry the endpoint sea's name.  This is the same
    # name tools/ou_validation.py gives it.
    path = (tmpdir / f"wave{wave_seed}" /
            "wave_data_jonswap_H4.000_L202.839_A-30.00_P120.00.csv")
    path.parent.mkdir(parents=True, exist_ok=True)
    ouv.write_wave_csv(path, columns, generated)
    return path


def run_one(family: str, record_path: Path, env_var: str, value: str,
            seed: int | None, segments: bool) -> dict[str, float]:
    """Replay one record under one setting and parse its metrics."""
    spec = FAMILIES[family]
    workdir = (ROOT / "tests" / spec["subdir"]).resolve()
    binary = workdir / spec["binary"]
    if not binary.exists():
        raise SystemExit(f"missing {binary}; run `make -C {workdir} build` first")
    if not record_path.exists():
        raise SystemExit(f"missing record {record_path}; run `make fetch-sim-data` first")

    env = dict(os.environ)
    env["W3D_WRITE_TIMESERIES"] = "0"
    # Every record has to be scored, including any that trips a historical gate;
    # this reports what the filter did, not whether it passed.
    env["W3D_COLLECT_ALL_GATES"] = "1"
    env["W3D_VALIDATION_WINDOW_SEC"] = str(WINDOW_SEC)
    env[env_var] = value
    if segments:
        env["W3D_VALIDATION_SEGMENTS"] = ",".join(
            f"{name}:{lo:.9g}:{hi:.9g}" for name, lo, hi in SEGMENTS
        )
    if seed is not None:
        env["W3D_IMU_SEED"] = str(seed)
        env["W3D_INIT_SEED"] = str(1000 + seed)

    completed = subprocess.run(
        [str(binary), "--input", str(record_path.resolve())],
        cwd=workdir,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    try:
        parsed = ouv.parse_validation_metrics(completed.stdout)
    except ValueError:
        sys.stderr.write(completed.stdout[-2000:] + completed.stderr[-2000:])
        raise SystemExit(
            f"no VALIDATION_METRICS for {family} {record_path.name} "
            f"{env_var}={value} seed={seed}")

    return {field: float(parsed.get(field, float("nan"))) for field in FIELDS}


def pct(new: float, old: float) -> str:
    if not (old > 0.0) or not (new == new):
        return "   n/a"
    return f"{100.0 * (new / old - 1.0):+6.2f}%"


def summarize_ratio(logs: list[float]) -> tuple[float, float, float]:
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


def paired_logs(raw, family, metric, value, baseline, keys, seeds) -> list[float]:
    out = []
    for key in keys:
        for seed in seeds:
            num = raw[(family, key, value, seed)][metric]
            den = raw[(family, key, baseline, seed)][metric]
            if num > 0.0 and den > 0.0 and num == num and den == den:
                out.append(math.log(num / den))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--axis", choices=sorted(AXES), default="k_periods")
    ap.add_argument("--values", help="comma-separated override of the axis values")
    ap.add_argument("--baseline", help="override the ratio denominator")
    ap.add_argument("--family", choices=sorted(FAMILIES), action="append")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--seeds", type=int, default=1)
    ap.add_argument("--wave-seed", type=int, default=11,
                    help="phase seed of the synthesized transition record")
    ap.add_argument("--skip-transition", action="store_true")
    ap.add_argument("--data-dir", type=Path,
                    default=ROOT / "tests" / "kalman_ou_iii")
    ap.add_argument("--csv", type=Path, help="write the raw rows here")
    args = ap.parse_args()

    families = args.family or ["OU_II", "OU_III"]
    env_var, values, baseline = AXES[args.axis]
    if args.values:
        values = tuple(v.strip() for v in args.values.split(",") if v.strip())
    if args.baseline:
        baseline = args.baseline
    if baseline not in values:
        raise SystemExit(f"baseline {baseline!r} is not among {values}")
    others = [v for v in values if v != baseline]

    # seed=None reproduces the simulator's default seeds exactly, so --seeds 1
    # stays bit-identical to the deterministic protocol.
    seeds = [None] if args.seeds <= 1 else list(range(1, args.seeds + 1))

    tmpdir_holder = tempfile.TemporaryDirectory(prefix="ou_sigma_horizon_")
    tmpdir = Path(tmpdir_holder.name)

    scenarios: list[tuple[str, str, float, bool]] = [
        (record, spectrum, hs, False) for spectrum, hs, record in STATIONARY
    ]
    # The stationary records are replayed as released, so a seed varies only the
    # IMU noise; the transition record is synthesized, so its wave phases are
    # re-drawn per seed as well.  One crossfade realization is not enough to
    # read a blend-interval effect off, and the pairing is unaffected: both
    # settings see the identical record at a given seed.
    paths: dict[tuple[str, int | None], Path] = {}
    for record, _, _, _ in scenarios:
        for seed in seeds:
            paths[(record, seed)] = args.data_dir / record
    if not args.skip_transition:
        scenarios.append((TRANSITION_KEY, TRANSITION_LABEL, 4.0, True))
        for seed in seeds:
            wave_seed = args.wave_seed if seed is None else args.wave_seed + seed
            paths[(TRANSITION_KEY, seed)] = transition_record(
                wave_seed, args.data_dir, tmpdir)

    keys = [key for key, _, _, _ in scenarios]
    stationary_keys = [key for key, _, _, seg in scenarios if not seg]
    by_key = {key: (spectrum, hs, seg) for key, spectrum, hs, seg in scenarios}

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {
            (family, key, value, seed): pool.submit(
                run_one, family, paths[(key, seed)], env_var, value, seed,
                by_key[key][2])
            for family in families
            for key in keys
            for value in values
            for seed in seeds
        }
        raw = {k: f.result() for k, f in futures.items()}

    results = {
        (family, key, value): {
            field: sum(raw[(family, key, value, s)][field] for s in seeds) / len(seeds)
            for field in FIELDS
        }
        for family in families
        for key in keys
        for value in values
    }

    rows = []
    for family in families:
        metrics = list(RATIO_METRICS)
        if not args.skip_transition:
            metrics += list(TRANSITION_RATIO_METRICS)

        for metric, label in metrics:
            transition_only = metric.startswith("seg_")
            shown = [TRANSITION_KEY] if transition_only else keys
            pooled_over = [TRANSITION_KEY] if transition_only else stationary_keys
            print()
            print(f"=== {family}: {label}, last {int(WINDOW_SEC)} s, "
                  f"{env_var} vs {baseline} ===")
            head = f"{'record':<12}{'Hs':>6}{baseline:>13}"
            for value in others:
                head += f"{value:>13}{'delta':>9}"
            print(head)
            for key in shown:
                spectrum, hs, _ = by_key[key]
                base = results[(family, key, baseline)][metric]
                line = f"{spectrum:<12}{hs:>6.2f}{base:>13.4f}"
                for value in others:
                    got = results[(family, key, value)][metric]
                    line += f"{got:>13.4f}{pct(got, base):>9}"
                print(line)

            tail = f"{'pooled ratio':<18}{'':>13}"
            for value in others:
                ratio, _, _ = summarize_ratio(
                    paired_logs(raw, family, metric, value, baseline,
                                pooled_over, seeds))
                tail += f"{ratio:>12.4f}x{'':>9}"
            print(tail)
            if len(seeds) > 1:
                ci = f"{'  paired 95% CI':<18}{'':>13}"
                for value in others:
                    _, lo, hi = summarize_ratio(
                        paired_logs(raw, family, metric, value, baseline,
                                    pooled_over, seeds))
                    ci += f"{'[' + f'{lo:.4f}, {hi:.4f}' + ']':>22}"
                print(ci)

        for key in keys:
            spectrum, hs, _ = by_key[key]
            for value in values:
                row = dict(results[(family, key, value)])
                row.update(family=family, record=key, sea_state=spectrum,
                           hs_m=hs, axis=args.axis, setting=value,
                           seeds=len(seeds))
                rows.append(row)

    if args.csv:
        with args.csv.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["family", "sea_state", "hs_m", "record", "axis",
                            "setting", "seeds", *FIELDS])
            writer.writeheader()
            writer.writerows(rows)
        print(f"\nwrote {args.csv}")

    tmpdir_holder.cleanup()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
