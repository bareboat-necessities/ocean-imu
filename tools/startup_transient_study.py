#!/usr/bin/env python3
"""Startup transient of NLO, TFG, PII, OU-II, and OU-III under engine noise.

The truth is vessel motion, not the sea surface: every record is the pinned
oceanography-waves-lib v1.2.1 wave field passed through the estimated 28 ft
fin-keel sailboat RAO at the CG (``tools/sim_dataset.py``), and that motion is
what the simulated IMU measures.

Replays the first three minutes of every versioned wave record (Gerstner,
cnoidal, Fenton, JONSWAP, and PM-Stokes at the four pinned heights) through
each filter family with the nominal inboard-diesel vibration of
``W3dSimCommon`` switched on (2400 rpm, 0.60 m/s^2 hull RMS, 80 Hz sensor
bandwidth).  Everything else is the deployed configuration of each simulator,
including the front-end vibration guards where a family has one.

Each record is replayed in full (20 min) so that every estimator's own
steady-state error is known; the charts and transient metrics cover the first
three minutes.  An engine-off arm of the same replay separates what the engine
adds from the estimator's intrinsic startup.  Per family, record, and axis:

``settle_own_s``
    The last time in the first 180 s at which the 10 s centred rolling RMS of
    the displacement error exceeds ``OWN_FACTOR`` times the estimator's own
    steady-state error RMS (300-1200 s).  This is how long the startup
    oscillation lasts, whatever the steady-state accuracy is.
``settle_usable_s``
    The same, against ``USABLE_FRACTION`` of the rolling RMS of the true
    displacement: how long until the output is a usable displacement.
    ``180`` means not within the window (possibly never, see ``ss_ratio``).
``peak_err_m`` / ``peak_t_s``
    The largest absolute error in the first 180 s and when it happens.
``ss_err_rms_m`` / ``ss_ratio``
    Steady-state error RMS and its ratio to the true displacement RMS.

PII is a vertical-only observer, so its horizontal channels are reported as
not estimated rather than scored as a zero estimate.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = ROOT / "reports" / "results" / "startup_transient"

WINDOW_S = 180.0
DT = 0.005
ROLL_S = 10.0
OWN_FACTOR = 2.0
USABLE_FRACTION = 0.5
STEADY_S = (300.0, 1200.0)

ARMS = {
    "engine": {
        "W3D_ENGINE_RPM": "2400",
        "W3D_ENGINE_LEVEL_MPS2": "0.6",
        "W3D_ENGINE_BANDWIDTH_HZ": "80",
    },
    "no-engine": {},
}

# name -> (binary, args, vertical_only)
FAMILIES = {
    "NLO": (ROOT / "tests/nlo/nlo-sim", [], False),
    "TFG": (ROOT / "tests/kalman_tfg/kalman_tfg-sim", [], False),
    "PII": (ROOT / "tests/pii_observer/pii_observer-adaptive", [], True),
    "OU-II": (ROOT / "tests/kalman_ou_ii/kalman_ou_ii-sim", [], False),
    "OU-III": (ROOT / "tests/kalman_ou_iii/kalman_ou_iii-sim", [], False),
}
COLORS = {
    "NLO": "#1f77b4",
    "TFG": "#ff7f0e",
    "PII": "#2ca02c",
    "OU-II": "#d62728",
    "OU-III": "#9467bd",
}
WAVES = ("gerstner", "cnoidal", "fenton", "jonswap", "pmstokes")
HEIGHTS = ("0.270", "1.500", "4.000", "8.500")
AXES = ("x", "y", "z")


def find_records(data_dir: Path) -> list[Path]:
    out = []
    for wave in WAVES:
        for h in HEIGHTS:
            hits = sorted(data_dir.rglob(f"wave_data_{wave}_H{h}_*.csv"))
            if len(hits) != 1:
                raise FileNotFoundError(f"{wave} H{h}: found {len(hits)} records")
            out.append(hits[0])
    return out


def run_one(family: str, arm: str, record: Path, work: Path) -> dict:
    """Replay one record and reduce its time series before deleting it."""

    binary, args, vertical_only = FAMILIES[family]
    run_dir = work / arm / family / record.stem
    run_dir.mkdir(parents=True, exist_ok=True)
    for old in run_dir.iterdir():
        old.unlink()
    (run_dir / record.name).symlink_to(record.resolve())
    env = dict(os.environ, W3D_COLLECT_ALL_GATES="1", NLO_NOGATE="1",
               W3D_ALL_WAVE_TYPES="1", **ARMS[arm])
    done = subprocess.run([str(binary), *args], cwd=run_dir, env=env,
                          capture_output=True, text=True)
    engaged = "ENGINE_VIBRATION" in done.stdout
    if done.returncode not in (0, 1) or engaged != bool(ARMS[arm]):
        raise RuntimeError(f"{family}/{arm} {record.name}: exit {done.returncode}\n"
                           + done.stdout[-1500:] + done.stderr[-1500:])
    outs = list(run_dir.glob("w3d_*.csv"))
    if len(outs) != 1:
        raise RuntimeError(f"{family}/{arm} {record.name}: {len(outs)} outputs")
    names = ["time"] + [f"disp_ref_{a}" for a in AXES] + [f"disp_est_{a}" for a in AXES]
    frame = pd.read_csv(outs[0], usecols=names)
    for f in run_dir.iterdir():
        f.unlink()
    t = frame["time"].to_numpy()
    steady = (t >= STEADY_S[0]) & (t <= STEADY_S[1])
    head = t <= WINDOW_S
    out = {"head": {n: frame[n].to_numpy()[head] for n in names}, "steady": {}}
    for axis in AXES:
        if vertical_only and axis != "z":
            continue
        err = (frame[f"disp_est_{axis}"] - frame[f"disp_ref_{axis}"]).to_numpy()[steady]
        ref = frame[f"disp_ref_{axis}"].to_numpy()[steady]
        out["steady"][axis] = (float(np.sqrt(np.nanmean(err ** 2))),
                               float(np.sqrt(np.nanmean(ref ** 2))))
    return out


def rolling_rms(x: np.ndarray, n: int) -> np.ndarray:
    x = np.nan_to_num(x)
    c = np.convolve(x * x, np.ones(n) / n, mode="same")
    return np.sqrt(np.maximum(c, 0.0))


def settle_time(t: np.ndarray, bad: np.ndarray) -> float:
    # Ignore the last half-window, where the centred window is truncated.
    valid = t <= t[-1] - ROLL_S / 2
    bad = bad & valid
    if not bad.any():
        return 0.0
    last = int(np.nonzero(bad)[0][-1])
    return WINDOW_S if last >= int(np.nonzero(valid)[0][-1]) else float(t[last])


def score(run: dict, axis: str) -> dict[str, float]:
    ts = run["head"]
    t = ts["time"]
    err = ts[f"disp_est_{axis}"] - ts[f"disp_ref_{axis}"]
    n = int(round(ROLL_S / DT))
    err_rms = rolling_rms(err, n)
    ss_err, ss_ref = run["steady"][axis]
    k = int(np.nanargmax(np.abs(err)))
    return {
        "settle_own_s": settle_time(t, err_rms > OWN_FACTOR * ss_err),
        "settle_usable_s": settle_time(
            t, err_rms > USABLE_FRACTION * rolling_rms(ts[f"disp_ref_{axis}"], n)),
        "peak_err_m": float(abs(err[k])),
        "peak_t_s": float(t[k]),
        "ss_err_rms_m": ss_err,
        "ss_ref_rms_m": ss_ref,
        "ss_ratio": ss_err / ss_ref if ss_ref > 0 else float("nan"),
    }


def plot_record(record: str, series: dict[str, dict[str, np.ndarray]], out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(3, 2, figsize=(16, 10), sharex=True)
    ref = next(iter(series.values()))
    for r, axis in enumerate(AXES):
        a_disp, a_err = axs[r, 0], axs[r, 1]
        a_disp.plot(ref["time"], ref[f"disp_ref_{axis}"], color="black", lw=1.4,
                    label="truth")
        # Wide enough for the truth and for every settled error, but clipping
        # the startup excursions, which can be orders of magnitude larger.
        tail = ref["time"] >= WINDOW_S - 60.0
        settled = [np.nanmax(np.abs(ts[f"disp_est_{axis}"] - ts[f"disp_ref_{axis}"])[tail])
                   for fam, ts in series.items() if not (FAMILIES[fam][2] and axis != "z")]
        ylim = max(3.0 * np.nanmax(np.abs(ref[f"disp_ref_{axis}"])),
                   1.5 * max(settled),
                   0.5 * np.nanmax(np.abs(ref["disp_ref_z"])))
        for fam, ts in series.items():
            if FAMILIES[fam][2] and axis != "z":
                continue
            est = ts[f"disp_est_{axis}"]
            a_disp.plot(ts["time"], est, color=COLORS[fam], lw=0.8, label=fam)
            a_err.plot(ts["time"], est - ts[f"disp_ref_{axis}"], color=COLORS[fam],
                       lw=0.8, label=fam)
        for a in (a_disp, a_err):
            a.set_ylim(-ylim, ylim)
            a.grid(alpha=0.3)
            a.axhline(0, color="gray", lw=0.5)
        a_disp.set_ylabel(f"{axis} displacement (m)")
        a_err.set_ylabel(f"{axis} error (m)")
    axs[0, 0].set_title("Estimate vs truth")
    axs[0, 1].set_title("Error (estimate - truth)")
    axs[0, 0].legend(loc="upper right", ncol=3, fontsize=8)
    axs[2, 0].set_xlabel("time (s)")
    axs[2, 1].set_xlabel("time (s)")
    axs[2, 0].set_xlim(0, WINDOW_S)
    fig.suptitle(f"28 ft sailboat RAO at CG, {record.replace('wave_data_', '')}: "
                 f"first {WINDOW_S:.0f} s with engine noise (2400 rpm, 0.60 m/s^2)\n"
                 f"truth = vessel motion; startup excursions are clipped; "
                 f"PII is vertical only")
    fig.tight_layout()
    fig.savefig(out, dpi=90)
    plt.close(fig)


def plot_settle(rows: list[dict], metric: str, title: str, out: Path) -> None:
    """Bars are the engine arm; black ticks mark the engine-off arm."""

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    records = list(dict.fromkeys(r["record"] for r in rows))
    fams = list(FAMILIES)
    index = {(r["arm"], r["family"], r["record"], r["axis"]): r[metric] for r in rows}
    fig, axs = plt.subplots(3, 1, figsize=(16, 11), sharex=True)
    width = 0.8 / len(fams)
    x = np.arange(len(records))
    for ax, axis in zip(axs, AXES):
        for j, fam in enumerate(fams):
            xs = x + (j - (len(fams) - 1) / 2) * width
            on = [index.get(("engine", fam, rec, axis), np.nan) for rec in records]
            off = [index.get(("no-engine", fam, rec, axis), np.nan) for rec in records]
            ax.bar(xs, on, width, color=COLORS[fam], label=fam)
            ax.scatter(xs, off, marker="_", s=90, color="black", lw=1.5,
                       label="engine off" if j == 0 else None, zorder=3)
        ax.axhline(WINDOW_S, color="gray", ls="--", lw=0.8)
        ax.set_ylabel(f"{axis} settle time (s)")
        ax.set_ylim(0, WINDOW_S * 1.05)
        ax.grid(axis="y", alpha=0.3)
    axs[0].legend(ncol=6, loc="upper left", fontsize=9)
    axs[0].set_title(f"28 ft sailboat RAO at CG. {title}")
    axs[-1].set_xticks(x)
    axs[-1].set_xticklabels([r.replace("wave_data_", "").split("_L")[0] for r in records],
                            rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(out, dpi=90)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, required=True)
    ap.add_argument("--work-dir", type=Path, required=True,
                    help="scratch directory for per-run inputs and time series")
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    ap.add_argument("--jobs", type=int, default=4)
    args = ap.parse_args()

    records = find_records(args.data_dir)
    tasks = [(arm, fam, rec) for arm in ARMS for rec in records for fam in FAMILIES]
    with ThreadPoolExecutor(args.jobs) as pool:
        runs = list(pool.map(lambda t: run_one(t[1], t[0], t[2], args.work_dir), tasks))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    engine_series: dict[str, dict[str, dict[str, np.ndarray]]] = {}
    for (arm, fam, rec), run in zip(tasks, runs):
        if arm == "engine":
            engine_series.setdefault(rec.stem, {})[fam] = run["head"]
        for axis in run["steady"]:
            rows.append({"arm": arm, "family": fam, "record": rec.stem, "axis": axis,
                         **score(run, axis)})

    for rec, series in engine_series.items():
        plot_record(rec, series, args.output_dir / f"{rec.replace('wave_data_', '')}.png")
    plot_settle(rows, "settle_own_s",
                f"Startup oscillation: time until 10 s rolling error RMS stays below "
                f"{OWN_FACTOR:g}x the filter's own steady-state error (180 = not within 3 min)",
                args.output_dir / "settle_own.png")
    plot_settle(rows, "settle_usable_s",
                f"Usable output: time until 10 s rolling error RMS stays below "
                f"{USABLE_FRACTION:.0%} of truth RMS (180 = not within 3 min)",
                args.output_dir / "settle_usable.png")

    with (args.output_dir / "startup_transient.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"wrote {len(rows)} rows and {len(engine_series) + 2} charts to {args.output_dir}")


if __name__ == "__main__":
    main()
