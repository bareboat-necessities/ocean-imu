#!/usr/bin/env python3
"""Charts for the two-frame Lie-group filter's deterministic replays.

Reads every ``*_fusion_tfg.csv`` the simulator wrote in this directory and
emits one PGF and one SVG per figure, which the CI step then moves into
``doc/kalman_tfg/`` for the manuscript to include.

Deliberately narrower than ``kalman_ou_iii-plots.py``: the TFG orchestrator
does not carry the wave-direction estimator or the alternative frequency
trackers, so the simulator writes NaN into those columns.  Plotting them would
render the absence of a feature as a flat line and invite it to be read as a
measurement.  The figures here cover what this filter actually estimates --
displacement, attitude, and the adaptation schedule.
"""

from __future__ import annotations

import glob
import os
import re

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "pgf.texsystem": "pdflatex",
        "font.family": "serif",
        "text.usetex": False,
        "pgf.rcfonts": False,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "figure.autolayout": True,
    }
)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def latex_safe(s: str) -> str:
    return s.replace("_", r"\_")


def make_subplots(nrows: int, title: str, width: float = 10.0, row_height: float = 2.4):
    fig, axes = plt.subplots(nrows, 1, figsize=(width, row_height * nrows), sharex=True)
    if nrows == 1:
        axes = [axes]
    fig.suptitle(title)
    return fig, axes


_PGF_OK = True


def finalize(fig, outbase: str, suffix: str = "", exts=("pgf", "svg")) -> None:
    """Write each requested format, tolerating a missing TeX installation.

    CI installs TeX Live and gets the PGF the manuscript includes. A developer
    without it should still get the SVGs rather than a traceback, so a PGF
    failure is reported once and then skipped.
    """
    global _PGF_OK
    for ext in exts:
        if ext == "pgf" and not _PGF_OK:
            continue
        try:
            fig.savefig(f"{outbase}{suffix}.{ext}", format=ext, bbox_inches="tight")
        except (RuntimeError, FileNotFoundError) as exc:
            if ext != "pgf":
                raise
            _PGF_OK = False
            print(f"PGF export unavailable ({exc}); writing SVG only")
    plt.close(fig)


# Charts show the last PLOT_WINDOW_S of the replay, matching the OU-III charts
# (which fix the window at 1140-1200 s of a 1200 s record). Taken from the
# record's own end rather than hardcoded, so a record of a different length
# still plots its tail instead of selecting nothing.
PLOT_WINDOW_S = 60.0

# The results table scores a longer window. RMS annotations quote that scored
# figure rather than the RMS of the 60 s on screen, so a chart and the table
# cannot appear to disagree about the same record.
SCORED_WINDOW_S = 900.0


def tail_window(df: pd.DataFrame, seconds: float) -> pd.DataFrame:
    """The last `seconds` of the replay."""
    if "time" not in df.columns or df.empty:
        return df
    t_end = float(df["time"].iloc[-1])
    return df[df["time"] >= t_end - seconds]


def plot_window(df: pd.DataFrame) -> pd.DataFrame:
    return tail_window(df, PLOT_WINDOW_S)


def scored_window(df: pd.DataFrame) -> pd.DataFrame:
    return tail_window(df, SCORED_WINDOW_S)


def rms(series: pd.Series) -> float:
    v = series.to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    return float(np.sqrt(np.mean(v * v))) if v.size else float("nan")


def plot_displacement(df: pd.DataFrame, base: str, label: str) -> None:
    win = plot_window(df)
    scored = scored_window(df)
    fig, axes = make_subplots(3, f"{label}: displacement (last {PLOT_WINDOW_S:.0f} s)")
    for ax, axis in zip(axes, ("x", "y", "z")):
        ax.plot(win["time"], win[f"disp_ref_{axis}"], lw=0.9, label="reference")
        ax.plot(win["time"], win[f"disp_est_{axis}"], lw=0.9, label="estimate")
        err = rms(scored[f"disp_est_{axis}"] - scored[f"disp_ref_{axis}"])
        ax.set_ylabel(f"{axis.upper()} [m]")
        ax.legend(loc="upper right", fontsize="small",
                  title=f"RMS err {err:.4f} m ({SCORED_WINDOW_S:.0f} s)",
                  title_fontsize="small")
    axes[-1].set_xlabel("time [s]")
    finalize(fig, base, "_disp")


def plot_heave(df: pd.DataFrame, base: str, label: str) -> None:
    """Vertical channel on its own -- it is the headline number."""
    win = plot_window(df)
    scored = scored_window(df)
    fig, axes = make_subplots(
        2, f"{label}: vertical channel (last {PLOT_WINDOW_S:.0f} s)", row_height=2.6)
    axes[0].plot(win["time"], win["disp_ref_z"], lw=0.9, label="reference")
    axes[0].plot(win["time"], win["disp_est_z"], lw=0.9, label="estimate")
    axes[0].set_ylabel("heave [m]")
    axes[0].legend(loc="upper right", fontsize="small")
    err = win["disp_est_z"] - win["disp_ref_z"]
    axes[1].plot(win["time"], err, lw=0.8, color="tab:red")
    axes[1].set_ylabel("error [m]")
    axes[1].set_xlabel("time [s]")
    scored_err = scored["disp_est_z"] - scored["disp_ref_z"]
    axes[1].legend([f"RMS {rms(scored_err):.4f} m over {SCORED_WINDOW_S:.0f} s"],
                   loc="upper right", fontsize="small")
    finalize(fig, base, "_heave")


def plot_attitude(df: pd.DataFrame, base: str, label: str) -> None:
    win = plot_window(df)
    scored = scored_window(df)
    fig, axes = make_subplots(3, f"{label}: attitude (last {PLOT_WINDOW_S:.0f} s)")
    for ax, name in zip(axes, ("roll", "pitch", "yaw")):
        ax.plot(win["time"], win[f"{name}_ref"], lw=0.9, label="reference")
        ax.plot(win["time"], win[f"{name}_est"], lw=0.9, label="estimate")
        # Yaw wraps; compare on the shortest arc so the RMS is not dominated by
        # a 360 degree jump that is not an error.
        d = scored[f"{name}_est"] - scored[f"{name}_ref"]
        if name == "yaw":
            d = (d + 180.0) % 360.0 - 180.0
        ax.set_ylabel(f"{name} [deg]")
        ax.legend(loc="upper right", fontsize="small",
                  title=f"RMS err {rms(d):.3f} deg ({SCORED_WINDOW_S:.0f} s)",
                  title_fontsize="small")
    axes[-1].set_xlabel("time [s]")
    finalize(fig, base, "_att")


def plot_tuner(df: pd.DataFrame, base: str, label: str) -> None:
    """The adaptation schedule the filter actually ran with.

    Full run, not the 60 s tail: this chart is about convergence and the
    startup transient, neither of which is visible in the last minute.
    """
    panels = [
        ("tau_applied", r"$\tau$ applied [s]"),
        ("sigma_a_applied", r"$\sigma_{a_w}$ applied [m/s$^2$]"),
        ("R_p0_applied", r"$r_S$ applied [m$\cdot$s]"),
        ("accel_var_tuner", r"accel variance [m$^2$/s$^4$]"),
    ]
    panels = [(k, y) for k, y in panels if k in df.columns]
    if not panels:
        return
    fig, axes = make_subplots(len(panels), f"{label}: adaptation schedule")
    for ax, (key, ylabel) in zip(axes, panels):
        ax.plot(df["time"], df[key], lw=0.9)
        ax.set_ylabel(ylabel)
    axes[-1].set_xlabel("time [s]")
    finalize(fig, base, "_tuner")


def plot_bias(df: pd.DataFrame, base: str, label: str) -> None:
    """Accelerometer bias, which is where a seeded tilt error goes to hide.

    Full run for the same reason as the tuner chart: a bias that walks away
    during startup and never comes back is the failure mode worth seeing.
    """
    cols = [f"acc_bias_{a}" for a in ("x", "y", "z")]
    if not all(c in df.columns for c in cols):
        return
    fig, axes = make_subplots(1, f"{label}: estimated accelerometer bias", row_height=2.8)
    for a in ("x", "y", "z"):
        axes[0].plot(df["time"], df[f"acc_bias_{a}"], lw=0.9, label=a.upper())
    axes[0].set_ylabel(r"bias [m/s$^2$]")
    axes[0].set_xlabel("time [s]")
    axes[0].legend(loc="upper right", fontsize="small", ncol=3)
    finalize(fig, base, "_bias")


def main() -> None:
    files = sorted(glob.glob("*_fusion_tfg.csv"))
    if not files:
        print("no *_fusion_tfg.csv in this directory; nothing to plot")
        return

    for path in files:
        df = pd.read_csv(path)
        if df.empty or "time" not in df.columns:
            print(f"skipping {path}: no usable rows")
            continue
        # w3d_<case>_fusion_tfg.csv -> tfg_<case>
        stem = os.path.splitext(os.path.basename(path))[0]
        case = re.sub(r"^w3d_", "", re.sub(r"_fusion_tfg$", "", stem))
        base = f"tfg_{case}"
        label = latex_safe(case)

        plot_displacement(df, base, label)
        plot_heave(df, base, label)
        plot_attitude(df, base, label)
        plot_tuner(df, base, label)
        plot_bias(df, base, label)
        print(f"plotted {path} -> {base}_*.pgf/.svg")


if __name__ == "__main__":
    main()
