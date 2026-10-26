#!/usr/bin/env python3
"""Charts for the three-axis CARMA filter's deterministic dataset replays.

Reads the ``*_fusion_CARMA3_A41.csv`` replays that ``tests/kalman_carma/carma3-sim``
writes (default variant A41: CARMA(4,1) on x, y and z) and emits one PGF and
one SVG per figure.  CI moves them into ``doc/kalman_carma/`` for
``kalman_carma-charts.tex`` and ``kalman_carma.tex``.

The figure set is OU-III's (``plots/kalman_ou_iii/kalman_ou_iii-plots.py``) for
the same low/medium/high records, so the two galleries compare page by page.
The tracking charts show the last 60 s; the adaptation, bias and direction
charts show the whole replay because they are about convergence.  RMS
annotations are taken over the scored 900 s window, so they agree with the
results table rather than describing only the minute on screen.
"""

from __future__ import annotations

import glob
import os
import re
import sys
from pathlib import Path

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
    }
)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.append(str(Path(__file__).resolve().parents[1]))
from plot_sampling import BASE_SAMPLE_RATE_HZ, get_decimation_step  # noqa: E402
from svg_determinism import configure_svg, save_svg  # noqa: E402

configure_svg(matplotlib)

FAMILY = "CARMA3_A41"
HEIGHT_GROUPS = {"low": 0.27, "medium": 1.50, "high": 8.50}
ALLOWED_WAVES = {"jonswap", "pmstokes"}
PLOT_WINDOW_S = 60.0
SCORED_WINDOW_S = 900.0
# Whole-run charts are thinned to this rate; at the 20 Hz tracking-chart rate
# a 1200 s trace is too many PGF points for TeX to typeset comfortably.
FULL_RUN_RATE_HZ = 2.0

PATTERN = re.compile(
    r"^w3d_(?P<wave>[a-z0-9]+)_H(?P<height>[0-9.]+)_.*_fusion_" + FAMILY + r"\.csv$"
)

_PGF_OK = True


def make_subplots(nrows: int, width: float = 10.0, row_height: float = 2.4):
    fig, axes = plt.subplots(nrows, 1, figsize=(width, row_height * nrows), sharex=True)
    if nrows == 1:
        axes = [axes]
    return fig, axes


def finalize(fig, outbase: str, suffix: str) -> None:
    """Write PGF and SVG, tolerating a missing TeX installation.

    CI installs TeX Live and gets the PGF the documents include.  Without it a
    developer still gets the SVGs rather than a traceback.
    """
    global _PGF_OK
    fig.tight_layout()
    if _PGF_OK:
        try:
            fig.savefig(f"{outbase}{suffix}.pgf", format="pgf", bbox_inches="tight")
        except (RuntimeError, FileNotFoundError) as exc:
            _PGF_OK = False
            print(f"PGF export unavailable ({exc}); writing SVG only")
    save_svg(fig, f"{outbase}{suffix}.svg", bbox_inches="tight")
    plt.close(fig)


def rms(values) -> float:
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    return float(np.sqrt(np.mean(v * v))) if v.size else float("nan")


def tail(df: pd.DataFrame, seconds: float) -> pd.DataFrame:
    t_end = float(df["time"].iloc[-1])
    return df[df["time"] >= t_end - seconds]


def wrap_deg(d):
    return (d + 180.0) % 360.0 - 180.0


def legend_rms(ax, text: str) -> None:
    ax.legend(loc="upper right", fontsize="small", ncol=2, title=text, title_fontsize="small")


def plot_attitude(win, scored, base):
    fig, axes = make_subplots(3)
    for ax, name in zip(axes, ("roll", "pitch", "yaw")):
        ax.plot(win["time"], win[f"{name}_ref"], lw=1.2, label="reference")
        ax.plot(win["time"], win[f"{name}_est"], lw=0.9, ls="--", label="estimate")
        err = scored[f"{name}_est"] - scored[f"{name}_ref"]
        if name == "yaw":
            err = wrap_deg(err)
        ax.set_ylabel(f"{name} [deg]")
        legend_rms(ax, f"RMS err {rms(err):.3f} deg ({SCORED_WINDOW_S:.0f} s)")
    axes[-1].set_xlabel("time [s]")
    finalize(fig, base, "_att")


def plot_vertical(win, scored, base):
    fig, axes = make_subplots(3)
    units = {"disp": "m", "vel": "m/s", "acc": "m/s$^2$"}
    for ax, prefix in zip(axes, ("disp", "vel", "acc")):
        ax.plot(win["time"], win[f"{prefix}_ref_z"], lw=1.2, label="reference")
        ax.plot(win["time"], win[f"{prefix}_est_z"], lw=0.9, ls="--", label="estimate")
        env_col = {"disp": "disp_scale_m", "vel": "vel_scale_mps"}.get(prefix)
        if env_col and env_col in win.columns:
            env = win[env_col].to_numpy(dtype=float)
            ax.plot(win["time"], env, ls=":", lw=0.9, color="0.4", label="$\\pm$ envelope")
            ax.plot(win["time"], -env, ls=":", lw=0.9, color="0.4")
        err = scored[f"{prefix}_est_z"] - scored[f"{prefix}_ref_z"]
        ax.set_ylabel(f"{prefix} Z [{units[prefix]}]")
        legend_rms(ax, f"RMS err {rms(err):.4f} {units[prefix]} ({SCORED_WINDOW_S:.0f} s)")
    axes[-1].set_xlabel("time [s]")
    finalize(fig, base, "_zkin")


def plot_horizontal(win, scored, base):
    fig, axes = make_subplots(3)
    units = {"disp": "m", "vel": "m/s", "acc": "m/s$^2$"}
    for ax, prefix in zip(axes, ("disp", "vel", "acc")):
        for axis, color in (("x", "tab:blue"), ("y", "tab:orange")):
            ax.plot(win["time"], win[f"{prefix}_ref_{axis}"], lw=1.2, color=color,
                    label=f"ref {axis.upper()}")
            ax.plot(win["time"], win[f"{prefix}_est_{axis}"], lw=0.9, ls="--", color=color,
                    label=f"est {axis.upper()}")
        ex = rms(scored[f"{prefix}_est_x"] - scored[f"{prefix}_ref_x"])
        ey = rms(scored[f"{prefix}_est_y"] - scored[f"{prefix}_ref_y"])
        ax.set_ylabel(f"{prefix} XY [{units[prefix]}]")
        legend_rms(ax, f"RMS err X {ex:.4f}, Y {ey:.4f} {units[prefix]} ({SCORED_WINDOW_S:.0f} s)")
    axes[-1].set_xlabel("time [s]")
    finalize(fig, base, "_xykin")


def plot_tuner(full, base):
    """The operating point the CARMA priors actually ran with.

    tau_applied is the smoothed tuner time constant; the filter places both
    CARMA shapes at f = 1/(2 tau).  R_p0_applied is the shared harness column
    that carries the vertical frequency scale kappa_z = f / f_z,base.
    """
    t = full["time"]
    fig, axes = make_subplots(4)
    axes[0].plot(t, full["freq_tracker_hz"], lw=1.0, label="frequency tracker")
    tau = full["tau_applied"].to_numpy(dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        f_op = np.where(tau > 0.0, 0.5 / tau, np.nan)
    axes[0].plot(t, f_op, lw=1.0, ls="--", label=r"$1/(2\tau_{\mathrm{applied}})$")
    axes[0].set_ylabel("frequency [Hz]")
    axes[0].legend(loc="upper right", fontsize="small")
    std = np.sqrt(np.clip(full["accel_var_tuner"].to_numpy(dtype=float), 0.0, None))
    axes[1].plot(t, std, lw=1.0, label="band accel std")
    axes[1].plot(t, full["sigma_a_applied"], lw=1.0, ls="--", label=r"$\sigma_a$ applied")
    axes[1].set_ylabel(r"accel [m/s$^2$]")
    axes[1].legend(loc="upper right", fontsize="small")
    axes[2].plot(t, full["tau_applied"], lw=1.0)
    axes[2].set_ylabel(r"$\tau$ applied [s]")
    axes[3].plot(t, full["R_p0_applied"], lw=1.0)
    axes[3].set_ylabel(r"$\kappa_z$ applied")
    axes[-1].set_xlabel("time [s]")
    finalize(fig, base, "_tuner")


def plot_bias(full, base, kind):
    unit = "m/s$^2$" if kind == "acc" else "rad/s"
    fig, axes = make_subplots(3)
    for ax, axis in zip(axes, ("x", "y", "z")):
        ax.plot(full["time"], full[f"{kind}_bias_{axis}"], lw=1.2, label="true")
        ax.plot(full["time"], full[f"{kind}_bias_est_{axis}"], lw=0.9, ls="--", label="estimate")
        ax.set_ylabel(f"{kind} bias {axis.upper()} [{unit}]")
        ax.legend(loc="upper right", fontsize="small")
    axes[-1].set_xlabel("time [s]")
    finalize(fig, base, f"_{kind}_bias")


def plot_direction(full, base):
    panels = [
        ("dir_axis_deg", "direction [deg, axial]"),
        ("dir_uncert_deg", "uncertainty [deg]"),
        ("dir_conf", "confidence"),
        ("dir_sign_num", "sense (+1/-1/0)"),
    ]
    panels = [(c, y) for c, y in panels if c in full.columns]
    if not panels:
        return
    fig, axes = make_subplots(len(panels))
    for ax, (col, ylabel) in zip(axes, panels):
        if col == "dir_sign_num":
            ax.step(full["time"], full[col], where="post", lw=1.0)
            ax.set_yticks([-1, 0, 1])
        else:
            ax.plot(full["time"], full[col], lw=1.0)
        ax.set_ylabel(ylabel)
    axes[-1].set_xlabel("time [s]")
    finalize(fig, base, "_dir")


def main() -> None:
    files = sorted(glob.glob(f"*_fusion_{FAMILY}.csv"))
    if not files:
        print(f"no *_fusion_{FAMILY}.csv in this directory; nothing to plot")
        return
    step_tracking = get_decimation_step(BASE_SAMPLE_RATE_HZ)
    step_full = get_decimation_step(BASE_SAMPLE_RATE_HZ, FULL_RUN_RATE_HZ)
    for path in files:
        m = PATTERN.match(os.path.basename(path))
        if not m or m.group("wave") not in ALLOWED_WAVES:
            continue
        height = float(m.group("height"))
        group = next((g for g, h in HEIGHT_GROUPS.items() if abs(height - h) < 1e-6), None)
        if group is None:
            continue
        df = pd.read_csv(path)
        if df.empty:
            print(f"skipping {path}: no rows")
            continue
        base = f"w3d_carma3_{m.group('wave')}_{group}"
        scored = tail(df, SCORED_WINDOW_S)
        win = tail(df, PLOT_WINDOW_S).iloc[::step_tracking]
        full = df.iloc[::step_full]
        plot_attitude(win, scored, base)
        plot_vertical(win, scored, base)
        plot_horizontal(win, scored, base)
        plot_tuner(full, base)
        plot_bias(full, base, "acc")
        plot_bias(full, base, "gyro")
        plot_direction(full, base)
        print(f"plotted {path} -> {base}_*.pgf/.svg")


if __name__ == "__main__":
    main()
