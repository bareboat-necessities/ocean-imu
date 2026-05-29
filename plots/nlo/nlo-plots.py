#!/usr/bin/env python3
import glob
import os
import re
import sys
from pathlib import Path

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))
from plot_sampling import BASE_SAMPLE_RATE_HZ, get_decimation_step

# === Matplotlib PGF/LaTeX config ===
mpl.use("pgf")
plt.rcParams.update({
    "pgf.texsystem": "xelatex",
    "font.family": "serif",
    "text.usetex": True,
    "pgf.rcfonts": False,
    "pgf.preamble": "\n".join([
        r"\usepackage{fontspec}",
        r"\usepackage{unicode-math}",
        r"\usepackage{amsmath}",
        r"\setmainfont{DejaVu Serif}",
        r"\setmathfont{Latin Modern Math}",
        r"\providecommand{\mathdefault}[1]{#1}"
    ])
})

# === Config ===
DATA_DIR = "./"
SAMPLE_RATE_HZ = BASE_SAMPLE_RATE_HZ
SKIP_TIME_S = 1140.0
PLOT_TIME_S = 60.0
MAX_TIME_S = SKIP_TIME_S + PLOT_TIME_S
MAX_ROWS = int(SAMPLE_RATE_HZ * MAX_TIME_S)

PLOT_ERRORS = False
PLOT_YAW = False   # TVG no-mag/no-GNSS has free yaw; keep false by default.

height_groups = {
    "low": 0.27,
    "medium": 1.50,
    "high": 8.50,
}

ALLOWED_WAVES = {"jonswap", "pmstokes"}

SUFFIX = "tvg_nlo_nomag_nognss"
FILE_GLOB = f"*_{SUFFIX}.csv"

pattern = re.compile(
    rf".*?_(?P<wave>[a-zA-Z0-9]+)_H(?P<height>[-0-9\.]+).*?_{SUFFIX}\.csv$"
)

VARIANT_LABEL = "TimeVarGain NLO, VVR-only, no GNSS, no mag"


def latex_safe(s: str) -> str:
    return s.replace("_", r"\_")


def col(df: pd.DataFrame, name: str):
    return df[name] if name in df.columns else None


def err_series(df: pd.DataFrame, est_col: str, ref_col: str):
    if est_col not in df.columns or ref_col not in df.columns:
        return None
    return df[est_col] - df[ref_col]


def make_subplots(nrows: int, title: str, width: float = 10.0,
                  row_height: float = 2.5, sharex: bool = True):
    fig_height = row_height * nrows
    fig, axes = plt.subplots(nrows, 1, figsize=(width, fig_height), sharex=sharex)
    fig.suptitle(title)
    if nrows == 1:
        axes = [axes]
    return fig, axes


def finalize_plot(fig, outbase: str, suffix: str = "", exts=("pgf", "svg")):
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    for ext in exts:
        fig.savefig(f"{outbase}{suffix}.{ext}", format=ext, bbox_inches="tight")
    plt.close(fig)


def plot_series_or_missing(ax, time, df, key, label=None, linestyle="-", linewidth=1.2):
    if key not in df.columns:
        ax.text(0.01, 0.50, "Missing: " + latex_safe(key), transform=ax.transAxes)
        return False
    ax.plot(time, df[key], label=label or latex_safe(key), linestyle=linestyle, linewidth=linewidth)
    return True


def plot_vec3_panel(ax, time, df, prefix, ylabel, norm_col=None):
    keys = [f"{prefix}_x", f"{prefix}_y", f"{prefix}_z"]
    labels = ["x", "y", "z"]

    have_any = False
    for key, lab in zip(keys, labels):
        if key in df.columns:
            ax.plot(time, df[key], linewidth=1.0, label=lab)
            have_any = True

    if norm_col and norm_col in df.columns:
        ax.plot(time, df[norm_col], linewidth=1.1, linestyle="--", label="norm")
        have_any = True

    if not have_any:
        ax.text(0.01, 0.50, "Missing: " + latex_safe(prefix + "_*"), transform=ax.transAxes)

    ax.set_ylabel(ylabel)
    ax.grid(True)
    if have_any:
        ax.legend(loc="upper right", fontsize=8)


files = glob.glob(os.path.join(DATA_DIR, FILE_GLOB))
if not files:
    print(f"No {FILE_GLOB} files found in", DATA_DIR)
    sys.exit(0)

for fname in files:
    basename = os.path.basename(fname)

    m = pattern.match(basename)
    if not m:
        continue

    wave_type = m.group("wave").lower()
    if wave_type not in ALLOWED_WAVES:
        continue

    try:
        height_val = float(m.group("height"))
    except (TypeError, ValueError):
        continue

    group_name = None
    for name, h in height_groups.items():
        if abs(height_val - h) < 1e-6:
            group_name = name
            break

    if group_name is None:
        continue

    outbase = f"w3d_tvg_nlo_{wave_type}_{group_name}"
    outbase = re.sub(r"[^A-Za-z0-9_\-]", "_", outbase)
    outbase = os.path.join(DATA_DIR, outbase)

    print(f"Plotting {fname} → {outbase} ({VARIANT_LABEL}) ...")

    df = pd.read_csv(fname, nrows=MAX_ROWS)
    df = df[(df["time"] >= SKIP_TIME_S) & (df["time"] <= MAX_TIME_S)].reset_index(drop=True)
    df = df.iloc[::get_decimation_step(SAMPLE_RATE_HZ)].reset_index(drop=True)

    if df.empty:
        print(f"Skipping empty time window for {fname}")
        continue

    time = df["time"]
    title = latex_safe(f"{VARIANT_LABEL}: {basename}")

    # === Attitude: roll/pitch only by default ===
    angle_specs = [
        ("roll_ref", "roll_est", "Roll (deg)"),
        ("pitch_ref", "pitch_est", "Pitch (deg)"),
    ]
    if PLOT_YAW:
        angle_specs.append(("yaw_ref", "yaw_est", "Yaw free (deg)"))

    nrows = len(angle_specs) if not PLOT_ERRORS else 2 * len(angle_specs)
    fig, axes = make_subplots(nrows, title)

    for i, (ref_col, est_col, label) in enumerate(angle_specs):
        if PLOT_ERRORS:
            ax_val = axes[2 * i]
            ax_err = axes[2 * i + 1]
        else:
            ax_val = axes[i]
            ax_err = None

        if ref_col in df.columns:
            ax_val.plot(time, df[ref_col], label="Reference", linewidth=1.5)
        else:
            ax_val.text(0.01, 0.70, "Missing: " + latex_safe(ref_col), transform=ax_val.transAxes)

        if est_col in df.columns:
            ax_val.plot(time, df[est_col], label="Estimated", linewidth=1.0, linestyle="--")
        else:
            ax_val.text(0.01, 0.50, "Missing: " + latex_safe(est_col), transform=ax_val.transAxes)

        ax_val.set_ylabel(latex_safe(label))
        ax_val.grid(True)
        ax_val.legend(loc="upper right", fontsize=8)

        if PLOT_ERRORS:
            e = err_series(df, est_col, ref_col)
            if e is not None:
                ax_err.plot(time, e)
            else:
                ax_err.text(0.01, 0.50, "Missing ref/est for error", transform=ax_err.transAxes)
            ax_err.set_ylabel("Error [deg]")
            ax_err.grid(True)

    axes[-1].set_xlabel("Time (s)")
    finalize_plot(fig, outbase, "_att")

    # === Z-axis kinematics ===
    nrows = 6 if PLOT_ERRORS else 3
    fig, axes = make_subplots(nrows, title + " (Z-axis)")

    for i, prefix in enumerate(["disp", "vel", "acc"]):
        if PLOT_ERRORS:
            ax_val = axes[2 * i]
            ax_err = axes[2 * i + 1]
        else:
            ax_val = axes[i]
            ax_err = None

        ref_col = f"{prefix}_ref_z"
        est_col = f"{prefix}_est_z"

        if ref_col in df.columns:
            ax_val.plot(time, df[ref_col], label="Ref", linewidth=1.4)
        else:
            ax_val.text(0.01, 0.70, "Missing: " + latex_safe(ref_col), transform=ax_val.transAxes)

        if est_col in df.columns:
            ax_val.plot(time, df[est_col], label="Est", linewidth=1.1, linestyle="--")
        else:
            ax_val.text(0.01, 0.50, "Missing: " + latex_safe(est_col), transform=ax_val.transAxes)

        ax_val.set_ylabel(f"{prefix.capitalize()} Z")
        ax_val.grid(True)
        ax_val.legend(loc="upper right", fontsize=8)

        if PLOT_ERRORS:
            e = err_series(df, est_col, ref_col)
            if e is not None:
                ax_err.plot(time, e)
            else:
                ax_err.text(0.01, 0.50, "Missing ref/est for error", transform=ax_err.transAxes)
            ax_err.set_ylabel("Error")
            ax_err.grid(True)

    axes[-1].set_xlabel("Time (s)")
    finalize_plot(fig, outbase, "_zkin")

    # === TVG scalar gains/states ===
    scalar_panels = [
        ("tvg_k1", r"$k_1$"),
        ("tvg_kI", r"$k_I$"),
        ("tvg_vartheta", r"$\vartheta$"),
        ("tvg_p0z_hat", r"$\hat p^0_z$"),
    ]

    fig, axes = make_subplots(len(scalar_panels), title + " (TVG scalar diagnostics)")
    for ax, (key, ylabel) in zip(axes, scalar_panels):
        ok = plot_series_or_missing(ax, time, df, key, label=latex_safe(key), linewidth=1.2)
        ax.set_ylabel(ylabel)
        ax.grid(True)
        if ok:
            ax.legend(loc="upper right", fontsize=8)

    axes[-1].set_xlabel("Time (s)")
    finalize_plot(fig, outbase, "_tvg_scalars")

    # === TVG vectors ===
    vector_panels = [
        ("tvg_xi_n", r"$\xi^n$", "tvg_xi_norm"),
        ("tvg_fhat_n", r"$\hat f^n$", "tvg_fhat_norm"),
        ("tvg_sigma_b", r"$\sigma^b$", "tvg_sigma_norm"),
        ("tvg_gyro_bias_b", r"$\hat b_g^b$", "tvg_gyro_bias_norm"),
    ]

    fig, axes = make_subplots(len(vector_panels), title + " (TVG vector diagnostics)")
    for ax, (prefix, ylabel, norm_col) in zip(axes, vector_panels):
        plot_vec3_panel(ax, time, df, prefix, ylabel, norm_col=norm_col)

    axes[-1].set_xlabel("Time (s)")
    finalize_plot(fig, outbase, "_tvg_vectors")

    # === Compact vertical observer internals ===
    compact_cols = [
        ("tvg_fhat_n_z", r"$\hat f^n_z$"),
        ("tvg_xi_n_z", r"$\xi^n_z$"),
        ("tvg_sigma_b_x", r"$\sigma^b_x$"),
        ("tvg_sigma_b_y", r"$\sigma^b_y$"),
    ]

    fig, axes = make_subplots(len(compact_cols), title + " (Vertical internals)")
    for ax, (key, ylabel) in zip(axes, compact_cols):
        ok = plot_series_or_missing(ax, time, df, key, label=latex_safe(key), linewidth=1.2)
        ax.set_ylabel(ylabel)
        ax.grid(True)
        if ok:
            ax.legend(loc="upper right", fontsize=8)

    axes[-1].set_xlabel("Time (s)")
    finalize_plot(fig, outbase, "_tvg_vertical")
