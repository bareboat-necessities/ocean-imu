#!/usr/bin/env python3
import glob
import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

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
        r"\providecommand{\mathdefault}[1]{#1}",
    ]),
})

# === Height classification ===
height_groups = {"low": 0.27, "medium": 1.50, "high": 8.50}


def classify_height(h):
    diffs = {lvl: abs(h - v) for lvl, v in height_groups.items()}
    return min(diffs, key=diffs.get)


# === Filename pattern ===
pattern = re.compile(
    r"reg_spectrum_"
    r"(?P<tracker>[^_]+)_"
    r"(?P<wave>[A-Za-z0-9]+)_"
    r"H(?P<height>[-0-9\.]+)"
    r"(?:_L(?P<length>[-0-9\.]+))?"
    r"(?:_A(?P<azimuth>[-0-9\.]+))?"
    r"(?:_P(?P<phase>[-0-9\.]+))?"
    r"_N(?P<noise>[-0-9\.]+)"
    r"_B(?P<bias>[-0-9\.]+)"
    r"\.csv"
)


def pick_column(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


# === Load all files and group by (wave, height) ===
groups = {}
for f in sorted(glob.glob("reg_spectrum_*.csv")):
    m = pattern.search(os.path.basename(f))
    if not m:
        print(f"Skipping unrecognized: {f}")
        continue
    wave = m.group("wave")
    height = float(m.group("height"))
    tracker = m.group("tracker")
    key = (wave, height)
    groups.setdefault(key, []).append((tracker, f))

if not groups:
    print("No recognized reg_spectrum_*.csv files found.")
    raise SystemExit(0)

# === Plot each wave × height group with all trackers ===
for (wave, height), tracker_files in groups.items():
    level = classify_height(height)
    print(f"\nPlotting {wave} {level} sea with {len(tracker_files)} trackers")

    fig, axes = plt.subplots(3, 1, figsize=(6.5, 8.5), sharex=True)

    for idx, (tracker, f) in enumerate(tracker_files):
        df = pd.read_csv(f, comment="#")
        if "freq_hz" not in df.columns:
            print(f"  skipping {f} (no freq_hz)")
            continue

        label = tracker.capitalize()
        lw = 1.8

        # 1. Amplitude spectra
        if "A_eta_est" in df.columns:
            axes[0].semilogx(df["freq_hz"], df["A_eta_est"], label=f"{label} est", lw=lw)
        if "A_eta_ref" in df.columns and idx == 0:
            axes[0].semilogx(df["freq_hz"], df["A_eta_ref"], "--", color="gray", label="Reference", lw=1.2)

        # 2. Energy density
        if "E_eta_est" in df.columns:
            axes[1].semilogx(df["freq_hz"], df["E_eta_est"], label=f"{label} est", lw=lw)
        if "E_eta_ref" in df.columns and idx == 0:
            axes[1].semilogx(df["freq_hz"], df["E_eta_ref"], "--", color="gray", label="Reference", lw=1.2)

        # 3. Cumulative variance
        if "CumVar_est" in df.columns:
            est_norm = df["CumVar_est"] / max(df["CumVar_est"].iloc[-1], 1e-12)
            axes[2].semilogx(df["freq_hz"], est_norm, label=f"{label} est", lw=lw)
        if "CumVar_ref" in df.columns and idx == 0:
            ref_norm = df["CumVar_ref"] / max(df["CumVar_ref"].iloc[-1], 1e-12)
            axes[2].semilogx(df["freq_hz"], ref_norm, "--", color="gray", label="Reference", lw=1.2)

    axes[0].set_ylabel(r"$A_\eta(f)$ [m]")
    axes[1].set_ylabel(r"$E_\eta(f)=fS_\eta(f)$ [m$^2$]")
    axes[2].set_ylabel("Cumulative variance")
    axes[2].set_xlabel("Frequency [Hz]")
    for ax in axes:
        ax.grid(True, which="both", lw=0.3, ls=":")
        ax.legend(fontsize=7)
        ax.set_xlim(left=0.02, right=2.0)

    fig.suptitle(f"{wave.upper()} — {level.capitalize()} sea ($H_s={height:.2f}$ m)", fontsize=11, y=0.97)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    out_spectrum = f"spectrum_estimates_{wave}_{level}_spectrum.pgf"
    fig.savefig(out_spectrum, bbox_inches="tight")
    print(f"  Saved → {out_spectrum}")
    plt.close(fig)

    # --- Time-series plot expected by LaTeX doc ---
    fig2, ax2 = plt.subplots(figsize=(6.5, 3.2))
    plotted_any = False

    for tracker, f in tracker_files:
        df = pd.read_csv(f, comment="#")

        xcol = pick_column(df, ["time", "time_s", "t", "block_time_s", "block"])
        if xcol is None:
            xcol = df.columns[0]

        hs_col = pick_column(df, ["Hs", "Hs_est", "hs", "hs_est"])
        fp_col = pick_column(df, ["Fp", "Fp_est", "fp", "fp_est"])
        tp_col = pick_column(df, ["Tp", "Tp_est", "tp", "tp_est"])

        label = tracker.capitalize()
        if hs_col is not None:
            ax2.plot(df[xcol], df[hs_col], lw=1.6, label=f"{label} $H_s$")
            plotted_any = True

        if tp_col is not None:
            ax2.plot(df[xcol], df[tp_col], lw=1.6, linestyle="--", label=f"{label} $T_p$")
            plotted_any = True
        elif fp_col is not None:
            safe_fp = df[fp_col].clip(lower=1e-9)
            ax2.plot(df[xcol], 1.0 / safe_fp, lw=1.6, linestyle="--", label=f"{label} $T_p$")
            plotted_any = True

    if not plotted_any:
        ax2.text(0.5, 0.5, "No block time-series columns found", ha="center", va="center", transform=ax2.transAxes)

    ax2.set_xlabel("Block index / time")
    ax2.set_ylabel("Estimate")
    ax2.grid(True, alpha=0.3)
    if plotted_any:
        ax2.legend(fontsize=7)
    ax2.set_title(f"{wave.upper()} — {level.capitalize()} sea")

    out_timeseries = f"spectrum_estimates_{wave}_{level}_timeseries.pgf"
    fig2.tight_layout()
    fig2.savefig(out_timeseries, bbox_inches="tight")
    print(f"  Saved → {out_timeseries}")
    plt.close(fig2)

print("\nAll spectrum estimator plots generated successfully.")
