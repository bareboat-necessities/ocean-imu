#!/usr/bin/env python3
import argparse
import glob
import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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

height_groups = {"low": 0.27, "medium": 1.50, "high": 8.50}

classic_pattern = re.compile(
    r"spectrum_"
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

wavelets_pattern = re.compile(
    r"spectrum_wavelets_"
    r"(?P<wave>[A-Za-z0-9]+)_"
    r"H(?P<height>[-0-9\.]+)"
    r"(?:_L(?P<length>[-0-9\.]+))?"
    r"(?:_A(?P<azimuth>[-0-9\.]+))?"
    r"(?:_P(?P<phase>[-0-9\.]+))?"
    r"_N(?P<noise>[-0-9\.]+)"
    r"_B(?P<bias>[-0-9\.]+)"
    r"\.csv"
)


def classify_height(h):
    diffs = {level: abs(h - value) for level, value in height_groups.items()}
    return min(diffs, key=diffs.get)


def pick_column(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def collect_groups(mode):
    groups = {}
    if mode == "classic":
        files = sorted(glob.glob("spectrum_*.csv"))
        for f in files:
            m = classic_pattern.search(os.path.basename(f))
            if not m:
                continue
            tracker = m.group("tracker")
            wave = m.group("wave")
            height = float(m.group("height"))
            groups.setdefault((wave, height), []).append((tracker, f))
    else:
        files = sorted(glob.glob("spectrum_wavelets_*.csv"))
        for f in files:
            m = wavelets_pattern.search(os.path.basename(f))
            if not m:
                continue
            wave = m.group("wave")
            height = float(m.group("height"))
            groups.setdefault((wave, height), []).append(("wavelets", f))
    return groups


def output_prefix(mode):
    return "spectrum_estimates" if mode == "classic" else "spectrum_wavelets_estimates"


def estimate_bulk_from_spectrum(df):
    fcol = pick_column(df, ["freq_hz", "f_hz", "freq"])
    scol = pick_column(df, ["S_eta_hz", "S_eta_est", "S_eta", "S_est"])
    if fcol is None or scol is None:
        return None, None

    f = df[fcol]
    s = df[scol].clip(lower=0.0)
    if len(f) < 2 or len(s) < 2:
        return None, None

    trapz = getattr(np, "trapezoid", np.trapz)
    m0 = trapz(s.to_numpy(), f.to_numpy())
    hs = 4.0 * (m0 ** 0.5) if m0 > 0 else None

    fp_idx = s.idxmax() if len(s) > 0 else None
    tp = None
    if fp_idx is not None:
        fp = float(f.loc[fp_idx])
        if fp > 0.0:
            tp = 1.0 / fp

    return hs, tp


def estimate_tp_from_spectrum(df, min_freq_hz=0.08):
    fcol = pick_column(df, ["freq_hz", "f_hz", "freq"])
    scol = pick_column(df, ["S_eta_hz", "S_eta_est", "S_eta", "S_est"])
    if fcol is None or scol is None:
        return None

    work = df[[fcol, scol]].copy()
    work = work[np.isfinite(work[fcol]) & np.isfinite(work[scol])]
    work = work[work[fcol] > 0.0]
    if work.empty:
        return None

    work[scol] = work[scol].clip(lower=0.0)
    filtered = work[work[fcol] >= min_freq_hz]
    if filtered.empty:
        filtered = work

    if filtered[scol].max() <= 0.0:
        return None

    peak_row = filtered.loc[filtered[scol].idxmax()]
    fp = float(peak_row[fcol])
    if fp <= 0.0:
        return None
    return 1.0 / fp


def plot_mode(mode):
    groups = collect_groups(mode)
    if not groups:
        print(f"No recognized CSV files for mode={mode}.")
        return

    out_prefix = output_prefix(mode)

    for (wave, height), tracker_files in groups.items():
        level = classify_height(height)
        print(f"Plotting {mode}: {wave} {level} sea ({len(tracker_files)} file(s))")

        fig, axes = plt.subplots(3, 1, figsize=(6.5, 8.5), sharex=True)

        for idx, (tracker, f) in enumerate(tracker_files):
            df = pd.read_csv(f, comment="#")
            if "freq_hz" not in df.columns:
                continue

            label = tracker.capitalize()
            lw = 1.8

            if "S_eta_hz" in df.columns:
                axes[0].loglog(df["freq_hz"], df["S_eta_hz"], label=f"{label} est", lw=lw)
            if "S_ref_interp" in df.columns and idx == 0:
                axes[0].loglog(df["freq_hz"], df["S_ref_interp"], "--", color="gray", label="Reference", lw=1.2)

            if "S_ratio" in df.columns:
                axes[1].semilogx(df["freq_hz"], df["S_ratio"], label=f"{label} est/ref", lw=lw)
            axes[1].axhline(1.0, color="gray", ls="--", lw=1.0)

            if "CumVar_est" in df.columns:
                est_norm = df["CumVar_est"] / max(df["CumVar_est"].iloc[-1], 1e-12)
                axes[2].semilogx(df["freq_hz"], est_norm, label=f"{label} est", lw=lw)
            if "CumVar_ref" in df.columns and idx == 0:
                ref_norm = df["CumVar_ref"] / max(df["CumVar_ref"].iloc[-1], 1e-12)
                axes[2].semilogx(df["freq_hz"], ref_norm, "--", color="gray", label="Reference", lw=1.2)

        axes[0].set_ylabel(r"$S_\eta(f)$ [m$^2$/Hz]")
        axes[1].set_ylabel(r"$S_{\eta,\mathrm{est}}/S_{\eta,\mathrm{ref}}$")
        axes[2].set_ylabel("Cumulative variance")
        axes[2].set_xlabel("Frequency [Hz]")

        for ax in axes:
            ax.grid(True, which="both", lw=0.3, ls=":")
            ax.legend(fontsize=7)
            ax.set_xlim(left=0.02, right=2.0)

        fig.suptitle(f"{wave.upper()} — {level.capitalize()} sea ($H_s={height:.2f}$ m)", fontsize=11, y=0.97)
        fig.tight_layout(rect=[0, 0, 1, 0.96])

        out_spectrum = f"{out_prefix}_{wave}_{level}_spectrum.pgf"
        fig.savefig(out_spectrum, bbox_inches="tight")
        plt.close(fig)

        hs_values = []
        tp_values = []

        for tracker, f in tracker_files:
            df = pd.read_csv(f, comment="#")
            label = tracker.capitalize()

            hs_col = pick_column(df, ["Hs", "Hs_est", "hs", "hs_est"])
            tp_col = pick_column(df, ["Tp", "Tp_est", "tp", "tp_est"])
            fp_col = pick_column(df, ["Fp", "Fp_est", "fp", "fp_est"])

            if hs_col is not None:
                hs_values.append((label, float(df[hs_col].iloc[-1])))
            tp_val = None
            if tp_col is not None:
                tp_val = float(df[tp_col].iloc[-1])
            elif fp_col is not None:
                safe_fp = max(float(df[fp_col].iloc[-1]), 1e-9)
                tp_val = 1.0 / safe_fp

            tp_from_spec = estimate_tp_from_spectrum(df)
            if tp_from_spec is not None:
                if tp_val is None or not np.isfinite(tp_val) or tp_val > 20.0 or tp_val < 0.5:
                    tp_val = tp_from_spec
                else:
                    fmin = float(df["freq_hz"].iloc[0]) if "freq_hz" in df.columns else None
                    if fmin is not None and fmin > 0.0:
                        edge_tp = 1.0 / fmin
                        if abs(tp_val - edge_tp) / edge_tp < 0.05:
                            tp_val = tp_from_spec
            if tp_val is not None:
                tp_values.append((label, tp_val))

            if hs_col is None or tp_val is None:
                hs_est, tp_est = estimate_bulk_from_spectrum(df)
                if hs_col is None and hs_est is not None:
                    hs_values.append((label, hs_est))
                if tp_val is None and tp_est is not None:
                    tp_values.append((label, tp_est))

        fig_hs, ax_hs = plt.subplots(figsize=(6.5, 2.8))
        if hs_values:
            labels = [x[0] for x in hs_values]
            values = [x[1] for x in hs_values]
            ax_hs.bar(labels, values, color="#4C78A8")
        else:
            ax_hs.text(0.5, 0.5, "No $H_s$ estimates found", ha="center", va="center", transform=ax_hs.transAxes)
        ax_hs.set_ylabel(r"$H_s$ [m]")
        ax_hs.set_title(f"{wave.upper()} — {level.capitalize()} sea")
        ax_hs.grid(True, axis="y", alpha=0.3)
        fig_hs.tight_layout()
        out_hs = f"{out_prefix}_{wave}_{level}_hs.pgf"
        fig_hs.savefig(out_hs, bbox_inches="tight")
        plt.close(fig_hs)

        fig_tp, ax_tp = plt.subplots(figsize=(6.5, 2.8))
        if tp_values:
            labels = [x[0] for x in tp_values]
            values = [x[1] for x in tp_values]
            ax_tp.bar(labels, values, color="#F58518")
        else:
            ax_tp.text(0.5, 0.5, "No $T_p$ estimates found", ha="center", va="center", transform=ax_tp.transAxes)
        ax_tp.set_ylabel(r"$T_p$ [s]")
        ax_tp.set_title(f"{wave.upper()} — {level.capitalize()} sea")
        ax_tp.grid(True, axis="y", alpha=0.3)
        fig_tp.tight_layout()
        out_tp = f"{out_prefix}_{wave}_{level}_tp.pgf"
        fig_tp.savefig(out_tp, bbox_inches="tight")
        plt.close(fig_tp)

        print(f"  Saved → {out_spectrum}")
        print(f"  Saved → {out_hs}")
        print(f"  Saved → {out_tp}")


def main():
    parser = argparse.ArgumentParser(description="Generate spectrum estimator PGF plots.")
    parser.add_argument(
        "--mode",
        choices=["classic", "wavelets"],
        default="classic",
        help="Use 'classic' for spectrum_estimator_*.csv or 'wavelets' for spectrum_wavelets_*.csv",
    )
    args = parser.parse_args()
    plot_mode(args.mode)


if __name__ == "__main__":
    main()
