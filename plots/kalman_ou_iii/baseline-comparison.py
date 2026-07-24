#!/usr/bin/env python3
import glob
import math
import os
import re
from pathlib import Path

import matplotlib as mpl
mpl.use("pgf")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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

OU_RE = re.compile(r".*?_(?P<wave>jonswap|pmstokes)_H(?P<height>[0-9.]+).*?_fusion_ou3\.csv$")
NLO_RE = re.compile(r".*?_(?P<wave>jonswap|pmstokes)_H(?P<height>[0-9.]+).*?_tvg_nlo_nomag_nognss\.csv$")
HEIGHTS = (0.27, 1.50, 4.00, 8.50)
WINDOW_S = 60.0


def index_files(pattern, regex):
    out = {}
    for name in glob.glob(pattern):
        match = regex.match(os.path.basename(name))
        if not match:
            continue
        key = (match.group("wave"), round(float(match.group("height")), 2))
        out[key] = name
    return out


def tail(df):
    if "time" not in df.columns or df.empty:
        raise ValueError("comparison CSV must contain a non-empty time column")
    end = float(df["time"].max())
    return df[df["time"] >= end - WINDOW_S].copy()


def rms(series):
    values = np.asarray(series, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return math.nan
    return float(np.sqrt(np.mean(values * values)))


def metrics(path, hs):
    df = tail(pd.read_csv(path))
    err_z = df["disp_est_z"] - df["disp_ref_z"]
    mean_z = float(np.nanmean(err_z))
    result = {
        "z_rms": rms(err_z),
        "z_demeaned_rms": rms(err_z - mean_z),
        "z_mean": mean_z,
        "roll_rms": rms(df["roll_est"] - df["roll_ref"]),
        "pitch_rms": rms(df["pitch_est"] - df["pitch_ref"]),
    }
    result["z_pct"] = 100.0 * result["z_rms"] / hs
    return result, df


def tex_num(x, digits=3):
    return "--" if not math.isfinite(x) else f"{x:.{digits}f}"


def write_table(rows, output):
    lines = [
        r"\begin{table*}[t]",
        r"  \centering",
        r"  \caption{Baseline comparison against the adapted time-varying-gain nonlinear observer (TVG--NLO), using the same reference records and final \SI{60}{s} scoring window.}",
        r"  \label{tab:nlo_baseline_comparison}",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{3.0pt}",
        r"  \begin{tabular}{@{}llrrrrrrrr@{}}",
        r"    \toprule",
        r"    Spectrum & $H_s$ & \multicolumn{2}{c}{Z RMS [m]} & \multicolumn{2}{c}{Z RMS [$\%H_s$]} & \multicolumn{2}{c}{Roll RMS [$^\circ$]} & \multicolumn{2}{c}{Pitch RMS [$^\circ$]} \\",
        r"    \cmidrule(lr){3-4}\cmidrule(lr){5-6}\cmidrule(lr){7-8}\cmidrule(lr){9-10}",
        r"    & & OU--III & TVG--NLO & OU--III & TVG--NLO & OU--III & TVG--NLO & OU--III & TVG--NLO \\",
        r"    \midrule",
    ]
    for row in rows:
        spectrum = "JONSWAP" if row["wave"] == "jonswap" else "PM--Stokes"
        ou, nlo = row["ou"], row["nlo"]
        lines.append(
            f"    {spectrum} & {row['hs']:.2f} & "
            f"{tex_num(ou['z_rms'])} & {tex_num(nlo['z_rms'])} & "
            f"{tex_num(ou['z_pct'], 1)} & {tex_num(nlo['z_pct'], 1)} & "
            f"{tex_num(ou['roll_rms'], 2)} & {tex_num(nlo['roll_rms'], 2)} & "
            f"{tex_num(ou['pitch_rms'], 2)} & {tex_num(nlo['pitch_rms'], 2)} \\\\" 
        )
    lines += [
        r"    \bottomrule",
        r"  \end{tabular}",
        r"\end{table*}",
    ]
    Path(output).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(ou_df, nlo_df, output_base):
    ou = tail(ou_df)
    nlo = tail(nlo_df)
    common_start = max(float(ou["time"].min()), float(nlo["time"].min()))
    common_end = min(float(ou["time"].max()), float(nlo["time"].max()))
    ou = ou[(ou["time"] >= common_start) & (ou["time"] <= common_end)]
    nlo = nlo[(nlo["time"] >= common_start) & (nlo["time"] <= common_end)]

    # Both runners use the same 200 Hz reference record, so nearest-time interpolation
    # is only a safeguard against a one-sample output offset.
    t = ou["time"].to_numpy()
    ref = ou["disp_ref_z"].to_numpy()
    ou_est = ou["disp_est_z"].to_numpy()
    nlo_est = np.interp(t, nlo["time"], nlo["disp_est_z"])

    fig, axes = plt.subplots(2, 1, figsize=(7.0, 4.8), sharex=True)
    axes[0].plot(t, ref, label="Reference", linewidth=1.3)
    axes[0].plot(t, ou_est, label="OU--III", linewidth=1.0)
    axes[0].plot(t, nlo_est, label="TVG--NLO", linewidth=1.0, linestyle="--")
    axes[0].set_ylabel("Vertical displacement [m]")
    axes[0].grid(True)
    axes[0].legend(loc="upper right", fontsize=8)

    axes[1].plot(t, ou_est - ref, label="OU--III error", linewidth=1.0)
    axes[1].plot(t, nlo_est - ref, label="TVG--NLO error", linewidth=1.0, linestyle="--")
    axes[1].set_ylabel("Error [m]")
    axes[1].set_xlabel("Time [s]")
    axes[1].grid(True)
    axes[1].legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    for ext in ("pgf", "svg"):
        fig.savefig(f"{output_base}.{ext}", format=ext, bbox_inches="tight")
    plt.close(fig)


def main():
    ou_files = index_files("*_fusion_ou3.csv", OU_RE)
    nlo_files = index_files("*_tvg_nlo_nomag_nognss.csv", NLO_RE)
    expected = {(wave, round(hs, 2)) for wave in ("jonswap", "pmstokes") for hs in HEIGHTS}
    missing = sorted(expected - set(ou_files) | expected - set(nlo_files))
    if missing:
        raise RuntimeError(f"missing OU/NLO comparison cases: {missing}")

    rows = []
    representative = None
    for wave in ("jonswap", "pmstokes"):
        for hs in HEIGHTS:
            key = (wave, round(hs, 2))
            ou_metrics, ou_df = metrics(ou_files[key], hs)
            nlo_metrics, nlo_df = metrics(nlo_files[key], hs)
            rows.append({"wave": wave, "hs": hs, "ou": ou_metrics, "nlo": nlo_metrics})
            if key == ("jonswap", 1.50):
                representative = (ou_df, nlo_df)

    write_table(rows, "w3d-baseline-results-generated.tex-part")
    if representative is None:
        raise RuntimeError("representative JONSWAP Hs=1.50 case not found")
    write_plot(*representative, "w3d_ou3_vs_tvg_nlo_jonswap_medium")


if __name__ == "__main__":
    main()
