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

PATTERNS = {
    "ou3": (re.compile(r".*?_(?P<wave>jonswap|pmstokes)_H(?P<height>[0-9.]+).*?_fusion_ou3\.csv$"), "*_fusion_ou3.csv"),
    "ou2": (re.compile(r".*?_(?P<wave>jonswap|pmstokes)_H(?P<height>[0-9.]+).*?_fusion_ou2\.csv$"), "*_fusion_ou2.csv"),
    "pii": (re.compile(r".*?_(?P<wave>jonswap|pmstokes)_H(?P<height>[0-9.]+).*?_nonkalman_fusion\.csv$"), "*_nonkalman_fusion.csv"),
    "nlo": (re.compile(r".*?_(?P<wave>jonswap|pmstokes)_H(?P<height>[0-9.]+).*?_tvg_nlo_nomag_nognss\.csv$"), "*_tvg_nlo_nomag_nognss.csv"),
}
LABELS = {"ou3": "OU--III", "ou2": "OU--II", "pii": "Adaptive PII", "nlo": "TVG--NLO"}
HEIGHTS = (0.27, 1.50, 4.00, 8.50)
WINDOW_S = 60.0


def index_files(pattern, regex):
    out = {}
    for name in glob.glob(pattern):
        match = regex.match(os.path.basename(name))
        if match:
            out[(match.group("wave"), round(float(match.group("height")), 2))] = name
    return out


def tail(df):
    if "time" not in df.columns or df.empty:
        raise ValueError("comparison CSV must contain a non-empty time column")
    end = float(df["time"].max())
    return df[df["time"] >= end - WINDOW_S].copy()


def rms(series):
    values = np.asarray(series, dtype=float)
    values = values[np.isfinite(values)]
    return math.nan if values.size == 0 else float(np.sqrt(np.mean(values * values)))


def metrics(path, hs):
    df = tail(pd.read_csv(path))
    err_z = df["disp_est_z"] - df["disp_ref_z"]
    result = {
        "z_rms": rms(err_z),
        "roll_rms": rms(df["roll_est"] - df["roll_ref"]),
        "pitch_rms": rms(df["pitch_est"] - df["pitch_ref"]),
    }
    result["z_pct"] = 100.0 * result["z_rms"] / hs
    return result, df


def tex_num(x, digits=2):
    return "--" if not math.isfinite(x) else f"{x:.{digits}f}"


def write_tables(rows, output):
    methods = ("pii", "nlo", "ou2", "ou3")
    lines = [
        r"\begin{table*}[t]",
        r"  \centering",
        r"  \caption{Paired vertical-displacement RMS error for the proposed OU--III estimator, its OU--II predecessor, the adaptive PII observer, and the published TVG--NLO baseline. All values use the same reference records and final \SI{60}{s} scoring window.}",
        r"  \label{tab:multi_observer_scenario_comparison}",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{4.0pt}",
        r"  \begin{tabular}{@{}llrrrrrrrr@{}}",
        r"    \toprule",
        r"    Spectrum & $H_s$ & \multicolumn{4}{c}{Z RMS [m]} & \multicolumn{4}{c}{Z RMS [$\%H_s$]} \\",
        r"    \cmidrule(lr){3-6}\cmidrule(lr){7-10}",
        r"    & & Adaptive PII & TVG--NLO & OU--II & OU--III & Adaptive PII & TVG--NLO & OU--II & OU--III \\",
        r"    \midrule",
    ]
    for row in rows:
        spectrum = "JONSWAP" if row["wave"] == "jonswap" else "PM--Stokes"
        lines.append(
            f"    {spectrum} & {row['hs']:.2f} & " +
            " & ".join(tex_num(row[m]["z_rms"], 3) for m in methods) + " & " +
            " & ".join(tex_num(row[m]["z_pct"], 1) for m in methods) + r" \\"
        )
    lines += [r"    \bottomrule", r"  \end{tabular}", r"\end{table*}", ""]

    aggregates = {}
    for method in methods:
        vals = [r[method] for r in rows]
        aggregates[method] = {
            "mean_z_pct": float(np.mean([v["z_pct"] for v in vals])),
            "best_z_pct": float(np.min([v["z_pct"] for v in vals])),
            "worst_z_pct": float(np.max([v["z_pct"] for v in vals])),
            "mean_roll": float(np.mean([v["roll_rms"] for v in vals])),
            "mean_pitch": float(np.mean([v["pitch_rms"] for v in vals])),
        }

    lines += [
        r"\begin{table}[t]",
        r"  \centering",
        r"  \caption{Aggregate comparison over all eight paired wave cases.}",
        r"  \label{tab:multi_observer_aggregate_comparison}",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{3.0pt}",
        r"  \begin{tabular}{@{}lrrrrr@{}}",
        r"    \toprule",
        r"    Method & Mean Z\%$H_s$ & Best & Worst & Mean roll & Mean pitch \\",
        r"    \midrule",
    ]
    for method in methods:
        a = aggregates[method]
        label = r"\textbf{OU--III}" if method == "ou3" else LABELS[method]
        lines.append(
            f"    {label} & {tex_num(a['mean_z_pct'],1)} & {tex_num(a['best_z_pct'],1)} & "
            f"{tex_num(a['worst_z_pct'],1)} & {tex_num(a['mean_roll'],2)} & {tex_num(a['mean_pitch'],2)} \\\\"
        )
    lines += [r"    \bottomrule", r"  \end{tabular}", r"\end{table}"]
    Path(output).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(dfs, output_base):
    trimmed = {name: tail(df) for name, df in dfs.items()}
    start = max(float(df["time"].min()) for df in trimmed.values())
    end = min(float(df["time"].max()) for df in trimmed.values())
    base = trimmed["ou3"]
    base = base[(base["time"] >= start) & (base["time"] <= end)]
    t = base["time"].to_numpy()
    ref = base["disp_ref_z"].to_numpy()

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.0), sharex=True)
    axes[0].plot(t, ref, label="Reference", linewidth=1.3)
    styles = {"ou3": "-", "ou2": "-.", "pii": ":", "nlo": "--"}
    for method in ("ou3", "ou2", "pii", "nlo"):
        df = trimmed[method]
        est = np.interp(t, df["time"], df["disp_est_z"])
        axes[0].plot(t, est, label=LABELS[method], linewidth=1.0, linestyle=styles[method])
        axes[1].plot(t, est - ref, label=f"{LABELS[method]} error", linewidth=0.9, linestyle=styles[method])
    axes[0].set_ylabel("Vertical displacement [m]")
    axes[0].grid(True)
    axes[0].legend(loc="upper right", fontsize=7, ncol=2)
    axes[1].set_ylabel("Error [m]")
    axes[1].set_xlabel("Time [s]")
    axes[1].grid(True)
    axes[1].legend(loc="upper right", fontsize=7, ncol=2)
    fig.tight_layout()
    for ext in ("pgf", "svg"):
        fig.savefig(f"{output_base}.{ext}", format=ext, bbox_inches="tight")
    plt.close(fig)


def main():
    files = {name: index_files(pattern, regex) for name, (regex, pattern) in PATTERNS.items()}
    expected = {(wave, round(hs, 2)) for wave in ("jonswap", "pmstokes") for hs in HEIGHTS}
    missing = {name: sorted(expected - set(index)) for name, index in files.items() if expected - set(index)}
    if missing:
        raise RuntimeError(f"missing comparison cases: {missing}")

    rows = []
    representative = {}
    for wave in ("jonswap", "pmstokes"):
        for hs in HEIGHTS:
            key = (wave, round(hs, 2))
            row = {"wave": wave, "hs": hs}
            for method in ("ou3", "ou2", "pii", "nlo"):
                row[method], df = metrics(files[method][key], hs)
                if key == ("jonswap", 1.50):
                    representative[method] = df
            rows.append(row)

    write_tables(rows, "w3d-baseline-results-generated.tex-part")
    if len(representative) != 4:
        raise RuntimeError("representative JONSWAP Hs=1.50 cases not found")
    write_plot(representative, "w3d_multi_observer_jonswap_medium")


if __name__ == "__main__":
    main()
