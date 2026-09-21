#!/usr/bin/env python3
import glob
import hashlib
import json
import subprocess
import sys
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
    "tfg": (re.compile(r".*?_(?P<wave>jonswap|pmstokes)_H(?P<height>[0-9.]+).*?_fusion_tfg\.csv$"), "*_fusion_tfg.csv"),
    "ou2": (re.compile(r".*?_(?P<wave>jonswap|pmstokes)_H(?P<height>[0-9.]+).*?_fusion_ou2\.csv$"), "*_fusion_ou2.csv"),
    "pii": (re.compile(r".*?_(?P<wave>jonswap|pmstokes)_H(?P<height>[0-9.]+).*?_nonkalman_fusion\.csv$"), "*_nonkalman_fusion.csv"),
    "nlo": (re.compile(r".*?_(?P<wave>jonswap|pmstokes)_H(?P<height>[0-9.]+).*?_tvg_nlo_nomag_nognss\.csv$"), "*_tvg_nlo_nomag_nognss.csv"),
}
LABELS = {"ou3": "OU--III", "ou2": "OU--II", "tfg": "TFG", "pii": "Adaptive PII", "nlo": "TVG--NLO"}
HEIGHTS = (0.27, 1.50, 4.00, 8.50)
REPO_ROOT = Path(__file__).resolve().parents[2]
# Scoring window for the tables: the same trailing window the simulators gate
# on and the ten-seed study scores.
WINDOW_S = 900.0
# The time-series figure stays at a minute.  Fifteen minutes of wave traces in
# a two-panel figure is a solid band, not a reconstruction anyone can read.
PLOT_WINDOW_S = 60.0


def index_files(pattern, regex):
    out = {}
    for name in glob.glob(pattern):
        match = regex.match(os.path.basename(name))
        if match:
            out[(match.group("wave"), round(float(match.group("height")), 2))] = name
    return out


def tail(df, window_s=WINDOW_S):
    if "time" not in df.columns or df.empty:
        raise ValueError("comparison CSV must contain a non-empty time column")
    end = float(df["time"].max())
    return df[df["time"] >= end - window_s].copy()


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
    methods = ("pii", "nlo", "ou2", "tfg", "ou3")
    lines = [
        r"\begin{table*}[t]",
        r"  \centering",
        r"  \caption{Common-record deterministic vertical-displacement RMS error for the proposed OU--III estimator, its OU--II predecessor, the adaptive PII observer, the adapted TVG--NLO baseline, and the two-frame Lie-group filter (TFG). All values use the same reference records and final \SI{900}{s} scoring window; this is not the ten-seed paired inferential comparison.}",
        r"  \label{tab:multi_observer_scenario_comparison}",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{3.0pt}",
        r"  \begin{tabular}{@{}llrrrrrrrrrr@{}}",
        r"    \toprule",
        r"    Spectrum & $H_s$ & \multicolumn{5}{c}{Z RMS [m]} & \multicolumn{5}{c}{Z RMS [$\%H_s$]} \\",
        r"    \cmidrule(lr){3-7}\cmidrule(lr){8-12}",
        r"    & & PII & TVG--NLO & OU--II & TFG & OU--III & PII & TVG--NLO & OU--II & TFG & OU--III \\",
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

    macro_names = {
        "pii": "OUBaselinePIIMeanVerticalPercent",
        "nlo": "OUBaselineTVGNLOMeanVerticalPercent",
        "ou2": "OUBaselineOUIIMeanVerticalPercent",
        "tfg": "OUBaselineTFGMeanVerticalPercent",
        "ou3": "OUBaselineOUIIIMeanVerticalPercent",
    }
    for method in methods:
        macro = macro_names[method]
        value = tex_num(aggregates[method]["mean_z_pct"], 1)
        # provide+renew works both when the article installed a fallback and
        # when this generated file is compiled on its own.
        lines += [
            rf"\providecommand{{\{macro}}}{{}}",
            rf"\renewcommand{{\{macro}}}{{{value}}}",
        ]
    lines.append("")

    lines += [
        r"\begin{table}[t]",
        r"  \centering",
        r"  \caption{Aggregate deterministic comparison over all eight common-record wave cases.}",
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
    trimmed = {name: tail(df, PLOT_WINDOW_S) for name, df in dfs.items()}
    start = max(float(df["time"].min()) for df in trimmed.values())
    end = min(float(df["time"].max()) for df in trimmed.values())
    base = trimmed["ou3"]
    base = base[(base["time"] >= start) & (base["time"] <= end)]
    t = base["time"].to_numpy()
    ref = base["disp_ref_z"].to_numpy()

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.0), sharex=True)
    axes[0].plot(t, ref, label="Reference", linewidth=1.3)
    styles = {"ou3": "-", "ou2": "-.", "pii": ":", "nlo": "--", "tfg": (0, (3, 1, 1, 1))}
    for method in ("ou3", "ou2", "pii", "nlo", "tfg"):
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


def write_provenance(files, rows):
    """Retain source and time-series identities alongside the generated tables."""
    sys.path.insert(0, str(REPO_ROOT / "tools"))
    import ou_evidence_provenance as provenance
    import sim_dataset

    def digest(path):
        h = hashlib.sha256()
        with Path(path).open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                h.update(block)
        return h.hexdigest()

    roots = [REPO_ROOT / "src/util/W3dSimCommon.cpp"] + [
        REPO_ROOT / "tests" / family / name for family, name in (
            ("kalman_ou_iii", "kalman_ou_iii-sim.cpp"),
            ("kalman_ou_ii", "kalman_ou_ii-sim.cpp"),
            ("kalman_tfg", "kalman_tfg-sim.cpp"),
            ("pii_observer", "pii_observer-adaptive.cpp"), ("nlo", "nlo-sim.cpp"))]
    sources = set(provenance.implementation_closure(roots))
    sources.update(p.parent / "Makefile" for p in roots[1:])
    sources.update([Path(__file__).resolve(), Path(__file__).with_name("draw_plots.sh").resolve()])
    output = REPO_ROOT / "reports/results/baseline_comparison"
    output.mkdir(parents=True, exist_ok=True)
    result = {"schema_version": 1, "protocol": "deterministic common-record, final 900 s; not Monte Carlo",
        "git_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True).strip(),
        "workflow": provenance.workflow_metadata(),
        "build_environment": provenance.environment_metadata(),
        "dataset": sim_dataset.input_provenance([Path(n) for n in sim_dataset.REFERENCE_CSV_SHA256]),
        "source_sha256": {str(p.relative_to(REPO_ROOT)): digest(p) for p in sorted(sources)},
        "timeseries_sha256": {m: {Path(p).name: digest(p) for p in index.values()} for m, index in files.items()},
        "rows": rows, "table_sha256": digest("w3d-baseline-results-generated.tex-part")}
    (output / "baseline_comparison.json").write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")


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
            for method in ("ou3", "ou2", "pii", "nlo", "tfg"):
                row[method], df = metrics(files[method][key], hs)
                if key == ("jonswap", 1.50):
                    representative[method] = df
            rows.append(row)

    write_tables(rows, "w3d-baseline-results-generated.tex-part")
    if len(representative) != len(PATTERNS):
        raise RuntimeError("representative JONSWAP Hs=1.50 cases not found")
    write_plot(representative, "w3d_multi_observer_jonswap_medium")
    write_provenance(files, rows)


if __name__ == "__main__":
    main()
