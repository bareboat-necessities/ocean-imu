#!/usr/bin/env python3
"""Summarize paired comparator sweeps without selecting or changing defaults."""
from __future__ import annotations

import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports/results/comparator_rao_retuning"


def summarize(rows):
    """Equal record weights; unsupported channels remain absent."""
    if "metrics" in rows[0]:
        names = ("disp_z_pct_hs", "disp_3d_rms_m", "roll_rms_deg",
                 "pitch_rms_deg", "yaw_rms_deg")
        result = {name: statistics.fmean(r["metrics"][name] for r in rows)
                  for name in names}
        result.update(records=len(rows), violations=sum(len(r["violations"]) for r in rows))
    else:
        result = {"disp_z_pct_hs": statistics.fmean(v for r in rows for v in r["z_pct_hs"])}
        for i, name in enumerate(("roll_rms_deg", "pitch_rms_deg", "yaw_rms_deg")):
            result[name] = statistics.fmean(a[i] for r in rows for a in r["angles_deg"])
        result.update(records=sum(len(r["z_pct_hs"]) for r in rows),
                      violations=sum(len(r["failures"]) for r in rows))
    return result


def main():
    summaries = []
    for path in sorted(OUT.glob("*/runs.json")):
        rows = json.loads(path.read_text())
        for config in sorted({r["config"] for r in rows}):
            arm = [r for r in rows if r["config"] == config]
            summaries.append(dict(study=path.parent.name, config=config, **summarize(arm)))
    (OUT / "sweep-summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
    lines = ["# Complete sweep results", "",
             "Equal-weight means over complete records; Z is percent of incident Hs. "
             "Violations count unchanged executable gates. NLO yaw is unscored. "
             "Training, refinement, and held-out stages remain separate.", "",
             "| Stage | Configuration | Records | Violations | Z, %Hs | 3-D, m | Roll, deg | Pitch, deg | Yaw, deg |",
             "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for row in summaries:
        def value(name, row=row):
            if name not in row or (row["study"].startswith("nlo-") and name == "yaw_rms_deg"):
                return "—"
            return f"{row[name]:.6f}"
        cells = [row["study"], row["config"], str(row["records"]), str(row["violations"])]
        cells += [value(k) for k in ("disp_z_pct_hs", "disp_3d_rms_m", "roll_rms_deg",
                                     "pitch_rms_deg", "yaw_rms_deg")]
        lines.append("| " + " | ".join(cells) + " |")
    (OUT / "sweep-results.md").write_text("\n".join(lines) + "\n")

    # Publication is allowed only after every family's paired holdout exists.
    pairs = []
    for family in ("ou2", "tfg", "nlo", "pii"):
        path = OUT / f"{family}-holdout/runs.json"
        if not path.exists():
            return
        rows = json.loads(path.read_text())
        arms = [[r for r in rows if r["config"] == c] for c in ("baseline", "candidate")]
        keys = [{(r["seed"], r.get("input", "all-eight")) for r in arm} for arm in arms]
        if not keys[0] or keys[0] != keys[1]:
            raise ValueError(f"{family}: incomplete paired holdout")
        pairs.append((family, *(summarize(arm) for arm in arms)))
    tex = [r"\begin{table}[t]", r"\centering\footnotesize",
           r"\caption{Comparator retuning on two additional paired sensor and initialization draws over eight fixed vessel-response records. Values are mean vertical RMS in percent of incident $H_s$; violations count all scored executable gates, whose sets differ between methods.}",
           r"\label{tab:comparator-retuning}", r"\begin{tabular}{@{}lrrrr@{}}\toprule",
           r"& \multicolumn{2}{c}{Z RMS, \%$H_s$} & \multicolumn{2}{c}{Violations} \\",
           r"Method & Before & After & Before & After \\\midrule"]
    labels = {"ou2": "OU--II", "tfg": "TFG", "nlo": "TVG--NLO", "pii": "PII"}
    for family, before, after in pairs:
        tex.append(f"{labels[family]} & {before['disp_z_pct_hs']:.3f} & "
                   f"{after['disp_z_pct_hs']:.3f} & {before['violations']} & {after['violations']} " + r"\\")
    tex += [r"\bottomrule\end{tabular}\end{table}"]
    (ROOT / "doc/kalman_ou_iii/w3d-comparator-retuning-generated.tex-part").write_text("\n".join(tex) + "\n")


if __name__ == "__main__":
    main()
