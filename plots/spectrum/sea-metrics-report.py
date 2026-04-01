#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def tex_escape(text: str) -> str:
    return text.replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")


def pass_cell(v: str) -> str:
    return r"\textbf{PASS}" if str(v).strip() in ("1", "true", "True") else r"\textbf{FAIL}"


def build_table(rows) -> str:
    lines = [
        r"\begin{table}[H]",
        r"\centering",
        r"\caption{SeaMetricsFromSpectrum validation using two synthetic spectra derived from spectrum-estimator scenarios.}",
        r"\label{tab:sea_metrics_from_spectrum}",
        r"\begin{tabular}{lrrrrrcl}",
        r"\toprule",
        r"Case & $H_{s,target}$ & $H_{s,est}$ & $H_s$ err [\%] & $T_{p,target}$ & $T_{p,est}$ & Type & Gate \\",
        r"\midrule",
    ]

    for row in rows:
        line = (
            f"{tex_escape(row['case_name'])} & "
            f"{float(row['hs_target_m']):.3f} & "
            f"{float(row['hs_est_m']):.3f} & "
            f"{float(row['hs_error_pct']):.2f} & "
            f"{float(row['tp_target_s']):.3f} & "
            f"{float(row['tp_est_s']):.3f} & "
            f"{tex_escape(row['spectrum_type'])} & "
            f"{pass_cell(row['pass'])} "
            r"\\"
        )
        lines.append(line)

    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])
    return "\n".join(lines) + "\n"


def build_fragment(rows) -> str:
    table_tex = build_table(rows).strip()
    return "\n".join([
        r"\documentclass[11pt,letterpaper]{article}",
        r"\usepackage{booktabs}",
        r"\usepackage{float}",
        r"\begin{document}",
        table_tex,
        r"\end{document}",
        "",
    ])


def build_main(rows) -> str:
    table_tex = build_table(rows).strip()
    return "\n".join([
        r"\documentclass[11pt,letterpaper]{article}",
        r"\usepackage{fullpage}",
        r"\usepackage{booktabs}",
        r"\usepackage{float}",
        r"\title{SeaMetricsFromSpectrum Test Report}",
        r"\author{Auto-generated}",
        r"\date{\today}",
        r"\begin{document}",
        r"\maketitle",
        table_tex,
        r"\end{document}",
        "",
    ])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", default="sea_metrics_from_spectrum_report.csv")
    parser.add_argument("--out-fragment", default="../../doc/spectrum/sea_metrics_from_spectrum_table_fragment.tex")
    parser.add_argument("--out-main", default="../../doc/spectrum/sea_metrics_from_spectrum_table.tex")
    args = parser.parse_args()

    csv_path = Path(args.csv)
    fragment_path = Path(args.out_fragment)
    main_path = Path(args.out_main)

    with csv_path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))

    fragment_path.parent.mkdir(parents=True, exist_ok=True)
    main_path.parent.mkdir(parents=True, exist_ok=True)
    fragment_path.write_text(build_fragment(rows), encoding="utf-8")
    main_path.write_text(build_main(rows), encoding="utf-8")

    print(f"Loaded CSV rows: {len(rows)} from {csv_path}")
    print(f"Wrote LaTeX fragment: {fragment_path}")
    print(f"Wrote LaTeX main: {main_path}")


if __name__ == "__main__":
    main()
