#!/usr/bin/env python3
"""Generate the compact OU-III paper fragment from a lever-arm summary CSV."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def load(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def value(row: dict[str, str], name: str) -> float:
    result = float(row[name])
    if not math.isfinite(result):
        raise ValueError(f"{name} is not finite for {row}")
    return result


def axis_tex(axis: str) -> str:
    return {
        "x-athwartships": "athwartships",
        "y-fore-aft": "fore--aft",
        "z-vertical": "vertical",
    }[axis]


def generate(rows: list[dict[str, str]]) -> str:
    baseline = next(
        row for row in rows
        if row["mode"] == "baseline" and row["axis"] == "cg"
    )
    distances = sorted(
        {
            value(row, "distance_m")
            for row in rows
            if row["mode"] == "unmodeled"
        }
    )
    table_rows: list[str] = []
    worst_at_max: tuple[dict[str, str], dict[str, str], float] | None = None

    for distance in distances:
        unmodeled = [
            row for row in rows
            if row["mode"] == "unmodeled"
            and math.isclose(value(row, "distance_m"), distance, abs_tol=1e-9)
        ]
        exact = [
            row for row in rows
            if row["mode"] == "exact"
            and math.isclose(value(row, "distance_m"), distance, abs_tol=1e-9)
        ]
        worst_3d = max(unmodeled, key=lambda row: value(row, "disp_3d_ratio_to_baseline"))
        worst_tilt = max(unmodeled, key=lambda row: value(row, "tilt_ratio_to_baseline"))
        exact_max = max(value(row, "disp_3d_ratio_to_baseline") for row in exact)
        table_rows.append(
            "    {cm:.0f} & {d3:.3f} ({daxis}) & {tilt:.3f} ({taxis}) & {exact:.3f} \\\\".format(
                cm=100.0 * distance,
                d3=value(worst_3d, "disp_3d_ratio_to_baseline"),
                daxis=axis_tex(worst_3d["axis"]),
                tilt=value(worst_tilt, "tilt_ratio_to_baseline"),
                taxis=axis_tex(worst_tilt["axis"]),
                exact=exact_max,
            )
        )
        if distance == distances[-1]:
            worst_at_max = (worst_3d, worst_tilt, exact_max)

    assert worst_at_max is not None
    worst_3d, worst_tilt, exact_max = worst_at_max
    return "\n".join(
        [
            "The CG-mounted baseline has pooled 3-D displacement RMS of",
            f"\\SI{{{value(baseline, 'disp_3d_rms_m'):.3f}}}{{m}} and maximum pooled",
            f"roll/pitch RMS of \\SI{{{value(baseline, 'max_tilt_rms_deg'):.3f}}}{{\\degree}}.",
            "Table~\\ref{tab:imu-lever-arm} reports the worst canonical direction at",
            "each installation distance.  At \\SI{30}{cm}, leaving the lever arm",
            f"unmodeled raises 3-D displacement to \\num{{{value(worst_3d, 'disp_3d_ratio_to_baseline'):.3f}}}",
            f"times the CG baseline in the {axis_tex(worst_3d['axis'])} direction, while",
            f"the largest tilt ratio is \\num{{{value(worst_tilt, 'tilt_ratio_to_baseline'):.3f}}}",
            f"in the {axis_tex(worst_tilt['axis'])} direction.  With the exact lever-arm",
            f"model, the largest 3-D ratio over the same three directions is",
            f"\\num{{{exact_max:.3f}}}, showing how much of the installation penalty is",
            "deterministic and recoverable rather than intrinsic OU--III error.",
            "",
            "\\begin{table}[t]",
            "  \\centering",
            "  \\caption{Worst pooled degradation over the three canonical IMU lever-arm directions.  Ratios are relative to the CG-mounted OU--III baseline; the parenthesized direction is the maximizing unmodeled case.}",
            "  \\label{tab:imu-lever-arm}",
            "  \\footnotesize",
            "  \\setlength{\\tabcolsep}{3.0pt}",
            "  \\begin{tabular}{@{}rrrr@{}}",
            "    \\toprule",
            "    Offset [cm] & Max 3-D / CG & Max tilt / CG & Exact max 3-D / CG \\\\",
            "    \\midrule",
            *table_rows,
            "    \\bottomrule",
            "  \\end{tabular}",
            "\\end{table}",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = load(args.summary)
    if not rows:
        raise SystemExit("summary is empty")
    text = generate(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
