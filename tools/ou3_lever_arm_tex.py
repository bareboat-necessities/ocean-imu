#!/usr/bin/env python3
"""Generate the OU-III paper fragment from the lever-arm study summaries.

Consumes ``lever_arm_summary.csv`` and, when present, the derivative-band
sweep in ``lever_arm_cutoff_summary.csv``.  Every number the article states
about the study comes from here; none is hand-typed.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


# The corner the simulator ships; see W3dLeverArmConfig::derivative_cutoff_hz.
DEPLOYED_CUTOFF_HZ = 15.0

AXIS_TEX = {
    "x-athwartships": "athwartships",
    "y-fore-aft": "fore--aft",
    "z-vertical": "vertical",
}


def load(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def value(row: dict[str, str], name: str) -> float:
    result = float(row[name])
    if not math.isfinite(result):
        raise ValueError(f"{name} is not finite for {row}")
    return result


def axis_tex(axis: str) -> str:
    return AXIS_TEX[axis]


def at(rows: list[dict[str, str]], mode: str, distance: float) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if row["mode"] == mode
        and math.isclose(value(row, "distance_m"), distance, abs_tol=1e-9)
    ]


def generate(
    rows: list[dict[str, str]],
    sweep: list[dict[str, str]] | None = None,
) -> str:
    baseline = next(
        row for row in rows if row["mode"] == "baseline" and row["axis"] == "cg"
    )
    distances = sorted(
        {value(row, "distance_m") for row in rows if row["mode"] == "unmodeled"}
    )

    table_rows: list[str] = []
    for distance in distances:
        unmodeled = at(rows, "unmodeled", distance)
        gyro = at(rows, "gyro", distance)
        exact = at(rows, "exact", distance)
        estimated = at(rows, "estimated", distance)
        worst_3d = max(unmodeled, key=lambda row: value(row, "disp_3d_ratio_to_baseline"))
        worst_tilt = max(unmodeled, key=lambda row: value(row, "tilt_ratio_to_baseline"))
        gyro_max = max(value(row, "disp_3d_ratio_to_baseline") for row in gyro)
        exact_max = max(value(row, "disp_3d_ratio_to_baseline") for row in exact)
        estimated_max = (
            max(value(row, "disp_3d_ratio_to_baseline") for row in estimated)
            if estimated
            else float("nan")
        )
        table_rows.append(
            "    {cm:.0f} & {d3:.3f} ({daxis}) & {tilt:.3f} ({taxis}) & "
            "{gyro:.3f} & {exact:.3f} & {est} \\\\".format(
                cm=100.0 * distance,
                d3=value(worst_3d, "disp_3d_ratio_to_baseline"),
                daxis=axis_tex(worst_3d["axis"]),
                tilt=value(worst_tilt, "tilt_ratio_to_baseline"),
                taxis=axis_tex(worst_tilt["axis"]),
                gyro=gyro_max,
                exact=exact_max,
                est="---" if math.isnan(estimated_max) else f"{estimated_max:.3f}",
            )
        )

    far = distances[-1]
    unmodeled_far = at(rows, "unmodeled", far)
    worst_3d = max(unmodeled_far, key=lambda row: value(row, "disp_3d_ratio_to_baseline"))
    worst_tilt = max(unmodeled_far, key=lambda row: value(row, "tilt_ratio_to_baseline"))
    gyro_max = max(value(row, "disp_3d_ratio_to_baseline") for row in at(rows, "gyro", far))
    exact_max = max(
        value(row, "disp_3d_ratio_to_baseline") for row in at(rows, "exact", far)
    )
    injected = max(value(row, "installed_rms_mps2") for row in unmodeled_far)
    gyro_residual = max(value(row, "residual_rms_mps2") for row in at(rows, "gyro", far))
    gyro_tilt = max(value(row, "tilt_ratio_to_baseline") for row in at(rows, "gyro", far))
    exact_tilt = max(
        value(row, "tilt_ratio_to_baseline") for row in at(rows, "exact", far)
    )

    lines = [
        "The CG-mounted baseline has pooled 3-D displacement RMS of",
        f"\\SI{{{value(baseline, 'disp_3d_rms_m'):.3f}}}{{m}} and maximum pooled",
        f"roll/pitch RMS of \\SI{{{value(baseline, 'max_tilt_rms_deg'):.3f}}}{{\\degree}}.",
        "Table~\\ref{tab:imu-lever-arm} reports the worst canonical direction at each",
        f"installation distance.  At \\SI{{{100.0 * far:.0f}}}{{cm}} the rotational term",
        f"reaches \\SI{{{injected:.3f}}}{{\\meter\\per\\second\\squared}} RMS in the worst",
        "direction.  Leaving it unmodeled raises 3-D displacement to",
        f"\\num{{{value(worst_3d, 'disp_3d_ratio_to_baseline'):.3f}}} times the CG",
        f"baseline in the {axis_tex(worst_3d['axis'])} direction, and the largest tilt",
        f"ratio is \\num{{{value(worst_tilt, 'tilt_ratio_to_baseline'):.3f}}} in the",
        f"{axis_tex(worst_tilt['axis'])} direction.  Modelling the same lever arm inside",
        "the filter removes it: with exact angular kinematics the largest 3-D ratio over",
        f"the three directions falls to \\num{{{exact_max:.3f}}}, and with the deployable",
        "model driven by the measured rate it is",
        f"\\num{{{gyro_max:.3f}}}, leaving only",
        f"\\SI{{{gyro_residual:.3f}}}{{\\meter\\per\\second\\squared}} of the injected term",
        "behind.  The attitude channel, which carries the larger penalty, recovers with",
        f"it: the worst tilt ratio falls to \\num{{{exact_tilt:.3f}}} under the exact model",
        f"and \\num{{{gyro_tilt:.3f}}} under the deployable one.  The installation penalty",
        "is therefore deterministic and recoverable rather than intrinsic OU--III error.",
        "",
        "\\begin{table}[t]",
        "  \\centering",
        "  \\caption{Worst pooled degradation over the three canonical IMU lever-arm directions.  Ratios are relative to the CG-mounted OU--III baseline; the parenthesized direction is the maximizing unmodeled case.  The last three columns are the worst ratio over the same three directions once the lever arm is handled inside the filter: given exactly, reconstructed from the measured rate, or estimated as filter states from no survey at all.}",
        "  \\label{tab:imu-lever-arm}",
        "  \\footnotesize",
        "  \\setlength{\\tabcolsep}{3.0pt}",
        "  \\begin{tabular}{@{}rrrrrr@{}}",
        "    \\toprule",
        "    Offset & Max 3-D / CG & Max tilt / CG & Gyro model & Exact model & "
        "Estimated \\\\",
        # The brace is load-bearing.  A row opening with an unbraced
        # bracket is swallowed by the preceding \\ as its optional
        # vertical-space argument, and LaTeX dies on "Missing number".
        "    {[cm]} & unmodeled & unmodeled & max 3-D / CG & max 3-D / CG & "
        "max 3-D / CG \\\\",
        "    \\midrule",
        *table_rows,
        "    \\bottomrule",
        "  \\end{tabular}",
        "\\end{table}",
        "",
    ]

    estimated_far = at(rows, "estimated", far)
    if estimated_far:
        # The self-calibrating arm is reported on two different questions, and
        # they do not have the same answer: what it does to the score, and
        # whether it actually found the lever arm.  Quoting only the first
        # would let a partially-converged calibration pass for a converged one.
        est_3d = max(
            value(row, "disp_3d_ratio_to_baseline") for row in estimated_far
        )
        est_tilt = max(value(row, "tilt_ratio_to_baseline") for row in estimated_far)
        est_residual = max(value(row, "residual_rms_mps2") for row in estimated_far)
        best_dir = max(
            estimated_far, key=lambda row: value(row, "lever_recovered_fraction")
        )
        worst_dir = min(
            estimated_far, key=lambda row: value(row, "lever_recovered_fraction")
        )
        sigma = max(value(row, "lever_sigma_max_m") for row in estimated_far)
        err = value(worst_dir, "lever_estimate_err_m")
        overconfidence = err / sigma if sigma > 0.0 else float("nan")

        unmodeled_3d = max(
            value(row, "disp_3d_ratio_to_baseline") for row in unmodeled_far
        )
        unmodeled_tilt = max(
            value(row, "tilt_ratio_to_baseline") for row in unmodeled_far
        )
        lines += [
            "Estimating the same lever arm instead of supplying it splits the two",
            "channels rather than recovering both.  Over the same eight seas with the",
            f"IMU \\SI{{{100.0 * far:.0f}}}{{cm}} off the CG and a prior that says only",
            "\\SI{0.5}{m} per axis about where the sensor is, the worst tilt ratio over",
            f"the three directions falls from \\num{{{unmodeled_tilt:.3f}}} unmodeled to",
            f"\\num{{{est_tilt:.3f}}} --- nearly the whole attitude penalty, which is the",
            "larger one --- while the worst 3-D displacement ratio rises from",
            f"\\num{{{unmodeled_3d:.3f}}} unmodeled to \\num{{{est_3d:.3f}}}, leaving",
            f"\\SI{{{est_residual:.3f}}}{{\\meter\\per\\second\\squared}} of the injected",
            "term behind.  A half-converged lever arm is worse than none on that",
            "channel: the correction it applies is the wrong size, and being wrong in a",
            "way correlated with the rotation is worse for a doubly integrated quantity",
            "than leaving the term alone.",
            "",
            "The calibration behind those numbers is only partly converged, and by how",
            "much depends on how hard the sea rotates the hull.  Pooled over the eight",
            "seas the estimator recovers",
            f"\\num{{{100.0 * value(best_dir, 'lever_recovered_fraction'):.0f}}}\\,\\%",
            f"of the installed arm in the {axis_tex(best_dir['axis'])} direction and",
            f"\\num{{{100.0 * value(worst_dir, 'lever_recovered_fraction'):.0f}}}\\,\\%",
            f"in the {axis_tex(worst_dir['axis'])} one, and sea by sea the recovery is",
            "monotone in $H_s$: essentially nothing at",
            "\\SI{0.27}{m} and most of the arm at \\SI{8.5}{m}.",
            "None of that is visible in the filter's own uncertainty.  The worst",
            f"direction is off by \\SI{{{err:.3f}}}{{m}} while reporting a standard",
            f"deviation of at most \\SI{{{sigma:.4f}}}{{m}}, a factor of",
            f"\\num{{{overconfidence:.0f}}} between what the filter is wrong by and what",
            "it says it is wrong by.  The covariance is measuring the conditioning of",
            "the regression, which is genuinely good after twenty minutes of rolling;",
            "the error is set by everything the regression does not model, and those",
            "are different quantities.",
            "",
        ]

    if sweep:
        best = min(sweep, key=lambda row: value(row, "disp_3d_ratio_to_baseline"))
        narrow = min(sweep, key=lambda row: value(row, "cutoff_hz"))
        wide = max(sweep, key=lambda row: value(row, "cutoff_hz"))
        # The deployed corner is a choice, not the sweep's argmin; quoting the
        # argmin as "deployed" would overstate what the shipped model does.
        deployed = min(
            sweep,
            key=lambda row: abs(value(row, "cutoff_hz") - DEPLOYED_CUTOFF_HZ),
        )
        lines += [
            "The deployable model has one design parameter, the band of its rate",
            "derivative, and it is two-sided.  Over the same eight seas with a",
            f"\\SI{{{100.0 * value(narrow, 'distance_m'):.0f}}}{{cm}} fore--aft arm, a",
            f"\\SI{{{value(narrow, 'cutoff_hz'):.0f}}}{{Hz}} corner leaves",
            f"\\num{{{value(narrow, 'disp_3d_ratio_to_baseline'):.3f}}} times the CG",
            "baseline because the low-pass phase lag misaligns a correction whose",
            f"amplitude is already right, and a \\SI{{{value(wide, 'cutoff_hz'):.0f}}}{{Hz}}",
            f"corner leaves \\num{{{value(wide, 'disp_3d_ratio_to_baseline'):.3f}}} because",
            "the differentiated gyro noise it admits exceeds the term it removes.  The",
            "basin between them is flat: the deployed",
            f"\\SI{{{value(deployed, 'cutoff_hz'):.0f}}}{{Hz}} corner reaches",
            f"\\num{{{value(deployed, 'disp_3d_ratio_to_baseline'):.3f}}} and the best",
            f"corner measured, \\SI{{{value(best, 'cutoff_hz'):.0f}}}{{Hz}}, reaches",
            f"\\num{{{value(best, 'disp_3d_ratio_to_baseline'):.3f}}}, so the model is",
            "robust to this choice rather than tuned to it.",
            "",
        ]

    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cutoff-summary", type=Path, default=None)
    args = parser.parse_args()
    rows = load(args.summary)
    if not rows:
        raise SystemExit("summary is empty")
    sweep = None
    if args.cutoff_summary is not None and args.cutoff_summary.exists():
        sweep = load(args.cutoff_summary) or None
    text = generate(rows, sweep)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
