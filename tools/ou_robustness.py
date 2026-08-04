#!/usr/bin/env python3
"""Ten-seed OU-III sensitivity and degradation-case evaluation.

The study complements ``ou_validation.py`` without changing its confirmatory
OU-II/OU-III comparison.  It varies one fixed OU-III parameter at a time around
the nominal noise-free operating point and separately evaluates a low-motion
noise-floor case and a rapid sea-state transition.  Every comparison is paired
by wave-realization, IMU-noise, and initialization seed.
"""

from __future__ import annotations

import argparse
import math
import os
import platform
import re
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any, Mapping, Sequence

os.environ.setdefault("MPLCONFIGDIR", "/tmp/ocean-imu-matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import ou_validation as core  # noqa: E402

matplotlib.rcParams["svg.hashsalt"] = "ocean-imu-ou-robustness"


REPO_ROOT = Path(__file__).resolve().parents[1]
# One-factor-at-a-time parameters, followed by the two coupled sweeps that
# move r_S the way the deployed adaptation law does
# (r_S = clip(0.35 sigma_aw tau^3, 0.4, 400)).  Freezing r_S while sweeping
# tau or sigma_aw measures only the direct OU process-covariance effect, which
# is not what the online tuner ever does.
OFAT_PARAMETERS = ("tau", "sigma_aw", "r_s")
COUPLED_PARAMETERS = ("sigma_aw_rs", "tau_rs")
SENSITIVITY_PARAMETERS = OFAT_PARAMETERS + COUPLED_PARAMETERS
FULL_SCALES = (0.5, 0.75, 1.0, 1.25, 1.5)
SMOKE_SCALES = (0.75, 1.0, 1.25)
LOW_REFERENCE_HS_M = 0.27
LOW_STRESS_HS_M = 0.05
# Mirrors the implementation bounds in SeaStateFusionFilter_OU_III.h.  tau is
# now tied to the zero-crossing wave period rather than to an acceleration-band
# frequency, so both it and the r_S it cubes reach much further.
TAU_BOUNDS_S = (0.02, 12.0)
SIGMA_AW_BOUNDS_MPS2 = (1e-9, 6.0)
R_S_BOUNDS_MS = (0.4, 400.0)


@dataclass(frozen=True)
class TransitionCase:
    name: str
    start_sec: float
    end_sec: float

    @property
    def duration_sec(self) -> float:
        return self.end_sec - self.start_sec


def parse_float_list(text: str) -> list[float]:
    values = [float(token.strip()) for token in text.split(",") if token.strip()]
    if not values or any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise argparse.ArgumentTypeError("scale lists require positive finite values")
    if not any(math.isclose(value, 1.0, rel_tol=0.0, abs_tol=1e-12) for value in values):
        raise argparse.ArgumentTypeError("scale lists must include the 1.0 reference")
    return values


def scaled_tuning_point(
    baseline: core.TuningPoint,
    parameter: str,
    scale: float,
) -> core.TuningPoint:
    if parameter not in SENSITIVITY_PARAMETERS:
        raise ValueError(f"unknown sensitivity parameter: {parameter}")
    if not (math.isfinite(scale) and scale > 0.0):
        raise ValueError("sensitivity scale must be positive and finite")
    if baseline.RS_ms is None:
        raise ValueError("OU-III sensitivity requires an r_S tuning value")
    if parameter == "tau":
        return replace(baseline, tau_s=baseline.tau_s * scale)
    if parameter == "sigma_aw":
        return replace(baseline, sigma_a_mps2=baseline.sigma_a_mps2 * scale)
    if parameter == "r_s":
        return replace(baseline, RS_ms=baseline.RS_ms * scale)
    if parameter == "sigma_aw_rs":
        # r_S is linear in sigma_aw under the deployed law.
        return replace(
            baseline,
            sigma_a_mps2=baseline.sigma_a_mps2 * scale,
            RS_ms=baseline.RS_ms * scale,
        )
    # tau_rs: r_S is cubic in tau under the deployed law.
    return replace(
        baseline,
        tau_s=baseline.tau_s * scale,
        RS_ms=baseline.RS_ms * scale**3,
    )


def validate_tuning_point(point: core.TuningPoint) -> None:
    if point.RS_ms is None:
        raise ValueError("OU-III tuning point lacks r_S")
    for label, value, bounds in (
        ("tau", point.tau_s, TAU_BOUNDS_S),
        ("sigma_aw", point.sigma_a_mps2, SIGMA_AW_BOUNDS_MPS2),
        ("r_s", point.RS_ms, R_S_BOUNDS_MS),
    ):
        if not bounds[0] <= value <= bounds[1]:
            raise ValueError(
                f"{label}={value:.6g} lies outside implementation bounds "
                f"[{bounds[0]:.6g}, {bounds[1]:.6g}]"
            )


def summarize_rows(
    rows: Sequence[Mapping[str, Any]],
    bootstrap_resamples: int,
    stats_seed: int,
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, str, str, str], list[Mapping[str, Any]]] = {}
    for row in rows:
        key = tuple(
            str(row[name])
            for name in ("experiment", "case", "parameter", "scale_label", "mode")
        )
        grouped.setdefault(key, []).append(row)

    rng = np.random.default_rng(stats_seed)
    summary: list[dict[str, Any]] = []
    for key, group in sorted(grouped.items()):
        experiment, case, parameter, scale_label, mode = key
        for metric in core.NON_SEGMENT_METRIC_NAMES:
            values = core._finite_values(group, metric)
            n = int(values.size)
            mean = float(np.mean(values)) if n else math.nan
            std = (
                float(np.std(values, ddof=1))
                if n > 1
                else 0.0 if n == 1 else math.nan
            )
            bootstrap_low, bootstrap_high = core._bootstrap_mean_ci(
                values, bootstrap_resamples, rng
            )
            summary.append(
                {
                    "experiment": experiment,
                    "case": case,
                    "parameter": parameter,
                    "scale_label": scale_label,
                    "mode": mode,
                    "metric": metric,
                    "n": n,
                    "mean": mean,
                    "std": std,
                    "bootstrap_ci95_low": bootstrap_low,
                    "bootstrap_ci95_high": bootstrap_high,
                }
            )
    return summary


def paired_effect_rows(
    rows: Sequence[Mapping[str, Any]],
    bootstrap_resamples: int,
    stats_seed: int,
) -> list[dict[str, Any]]:
    pair_keys = ("wave_phase_seed", "imu_noise_seed", "initialization_seed")
    rng = np.random.default_rng(stats_seed + 1)
    effects: list[dict[str, Any]] = []

    def selected(**conditions: Any) -> list[Mapping[str, Any]]:
        return [
            row for row in rows
            if all(row.get(name) == value for name, value in conditions.items())
        ]

    def append_comparison(
        comparison: str,
        left: Sequence[Mapping[str, Any]],
        right: Sequence[Mapping[str, Any]],
        context: Mapping[str, Any],
    ) -> None:
        for metric in core.NON_SEGMENT_METRIC_NAMES:
            effect = core._paired_effect(
                left,
                right,
                metric,
                pair_keys,
                bootstrap_resamples,
                rng,
            )
            if effect is not None:
                effects.append(
                    {
                        **context,
                        "comparison": comparison,
                        "metric": metric,
                        **effect,
                    }
                )

    sensitivity = selected(experiment="sensitivity")
    scales = sorted({float(row["scale_multiplier"]) for row in sensitivity})
    for parameter in SENSITIVITY_PARAMETERS:
        reference = selected(
            experiment="sensitivity", parameter=parameter, scale_multiplier=1.0
        )
        for scale in scales:
            if math.isclose(scale, 1.0, rel_tol=0.0, abs_tol=1e-12):
                continue
            candidate = selected(
                experiment="sensitivity", parameter=parameter, scale_multiplier=scale
            )
            append_comparison(
                f"{parameter}_x{scale:g}_minus_x1",
                candidate,
                reference,
                {
                    "experiment": "sensitivity",
                    "parameter": parameter,
                    "left_case": f"x{scale:g}",
                    "right_case": "x1",
                },
            )

    append_comparison(
        "Hs0.05_minus_Hs0.27",
        selected(experiment="low_motion", case="Hs0.05"),
        selected(experiment="low_motion", case="Hs0.27"),
        {
            "experiment": "low_motion",
            "parameter": "",
            "left_case": "Hs0.05",
            "right_case": "Hs0.27",
        },
    )

    for mode in ("Adaptive", "FixedNominal"):
        append_comparison(
            f"rapid_minus_controlled_{mode}",
            selected(experiment="transition_rate", case="rapid", mode=mode),
            selected(experiment="transition_rate", case="controlled", mode=mode),
            {
                "experiment": "transition_rate",
                "parameter": "",
                "left_case": f"rapid/{mode}",
                "right_case": f"controlled/{mode}",
            },
        )
    for case in ("controlled", "rapid"):
        append_comparison(
            f"Adaptive_minus_FixedNominal_{case}",
            selected(experiment="transition_rate", case=case, mode="Adaptive"),
            selected(experiment="transition_rate", case=case, mode="FixedNominal"),
            {
                "experiment": "transition_rate",
                "parameter": "",
                "left_case": f"{case}/Adaptive",
                "right_case": f"{case}/FixedNominal",
            },
        )
    return effects


def _summary_index(
    summary: Sequence[Mapping[str, Any]],
) -> dict[tuple[str, str, str, str], Mapping[str, Any]]:
    return {
        (
            str(row["experiment"]),
            str(row["case"]),
            str(row["mode"]),
            str(row["metric"]),
        ): row
        for row in summary
    }


def _effect_index(
    effects: Sequence[Mapping[str, Any]],
) -> dict[tuple[str, str, str], Mapping[str, Any]]:
    return {
        (str(row["experiment"]), str(row["comparison"]), str(row["metric"])): row
        for row in effects
    }


def write_sensitivity_plot(
    path: Path,
    summary: Sequence[Mapping[str, Any]],
) -> None:
    labels = {
        "tau": r"$\tau$ only",
        "sigma_aw": r"$\sigma_{aw}$ only",
        "r_s": r"$r_S$ only ($R_S=\mathrm{diag}(r_S^2)$)",
        "sigma_aw_rs": r"$\sigma_{aw}$, $r_S$ coupled",
        "tau_rs": r"$\tau$, $r_S\propto\tau^3$ coupled",
    }
    colors = {
        "tau": "#0072B2",
        "sigma_aw": "#D55E00",
        "r_s": "#009E73",
        "sigma_aw_rs": "#CC79A7",
        "tau_rs": "#56B4E9",
    }
    styles = {
        "tau": "-",
        "sigma_aw": "-",
        "r_s": "-",
        "sigma_aw_rs": "--",
        "tau_rs": "--",
    }
    metrics = (
        ("disp_z_pct_hs", "Vertical reconstruction", r"RMS error (\% $H_s$)"),
        ("disp_3d_rms_m", "Full 3D displacement", "RMS error (m)"),
        ("pitch_rms_deg", "Pitch", "RMS error (deg)"),
    )
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.35))
    handles = []
    for axis, (metric, title, ylabel) in zip(axes, metrics):
        for parameter in SENSITIVITY_PARAMETERS:
            rows = sorted(
                (
                    row for row in summary
                    if row["experiment"] == "sensitivity"
                    and row["parameter"] == parameter
                    and row["metric"] == metric
                ),
                key=lambda row: float(row["scale_label"]),
            )
            x = np.asarray([float(row["scale_label"]) for row in rows])
            y = np.asarray([float(row["mean"]) for row in rows])
            low = y - np.asarray(
                [float(row["bootstrap_ci95_low"]) for row in rows]
            )
            high = np.asarray(
                [float(row["bootstrap_ci95_high"]) for row in rows]
            ) - y
            handle = axis.errorbar(
                x,
                y,
                yerr=np.vstack((low, high)),
                marker="o",
                linewidth=1.6,
                capsize=2.5,
                color=colors[parameter],
                linestyle=styles[parameter],
                label=labels[parameter],
            )
            if axis is axes[0]:
                handles.append(handle)
        axis.axvline(1.0, color="#555555", linewidth=0.9, linestyle="--")
        axis.set_xlabel("Nominal multiplier")
        axis.set_ylabel(ylabel)
        axis.set_title(title)
        axis.grid(True, alpha=0.25)
    fig.legend(
        handles,
        [labels[parameter] for parameter in SENSITIVITY_PARAMETERS],
        frameon=False,
        ncol=3,
        fontsize=8,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.01),
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.89))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format="svg", metadata={"Date": None})
    plt.close(fig)


def write_stress_plot(
    path: Path,
    summary: Sequence[Mapping[str, Any]],
) -> None:
    index = _summary_index(summary)
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.5))

    low_cases = ("Hs0.27", "Hs0.05")
    low_values = [
        float(index[("low_motion", case, "Adaptive", "disp_z_pct_hs")]["mean"])
        for case in low_cases
    ]
    low_error = [
        (
            value - float(index[("low_motion", case, "Adaptive", "disp_z_pct_hs")]["bootstrap_ci95_low"]),
            float(index[("low_motion", case, "Adaptive", "disp_z_pct_hs")]["bootstrap_ci95_high"]) - value,
        )
        for case, value in zip(low_cases, low_values)
    ]
    axes[0].bar(
        range(2),
        low_values,
        yerr=np.asarray(low_error).T,
        capsize=3,
        color=("#56B4E9", "#CC79A7"),
    )
    axes[0].set_xticks(range(2), (r"$H_s=0.27$ m", r"$H_s=0.05$ m"))
    axes[0].set_title("Low-motion noise floor")
    axes[0].set_ylabel(r"Vertical RMS error (\% $H_s$)")
    axes[0].grid(True, axis="y", alpha=0.25)

    x = np.arange(2, dtype=float)
    width = 0.34
    transition_upper = 0.0
    for offset, mode, color in (
        (-width / 2, "Adaptive", "#0072B2"),
        (width / 2, "FixedNominal", "#E69F00"),
    ):
        values = [
            float(index[("transition_rate", case, mode, "disp_z_pct_hs")]["mean"])
            for case in ("controlled", "rapid")
        ]
        errors = [
            (
                value - float(index[("transition_rate", case, mode, "disp_z_pct_hs")]["bootstrap_ci95_low"]),
                float(index[("transition_rate", case, mode, "disp_z_pct_hs")]["bootstrap_ci95_high"]) - value,
            )
            for case, value in zip(("controlled", "rapid"), values)
        ]
        transition_upper = max(
            transition_upper,
            max(value + high for value, (_, high) in zip(values, errors)),
        )
        axes[1].bar(
            x + offset,
            values,
            width,
            yerr=np.asarray(errors).T,
            capsize=3,
            color=color,
            label=mode,
        )
    axes[1].set_xticks(x, ("360 s ramp", "30 s ramp"))
    axes[1].set_title("Sea-state transition rate")
    axes[1].set_ylim(0.0, transition_upper * 1.22)
    axes[1].grid(True, axis="y", alpha=0.25)
    axes[1].legend(frameon=False, loc="upper center", ncol=2)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format="svg", metadata={"Date": None})
    plt.close(fig)


def write_publication_table(
    path: Path,
    summary: Sequence[Mapping[str, Any]],
    effects: Sequence[Mapping[str, Any]],
) -> None:
    index = _summary_index(summary)
    effect_index = _effect_index(effects)

    def summary_row(experiment: str, case: str, mode: str, metric: str) -> Mapping[str, Any]:
        return index[(experiment, case, mode, metric)]

    def mean_std(row: Mapping[str, Any], digits: int = 2) -> str:
        return f"{float(row['mean']):.{digits}f} $\\pm$ {float(row['std']):.{digits}f}"

    sensitivity_rows = [
        row for row in summary
        if row["experiment"] == "sensitivity" and row["metric"] == "disp_z_pct_hs"
    ]
    baseline_values = [
        float(row["mean"]) for row in sensitivity_rows
        if math.isclose(float(row["scale_label"]), 1.0, abs_tol=1e-12)
    ]
    baseline_mean = float(np.mean(baseline_values))
    worst = max(sensitivity_rows, key=lambda row: float(row["mean"]))
    if math.isclose(float(worst["scale_label"]), 1.0, abs_tol=1e-12):
        worst_effect = {
            "mean_paired_difference": 0.0,
            "bootstrap_ci95_low": 0.0,
            "bootstrap_ci95_high": 0.0,
        }
    else:
        worst_effect = effect_index[
            (
                "sensitivity",
                f"{worst['parameter']}_x{float(worst['scale_label']):g}_minus_x1",
                "disp_z_pct_hs",
            )
        ]
    low_effect = effect_index[("low_motion", "Hs0.05_minus_Hs0.27", "disp_z_pct_hs")]
    rapid_effect = effect_index[
        ("transition_rate", "rapid_minus_controlled_Adaptive", "disp_z_pct_hs")
    ]
    rapid_adaptation = effect_index[
        ("transition_rate", "Adaptive_minus_FixedNominal_rapid", "disp_z_pct_hs")
    ]
    r_s_half_3d = effect_index[
        ("sensitivity", "r_s_x0.5_minus_x1", "disp_3d_rms_m")
    ]
    sigma_half_pitch = effect_index[
        ("sensitivity", "sigma_aw_x0.5_minus_x1", "pitch_rms_deg")
    ]

    def vertical_span(parameter: str) -> float:
        """Peak-to-peak mean vertical error across a sweep, in points of H_s."""

        values = [
            float(row["mean"])
            for row in sensitivity_rows
            if row["parameter"] == parameter
        ]
        return max(values) - min(values) if values else math.nan

    parameter_word = {
        "tau": "tau",
        "sigma_aw": "sigma-aw",
        "r_s": "r-S",
        "sigma_aw_rs": "coupled sigma-aw and r-S",
        "tau_rs": "coupled tau and r-S",
    }
    lines = [
        r"% Generated by tools/ou_robustness.py from the committed full-study rows.",
        rf"\providecommand{{\OURobustnessPairs}}{{{int(worst['n'])}}}",
        rf"\providecommand{{\OURobustnessNominalMean}}{{{baseline_mean:.2f}}}",
        rf"\providecommand{{\OURobustnessWorstParameter}}{{{parameter_word[str(worst['parameter'])]}}}",
        rf"\providecommand{{\OURobustnessWorstScale}}{{{float(worst['scale_label']):g}}}",
        rf"\providecommand{{\OURobustnessWorstMean}}{{{float(worst['mean']):.2f}}}",
        rf"\providecommand{{\OURobustnessWorstDifference}}{{{float(worst_effect['mean_paired_difference']):+.2f}}}",
        rf"\providecommand{{\OURobustnessWorstDifferenceLow}}{{{float(worst_effect['bootstrap_ci95_low']):+.2f}}}",
        rf"\providecommand{{\OURobustnessWorstDifferenceHigh}}{{{float(worst_effect['bootstrap_ci95_high']):+.2f}}}",
        rf"\providecommand{{\OURobustnessRSHalfThreeDDifference}}{{{float(r_s_half_3d['mean_paired_difference']):+.3f}}}",
        rf"\providecommand{{\OURobustnessRSHalfThreeDLow}}{{{float(r_s_half_3d['bootstrap_ci95_low']):+.3f}}}",
        rf"\providecommand{{\OURobustnessRSHalfThreeDHigh}}{{{float(r_s_half_3d['bootstrap_ci95_high']):+.3f}}}",
        rf"\providecommand{{\OURobustnessSigmaHalfPitchDifference}}{{{float(sigma_half_pitch['mean_paired_difference']):+.3f}}}",
        rf"\providecommand{{\OURobustnessSigmaHalfPitchLow}}{{{float(sigma_half_pitch['bootstrap_ci95_low']):+.3f}}}",
        rf"\providecommand{{\OURobustnessSigmaHalfPitchHigh}}{{{float(sigma_half_pitch['bootstrap_ci95_high']):+.3f}}}",
        rf"\providecommand{{\OURobustnessTauSpan}}{{{vertical_span('tau'):.3f}}}",
        rf"\providecommand{{\OURobustnessSigmaSpan}}{{{vertical_span('sigma_aw'):.3f}}}",
        rf"\providecommand{{\OURobustnessRSSpan}}{{{vertical_span('r_s'):.3f}}}",
        rf"\providecommand{{\OURobustnessCoupledSigmaSpan}}{{{vertical_span('sigma_aw_rs'):.3f}}}",
        rf"\providecommand{{\OURobustnessCoupledTauSpan}}{{{vertical_span('tau_rs'):.3f}}}",
        rf"\providecommand{{\OURobustnessLowReferenceMean}}{{{float(summary_row('low_motion', 'Hs0.27', 'Adaptive', 'disp_z_pct_hs')['mean']):.2f}}}",
        rf"\providecommand{{\OURobustnessLowStressMean}}{{{float(summary_row('low_motion', 'Hs0.05', 'Adaptive', 'disp_z_pct_hs')['mean']):.2f}}}",
        rf"\providecommand{{\OURobustnessLowStressAbsolute}}{{{float(summary_row('low_motion', 'Hs0.05', 'Adaptive', 'disp_z_rms_m')['mean']):.4f}}}",
        rf"\providecommand{{\OURobustnessLowDifference}}{{{float(low_effect['mean_paired_difference']):+.2f}}}",
        rf"\providecommand{{\OURobustnessLowDifferenceLow}}{{{float(low_effect['bootstrap_ci95_low']):+.2f}}}",
        rf"\providecommand{{\OURobustnessLowDifferenceHigh}}{{{float(low_effect['bootstrap_ci95_high']):+.2f}}}",
        rf"\providecommand{{\OURobustnessControlledAdaptiveMean}}{{{float(summary_row('transition_rate', 'controlled', 'Adaptive', 'disp_z_pct_hs')['mean']):.2f}}}",
        rf"\providecommand{{\OURobustnessRapidAdaptiveMean}}{{{float(summary_row('transition_rate', 'rapid', 'Adaptive', 'disp_z_pct_hs')['mean']):.2f}}}",
        rf"\providecommand{{\OURobustnessRapidDifference}}{{{float(rapid_effect['mean_paired_difference']):+.2f}}}",
        rf"\providecommand{{\OURobustnessRapidDifferenceLow}}{{{float(rapid_effect['bootstrap_ci95_low']):+.2f}}}",
        rf"\providecommand{{\OURobustnessRapidDifferenceHigh}}{{{float(rapid_effect['bootstrap_ci95_high']):+.2f}}}",
        rf"\providecommand{{\OURobustnessRapidAdaptationDifference}}{{{float(rapid_adaptation['mean_paired_difference']):+.2f}}}",
        rf"\providecommand{{\OURobustnessRapidAdaptationLow}}{{{float(rapid_adaptation['bootstrap_ci95_low']):+.2f}}}",
        rf"\providecommand{{\OURobustnessRapidAdaptationHigh}}{{{float(rapid_adaptation['bootstrap_ci95_high']):+.2f}}}",
        "",
        r"\begin{table*}[t]",
        r"  \centering",
        r"  \caption{OU--III tuning sensitivity at the nominal $H_s=1.5$ m sea. Each entry is vertical-displacement RMS error in percent of $H_s$ (mean $\pm$ sample standard deviation, $n=\OURobustnessPairs$ paired seed triplets). The first three columns hold the other two parameters frozen; the last two move $r_S$ as the deployed law $r_S=\operatorname{clip}(0.35\,\sigma_{aw}\tau^3,0.4,400)$ does. The $r_S$ multiplier acts on the pseudo-measurement standard deviation, so the corresponding covariance scales quadratically.}",
        r"  \label{tab:ou_robustness_sensitivity}",
        r"  \footnotesize",
        r"  \setlength{\tabcolsep}{7pt}",
        r"  \begin{tabular}{@{}rrrrrr@{}}",
        r"    \toprule",
        r"    & \multicolumn{3}{c}{one factor at a time} & \multicolumn{2}{c}{coupled} \\",
        r"    \cmidrule(lr){2-4}\cmidrule(lr){5-6}",
        r"    Multiplier & $\tau$ & $\sigma_{aw}$ & $r_S$ & $\sigma_{aw},r_S$ & $\tau,r_S\propto\tau^3$ \\",
        r"    \midrule",
    ]
    scales = sorted({float(row["scale_label"]) for row in sensitivity_rows})
    for scale in scales:
        cells = []
        for parameter in SENSITIVITY_PARAMETERS:
            row = next(
                item for item in sensitivity_rows
                if item["parameter"] == parameter
                and math.isclose(float(item["scale_label"]), scale, abs_tol=1e-12)
            )
            cells.append(mean_std(row))
        lines.append(f"    {scale:g} & " + " & ".join(cells) + r" \\")
    lines.extend(
        (
            r"    \bottomrule",
            r"  \end{tabular}",
            r"\end{table*}",
            "",
            r"\begin{table*}[t]",
            r"  \centering",
            r"  \caption{OU--III degradation cases over the long scoring window (mean $\pm$ sample standard deviation, $n=\OURobustnessPairs$ paired seed triplets). Transition rows are normalized by the final $H_s=4.0$ m sea.}",
            r"  \label{tab:ou_robustness_stress}",
            r"  \footnotesize",
            r"  \setlength{\tabcolsep}{5pt}",
            r"  \begin{tabular}{@{}llrrrr@{}}",
            r"    \toprule",
            r"    Case & Mode & Z [m] & Z [\%$H_s$] & 3D [m] & Pitch [$^\circ$] \\",
            r"    \midrule",
        )
    )
    stress_rows = (
        ("Low motion, $H_s=0.27$ m", "low_motion", "Hs0.27", "Adaptive"),
        ("Low motion, $H_s=0.05$ m", "low_motion", "Hs0.05", "Adaptive"),
        ("Controlled 360 s ramp", "transition_rate", "controlled", "Adaptive"),
        ("Controlled 360 s ramp", "transition_rate", "controlled", "FixedNominal"),
        ("Rapid 30 s ramp", "transition_rate", "rapid", "Adaptive"),
        ("Rapid 30 s ramp", "transition_rate", "rapid", "FixedNominal"),
    )
    for label, experiment, case, mode in stress_rows:
        values = [
            mean_std(summary_row(experiment, case, mode, metric), digits)
            for metric, digits in (
                ("disp_z_rms_m", 3),
                ("disp_z_pct_hs", 2),
                ("disp_3d_rms_m", 3),
                ("pitch_rms_deg", 2),
            )
        ]
        lines.append(f"    {label} & {mode} & " + " & ".join(values) + r" \\")
    lines.extend((r"    \bottomrule", r"  \end{tabular}", r"\end{table*}"))
    core.assert_latex_macro_names(lines)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _run_unit(
    *,
    experiment: str,
    case: str,
    parameter: str,
    scale_label: str,
    scale_multiplier: float | str,
    mode: str,
    repetition: int,
    seed: core.SeedTriplet,
    true_hs_m: float,
    transition_duration_sec: float | str,
    input_path: Path,
    window_sec: float,
    tuning_mode: str,
    tuning_point: core.TuningPoint | None,
) -> dict[str, Any]:
    metrics, gate_pass, return_code = core.run_simulator(
        "OU_III",
        input_path,
        window_sec,
        seed.imu_noise_seed,
        seed.initialization_seed,
        tuning_mode=tuning_mode,
        tuning_point=tuning_point,
    )
    return (
        {
            "run_id": f"{experiment}:{case}:{mode}:{repetition}",
            "experiment": experiment,
            "case": case,
            "parameter": parameter,
            "scale_label": scale_label,
            "scale_multiplier": scale_multiplier,
            "family": "OU_III",
            "mode": mode,
            "repetition": repetition,
            "wave_phase_seed": seed.wave_phase_seed,
            "imu_noise_seed": seed.imu_noise_seed,
            "initialization_seed": seed.initialization_seed,
            "true_hs_m": true_hs_m,
            "transition_duration_sec": transition_duration_sec,
            "score_window_sec": window_sec,
            "configured_tau_s": tuning_point.tau_s if tuning_point else "adaptive",
            "configured_sigma_aw_mps2": (
                tuning_point.sigma_a_mps2 if tuning_point else "adaptive"
            ),
            "configured_r_s_ms": tuning_point.RS_ms if tuning_point else "adaptive",
            "historical_60s_gate_pass": int(gate_pass),
            "simulator_return_code": return_code,
            **{
                key: value for key, value in metrics.items()
                if key not in ("family", "tuning_mode", "input")
            },
        }
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("smoke", "full"), default="smoke")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=REPO_ROOT / "plots" / "kalman_ou_ii",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "reports" / "results" / "ou_robustness",
    )
    parser.add_argument("--duration-sec", type=float)
    parser.add_argument("--window-sec", type=float)
    parser.add_argument("--sensitivity-scales", type=parse_float_list)
    parser.add_argument("--wave-seeds", type=core.parse_int_list)
    parser.add_argument("--imu-seeds", type=core.parse_int_list)
    parser.add_argument("--initialization-seeds", type=core.parse_int_list)
    parser.add_argument("--bootstrap-resamples", type=int, default=10_000)
    parser.add_argument("--stats-seed", type=int, default=20260329)
    parser.add_argument(
        "--jobs",
        type=int,
        default=max(1, min(8, (os.cpu_count() or 1))),
        help="parallel simulator replays within a repetition; every run is "
        "fixed by its seed triplet and tuning point, so this changes runtime "
        "only",
    )
    parser.add_argument("--eigen-dir")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.data_dir = args.data_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    duration_sec = args.duration_sec or (180.0 if args.mode == "smoke" else 1200.0)
    window_sec = args.window_sec or (120.0 if args.mode == "smoke" else 900.0)
    if not duration_sec > window_sec > 0.0:
        raise ValueError("duration must be greater than the positive scoring window")
    scales = args.sensitivity_scales or (
        list(SMOKE_SCALES) if args.mode == "smoke" else list(FULL_SCALES)
    )

    if args.mode == "smoke":
        default_wave, default_imu, default_init = ([11], [101], [1009])
        transition_cases = (
            TransitionCase("controlled", 60.0, 120.0),
            TransitionCase("rapid", 85.0, 95.0),
        )
    else:
        default_wave = list(core.DEFAULT_FULL_WAVE_SEEDS)
        default_imu = list(core.DEFAULT_FULL_IMU_SEEDS)
        default_init = list(core.DEFAULT_FULL_INIT_SEEDS)
        transition_cases = (
            TransitionCase("controlled", 420.0, 780.0),
            TransitionCase("rapid", 585.0, 615.0),
        )
    if duration_sec < max(case.end_sec for case in transition_cases):
        raise ValueError("duration does not contain the configured transition cases")
    seeds = core.broadcast_seed_triplets(
        args.wave_seeds or default_wave,
        args.imu_seeds or default_imu,
        args.initialization_seeds or default_init,
    )

    nominal_input = core.find_default_input(args.data_dir, "1.500", "50.710").resolve()
    low_input = core.find_default_input(args.data_dir, "0.270", "14.047").resolve()
    endpoint_input = core.find_default_input(args.data_dir, "8.500", "202.839").resolve()
    if not args.skip_build:
        core.build_simulators(("OU_III",), args.eigen_dir)

    baseline = core.calibrate_tuning_point("OU_III", nominal_input, window_sec)
    validate_tuning_point(baseline)
    for parameter in SENSITIVITY_PARAMETERS:
        for scale in scales:
            validate_tuning_point(scaled_tuning_point(baseline, parameter, scale))
    nominal_columns, nominal_data = core.read_wave_csv(nominal_input, duration_sec)
    low_columns, low_data = core.read_wave_csv(low_input, duration_sec)
    endpoint_columns, endpoint_data = core.read_wave_csv(endpoint_input, duration_sec)
    if nominal_columns != endpoint_columns:
        raise ValueError("transition source columns do not match")

    rows: list[dict[str, Any]] = []
    with tempfile.TemporaryDirectory(prefix="ocean-imu-ou-robustness-") as temporary:
        work_dir = Path(temporary)
        for repetition, seed in enumerate(seeds, start=1):
            print(f"Repetition {repetition}/{len(seeds)}", flush=True)

            # Every unit in a repetition is an independent simulator replay
            # whose result is fixed by its seed triplet and tuning point, so
            # they are submitted together and drained in submission order.
            # Order is preserved because the raw CSV is hashed in the manifest;
            # concurrency must not reorder it.
            pending: list[Any] = []
            scratch: list[Path] = []
            pool = ThreadPoolExecutor(max_workers=max(1, args.jobs))

            sensitivity_path = work_dir / "sensitivity" / str(repetition) / nominal_input.name
            core.write_wave_csv(
                sensitivity_path,
                nominal_columns,
                core.phase_randomize_wave(nominal_columns, nominal_data, seed.wave_phase_seed),
            )
            scratch.append(sensitivity_path)
            for parameter in SENSITIVITY_PARAMETERS:
                for scale in scales:
                    point = scaled_tuning_point(baseline, parameter, scale)
                    pending.append(pool.submit(
                        _run_unit,
                        experiment="sensitivity",
                        case=f"{parameter}_x{scale:g}",
                        parameter=parameter,
                        scale_label=f"{scale:g}",
                        scale_multiplier=float(scale),
                        mode="FixedSensitivity",
                        repetition=repetition,
                        seed=seed,
                        true_hs_m=1.5,
                        transition_duration_sec="",
                        input_path=sensitivity_path,
                        window_sec=window_sec,
                        tuning_mode="fixed_sensitivity",
                        tuning_point=point,
                    ))

            low_randomized = core.phase_randomize_wave(
                low_columns, low_data, seed.wave_phase_seed
            )
            for height in (LOW_REFERENCE_HS_M, LOW_STRESS_HS_M):
                scale = height / LOW_REFERENCE_HS_M
                filename = re.sub(r"_H0\.270_", f"_H{height:.3f}_", low_input.name)
                low_path = work_dir / "low_motion" / str(repetition) / filename
                core.write_wave_csv(
                    low_path,
                    low_columns,
                    low_randomized if math.isclose(scale, 1.0) else core.scale_wave_motion(
                        low_columns, low_randomized, scale
                    ),
                )
                scratch.append(low_path)
                case = f"Hs{height:.2f}"
                pending.append(pool.submit(
                    _run_unit,
                    experiment="low_motion",
                    case=case,
                    parameter="",
                    scale_label="",
                    scale_multiplier="",
                    mode="Adaptive",
                    repetition=repetition,
                    seed=seed,
                    true_hs_m=height,
                    transition_duration_sec="",
                    input_path=low_path,
                    window_sec=window_sec,
                    tuning_mode="adaptive",
                    tuning_point=None,
                ))

            for transition in transition_cases:
                transition_path = (
                    work_dir
                    / "transition"
                    / transition.name
                    / str(repetition)
                    / "wave_data_jonswap_H4.000_L202.839_A-30.00_P120.00.csv"
                )
                generated = core.make_nonstationary_wave(
                    nominal_columns,
                    nominal_data,
                    endpoint_data,
                    seed.wave_phase_seed,
                    4.0 / 8.5,
                    transition.start_sec,
                    transition.end_sec,
                )
                core.write_wave_csv(transition_path, nominal_columns, generated)
                scratch.append(transition_path)
                for mode, tuning_mode, point in (
                    ("Adaptive", "adaptive", None),
                    ("FixedNominal", "fixed_nominal", baseline),
                ):
                    pending.append(pool.submit(
                        _run_unit,
                        experiment="transition_rate",
                        case=transition.name,
                        parameter="",
                        scale_label="",
                        scale_multiplier="",
                        mode=mode,
                        repetition=repetition,
                        seed=seed,
                        true_hs_m=4.0,
                        transition_duration_sec=transition.duration_sec,
                        input_path=transition_path,
                        window_sec=window_sec,
                        tuning_mode=tuning_mode,
                        tuning_point=point,
                    ))

            try:
                for index, future in enumerate(pending, start=1):
                    rows.append(future.result())
                    print(
                        f"  repetition {repetition}: {index}/{len(pending)} runs",
                        flush=True,
                    )
            finally:
                pool.shutdown()
                for path in scratch:
                    path.unlink(missing_ok=True)

    summary = summarize_rows(rows, args.bootstrap_resamples, args.stats_seed)
    effects = paired_effect_rows(rows, args.bootstrap_resamples, args.stats_seed)

    raw_path = args.output_dir / "ou_robustness_raw.csv"
    summary_path = args.output_dir / "ou_robustness_summary.csv"
    effects_path = args.output_dir / "ou_robustness_paired_effects.csv"
    json_path = args.output_dir / "ou_robustness.json"
    publication_path = args.output_dir / "ou_robustness_publication.tex"
    manifest_path = args.output_dir / "ou_robustness_manifest.json"
    sensitivity_plot = args.output_dir / "ou_robustness_sensitivity.svg"
    stress_plot = args.output_dir / "ou_robustness_stress.svg"

    core.write_csv(raw_path, rows)
    core.write_csv(summary_path, summary)
    core.write_csv(effects_path, effects)
    write_publication_table(publication_path, summary, effects)
    plot_paths: list[Path] = []
    if not args.no_plots:
        write_sensitivity_plot(sensitivity_plot, summary)
        write_stress_plot(stress_plot, summary)
        plot_paths = [sensitivity_plot, stress_plot]

    protocol = {
        "mode": args.mode,
        "duration_sec": duration_sec,
        "score_window_sec": window_sec,
        "seed_triplets": [asdict(seed) for seed in seeds],
        "sensitivity_parameters": list(SENSITIVITY_PARAMETERS),
        "sensitivity_scales": list(scales),
        "sensitivity_design": "one factor at a time around one frozen nominal OU-III point",
        "r_s_semantics": "standard deviation; R_S = diag(r_S^2)",
        "low_motion_heights_m": [LOW_REFERENCE_HS_M, LOW_STRESS_HS_M],
        "transition_cases": [asdict(case) for case in transition_cases],
        "wave_phase_method": core.WAVE_PHASE_METHOD,
        "transition_method": core.TRANSITION_METHOD,
        "bootstrap_resamples": args.bootstrap_resamples,
        "stats_seed": args.stats_seed,
        "pairing": "identical wave-realization, IMU-noise, and initialization seeds",
    }
    core.write_json(
        json_path,
        {
            "protocol": protocol,
            "nominal_tuning_point": asdict(baseline),
            "raw_runs": rows,
            "summary": summary,
            "paired_effects": effects,
        },
    )

    result_paths = [
        raw_path,
        summary_path,
        effects_path,
        json_path,
        publication_path,
        *plot_paths,
    ]
    source_paths = [
        nominal_input,
        low_input,
        endpoint_input,
        Path(__file__).resolve(),
        REPO_ROOT / "tools" / "ou_validation.py",
        REPO_ROOT / "src" / "util" / "W3dSimCommon.cpp",
        REPO_ROOT / "tests" / "kalman_ou_iii" / "kalman_ou_iii-sim.cpp",
    ]
    manifest = {
        "git_commit": core.git_output("rev-parse", "HEAD"),
        "git_branch": core.git_output("branch", "--show-current"),
        "git_diff_stat": core.git_output("diff", "--stat"),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
        "command": [sys.executable, *sys.argv],
        "protocol": protocol,
        "nominal_tuning_point": asdict(baseline),
        "source_files": {
            str(path.relative_to(REPO_ROOT)): {
                "sha256": core.sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for path in source_paths
        },
        "result_files": {
            path.name: {"sha256": core.sha256_file(path), "bytes": path.stat().st_size}
            for path in result_paths
        },
    }
    core.write_json(manifest_path, manifest)
    for path in (*result_paths, manifest_path):
        print(f"Wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
