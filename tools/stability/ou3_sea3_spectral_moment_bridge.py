#!/usr/bin/env python3
"""SEA0 spectral-moment subcertificate for the OU-III SEA3 theorem.

This module certifies only the *surface-elevation spectral* bridge between a
JONSWAP partition's peak period Tp and its zero-crossing/moment period Tz.  It
also records the exact multimodal moment-composition identity.  It does not
promote the full SEA0 certificate and it does not identify the deployed tuner's
Tz: the latter still requires the directional vessel/IMU response and the exact
WavePeriodEstimator front end.

No replay data are read.  The only numerical screen is a nondimensional JONSWAP
shape integral over the declared gamma interval.  Gamma cells exploit pointwise
monotonicity of gamma**r to enclose each continuum cell; fixed conservative
quadrature and tail pads are added.  Because the quadrature itself is not an
interval-arithmetic proof, the resulting gamma-continuum numbers remain a
non-promoting screening enclosure.  The PM baseline and multimodal composition
identities are analytical.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "OU3_SEA3_SPECTRAL_MOMENT_BRIDGE_V1"
M_MAX = 3
GAMMA_MIN = 1.0
GAMMA_MAX = 7.0
GAMMA_CELLS = 240
SIGMA_LOW = 0.07
SIGMA_HIGH = 0.09
PM_EXPONENT = 1.25
X_MIN = 0.1
X_MAX = 512.0
SIMPSON_PANELS = 16384
QUADRATURE_ABS_PAD = 5.0e-8


def _jonswap_peak_shape(x: float, gamma: float) -> float:
    sigma = SIGMA_LOW if x <= 1.0 else SIGMA_HIGH
    r = math.exp(-((x - 1.0) ** 2) / (2.0 * sigma * sigma))
    return gamma**r


def _dimensionless_integrand(x: float, gamma: float, moment_order: int) -> float:
    if not (x > 0.0):
        return 0.0
    base = x ** (moment_order - 5) * math.exp(-PM_EXPONENT * x**-4)
    return base * _jonswap_peak_shape(x, gamma)


def _simpson_log_integral(gamma: float, moment_order: int) -> float:
    """Composite Simpson integral on log(x), excluding analytical tail pads."""
    if SIMPSON_PANELS % 2:
        raise RuntimeError("SIMPSON_PANELS must be even")
    log_lo = math.log(X_MIN)
    log_hi = math.log(X_MAX)
    h = (log_hi - log_lo) / SIMPSON_PANELS
    total = 0.0
    for i in range(SIMPSON_PANELS + 1):
        t = log_lo + i * h
        x = math.exp(t)
        # dx = x dt.
        value = _dimensionless_integrand(x, gamma, moment_order) * x
        if i == 0 or i == SIMPSON_PANELS:
            weight = 1.0
        elif i % 2:
            weight = 4.0
        else:
            weight = 2.0
        total += weight * value
    return total * h / 3.0


def _lower_tail_upper(moment_order: int, gamma_upper: float) -> float:
    """Analytical upper bound for x in (0, X_MIN]."""
    exponent = math.exp(-PM_EXPONENT * X_MIN**-4)
    if moment_order == 0:
        # gamma/4 int y^0 exp(-a y) dy, y=x^-4.
        return gamma_upper * exponent / (4.0 * PM_EXPONENT)
    if moment_order == 2:
        # y^-1/2 <= y0^-1/2 on [y0, inf), giving a simple safe bound.
        return (
            gamma_upper
            * X_MIN**2
            * exponent
            / (4.0 * PM_EXPONENT)
        )
    raise ValueError("only m0 and m2 are certified")


def _upper_tail_upper(moment_order: int, gamma_upper: float) -> float:
    """Analytical upper bound for x in [X_MAX, inf)."""
    if moment_order == 0:
        return gamma_upper / (4.0 * X_MAX**4)
    if moment_order == 2:
        return gamma_upper / (2.0 * X_MAX**2)
    raise ValueError("only m0 and m2 are certified")


def _moment_interval(gamma: float, moment_order: int) -> tuple[float, float]:
    middle = _simpson_log_integral(gamma, moment_order)
    tail = _lower_tail_upper(moment_order, gamma) + _upper_tail_upper(
        moment_order, gamma
    )
    lo = max(0.0, middle - QUADRATURE_ABS_PAD)
    hi = middle + QUADRATURE_ABS_PAD + tail
    return lo, hi


def _ratio_interval_for_gamma_cell(
    gamma_lo: float,
    gamma_hi: float,
) -> tuple[float, float]:
    """Outer Tz/Tp interval over one gamma cell.

    The JONSWAP peak factor gamma**r is pointwise nondecreasing in gamma because
    r is nonnegative.  Therefore each positive spectral moment is bounded by
    its endpoint moments.  Cross-combining numerator/denominator endpoints is
    conservative and does not assume Tz/Tp itself is monotone in gamma.
    """
    m0_lo, _ = _moment_interval(gamma_lo, 0)
    _, m0_hi = _moment_interval(gamma_hi, 0)
    m2_lo, _ = _moment_interval(gamma_lo, 2)
    _, m2_hi = _moment_interval(gamma_hi, 2)
    if not (m0_lo > 0.0 and m2_lo > 0.0):
        raise RuntimeError("spectral moment enclosure lost positivity")
    return math.sqrt(m0_lo / m2_hi), math.sqrt(m0_hi / m2_lo)


def mixture_zero_crossing_period(
    heights: list[float],
    periods: list[float],
) -> float:
    """Exact Tz of a sum of uncorrelated normalized spectral partitions."""
    if len(heights) != len(periods) or not heights:
        raise ValueError("heights and periods must have the same nonzero length")
    energy = 0.0
    weighted_inverse_square = 0.0
    for height, period in zip(heights, periods):
        h = float(height)
        t = float(period)
        if not (math.isfinite(h) and h >= 0.0):
            raise ValueError("partition heights must be finite and nonnegative")
        if h == 0.0:
            continue
        if not (math.isfinite(t) and t > 0.0):
            raise ValueError("active partition periods must be positive and finite")
        e = h * h
        energy += e
        weighted_inverse_square += e / (t * t)
    if not (energy > 0.0 and weighted_inverse_square > 0.0):
        raise ValueError("at least one active partition is required")
    return math.sqrt(energy / weighted_inverse_square)


def _pm_exact_ratio() -> float:
    # For gamma=1, I0=1/(4a), I2=sqrt(pi)/(4 sqrt(a)), a=5/4.
    return (PM_EXPONENT * math.pi) ** (-0.25)


def build() -> dict[str, Any]:
    gamma_step = (GAMMA_MAX - GAMMA_MIN) / GAMMA_CELLS
    global_lo = math.inf
    global_hi = 0.0
    limiting_lo: list[float] | None = None
    limiting_hi: list[float] | None = None

    for cell in range(GAMMA_CELLS):
        gamma_lo = GAMMA_MIN + cell * gamma_step
        gamma_hi = GAMMA_MIN + (cell + 1) * gamma_step
        ratio_lo, ratio_hi = _ratio_interval_for_gamma_cell(gamma_lo, gamma_hi)
        if ratio_lo < global_lo:
            global_lo = ratio_lo
            limiting_lo = [gamma_lo, gamma_hi]
        if ratio_hi > global_hi:
            global_hi = ratio_hi
            limiting_hi = [gamma_lo, gamma_hi]

    pm_exact = _pm_exact_ratio()
    example_heights = [3.0, 4.0, 0.0]
    example_periods = [5.0, 10.0, 7.0]
    example_mix = mixture_zero_crossing_period(example_heights, example_periods)
    pm_screen = _ratio_interval_for_gamma_cell(GAMMA_MIN, GAMMA_MIN)
    gamma33_screen = _ratio_interval_for_gamma_cell(3.3, 3.3)
    gamma7_screen = _ratio_interval_for_gamma_cell(GAMMA_MAX, GAMMA_MAX)

    return {
        "schema_version": SCHEMA_VERSION,
        "qualification": "OU3_SEA0_SURFACE_SPECTRAL_MOMENT_SUBCERTIFICATE",
        "proof_status": "mixed_analytic_and_non_promoting_numerical_screen",
        "SEA0_full_certificate_promoted": False,
        "P2_promoted_from_this_artifact": False,
        "P3_promoted_from_this_artifact": False,
        "P4_promoted_from_this_artifact": False,
        "P5_promoted_from_this_artifact": False,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_operating_domain_shrunk": False,
        "sea_family": {
            "m_max": M_MAX,
            "partition_frequency_shape": "JONSWAP_with_PM_at_gamma_1",
            "declared_gamma_interval": [GAMMA_MIN, GAMMA_MAX],
            "jonswap_sigma_below_peak": SIGMA_LOW,
            "jonswap_sigma_above_peak": SIGMA_HIGH,
            "partition_energy_coordinate": "m0_r = H_r^2 / 16",
            "total_energy_coupling": "H_s^2 = sum_r H_r^2",
        },
        "analytical_lemmas": {
            "pm_gamma_1_exact_Tz_over_Tp": round(pm_exact, 12),
            "pm_identity": "Tz/Tp = ((5/4)*pi)^(-1/4)",
            "multimodal_zero_crossing_identity": (
                "1/Tz_mix^2 = sum_r w_r/Tz_r^2, "
                "w_r = H_r^2/sum_j H_j^2"
            ),
            "multimodal_period_between_component_extrema": True,
            "arithmetic_period_average_used": False,
            "multimodal_identity_example": {
                "H_r_m": example_heights,
                "Tz_r_s": example_periods,
                "Tz_mix_s": round(example_mix, 12),
                "lies_between_active_component_periods": (
                    min(example_periods[:2])
                    <= example_mix
                    <= max(example_periods[:2])
                ),
            },
            "unbanded_surface_acceleration_variance_finite": False,
            "unbanded_acceleration_tail_reason": (
                "JONSWAP S_eta ~ omega^-5, so omega^4*S_eta ~ omega^-1"
            ),
            "band_or_directional_response_required_for_acceleration_moments": True,
        },
        "gamma_continuum_screen": {
            "method": (
                "240 gamma cells; pointwise gamma monotonicity; log-space "
                "composite Simpson endpoint moments; conservative quadrature "
                "and analytical tail pads"
            ),
            "gamma_cells": GAMMA_CELLS,
            "quadrature": {
                "x_interval": [X_MIN, X_MAX],
                "simpson_panels": SIMPSON_PANELS,
                "absolute_pad_per_moment": QUADRATURE_ABS_PAD,
                "lower_and_upper_tails_bounded_analytically": True,
                "interval_arithmetic_used_for_quadrature": False,
            },
            "surface_elevation_Tz_over_Tp_outer": [
                round(global_lo, 12),
                round(global_hi, 12),
            ],
            "limiting_gamma_cell_for_lower": [
                round(x, 12) for x in limiting_lo or []
            ],
            "limiting_gamma_cell_for_upper": [
                round(x, 12) for x in limiting_hi or []
            ],
            "point_checks": {
                "gamma_1_PM": [round(x, 12) for x in pm_screen],
                "gamma_3p3": [round(x, 12) for x in gamma33_screen],
                "gamma_7": [round(x, 12) for x in gamma7_screen],
            },
            "promotion_use": False,
            "reason_not_promoted": (
                "floating quadrature is padded but not a validated interval "
                "integration proof"
            ),
        },
        "tuner_bridge_contract": {
            "surface_elevation_Tz_may_be_substituted_for_tuner_Tz": False,
            "directional_vessel_IMU_RAO_required": True,
            "exact_WavePeriodEstimator_front_end_required": True,
            "finite_EMA_and_log_period_state_required": True,
            "mixture_moments_must_be_formed_before_period_ratio": True,
        },
        "next_obligation": (
            "certify the directional vessel/IMU response and deployed band, "
            "then propagate the resulting response-weighted moments through "
            "the exact WavePeriodEstimator before using a Tz/Tp restriction "
            "to prune the P2 source language"
        ),
    }


def validate(payload: dict[str, Any]) -> list[str]:
    failures: list[str] = []
    if payload.get("schema_version") != SCHEMA_VERSION:
        failures.append("schema mismatch")
    if payload.get("trajectory_replay_used") is not False:
        failures.append("spectral bridge must be replay free")
    if payload.get("SEA0_full_certificate_promoted") is not False:
        failures.append("surface-only bridge must not promote SEA0")
    sea = payload.get("sea_family", {})
    if sea.get("m_max") != M_MAX:
        failures.append("SEA3 mode count changed")
    if sea.get("declared_gamma_interval") != [GAMMA_MIN, GAMMA_MAX]:
        failures.append("declared JONSWAP gamma interval changed")
    analytical = payload.get("analytical_lemmas", {})
    if analytical.get("multimodal_period_between_component_extrema") is not True:
        failures.append("multimodal period identity missing")
    if analytical.get("arithmetic_period_average_used") is not False:
        failures.append("arithmetic period averaging is forbidden")
    example = analytical.get("multimodal_identity_example", {})
    if example.get("lies_between_active_component_periods") is not True:
        failures.append("multimodal identity example escaped component extrema")
    if analytical.get("unbanded_surface_acceleration_variance_finite") is not False:
        failures.append(
            "unbanded JONSWAP acceleration variance must not be treated as finite"
        )
    screen = payload.get("gamma_continuum_screen", {})
    bounds = screen.get("surface_elevation_Tz_over_Tp_outer", [])
    if len(bounds) != 2 or not (
        0.0 < float(bounds[0]) < float(bounds[1]) < 1.0
    ):
        failures.append("invalid surface Tz/Tp screen")
    if not (float(bounds[0]) < _pm_exact_ratio() < float(bounds[1])):
        failures.append("PM exact ratio escaped the gamma screen")
    if screen.get("promotion_use") is not False:
        failures.append("floating spectral screen must remain non-promoting")
    bridge = payload.get("tuner_bridge_contract", {})
    if bridge.get("surface_elevation_Tz_may_be_substituted_for_tuner_Tz") is not False:
        failures.append("surface Tz must not be silently substituted for tuner Tz")
    if bridge.get("directional_vessel_IMU_RAO_required") is not True:
        failures.append("RAO obligation was dropped")
    return failures


def _render(payload: dict[str, Any]) -> str:
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="optional path for the generated machine-readable subcertificate",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="check the committed JSON artifact against the current implementation",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    args = parser.parse_args()

    payload = build()
    failures = validate(payload)
    payload["validation_pass"] = not failures
    payload["validation_failures"] = failures
    rendered = _render(payload)

    committed = args.repo_root / "tools" / "ou3_sea3_spectral_moment_bridge.json"
    if args.check:
        if not committed.exists():
            raise SystemExit(f"missing {committed}")
        if committed.read_text(encoding="utf-8") != rendered:
            raise SystemExit(
                "SEA3 spectral-moment bridge is stale; rerun "
                "python3 tools/stability/ou3_sea3_spectral_moment_bridge.py "
                "--output tools/ou3_sea3_spectral_moment_bridge.json"
            )
        if failures:
            raise SystemExit(f"SEA3 spectral bridge validation failed: {failures}")
        return 0

    if args.output is not None:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
