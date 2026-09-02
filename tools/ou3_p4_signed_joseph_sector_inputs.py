#!/usr/bin/env python3
"""Retained P3 -> broad-sector P4 composition contract for OU-III.

This producer deliberately stops one step before theorem promotion.  It binds
four already-certified ingredients that the complete H/A nonlinear word must
use together:

* the dependency-preserving source-uniform P3 information metric;
* the global 0.8-rad Cayley chart/geometry certificate;
* source-complete timing, with uncertain S=0 firings left inside linear P3;
* the exact effective-vector-input identities for accelerometer/magnetometer.

The active route does **not** subtract an independently boxed
``eta^T R^-1 eta`` from the translation-limited scalar P3 gap.  The magnetic
radial residual is in the exact Kalman-gain nullspace and the accelerometer
finite-angle remainder is exactly representable through the source a_w
measurement column.  Consequently the remaining numerical obligation is a
source-correlated word calculation carrying P,H,R,S,r, the effective vector
input, covariance recursion, quaternion injection and immediate reset in one
jointly reachable tuple.

Likewise this module does not insert a global covariance condition-number
conversion between P3 and P4.  P4 must use the same node metric

    M_i = s_m Sigma_i^-1

with all attitude-linear cross terms retained.  A complete P4 producer may
promote only after it proves a strict H/A generalized word margin and prefix
safety over the entire declared 0.8-rad entrance sector.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p3_source_uniform_certificate as P3
import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_effective_vector_inputs as EFFECTIVE
import ou3_p4_source_word_timing as TIMING

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2


def _finite_positive(x, label: str) -> float:
    y = float(x)
    if not math.isfinite(y) or y <= 0.0:
        raise RuntimeError(f"{label} must be finite positive, got {x!r}")
    return y


def _slice_positive(values, sl: slice, label: str) -> dict:
    xs = [float(x) for x in values[sl]]
    if not xs or any(not math.isfinite(x) or x <= 0.0 for x in xs):
        raise RuntimeError(f"{label} must contain finite positive entries")
    return {"min": min(xs), "max": max(xs), "values": xs}


def _mode(mode: str, p3: dict, effective: dict) -> dict:
    m = p3["modes"][mode]
    n = 18 if mode == "H" else 21
    delta = _finite_positive(
        m["relative_Riccati_injection_margin_lower"], f"{mode} P3 delta"
    )
    sigma_diag = list(m["Sigma_diagonal_upper"])
    scale2 = list(m["comparison_scale_diagonal_squared"])
    if len(sigma_diag) != n or len(scale2) != n:
        raise RuntimeError(f"{mode} P3 directional metric dimension mismatch")

    # These are reported as directional inputs only.  They are intentionally
    # not collapsed to max(Sigma)/min(Sigma), because that would reintroduce
    # the translation-dominated condition-number loss the broad-sector route
    # is designed to avoid.
    directional = {
        "attitude_comparison_scale_squared": _slice_positive(
            scale2, slice(0, 3), f"{mode} attitude comparison scale"
        ),
        "aw_comparison_scale_squared": _slice_positive(
            scale2, slice(15, 18), f"{mode} aw comparison scale"
        ),
        "attitude_Sigma_diagonal_upper": _slice_positive(
            sigma_diag, slice(0, 3), f"{mode} attitude Sigma upper"
        ),
        "aw_Sigma_diagonal_upper": _slice_positive(
            sigma_diag, slice(15, 18), f"{mode} aw Sigma upper"
        ),
    }
    if mode == "A":
        directional["accelerometer_bias_comparison_scale_squared"] = _slice_positive(
            scale2, slice(18, 21), "A accelerometer-bias comparison scale"
        )
        directional["accelerometer_bias_Sigma_diagonal_upper"] = _slice_positive(
            sigma_diag, slice(18, 21), "A accelerometer-bias Sigma upper"
        )

    return {
        "dimension": n,
        "P3_relative_Riccati_injection_margin_lower": delta,
        "P3_prefix_information_gain_upper": m["prefix_information_gain_upper"],
        "P3_Sigma_lambda_min_lower": m["Sigma_lambda_min_lower"],
        "P3_Sigma_lambda_max_upper": m["Sigma_lambda_max_upper"],
        "P3_word_noise_Omega_lambda_min_lower": m[
            "word_noise_Omega_lambda_min_lower"
        ],
        "directional_same_metric_inputs": directional,
        "mag_radial_residual_gain_null_exact": effective[
            "mag_radial_residual_gain_null_exact"
        ],
        "mag_effective_coordinate_nonexpansive_factor_upper": effective[
            "mag_effective_coordinate_nonexpansive_factor_upper"
        ],
        "mag_effective_coordinate_tangent_defect_factor_upper": effective[
            "mag_effective_coordinate_tangent_defect_factor_upper"
        ],
        "acc_eta_in_aw_measurement_range_exact": effective[
            "acc_eta_in_aw_measurement_range_exact"
        ],
        "acc_effective_aw_input_isometry_exact": effective[
            "acc_effective_aw_input_isometry_exact"
        ],
        "acc_force_remainder_factor_upper": effective[
            "acc_force_remainder_factor_upper"
        ],
        "acc_aw_rotation_factor_upper": effective[
            "acc_aw_rotation_factor_upper"
        ],
        "acc_effective_aw_input_norm_upper_mps2": effective[
            "acc_effective_aw_input_norm_upper_mps2"
        ],
        "accelerometer_bias_standalone_nonlinear_penalty": 0.0,
        "standalone_vector_eta_penalty_active": False,
        "condition_number_conversion_inserted_between_P3_and_P4": False,
        "positive_source_correlated_word_form_built_here": False,
        "signed_word_generalized_margin_lower": None,
        "rho_full_nonlinear_word_upper": None,
        "P4_PROMOTED": False,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    p3 = P3.build(path)
    p3f = P3.validate(p3)
    cayley = CAYLEY.build(path)
    cf = CAYLEY.validate(cayley)
    timing = TIMING.build(path)
    tf = TIMING.validate(timing)
    effective = EFFECTIVE.build(path, float(cayley["outer_angle_rad"]))
    ef = EFFECTIVE.validate(effective)

    failures = [f"P3: {x}" for x in p3f]
    failures += [f"Cayley: {x}" for x in cf]
    failures += [f"timing: {x}" for x in tf]
    failures += [f"effective-input: {x}" for x in ef]

    modes = {}
    if not failures:
        modes = {mode: _mode(mode, p3, effective) for mode in ("H", "A")}

    ready = not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_EFFECTIVE_INPUT_SIGNED_JOSEPH_COMPOSITION_CONTRACT",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "outer_angle_rad": cayley.get("outer_angle_rad"),
        "declared_filter_entrance_covered": cayley.get(
            "declared_filter_entrance_covered"
        ),
        "cayley_radius_upper": cayley.get("cayley_radius_upper"),
        "exact_vector_information_retention_factor_lower": cayley.get(
            "exact_vector_information_retention_factor_lower"
        ),
        "S_timing_consumed_by_linear_P3": timing.get(
            "S_timing_consumed_by_linear_P3_translation_UCO"
        ),
        "S_nonlinear_eta_identically_zero": timing.get(
            "S_nonlinear_eta_identically_zero"
        ),
        "effective_vector_input_certificate": {
            "qualification": effective.get("qualification"),
            "mag_radial_residual_gain_null_exact": effective.get(
                "mag_radial_residual_gain_null_exact"
            ),
            "acc_eta_in_aw_measurement_range_exact": effective.get(
                "acc_eta_in_aw_measurement_range_exact"
            ),
            "standalone_vector_eta_penalty_retired_from_active_word_route": effective.get(
                "standalone_vector_eta_penalty_retired_from_active_word_route"
            ),
        },
        "modes": modes,
        "P3_ESTABLISHED": not p3f,
        "P4_COMPOSITION_PREREQUISITES_ESTABLISHED": ready,
        "full_source_correlated_word_form_remaining": True,
        "same_information_metric_retained": True,
        "condition_number_conversion_inserted_between_P3_and_P4": False,
        "standalone_vector_eta_penalty_active": False,
        "fixed_terminal_schedule_used": False,
        "interval_AD_long_prefix_used": False,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED": False,
        "P5_FINITE_CAPTURE_ESTABLISHED": False,
        "next_obligation": (
            "propagate each source-complete H/A vector word with jointly reachable P,H,R,S,r, "
            "the exact magnetic effective coordinate and accelerometer a_w effective input, "
            "including immediate quaternion injection/reset and prefix safety; certify a strict "
            "generalized endpoint margin mu>0 (rho<1) in M_i=s_m Sigma_i^-1 over theta<=0.8 rad"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "declared_filter_entrance_covered",
        "S_timing_consumed_by_linear_P3",
        "S_nonlinear_eta_identically_zero",
        "P3_ESTABLISHED",
        "P4_COMPOSITION_PREREQUISITES_ESTABLISHED",
        "full_source_correlated_word_form_remaining",
        "same_information_metric_retained",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "condition_number_conversion_inserted_between_P3_and_P4",
        "standalone_vector_eta_penalty_active",
        "fixed_terminal_schedule_used",
        "interval_AD_long_prefix_used",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED",
        "P5_FINITE_CAPTURE_ESTABLISHED",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    ev = d.get("effective_vector_input_certificate", {})
    if ev.get("mag_radial_residual_gain_null_exact") is not True:
        failures.append("magnetic radial-null identity missing")
    if ev.get("acc_eta_in_aw_measurement_range_exact") is not True:
        failures.append("accelerometer effective a_w identity missing")
    if ev.get("standalone_vector_eta_penalty_retired_from_active_word_route") is not True:
        failures.append("standalone vector eta penalty was not retired")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {})
        if m.get("P3_prefix_information_gain_upper") != 1.0:
            failures.append(f"{mode} P3 prefix gain changed")
        if m.get("standalone_vector_eta_penalty_active") is not False:
            failures.append(f"{mode} standalone eta penalty became active")
        if m.get("condition_number_conversion_inserted_between_P3_and_P4") is not False:
            failures.append(f"{mode} condition-number conversion inserted")
        if m.get("positive_source_correlated_word_form_built_here") is not False:
            failures.append(f"{mode} falsely claims complete word form")
        if m.get("P4_PROMOTED") is not False:
            failures.append(f"{mode} falsely promotes P4")
        q = m.get("mag_effective_coordinate_tangent_defect_factor_upper")
        if not isinstance(q, (int, float)) or not math.isfinite(float(q)) or not (0.0 <= float(q) < 1.0):
            failures.append(f"{mode} magnetic tangent defect bound invalid")
    return list(dict.fromkeys(failures))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "P3": d["P3_ESTABLISHED"],
        "composition_prerequisites": d["P4_COMPOSITION_PREREQUISITES_ESTABLISHED"],
        "outer_angle_rad": d["outer_angle_rad"],
        "H": d.get("modes", {}).get("H"),
        "A": d.get("modes", {}).get("A"),
        "P4": d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"],
        "next_obligation": d["next_obligation"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
