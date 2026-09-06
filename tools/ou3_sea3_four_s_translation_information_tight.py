#!/usr/bin/env python3
"""Tight covariance read of the corrected complete-SEA3 four-S information lemma.

This is not a second source, word, estimator or proof architecture.  It invokes
``ou3_sea3_four_s_translation_information`` unchanged and tightens one PSD
covariance inequality that was deliberately loose there.

For the four selected actual S records, configured measurement noise is
independent between updates.  Therefore its 4x4 covariance is diagonal and

    lambda_max(Sigma_R) <= max_i R_{S,i}^2,

not four times that value.  Process nuisance may be correlated across records;
for it we retain the conservative PSD trace bound

    lambda_max(Sigma_Q) <= trace(Sigma_Q)
                        <= 4 q_S,max.

Hence

    Sigma_Y <= (R_max^2 + 4 q_S,max) I.

The corrected base lemma already includes the physical map from Newton divided
differences back to [S,gp,g^2v,g^3 a_w], including the lower third divided
difference of the OU response.  This module changes only the covariance upper
used with that physical matrix.  It cannot promote P3 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_sea3_four_s_translation_information as BASE

DEFAULT_DOMAIN = BASE.DEFAULT_DOMAIN
SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_FOUR_S_PHYSICAL_INFORMATION_TIGHT_COVARIANCE"


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    d = BASE.build(Path(domain_path).resolve())
    failures = BASE.validate(d)
    if failures:
        raise RuntimeError(f"base four-S physical information invalid: {failures}")

    noise = dict(d["selected_S_record_noise"])
    measurement_lambda = up(max(map(float, noise["measurement_variance_axis_upper"])))
    process_per_record = float(noise["process_S_variance_per_record_upper"])
    process_stack_lambda = up(4.0 * process_per_record)
    total_lambda = up(measurement_lambda + process_stack_lambda)
    inv_lower = down(1.0 / total_lambda)

    ni = dict(d["newton_coordinate_information"])
    physical = ni["physical_state_recovery"]
    mtm = float(physical["physical_observation_MtM_lambda_min_lower"])
    info = down(inv_lower * mtm)
    if not (math.isfinite(info) and info > 0.0):
        raise RuntimeError("tight four-S physical information lower is not strict")

    loose_lambda = float(noise["four_record_covariance_lambda_max_upper"])
    loose_info = float(ni["D_S_physical_lambda_min_lower"])
    noise.update({
        "legacy_loose_four_record_covariance_lambda_max_upper_diagnostic": loose_lambda,
        "measurement_covariance_is_diagonal_across_selected_updates": True,
        "measurement_covariance_lambda_max_upper": measurement_lambda,
        "process_covariance_may_be_correlated_across_selected_updates": True,
        "process_four_record_covariance_lambda_max_trace_upper": process_stack_lambda,
        "four_record_covariance_lambda_max_upper": total_lambda,
        "Sigma_Y_inverse_scalar_lower": inv_lower,
        "measurement_noise_not_multiplied_by_record_count": True,
        "process_cross_record_correlation_still_covered_by_trace_bound": True,
    })
    ni.update({
        "raw_record_inverse_covariance_scalar_lower": inv_lower,
        "D_S_physical_lambda_min_lower": info,
        "D_S_newton_lambda_min_lower": info,
        "D_S_newton_matrix_lower": f"D_S,z >= {info:.17g} * I_4",
        "tight_covariance_bound_consumed": True,
        "base_corrected_physical_information_lower_diagnostic": loose_info,
    })

    out = dict(d)
    out.update({
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "selected_S_record_noise": noise,
        "newton_coordinate_information": ni,
        "same_complete_SEA3_component_sharpened": True,
        "source_family_replaced": False,
        "P3_architecture_replaced": False,
        "P3_promoted": False,
    })
    return out


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "component_of_complete_SEA3_full_word",
        "actual_applied_SpectralMSE_R_S_consumed",
        "all_due_S_updates_remain_in_literal_word",
        "same_complete_SEA3_component_sharpened",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in ("source_family_replaced", "P3_architecture_replaced", "P3_promoted"):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    noise = d.get("selected_S_record_noise", {})
    for key in (
        "measurement_covariance_is_diagonal_across_selected_updates",
        "process_covariance_may_be_correlated_across_selected_updates",
        "measurement_noise_not_multiplied_by_record_count",
        "process_cross_record_correlation_still_covered_by_trace_bound",
    ):
        if noise.get(key) is not True:
            f.append(f"noise structure missing {key}")
    tight = float(noise.get("four_record_covariance_lambda_max_upper", math.inf))
    loose = float(noise.get("legacy_loose_four_record_covariance_lambda_max_upper_diagnostic", 0.0))
    if not (0.0 < tight < loose):
        f.append("tight covariance upper did not strictly improve the loose trace bound")
    ni = d.get("newton_coordinate_information", {})
    if ni.get("third_divided_difference_lower_enters_quantitative_bound") is not True:
        f.append("physical divided-difference correction was lost")
    if ni.get("tight_covariance_bound_consumed") is not True:
        f.append("tight covariance bound not consumed")
    info = ni.get("D_S_physical_lambda_min_lower")
    old = ni.get("base_corrected_physical_information_lower_diagnostic")
    if not isinstance(info, (int, float)) or not (math.isfinite(float(info)) and float(info) > 0.0):
        f.append("tight physical information lower is not strict")
    if not isinstance(old, (int, float)) or not float(info) > float(old):
        f.append("tight physical information did not improve corrected base")
    return list(dict.fromkeys(f))


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
    n = d["selected_S_record_noise"]
    ni = d["newton_coordinate_information"]
    print(json.dumps({
        "measurement_cov_lambda_max": n["measurement_covariance_lambda_max_upper"],
        "process_stack_cov_lambda_max": n["process_four_record_covariance_lambda_max_trace_upper"],
        "total_cov_lambda_max": n["four_record_covariance_lambda_max_upper"],
        "base_physical_information": ni["base_corrected_physical_information_lower_diagnostic"],
        "tight_physical_information": ni["D_S_physical_lambda_min_lower"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
