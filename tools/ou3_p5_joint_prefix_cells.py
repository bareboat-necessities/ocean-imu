#!/usr/bin/env python3
"""Active joint P5 prefix interface backed by the V3 full 18x18 H cell.

The former version of this producer stopped after a directional P envelope and
therefore could not form the signed correction direction required by the exact
Cayley denominator.  The active interface now delegates covariance and
correction propagation to :mod:`ou3_p5_full_h_prefix_cells_v3`, which carries
an outward full 18x18 H covariance cell, uses the dependency-preserving OU
kernel bounds, composes the deployed quaternion correction before forming the
resulting Cayley coordinate, and recomputes P,H,R,S,K,r,d_eff at every later
prediction/vector/S prefix.

The old scalar/directional calculation and the V1/V2 experimental backends are
intentionally not retained as fallback theorem routes.  The exact tangent-only
magnetic identity, effective accelerometer a_w input, signed group composition,
Joseph update and immediate reset congruence are mandatory properties of the
active backend.  If the broad source-complete cell cannot close q<=8, this
interface reports its first V3 full-matrix obstruction and remains
NOT_ESTABLISHED; it does not revert to an older norm-only or 3-rad route.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_effective_vector_input as VEFF
import ou3_p5_first_s_exact_prefix as FIRSTS
import ou3_p5_full_h_prefix_cells_v3 as FULL
import ou3_p5_mag_information_reduction as MAGINFO
import ou3_p5_signed_cayley_cell as SIGNED

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("joint full-matrix prefix domain must not be trajectory fitted")

    first = FIRSTS.build(domain_path)
    veff = VEFF.build(domain_path)
    mag = MAGINFO.build(domain_path)
    signed = SIGNED.build(domain_path)
    full = FULL.build(domain_path)

    failures = [f"first-S: {x}" for x in FIRSTS.validate(first)]
    failures += [f"effective-vector: {x}" for x in VEFF.validate(veff)]
    failures += [f"mag-information: {x}" for x in MAGINFO.validate(mag)]
    failures += [f"signed-Cayley: {x}" for x in SIGNED.validate(signed)]
    failures += [f"full-H-prefix-v3: {x}" for x in FULL.validate(full)]

    q8_closed = bool(full["complete_q_le_8_prefix_family_closed"])
    matrix_status = full["P5_FULL_H_PREFIX_MATRIX_CERTIFICATE"]
    promoted = not failures and q8_closed and matrix_status == "PASS"
    obstruction = full.get("first_failure")
    if promoted:
        first_unclosed = "NONE_AT_COMPLETE_GAUGED_H_PREFIX"
    elif obstruction:
        first_unclosed = "FULL_MATRIX_PREFIX_SUBDIVISION_REQUIRED"
    else:
        first_unclosed = "FULL_MATRIX_COMPLETE_WORD_NOT_CERTIFIED"

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_ACTIVE_JOINT_FULL_MATRIX_PREFIX_INTERFACE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "active_P_payload": "OUTWARD_FULL_18X18_H_COVARIANCE_CELL",
        "active_full_matrix_backend": full["active_full_matrix_backend"],
        "active_backend_is_v3_deployed_quaternion": True,
        "directional_P_payload_retained_as_active_backend": False,
        "old_directional_scalar_route_used_for_promotion": False,
        "independent_global_extrema_product_used": False,
        "full_signed_matrix_covariance_cells_available": full["full_18x18_covariance_propagated"],
        "P_H_R_K_S_r_d_eff_recomputed_in_same_prefix_cell": full["H_R_S_K_r_d_eff_recomputed_in_same_prefix_cell"],
        "shipping_Joseph_update_used": full["shipping_Joseph_update_used"],
        "immediate_left_error_reset_congruence_used": full["immediate_left_error_reset_congruence_used"],
        "physical_attitude_correction_is_minus_Etheta_Kr": full["physical_attitude_correction_is_minus_Etheta_Kr"],
        "signed_cayley_primitive_consumes_actual_interval_d": full["signed_cayley_primitive_consumes_actual_interval_d"],
        "signed_a_dot_c_replaced_by_independent_abs_product": full["signed_a_dot_c_replaced_by_independent_abs_product"],
        "magnetometer_radial_K_action_exact_zero": full["magnetometer_radial_K_action_exact_zero"],
        "magnetometer_radial_Joseph_information_exact_zero": True,
        "standalone_vector_eta_penalty_used": full["standalone_vector_eta_penalty_used"],
        "accelerometer_effective_aw_input_used": full["accelerometer_effective_aw_input_used"],
        "source_complete_rejection_identity_hulls": full["source_complete_rejection_identity_hulls"],
        "first_due_S_seed": {
            "status": first["P5_FIRST_DUE_S_EXACT_CAYLEY_PREFIX_CERTIFICATE"],
            "widened_cayley_norm_upper": float(first["widened_prefix_cayley_norm_upper"]),
            "required_post_cayley_norm_upper": float(first["required_first_S_post_cayley_norm_upper"]),
        },
        "full_matrix_prefix": {
            "status": matrix_status,
            "active_full_matrix_backend": full["active_full_matrix_backend"],
            "maximum_validated_deployed_correction_norm_rad": full["maximum_validated_deployed_correction_norm_rad"],
            "correction_norm_three_rad_is_promotion_gate": full["correction_norm_three_rad_is_promotion_gate"],
            "deployed_quaternion_composed_before_result_cayley": full["deployed_quaternion_composed_before_result_cayley"],
            "source_cell": full["source_cell"],
            "word_samples_upper": full["word_samples_upper"],
            "inverse_backend_counts": full["inverse_backend_counts"],
            "q_chart_upper": full["q_chart_upper"],
            "max_reached_cayley_norm_upper": full["max_reached_cayley_norm_upper"],
            "smaller_source_reachable_chart_upper": full["smaller_source_reachable_chart_upper"],
            "numerical_extrema": full["numerical_extrema"],
            "last_prefix_cells": full["last_prefix_cells"],
            "final_covariance": full["final_covariance"],
            "first_failure": obstruction,
        },
        "effective_vector_reductions": {
            "magnetometer_exact_state_correction_identity": veff["magnetometer"]["exact_state_correction_identity"],
            "magnetometer_radial_gain_action_exact_zero": veff["magnetometer"]["kalman_gain_radial_action_exact_zero"],
            "accelerometer_exact_state_correction_identity": veff["accelerometer"]["exact_state_correction_identity"],
            "standalone_vector_eta_penalty_retired": veff["standalone_vector_eta_penalty_retired_from_P5_numerical_route"],
            "magnetometer_information_reduction_status": mag["P5_MAGNETOMETER_INFORMATION_REDUCTION_CERTIFICATE"],
        },
        "signed_cayley": {
            "primitive_status": signed["P5_SIGNED_CAYLEY_CELL_PRIMITIVE"],
            "signed_a_dot_c_retained": signed["signed_a_dot_c_retained"],
            "independent_abs_a_abs_c_denominator_used": signed["independent_abs_a_abs_c_denominator_used"],
            "complete_q_le_8_prefix_family_closed": q8_closed,
        },
        "signed_cayley_prefix_composition_closed": q8_closed,
        "P5_JOINT_PREFIX_FULL_MATRIX_CERTIFICATE": "PASS" if promoted else "NOT_ESTABLISHED",
        "P5_JOINT_PREFIX_SCALAR_CELL_CERTIFICATE": "RETIRED_AS_ACTIVE_ROUTE",
        "P5_numerical_status_can_promote_from_this_stage": promoted,
        "N_H_words_set_here": False,
        "first_unclosed_numerical_obligation": first_unclosed,
        "next_obligation": (
            "compose the certified full-matrix H word with the P4 inner overlap and set N_H_words"
            if promoted else
            "subdivide the reported full-matrix prefix obstruction without dropping P direction, effective-vector identities, signed group composition, Joseph, or immediate reset congruence"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit", "active_backend_is_v3_deployed_quaternion",
        "full_signed_matrix_covariance_cells_available",
        "P_H_R_K_S_r_d_eff_recomputed_in_same_prefix_cell", "shipping_Joseph_update_used",
        "immediate_left_error_reset_congruence_used", "physical_attitude_correction_is_minus_Etheta_Kr",
        "signed_cayley_primitive_consumes_actual_interval_d", "magnetometer_radial_K_action_exact_zero",
        "magnetometer_radial_Joseph_information_exact_zero", "accelerometer_effective_aw_input_used",
        "source_complete_rejection_identity_hulls",
    ):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in (
        "source_replay_used", "filter_changed", "directional_P_payload_retained_as_active_backend",
        "old_directional_scalar_route_used_for_promotion", "independent_global_extrema_product_used",
        "signed_a_dot_c_replaced_by_independent_abs_product", "standalone_vector_eta_penalty_used",
        "N_H_words_set_here",
    ):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")
    if d.get("active_P_payload") != "OUTWARD_FULL_18X18_H_COVARIANCE_CELL":
        failures.append("active P payload is not the full matrix cell")
    if "DEPLOYED_QUATERNION_COMPOSITION" not in str(d.get("active_full_matrix_backend", "")):
        failures.append("active joint backend is not V3 deployed-quaternion backend")
    if d.get("P5_JOINT_PREFIX_SCALAR_CELL_CERTIFICATE") != "RETIRED_AS_ACTIVE_ROUTE":
        failures.append("old scalar route was not retired")
    sc = d.get("signed_cayley", {})
    if sc.get("signed_a_dot_c_retained") is not True or sc.get("independent_abs_a_abs_c_denominator_used") is not False:
        failures.append("signed Cayley denominator semantics changed")
    fm = d.get("full_matrix_prefix", {})
    if float(fm.get("maximum_validated_deployed_correction_norm_rad", 0.0)) < 6.0:
        failures.append("joint report regressed to pre-V3 correction range")
    if fm.get("correction_norm_three_rad_is_promotion_gate") is not False:
        failures.append("joint report restored retired 3-rad gate")
    if fm.get("deployed_quaternion_composed_before_result_cayley") is not True:
        failures.append("joint report is not using deployed quaternion composition")
    q = fm.get("max_reached_cayley_norm_upper")
    if not isinstance(q, (int, float)) or not math.isfinite(float(q)) or float(q) < 0.0:
        failures.append("full matrix prefix did not emit a finite reached Cayley bound")
    if d.get("P5_JOINT_PREFIX_FULL_MATRIX_CERTIFICATE") == "PASS":
        if d.get("signed_cayley_prefix_composition_closed") is not True:
            failures.append("full matrix prefix promoted without signed composition closure")
        if d.get("P5_numerical_status_can_promote_from_this_stage") is not True:
            failures.append("PASS full matrix prefix not marked promotable")
    else:
        if d.get("P5_numerical_status_can_promote_from_this_stage") is not False:
            failures.append("non-PASS full matrix prefix marked promotable")
        if d.get("first_unclosed_numerical_obligation") == "NONE_AT_COMPLETE_GAUGED_H_PREFIX":
            failures.append("non-PASS full matrix prefix has no obstruction")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_JOINT_PREFIX_FULL_MATRIX_CERTIFICATE"],
        "q8_closed": out["signed_cayley_prefix_composition_closed"],
        "full_matrix": out["full_matrix_prefix"],
        "next": out["first_unclosed_numerical_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
