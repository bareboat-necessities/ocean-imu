#!/usr/bin/env python3
"""Source-correlated exact transport contract for the usable OU-III P4 word.

This is the semantic/numerical boundary of the retained broad-sector route.
It deliberately does not restore any deleted fixed-schedule or long-prefix
interval-AD producer.  Instead it binds the current source-complete word
language to the exact operation classes that a future numerical word enclosure
must propagate.

The key structural point is that one accepted accelerometer+magnetometer packet
cannot have a positive scalar full-state information margin.  In the rotation
gauge the measurement-active H coordinates are (theta,a_w), with

    H_a = [-[f]x  I],       H_m = [-[m]x  0].

H_a has rank 3, H_m has rank 2, and the stacked packet has the exact nonzero
null family

    dtheta = alpha m,       da_w = [f]x dtheta.

Hence its exact rank is five on the 6-D H active block, and still five on the
9-D A active block after adding b_a.  P3 obtains full word detectability only
because prediction transports those directional nullspaces and recurrent
vector/S information accumulates.  P4 must preserve that structure.

For accepted vector updates this route uses the exact effective-input identities
certified by :mod:`ou3_p4_effective_vector_inputs`: magnetic radial residual is
annihilated by the Kalman gain, while accelerometer finite-angle nonlinearity is
represented in the a_w measurement range.  S=0 is exactly linear.  Joseph
covariance updates, immediate quaternion injection, and immediate left-error
reset remain sequential source operations.  Rejected/not-due updates are exact
identity corrections.  The pending a_w covariance synchronization is a PSD
increment and is information-nonexpansive at fixed physical error.

The numerical successor must therefore accumulate source-correlated directional
PSD forms through prediction/reset/effective coordinates over a complete H/A
word, and only then take a generalized scalar margin in the same P3 metric
M_i=s_m Sigma_i^-1.  This module fails closed until that margin is actually
computed and strictly positive.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_implementation_proof_manifest as MANIFEST
import ou3_implementation_word_language as WORDS
import ou3_p3_source_uniform_certificate as P3
import ou3_p4_effective_vector_inputs as EFFECTIVE
import ou3_p4_joint_joseph as JOSEPH
import ou3_p4_source_word_timing as TIMING

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1

EXPECTED_ORDER = [
    "commit_previous_tune",
    "vibration_guard_conditioning",
    "prediction",
    "apply_pending_aw_covariance_psd_increment",
    "periodic_S_zero_when_due_then_immediate_quaternion_injection_and_left_error_reset",
    "accelerometer_correction_or_rejection_then_immediate_quaternion_injection_and_left_error_reset_if_accepted",
    "source_tuner_evolution_and_stage_next_tune",
    "periodic_aw_covariance_sync_tick_stages_future_psd_increment",
    "asynchronous_magnetometer_correction_or_rejection_then_immediate_quaternion_injection_and_left_error_reset_if_accepted",
]


def _operation_calculus() -> list[dict]:
    return [
        {
            "operation": "commit_previous_tune",
            "state_map": "identity_on_error_state",
            "covariance_map": "source_parameter_commit_only",
            "nonlinear_budget": "none; committed parameters remain one jointly reachable source tuple",
        },
        {
            "operation": "vibration_guard_conditioning",
            "state_map": "identity_on_error_state",
            "measurement_map": "acc_in=acc",
            "covariance_map": "identity",
            "certified_branch": "dormant_zero_engagement_bit_exact_transparent_only",
            "active_or_transitioning_guard_covered": False,
            "active_or_transitioning_guard_requires_separate_source_certificate": True,
        },
        {
            "operation": "prediction",
            "state_map": "shipping fixed-mode tangent prediction plus exact finite-angle Cayley transport",
            "covariance_map": "F P F^T + Q",
            "directional_role": "transport packet nullspaces and accumulated PSD information",
        },
        {
            "operation": "apply_pending_aw_covariance_psd_increment",
            "state_map": "identity",
            "covariance_map": "P_awaw <- P_awaw + Delta_aw, Delta_aw >= 0",
            "information_role": "nonexpansive_at_fixed_physical_error",
        },
        {
            "operation": "periodic_S_zero_when_due_then_immediate_quaternion_injection_and_left_error_reset",
            "branch_family": ["not_due", "due"],
            "accepted_state_map": "full implemented K_S correction then immediate quaternion injection/reset",
            "covariance_map": "Joseph then immediate left-error reset when due",
            "nonlinear_measurement_eta": "IDENTICALLY_ZERO",
            "timing_role": "uncertain firing time remains in linear P3 translational UCO",
        },
        {
            "operation": "accelerometer_correction_or_rejection_then_immediate_quaternion_injection_and_left_error_reset_if_accepted",
            "branch_family": ["accepted", "rejected"],
            "accepted_state_map": "r_a=H_a(z+E_aw e_eta), full K_a, immediate deployed quaternion injection/reset",
            "effective_input": "e_eta=R_wb^T eta_a in the a_w coordinate",
            "effective_input_isometry": True,
            "standalone_eta_penalty_active": False,
            "covariance_map": "Joseph then immediate left-error reset when accepted",
        },
        {
            "operation": "source_tuner_evolution_and_stage_next_tune",
            "state_map": "no retroactive estimator-state correction",
            "covariance_map": "future source parameters only",
            "source_role": "joint source tuple evolves; Cartesian extrema products forbidden",
        },
        {
            "operation": "periodic_aw_covariance_sync_tick_stages_future_psd_increment",
            "state_map": "identity",
            "covariance_map": "future PSD increment only",
        },
        {
            "operation": "asynchronous_magnetometer_correction_or_rejection_then_immediate_quaternion_injection_and_left_error_reset_if_accepted",
            "branch_family": ["not_due", "accepted", "rejected"],
            "accepted_state_map": "useful residual H_theta d_eff, full K_m, immediate deployed quaternion injection/reset",
            "radial_residual_gain_action": "EXACTLY_ZERO",
            "effective_coordinate_nonexpansive": True,
            "standalone_eta_penalty_active": False,
            "covariance_map": "Joseph then immediate left-error reset when accepted",
        },
    ]


def _rank_certificate() -> dict:
    return {
        "rotation_gauge_accelerometer": "H_a=[-[f]_x, I_aw]",
        "rotation_gauge_magnetometer": "H_m=[-[m]_x, 0_aw]",
        "accelerometer_rank_exact": 3,
        "magnetometer_rank_exact": 2,
        "stacked_vector_packet_rank_exact": 5,
        "S_zero_rank_exact_when_due": 3,
        "stacked_vector_plus_due_S_rank_exact": 8,
        "exact_vector_packet_null_witness": {
            "parameter": "alpha != 0",
            "delta_theta": "alpha*m",
            "delta_a_w": "[f]_x*delta_theta",
            "magnetometer_residual": "-[m]_x*delta_theta=0",
            "accelerometer_residual": "-[f]_x*delta_theta+delta_a_w=0",
        },
        "H_vector_active_dimension": 6,
        "H_vector_packet_nullity_exact": 1,
        "A_vector_active_dimension": 9,
        "A_vector_packet_nullity_exact": 4,
        "instantaneous_positive_scalar_full_state_packet_margin_possible": False,
        "directional_PSD_word_accumulation_required": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("effective word transport domain must not be trajectory fitted")

    manifest = MANIFEST.build()
    words = WORDS.build(path)
    p3 = P3.build(path)
    effective = EFFECTIVE.build(path)
    timing = TIMING.build(path)

    failures = [f"manifest: {x}" for x in MANIFEST.validate(manifest)]
    failures += [f"word-language: {x}" for x in WORDS.validate(words)]
    failures += [f"P3: {x}" for x in P3.validate(p3)]
    failures += [f"effective-input: {x}" for x in EFFECTIVE.validate(effective)]
    failures += [f"timing: {x}" for x in TIMING.validate(timing)]

    order = list(manifest.get("normal_live_update_order", []))
    if order != EXPECTED_ORDER:
        failures.append("shipping normal-Live operation order changed")
    reset = manifest.get("same_sample_reset_policy", {})
    if reset.get("single_shared_end_of_sample_reset") is not False:
        failures.append("same-sample corrections were incorrectly merged into one reset")
    guard = manifest.get("vibration_guard", {})
    if guard.get("zero_engagement_is_bit_exact_transparent") is not True:
        failures.append("dormant vibration-guard branch is not source-certified transparent")
    if guard.get("active_guard_requires_separate_source_certificate") is not True:
        failures.append("active vibration-guard branch is not fail-closed")

    calculus = _operation_calculus()
    if [x["operation"] for x in calculus] != order:
        failures.append("operation calculus does not exactly match source operation order")

    rank = _rank_certificate()
    live = domain.get("normal_live", {})
    if not (
        float(live.get("specific_force_norm_lower_mps2", 0.0)) > 0.0
        and float(live.get("magnetic_vector_norm_lower_uT", 0.0)) > 0.0
        and 0.0 < float(live.get("vector_sine_separation_lower", 0.0)) < 1.0
    ):
        failures.append("declared vector geometry lost the nonzero/noncollinear rank premise")

    if effective.get("mag_radial_residual_gain_null_exact") is not True:
        failures.append("magnetic radial gain-null identity missing")
    if effective.get("acc_eta_in_aw_measurement_range_exact") is not True:
        failures.append("accelerometer effective a_w range identity missing")
    if timing.get("S_nonlinear_eta_identically_zero") is not True:
        failures.append("S=0 nonlinear remainder is not exact zero")

    word_contract = words["word_contract"]
    cw = word_contract["conditional_word_language"]
    dims = manifest["state_coordinates"]
    modes = {}
    for mode, n in (("H", 18), ("A", 21)):
        pm = p3["modes"][mode]
        modes[mode] = {
            "dimension": n,
            "fixed_dimension_inside_word": True,
            "P3_relative_Riccati_injection_margin_lower": pm[
                "relative_Riccati_injection_margin_lower"
            ],
            "P3_prefix_information_gain_upper": pm["prefix_information_gain_upper"],
            "same_information_metric": "M_i=s_m Sigma_i^-1",
            "full_attitude_linear_cross_terms_retained": True,
            "vector_packet_rank_exact": 5,
            "vector_active_nullity_exact": 1 if mode == "H" else 4,
            "word_directional_forms_numerically_accumulated_here": False,
            "strict_generalized_word_margin_lower": None,
            "rho_full_nonlinear_word_upper": None,
            "P4_PROMOTED": False,
        }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_SOURCE_CORRELATED_EFFECTIVE_WORD_TRANSPORT_CONTRACT",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "source_word_horizon_s": cw["word_horizon_lower_s"],
        "word_samples_upper": cw["word_samples_upper_at_configured_dt"],
        "shipping_operation_order": order,
        "operation_transport_calculus": calculus,
        "measurement_directional_rank": rank,
        "H_dimension": dims["H_dimension"],
        "A_dimension": dims["A_dimension"],
        "joint_source_reachability_required": word_contract["source_branch_language"][
            "joint_source_reachability_required"
        ],
        "cartesian_extrema_products_valid": False,
        "S_timing_consumed_by_linear_P3": timing[
            "S_timing_consumed_by_linear_P3_translation_UCO"
        ],
        "standalone_vector_eta_penalty_active": False,
        "condition_number_conversion_inserted_between_P3_and_P4": False,
        "reset_condition_number_multiplier_used": False,
        "per_packet_scalarization_allowed": False,
        "directional_PSD_word_accumulation_required": True,
        "modes": modes,
        "P4_WORD_TRANSPORT_SEMANTICS_ESTABLISHED": not failures,
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED": False,
        "P5_FINITE_CAPTURE_ESTABLISHED": False,
        "next_obligation": (
            "materialize the jointly reachable H/A source tuple at each operation; form the actual Joseph S^-1 directional PSD credit, "
            "transport it through prediction and exact immediate reset/effective-coordinate maps, accumulate all recurrent vector/S directions "
            "over the complete word, then compute the generalized endpoint margin in M_i=s_m Sigma_i^-1 and prefix overshoot before scalarization"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "joint_source_reachability_required",
        "S_timing_consumed_by_linear_P3",
        "directional_PSD_word_accumulation_required",
        "P4_WORD_TRANSPORT_SEMANTICS_ESTABLISHED",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "source_replay_used",
        "filter_changed",
        "cartesian_extrema_products_valid",
        "standalone_vector_eta_penalty_active",
        "condition_number_conversion_inserted_between_P3_and_P4",
        "reset_condition_number_multiplier_used",
        "per_packet_scalarization_allowed",
        "P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED",
        "P5_FINITE_CAPTURE_ESTABLISHED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("shipping_operation_order") != EXPECTED_ORDER:
        f.append("shipping operation order mismatch")
    rank = d.get("measurement_directional_rank", {})
    if rank.get("accelerometer_rank_exact") != 3:
        f.append("accelerometer rank is not exactly three")
    if rank.get("magnetometer_rank_exact") != 2:
        f.append("magnetometer rank is not exactly two")
    if rank.get("stacked_vector_packet_rank_exact") != 5:
        f.append("vector packet rank is not exactly five")
    if rank.get("instantaneous_positive_scalar_full_state_packet_margin_possible") is not False:
        f.append("impossible instantaneous scalar packet margin reintroduced")
    if rank.get("directional_PSD_word_accumulation_required") is not True:
        f.append("directional word accumulation requirement missing")
    for mode, n, nullity in (("H", 18, 1), ("A", 21, 4)):
        m = d.get("modes", {}).get(mode, {})
        if m.get("dimension") != n:
            f.append(f"{mode} dimension mismatch")
        if m.get("vector_packet_rank_exact") != 5:
            f.append(f"{mode} vector rank mismatch")
        if m.get("vector_active_nullity_exact") != nullity:
            f.append(f"{mode} vector nullity mismatch")
        if m.get("P3_prefix_information_gain_upper") != 1.0:
            f.append(f"{mode} P3 prefix gain changed")
        if m.get("word_directional_forms_numerically_accumulated_here") is not False:
            f.append(f"{mode} falsely claims numerical directional accumulation")
        if m.get("strict_generalized_word_margin_lower") is not None:
            f.append(f"{mode} invents an uncomputed word margin")
        if m.get("rho_full_nonlinear_word_upper") is not None:
            f.append(f"{mode} invents an uncomputed rho")
        if m.get("P4_PROMOTED") is not False:
            f.append(f"{mode} falsely promotes P4")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": "PASS" if not vf else "FAIL",
        "rank": d["measurement_directional_rank"],
        "H": d["modes"]["H"],
        "A": d["modes"]["A"],
        "next": d["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
