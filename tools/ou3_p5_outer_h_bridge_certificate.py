#!/usr/bin/env python3
"""Compose the staged OU-III P5 outer-H capture bridge.

The bridge is fail-closed and consumes the corrected outer routes together with
the exact transport algebra now available for the complete H word.

For gauged H, the false raw V_R sector is retired.  Exact finite-angle
Cayley/vector information is positive on the normal and timeout handoff nodes,
the first due S correction is source-staged through the deployed quaternion map
into a finite widened Cayley chart, and signed Cayley correction cells preserve
the actual ``a^T c`` denominator rather than replacing it by an independent
norm product.  The exact Cayley eta geometry is retained only as an algebraic
support lemma: the active numerical route now uses
``ou3_p5_effective_vector_input``.  For the configured magnetometer the radial
finite-angle residual is annihilated exactly by the Kalman gain and the useful
residual is ``H_theta d_eff``; for the accelerometer the orthogonal full-rank
``J_aw=R_wb`` absorbs the nonlinear residual as a source-correlated effective
``a_w`` tangent input.  No standalone vector-eta information penalty remains in
the P5 route.

For ungauged H, the false yaw-only/full-gyro-bias contraction is retired.  The
corrected gravity quotient has strict transverse attitude/bias detectability,
retains the complete four-S translation word, and carries gravity-parallel gyro
bias as a bounded neutral input.  The same effective accelerometer-input and
exact Joseph/quaternion/reset calculus applies after quotient projection.

What remains is numerical rather than semantic: source-correlated ``P,H,R,K,r``
and effective correction-input cells must be propagated through all later
prefixes of the 1 s word, together with exact reset and prediction budgets.  No
filter/tuning value, S-to-attitude gain, or theorem gate is changed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_nonlinear_word_certificate as P4
import ou3_p5_cayley_eta_geometry as ETA
import ou3_p5_complete_word_transport as TRANSPORT
import ou3_p5_effective_vector_input as VEFF
import ou3_p5_first_s_exact_prefix as FIRSTEX
import ou3_p5_first_s_gain_certificate as FIRSTS
import ou3_p5_first_s_state_prefix_certificate as SPREFIX
import ou3_p5_go_live_covariance_stage as GOLIVE
import ou3_p5_gravity_quotient_certificate as GQUOT
import ou3_p5_heading_handoff_contract as HEADING
import ou3_p5_large_angle_sector_certificate as SECTOR
import ou3_p5_outer_h_word_certificate as OUTER
import ou3_p5_outer_information_geometry as OUTINFO
import ou3_p5_signed_cayley_cell as SIGNED
import ou3_p5_yaw_quotient_word_certificate as OLDQUOT
import ou3_startup_stability_certificate as P1

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 7


def _gauged_node(name: str, heading_row: dict, outer_name: str,
                 outer: dict, p4: dict, first_s: dict, s_prefix: dict,
                 sector: dict, outinfo: dict, firstex: dict) -> dict:
    o = outer["node_word_tests"][outer_name]
    cfull = float(heading_row["full_attitude_cosine_lower"])
    qfull = float(heading_row["full_attitude_cayley_norm_upper"])
    V_R_upper = math.nextafter(1.0 - cfull, math.inf)
    if not (0.0 <= V_R_upper < 2.0 and 0.0 < qfull < 1.0):
        raise RuntimeError(f"{name}: invalid gauged full-attitude handoff node")

    kS = float(first_s["K_thetaS_operator_norm_upper_first_due"])
    Sbound = float(s_prefix["first_due_S_error_norm_upper_m_s"])
    Sinjection = float(s_prefix["first_due_S_induced_attitude_correction_norm_upper_rad"])
    helper = float(s_prefix["deployed_group_helper_correction_limit_rad"])
    if not Sinjection < helper:
        raise RuntimeError(f"{name}: staged first-S correction left validated group helper")

    raw_status = sector["P5_RAW_VR_LARGE_ANGLE_SECTOR"]
    if raw_status != "DISPROVED_ON_DECLARED_SOURCE_FAMILY":
        raise RuntimeError("raw V_R sector audit did not produce its source counterexample")
    ir = outinfo["nodes"][name]
    if ir.get("strict") is not True:
        raise RuntimeError(f"{name}: finite-angle source information prerequisite is not strict")
    if not math.isclose(float(ir["cayley_norm_upper"]), qfull, rel_tol=1e-12, abs_tol=1e-12):
        raise RuntimeError(f"{name}: finite-angle information node does not match handoff node")
    sx = firstex["nodes"][name]
    if sx.get("inside_widened_prefix_chart") is not True:
        raise RuntimeError(f"{name}: first-S exact image is outside widened prefix chart")

    return {
        "node": name,
        "heading_branch": "GAUGED_FULL_HEADING",
        "P1_gravity_tilt_cosine_not_used_as_full_attitude_cosine": True,
        "full_attitude_cosine_lower": cfull,
        "handoff_cayley_norm_upper": qfull,
        "handoff_group_energy_V_R_upper": V_R_upper,
        "inside_antipodal_exclusion": cfull > -1.0,
        "inside_candidate_outer_cayley_bootstrap": qfull < float(s_prefix["candidate_outer_bootstrap"]["cayley_norm_upper"]),
        "goLive_S_to_attitude_gain_exact_zero": True,
        "global_P3_S_prefix_bound_used_for_P5_gate": False,
        "first_due_S_gain_certificate": "PASS",
        "first_due_K_thetaS_operator_norm_upper": kS,
        "first_due_S_state_prefix_certificate": "PASS_CONDITIONAL",
        "first_due_S_error_norm_upper_m_s": Sbound,
        "first_due_S_induced_attitude_correction_norm_upper_rad": Sinjection,
        "S_induced_correction_inside_group_helper": True,
        "first_due_S_exact_prefix_certificate": firstex["P5_FIRST_DUE_S_EXACT_CAYLEY_PREFIX_CERTIFICATE"],
        "first_due_S_post_cayley_norm_upper": sx["post_injection_cayley_norm_upper"],
        "inside_widened_first_S_chart": sx["inside_widened_prefix_chart"],
        "retired_global_outer_diagnostic_ratio_upper": float(o["outer_vector_nonlinear_information_ratio_upper"]),
        "retired_P3_translation_limited_sqrt_gap_lower": float(o["P3_homogeneous_sqrt_decrease_lower"]),
        "raw_V_R_large_angle_sector_status": raw_status,
        "raw_V_R_sector_witness_inside_node": float(sector["validated_counterexample"]["cayley_norm_interval"][1]) < qfull,
        "raw_V_R_sector_is_P5_promotion_route": False,
        "finite_angle_information_geometry_status": outinfo["P5_FINITE_ANGLE_INFORMATION_GEOMETRY_CERTIFICATE"],
        "exact_cayley_residual_factor_lower": ir["exact_cayley_residual_factor_lower"],
        "exact_pair_residual_information_per_cayley_norm_sq_lower": ir["exact_pair_residual_information_per_cayley_norm_sq_lower"],
        "exact_pair_residual_information_vs_goLive_attitude_metric_lower": ir["exact_pair_residual_information_vs_goLive_attitude_metric_lower"],
        "source_correlated_Joseph_information_identity_retained": True,
        "source_shaped_Cayley_information_outer_sector_required": True,
        "source_shaped_Cayley_information_outer_sector_status": "GEOMETRY_FIRST_S_AND_EFFECTIVE_VECTOR_INPUT_PASS_LATER_CELL_PROPAGATION_PENDING",
        "exact_group_backend_retained": True,
        "outer_S_state_prefix_status": "PASS_CONDITIONAL_ON_OUTER_NODE",
        "outer_prefix_domain_bootstrap_status": "PASS_THROUGH_FIRST_S_IN_WIDENED_CHART_LATER_PREFIXES_PENDING",
        "P4_inner_overlap_target_W": float(p4["modes"]["H"]["certified_level_W"]),
        "pass": False,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    p1 = P1.build(domain_path)
    p4 = P4.build(domain_path)
    outer = OUTER.build(domain_path)
    stage = GOLIVE.build(domain_path)
    first_s = FIRSTS.build(domain_path)
    s_prefix = SPREFIX.build(domain_path)
    firstex = FIRSTEX.build(domain_path)
    heading = HEADING.build(domain_path)
    sector = SECTOR.build(domain_path)
    oldq = OLDQUOT.build(domain_path)
    outinfo = OUTINFO.build(domain_path)
    gquot = GQUOT.build(domain_path)
    eta = ETA.build(domain_path)
    veff = VEFF.build(domain_path)
    signed = SIGNED.build(domain_path)
    transport = TRANSPORT.build(domain_path)

    prereq = []
    prereq += [f"P1: {x}" for x in P1.validate(p1)]
    prereq += [f"P4: {x}" for x in P4.validate(p4)]
    prereq += [f"outer-diagnostic: {x}" for x in OUTER.validate(outer)]
    prereq += [f"goLive-stage: {x}" for x in GOLIVE.validate(stage)]
    prereq += [f"first-S-gain: {x}" for x in FIRSTS.validate(first_s)]
    prereq += [f"first-S-state-prefix: {x}" for x in SPREFIX.validate(s_prefix)]
    prereq += [f"first-S-exact: {x}" for x in FIRSTEX.validate(firstex)]
    prereq += [f"heading-handoff: {x}" for x in HEADING.validate(heading)]
    prereq += [f"raw-VR-audit: {x}" for x in SECTOR.validate(sector)]
    prereq += [f"yaw-only-audit: {x}" for x in OLDQUOT.validate(oldq)]
    prereq += [f"outer-information: {x}" for x in OUTINFO.validate(outinfo)]
    prereq += [f"gravity-quotient: {x}" for x in GQUOT.validate(gquot)]
    prereq += [f"Cayley-eta-geometry: {x}" for x in ETA.validate(eta)]
    prereq += [f"effective-vector-input: {x}" for x in VEFF.validate(veff)]
    prereq += [f"signed-Cayley-cell: {x}" for x in SIGNED.validate(signed)]
    prereq += [f"complete-word-transport: {x}" for x in TRANSPORT.validate(transport)]
    if prereq:
        return {
            "schema": SCHEMA,
            "qualification": "OU3_P5_STAGED_OUTER_H_BRIDGE",
            "source_generated_not_trajectory_fit": True,
            "source_replay_used": False,
            "P5_OUTER_H_BRIDGE_CERTIFICATE": "NOT_ESTABLISHED",
            "first_failure": "UPSTREAM_PREREQUISITE_FAILURE",
            "failures": prereq,
        }

    seed = stage["goLive_H_covariance_seed"]
    if seed["theta_S_cross_covariance_operator_norm_upper"] != 0.0:
        raise RuntimeError("goLive theta-S covariance is not exact zero")
    if seed["S_to_attitude_gain_at_goLive_exact_zero"] is not True:
        raise RuntimeError("goLive S-to-attitude gain is not exact zero")

    nodes = {
        "normal_gauged": _gauged_node(
            "normal_gauged", heading["gauged_quality_handoff"], "normal",
            outer, p4, first_s, s_prefix, sector, outinfo, firstex,
        ),
        "timeout_gauged": _gauged_node(
            "timeout_gauged", heading["gauged_timeout_subbranch"], "timeout",
            outer, p4, first_s, s_prefix, sector, outinfo, firstex,
        ),
    }

    gauged_failure = transport["gauged_H"]["first_unclosed_numerical_obligation"]
    quotient_failure = transport["gravity_quotient_H"]["first_unclosed_numerical_obligation"]
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_STAGED_OUTER_H_BRIDGE",
        "claim": "SOURCE_STAGED_OUTER_H_CAPTURE_WITH_CORRECTED_GAUGED_AND_DETECTABLE_QUOTIENT_ROUTES",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "local_P4_recurrence_used_as_outer_bridge": False,
        "goLive_covariance_stage": {
            "status": stage["P5_GOLIVE_COVARIANCE_STAGE_CERTIFICATE"],
            "theta_S_cross_covariance_operator_norm_upper": 0.0,
            "S_to_attitude_gain_exact_zero": True,
            "pseudo_elapsed_s": seed["pseudo_update_elapsed_s_at_goLive"],
        },
        "first_due_S_gain_stage": {
            "status": first_s["P5_FIRST_DUE_S_GAIN_CERTIFICATE"],
            "theta_S_canonical_correlation_upper": first_s["theta_S_canonical_correlation_upper"],
            "K_thetaS_operator_norm_upper": first_s["K_thetaS_operator_norm_upper_first_due"],
        },
        "first_due_S_state_prefix_stage": {
            "status": s_prefix["P5_FIRST_S_STATE_PREFIX_CERTIFICATE"],
            "conditional_on_outer_node_bootstrap": True,
            "first_due_S_error_norm_upper_m_s": s_prefix["first_due_S_error_norm_upper_m_s"],
            "first_due_S_induced_attitude_correction_norm_upper_rad": s_prefix["first_due_S_induced_attitude_correction_norm_upper_rad"],
            "group_helper_limit_rad": s_prefix["deployed_group_helper_correction_limit_rad"],
        },
        "first_due_S_exact_prefix_stage": {
            "status": firstex["P5_FIRST_DUE_S_EXACT_CAYLEY_PREFIX_CERTIFICATE"],
            "required_post_cayley_norm_upper": firstex["required_first_S_post_cayley_norm_upper"],
            "widened_cayley_norm_upper": firstex["widened_prefix_cayley_norm_upper"],
            "antipodal_margin_lower": firstex["widened_prefix_antipodal_one_plus_cosine_margin_lower"],
            "vector_information_vs_goLive_metric_lower": firstex["widened_prefix_pair_information_vs_goLive_attitude_metric_lower"],
            "diagnostic_q_lt_1_is_promotion_gate": firstex["diagnostic_q_lt_1_is_promotion_gate"],
        },
        "exact_eta_geometry_support_stage": {
            "status": eta["P5_CAYLEY_ETA_GEOMETRY_CERTIFICATE"],
            "exact_eta_identity": eta["exact_eta_identity"],
            "widened_cayley_norm_upper": eta["widened_cayley_norm_upper"],
            "annular_subdivision_cell_count": eta["subdivision_cell_count"],
            "global_packet_count_times_Lipschitz_defect_used": eta["global_packet_count_times_Lipschitz_defect_used"],
            "standalone_eta_penalty_is_active_P5_route": False,
        },
        "effective_vector_input_stage": {
            "status": veff["P5_EFFECTIVE_VECTOR_INPUT_CERTIFICATE"],
            "standalone_vector_eta_penalty_retired": veff["standalone_vector_eta_penalty_retired_from_P5_numerical_route"],
            "magnetometer_radial_gain_action_exact_zero": veff["magnetometer"]["kalman_gain_radial_action_exact_zero"],
            "magnetometer_exact_state_correction_identity": veff["magnetometer"]["exact_state_correction_identity"],
            "magnetometer_effective_coordinate_nonexpansive": veff["magnetometer"]["effective_coordinate_nonexpansive"],
            "accelerometer_exact_state_correction_identity": veff["accelerometer"]["exact_state_correction_identity"],
            "accelerometer_effective_aw_norm_preserved": veff["accelerometer"]["effective_aw_defect_norm_equals_eta_norm"],
            "gravity_quotient_uses_effective_accelerometer_input": veff["gravity_quotient"]["accelerometer_effective_aw_input_descends_to_quotient"],
            "subdivision_cell_count": veff["subdivision_cell_count"],
        },
        "signed_cayley_cell_stage": {
            "status": signed["P5_SIGNED_CAYLEY_CELL_PRIMITIVE"],
            "signed_a_dot_c_retained": signed["signed_a_dot_c_retained"],
            "independent_abs_a_abs_c_denominator_used": signed["independent_abs_a_abs_c_denominator_used"],
        },
        "complete_word_transport_stage": {
            "status": transport["P5_COMPLETE_WORD_TRANSPORT_ALGEBRA_CERTIFICATE"],
            "gauged_numerical_status": transport["P5_GAUGED_COMPLETE_WORD_NUMERICAL_CERTIFICATE"],
            "quotient_numerical_status": transport["P5_GRAVITY_QUOTIENT_COMPLETE_WORD_NUMERICAL_CERTIFICATE"],
            "source_operation_order": transport["source_operation_order"],
            "full_S_to_attitude_gain_retained": transport["full_S_to_attitude_gain_retained"],
            "sequential_immediate_quaternion_resets_retained": transport["sequential_immediate_quaternion_resets_retained"],
            "standalone_vector_eta_penalty_used": transport["standalone_vector_eta_penalty_used"],
        },
        "heading_handoff_contract": {
            "P1_gravity_cosines_are_tilt_only": True,
            "gauged_quality_full_cayley_norm_upper": heading["gauged_quality_handoff"]["full_attitude_cayley_norm_upper"],
            "gauged_timeout_full_cayley_norm_upper": heading["gauged_timeout_subbranch"]["full_attitude_cayley_norm_upper"],
            "ungauged_timeout_full_heading_node_available": False,
            "ungauged_timeout_required_route": heading["ungauged_timeout_subbranch"]["required_route"],
        },
        "raw_V_R_large_angle_sector_audit": {
            "status": sector["P5_RAW_VR_LARGE_ANGLE_SECTOR"],
            "counterexample": sector["validated_counterexample"],
            "beta_cannot_repair_xi_zero_counterexample": sector["beta_cannot_repair_xi_zero_counterexample"],
        },
        "finite_angle_information_geometry": {
            "status": outinfo["P5_FINITE_ANGLE_INFORMATION_GEOMETRY_CERTIFICATE"],
            "exact_joseph_tangent_information_identity": outinfo["exact_joseph_tangent_information_identity"],
            "packet_geometry": outinfo["packet_geometry"],
            "nodes": outinfo["nodes"],
            "complete_word_sector": outinfo["P5_GAUGED_OUTER_CAYLEY_INFORMATION_WORD_SECTOR"],
        },
        "yaw_only_quotient_audit": {
            "status": oldq["P5_UNGAUGED_TIMEOUT_QUOTIENT_WORD_CERTIFICATE"],
            "obstruction_identified": oldq["P5_YAW_ONLY_QUOTIENT_OBSTRUCTION_IDENTIFIED"],
            "required_correction": oldq["required_quotient_correction"],
        },
        "detectable_gravity_quotient": {
            "status": gquot["P5_GRAVITY_QUOTIENT_REDUCED_DETECTABILITY_CERTIFICATE"],
            "complete_word_status": gquot["P5_UNGAUGED_TIMEOUT_QUOTIENT_WORD_CERTIFICATE"],
            "axial_gyro_bias_role": gquot["axial_gyro_bias_role"],
            "reduced_attitude_bias_information": gquot["reduced_attitude_bias_information"],
            "translation_word": gquot["translation_word"],
        },
        "global_normal_live_P3_covariance_bound_retired_as_outer_S_prefix_gate": True,
        "finite_angle_perturbation_vs_tiny_P3_gap_retired_as_P5_promotion_route": True,
        "raw_V_R_large_angle_sector_retired_as_P5_promotion_route": True,
        "source_shaped_Cayley_information_is_primary_full_heading_outer_route": True,
        "standalone_vector_eta_penalty_retired_as_P5_promotion_route": True,
        "yaw_only_full_bias_quotient_retired": True,
        "detectable_gravity_quotient_is_primary_ungauged_route": True,
        "gauged_full_heading_nodes": nodes,
        "ungauged_timeout_route": {
            "status": "REDUCED_DETECTABILITY_EFFECTIVE_ACCEL_INPUT_AND_EXACT_TRANSPORT_PASS_CELL_PROPAGATION_PENDING",
            "full_heading_cayley_bound_available": False,
            "yaw_only_quotient_disproved": True,
            "reduced_detectability_certificate": "PASS",
            "required_route": "DETECTABLE_GRAVITY_ONLY_QUOTIENT_WITH_AXIAL_GYRO_BIAS_NEUTRAL_BOUNDED_INPUT",
            "promotion_to_full_heading_requires": "magnetic_lock_or_regauge_hybrid_gauge_event",
            "current_numerical_obligation": quotient_failure,
        },
        "gauged_full_heading_first_failure": gauged_failure,
        "first_failure": quotient_failure,
        "next_full_heading_numerical_certificate": (
            "outward-propagate source-correlated P,H,R,K,r,d_eff cells through every later prefix of the complete 1 s gauged H word; use the exact magnetometer radial-null/effective-tangent reduction, accelerometer effective a_w input, and signed a^T c quaternion denominator, retaining full S-to-attitude and all non-attitude coordinates"
        ),
        "next_complete_startup_family_certificate": (
            "perform the analogous outward P,H,R,K,r,e_aw_eff cell propagation in the detectable gravity quotient, charging b_g_parallel only as the explicit bounded input, then certify the magnetic-gauge jump into a full-heading node"
        ),
        "P5_OUTER_H_BRIDGE_CERTIFICATE": "NOT_ESTABLISHED",
        "N_H_words": None,
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("staged bridge is not source generated")
    if d.get("source_replay_used") is not False:
        failures.append("staged bridge uses replay")
    if d.get("filter_changed") is not False:
        failures.append("staged bridge changes filter")
    for flag in (
        "global_normal_live_P3_covariance_bound_retired_as_outer_S_prefix_gate",
        "finite_angle_perturbation_vs_tiny_P3_gap_retired_as_P5_promotion_route",
        "raw_V_R_large_angle_sector_retired_as_P5_promotion_route",
        "source_shaped_Cayley_information_is_primary_full_heading_outer_route",
        "standalone_vector_eta_penalty_retired_as_P5_promotion_route",
        "yaw_only_full_bias_quotient_retired",
        "detectable_gravity_quotient_is_primary_ungauged_route",
    ):
        if d.get(flag) is not True:
            failures.append(f"{flag} is not true")

    sector = d.get("raw_V_R_large_angle_sector_audit", {})
    if sector.get("status") != "DISPROVED_ON_DECLARED_SOURCE_FAMILY":
        failures.append("raw V_R sector disproof not consumed")
    if sector.get("beta_cannot_repair_xi_zero_counterexample") is not True:
        failures.append("bridge allows beta to hide xi=0 raw V_R failure")

    fi = d.get("finite_angle_information_geometry", {})
    if fi.get("status") != "PASS":
        failures.append("finite-angle source information geometry did not pass")
    if fi.get("complete_word_sector") != "NOT_ESTABLISHED":
        failures.append("gauged complete word sector promoted prematurely")
    for name in ("normal_gauged", "timeout_gauged"):
        row = fi.get("nodes", {}).get(name, {})
        if row.get("strict") is not True:
            failures.append(f"{name}: finite-angle information prerequisite not strict")

    sx = d.get("first_due_S_exact_prefix_stage", {})
    if sx.get("status") != "PASS_WIDENED_CHART":
        failures.append("first-S exact widened chart was not consumed")
    if sx.get("diagnostic_q_lt_1_is_promotion_gate") is not False:
        failures.append("bridge reinstated diagnostic q<1 gate")
    if not float(sx.get("antipodal_margin_lower", 0.0)) > 0.0:
        failures.append("first-S widened chart has no antipodal margin")

    eg = d.get("exact_eta_geometry_support_stage", {})
    if eg.get("status") != "PASS":
        failures.append("exact Cayley eta support geometry did not pass")
    if eg.get("global_packet_count_times_Lipschitz_defect_used") is not False:
        failures.append("bridge uses retired packet-count Lipschitz eta penalty")
    if eg.get("standalone_eta_penalty_is_active_P5_route") is not False:
        failures.append("eta support geometry remains an active standalone penalty route")
    if int(eg.get("annular_subdivision_cell_count", 0)) < 1:
        failures.append("bridge has no Cayley geometry subdivision")

    ev = d.get("effective_vector_input_stage", {})
    if ev.get("status") != "PASS":
        failures.append("effective vector-input reduction did not pass")
    if ev.get("standalone_vector_eta_penalty_retired") is not True:
        failures.append("standalone vector eta penalty was not retired")
    if ev.get("magnetometer_radial_gain_action_exact_zero") is not True:
        failures.append("magnetometer radial finite-angle residual is not gain-null")
    if ev.get("magnetometer_effective_coordinate_nonexpansive") is not True:
        failures.append("magnetometer effective tangent coordinate is not nonexpansive")
    if ev.get("accelerometer_effective_aw_norm_preserved") is not True:
        failures.append("accelerometer effective aw-input norm identity missing")
    if ev.get("gravity_quotient_uses_effective_accelerometer_input") is not True:
        failures.append("gravity quotient does not consume effective accelerometer input")

    sc = d.get("signed_cayley_cell_stage", {})
    if sc.get("status") != "PASS":
        failures.append("signed Cayley correction primitive did not pass")
    if sc.get("signed_a_dot_c_retained") is not True or sc.get("independent_abs_a_abs_c_denominator_used") is not False:
        failures.append("bridge loses signed Cayley correction/source correlation")

    tr = d.get("complete_word_transport_stage", {})
    if tr.get("status") != "PASS":
        failures.append("complete-word exact transport algebra did not pass")
    if tr.get("full_S_to_attitude_gain_retained") is not True:
        failures.append("complete-word transport drops S-to-attitude gain")
    if tr.get("sequential_immediate_quaternion_resets_retained") is not True:
        failures.append("complete-word transport drops sequential resets")
    if tr.get("standalone_vector_eta_penalty_used") is not False:
        failures.append("complete-word transport reintroduced standalone vector eta")
    if tr.get("gauged_numerical_status") != "NOT_ESTABLISHED" or tr.get("quotient_numerical_status") != "NOT_ESTABLISHED":
        failures.append("bridge consumed a numerical word promotion that is not established")

    oldq = d.get("yaw_only_quotient_audit", {})
    if oldq.get("obstruction_identified") != "PASS" or oldq.get("status") != "NOT_ESTABLISHED":
        failures.append("yaw-only quotient obstruction not consumed")
    dq = d.get("detectable_gravity_quotient", {})
    if dq.get("status") != "PASS":
        failures.append("corrected detectable gravity quotient prerequisite did not pass")
    if dq.get("complete_word_status") != "NOT_ESTABLISHED":
        failures.append("complete gravity quotient word promoted prematurely")
    if "NEUTRAL_BOUNDED" not in str(dq.get("axial_gyro_bias_role", "")):
        failures.append("axial gyro bias is not carried as neutral bounded input")

    g = d.get("goLive_covariance_stage", {})
    if g.get("status") != "PASS" or g.get("S_to_attitude_gain_exact_zero") is not True:
        failures.append("goLive covariance seed is not closed")
    fs = d.get("first_due_S_gain_stage", {})
    if fs.get("status") != "PASS":
        failures.append("first-due S gain stage did not pass")
    sp = d.get("first_due_S_state_prefix_stage", {})
    if sp.get("status") != "PASS_CONDITIONAL":
        failures.append("first-S state prefix stage did not pass conditionally")
    if not (float(sp.get("first_due_S_induced_attitude_correction_norm_upper_rad", math.inf))
            < float(sp.get("group_helper_limit_rad", -math.inf))):
        failures.append("staged first-S correction exceeds exact helper")

    h = d.get("heading_handoff_contract", {})
    if h.get("P1_gravity_cosines_are_tilt_only") is not True:
        failures.append("bridge treats gravity cosine as full attitude")
    if h.get("ungauged_timeout_full_heading_node_available") is not False:
        failures.append("bridge invented full-heading ungauged timeout node")

    nodes = d.get("gauged_full_heading_nodes", {})
    if set(nodes) != {"normal_gauged", "timeout_gauged"}:
        failures.append("gauged normal/timeout nodes missing")
    for name, n in nodes.items():
        if n.get("inside_antipodal_exclusion") is not True:
            failures.append(f"{name}: handoff reaches antipodal set")
        if n.get("S_induced_correction_inside_group_helper") is not True:
            failures.append(f"{name}: first-S correction outside group helper")
        if n.get("first_due_S_exact_prefix_certificate") != "PASS_WIDENED_CHART" or n.get("inside_widened_first_S_chart") is not True:
            failures.append(f"{name}: exact first-S prefix not closed")
        if n.get("raw_V_R_sector_is_P5_promotion_route") is not False:
            failures.append(f"{name}: disproved raw V_R route still promotes")
        if n.get("finite_angle_information_geometry_status") != "PASS":
            failures.append(f"{name}: finite-angle information prerequisite missing")
        if not (float(n.get("exact_pair_residual_information_vs_goLive_attitude_metric_lower", 0.0)) > 0.0):
            failures.append(f"{name}: finite-angle information coefficient is not positive")
        if n.get("source_correlated_Joseph_information_identity_retained") is not True:
            failures.append(f"{name}: source-correlated Joseph identity not retained")

    u = d.get("ungauged_timeout_route", {})
    if u.get("reduced_detectability_certificate") != "PASS":
        failures.append("ungauged reduced detectability prerequisite missing")
    if "AXIAL_GYRO_BIAS" not in str(u.get("required_route", "")):
        failures.append("ungauged route does not name axial gyro-bias neutral direction")

    if d.get("gauged_full_heading_first_failure") != "COMPLETE_WORD_EFFECTIVE_VECTOR_INPUT_RESET_PREFIX_BUDGET_NOT_CERTIFIED":
        failures.append("wrong gauged numerical obstruction")
    if d.get("first_failure") != "GRAVITY_QUOTIENT_EFFECTIVE_ACCEL_INPUT_RESET_PREFIX_BUDGET_NOT_CERTIFIED":
        failures.append("wrong complete-family numerical obstruction")
    if d.get("P5_OUTER_H_BRIDGE_CERTIFICATE") != "NOT_ESTABLISHED":
        failures.append("P5 bridge promoted before complete nonlinear words close")
    if d.get("N_H_words") is not None:
        failures.append("finite H word count set before P5 bridge closes")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    failures = validate(out)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out.get("P5_OUTER_H_BRIDGE_CERTIFICATE"),
        "first_failure": out.get("first_failure"),
        "gauged_first_failure": out.get("gauged_full_heading_first_failure"),
        "first_S_exact": out.get("first_due_S_exact_prefix_stage"),
        "effective_vector_input": out.get("effective_vector_input_stage"),
        "transport": out.get("complete_word_transport_stage"),
        "signed_cayley": out.get("signed_cayley_cell_stage"),
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
