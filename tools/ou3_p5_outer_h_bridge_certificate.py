#!/usr/bin/env python3
"""Compose the staged OU-III P5 outer-H capture bridge.

This composer consumes the source-staged goLive/first-S work, the corrected
heading handoff split, and two theorem audits. The requested raw trace-energy
large-angle sector is evaluated on an actual source-correlated goLive
covariance/gain tuple with validated interval Kalman arithmetic and the exact
deployed quaternion/Rodrigues backend. The ungauged timeout yaw-only quotient
is evaluated on an exact gravity-only axial-gyro-bias zero-dynamics word.

Both audits are fail-closed. If the raw V_R sector is disproved at xi=0, beta
cannot repair it and the bridge must stay on the source-shaped Cayley/information
metric. If the yaw-only quotient retains an unobservable axial gyro-bias
coordinate, it cannot be assigned lambda<1 merely by quotienting attitude yaw.
No filter behavior is changed and no older radius/block certificate is revived.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_nonlinear_word_certificate as P4
import ou3_p5_first_s_gain_certificate as FIRSTS
import ou3_p5_first_s_state_prefix_certificate as SPREFIX
import ou3_p5_go_live_covariance_stage as GOLIVE
import ou3_p5_heading_handoff_contract as HEADING
import ou3_p5_large_angle_sector_certificate as SECTOR
import ou3_p5_outer_h_word_certificate as OUTER
import ou3_p5_yaw_quotient_word_certificate as QUOTIENT
import ou3_startup_stability_certificate as P1

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 4


def _gauged_node(name: str, heading_row: dict, outer_name: str,
                 outer: dict, p4: dict, first_s: dict, s_prefix: dict,
                 sector: dict) -> dict:
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
        raise RuntimeError("raw V_R sector audit did not produce the validated source counterexample")

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
        "retired_global_outer_diagnostic_ratio_upper": float(o["outer_vector_nonlinear_information_ratio_upper"]),
        "retired_P3_translation_limited_sqrt_gap_lower": float(o["P3_homogeneous_sqrt_decrease_lower"]),
        "raw_V_R_large_angle_sector_status": raw_status,
        "raw_V_R_sector_witness_inside_node": float(sector["validated_counterexample"]["cayley_norm_interval"][1]) < qfull,
        "raw_V_R_sector_is_P5_promotion_route": False,
        "source_shaped_Cayley_information_outer_sector_required": True,
        "source_shaped_Cayley_information_outer_sector_status": "NOT_ESTABLISHED",
        "exact_group_backend_retained": True,
        "outer_S_state_prefix_status": "PASS_CONDITIONAL_ON_OUTER_NODE",
        "outer_prefix_domain_bootstrap_status": "NOT_ESTABLISHED",
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
    heading = HEADING.build(domain_path)
    sector = SECTOR.build(domain_path)
    quotient = QUOTIENT.build(domain_path)

    prereq = []
    prereq += [f"P1: {x}" for x in P1.validate(p1)]
    prereq += [f"P4: {x}" for x in P4.validate(p4)]
    prereq += [f"outer-diagnostic: {x}" for x in OUTER.validate(outer)]
    prereq += [f"goLive-stage: {x}" for x in GOLIVE.validate(stage)]
    prereq += [f"first-S-gain: {x}" for x in FIRSTS.validate(first_s)]
    prereq += [f"first-S-state-prefix: {x}" for x in SPREFIX.validate(s_prefix)]
    prereq += [f"heading-handoff: {x}" for x in HEADING.validate(heading)]
    prereq += [f"large-angle-sector-audit: {x}" for x in SECTOR.validate(sector)]
    prereq += [f"yaw-quotient-audit: {x}" for x in QUOTIENT.validate(quotient)]
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
            outer, p4, first_s, s_prefix, sector,
        ),
        "timeout_gauged": _gauged_node(
            "timeout_gauged", heading["gauged_timeout_subbranch"], "timeout",
            outer, p4, first_s, s_prefix, sector,
        ),
    }

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_STAGED_OUTER_H_BRIDGE",
        "claim": "SOURCE_STAGED_OUTER_H_CAPTURE_AFTER_RAW_VR_AND_YAW_ONLY_QUOTIENT_AUDITS",
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
            "replacement": sector["required_theorem_correction"],
        },
        "yaw_only_quotient_audit": {
            "status": quotient["P5_UNGAUGED_TIMEOUT_QUOTIENT_WORD_CERTIFICATE"],
            "obstruction_identified": quotient["P5_YAW_ONLY_QUOTIENT_OBSTRUCTION_IDENTIFIED"],
            "witness": quotient["witness"],
            "required_correction": quotient["required_quotient_correction"],
        },
        "global_normal_live_P3_covariance_bound_retired_as_outer_S_prefix_gate": True,
        "finite_angle_perturbation_vs_tiny_P3_gap_retired_as_P5_promotion_route": True,
        "raw_V_R_large_angle_sector_retired_as_P5_promotion_route": True,
        "source_shaped_Cayley_information_is_primary_full_heading_outer_route": True,
        "gauged_full_heading_nodes": nodes,
        "ungauged_timeout_route": {
            "status": "NOT_ESTABLISHED",
            "full_heading_cayley_bound_available": False,
            "yaw_only_quotient_disproved": True,
            "required_route": "OBSERVABLE_GRAVITY_ONLY_QUOTIENT_WITH_AXIAL_GYRO_BIAS_NEUTRAL_OR_BOUNDED_INPUT",
            "promotion_to_full_heading_requires": "magnetic_lock_or_regauge_hybrid_gauge_event",
        },
        "gauged_full_heading_first_failure": "SOURCE_SHAPED_CAYLEY_INFORMATION_OUTER_SECTOR_NOT_CERTIFIED",
        "first_failure": "OBSERVABLE_GRAVITY_ONLY_QUOTIENT_WORD_NOT_CERTIFIED",
        "next_full_heading_numerical_certificate": (
            "validated finite-angle source-correlated outer sector in the same Cayley/information path metric used by P3/P4, with exact quaternion/Rodrigues prefixes on the q~0.27 and q~0.60 gauged nodes"
        ),
        "next_complete_startup_family_certificate": (
            "observable/detectable gravity-only quotient word that removes yaw and does not demand strict contraction of the axial gyro-bias zero dynamics, followed by the magnetic-gauge jump"
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
    if d.get("local_P4_recurrence_used_as_outer_bridge") is not False:
        failures.append("staged bridge reuses local P4 recurrence")
    if d.get("global_normal_live_P3_covariance_bound_retired_as_outer_S_prefix_gate") is not True:
        failures.append("global P3 covariance still gates startup S prefix")
    if d.get("finite_angle_perturbation_vs_tiny_P3_gap_retired_as_P5_promotion_route") is not True:
        failures.append("P5 still compares outer physical dissipation to tiny P3 gap")
    if d.get("raw_V_R_large_angle_sector_retired_as_P5_promotion_route") is not True:
        failures.append("disproved raw V_R sector still promotes P5")
    if d.get("source_shaped_Cayley_information_is_primary_full_heading_outer_route") is not True:
        failures.append("source-shaped Cayley/information metric is not primary outer route")

    sector = d.get("raw_V_R_large_angle_sector_audit", {})
    if sector.get("status") != "DISPROVED_ON_DECLARED_SOURCE_FAMILY":
        failures.append("raw V_R sector disproof not consumed")
    if sector.get("beta_cannot_repair_xi_zero_counterexample") is not True:
        failures.append("bridge allows beta to hide xi=0 raw V_R failure")

    q = d.get("yaw_only_quotient_audit", {})
    if q.get("obstruction_identified") != "PASS":
        failures.append("yaw-only quotient zero dynamics not consumed")
    if q.get("status") != "NOT_ESTABLISHED":
        failures.append("yaw-only quotient promoted despite zero dynamics")

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
        if n.get("raw_V_R_large_angle_sector_status") != "DISPROVED_ON_DECLARED_SOURCE_FAMILY":
            failures.append(f"{name}: raw V_R disproof not attached")
        if n.get("raw_V_R_sector_witness_inside_node") is not True:
            failures.append(f"{name}: validated raw V_R witness not inside node")
        if n.get("raw_V_R_sector_is_P5_promotion_route") is not False:
            failures.append(f"{name}: disproved raw V_R route still promotes")
        if n.get("source_shaped_Cayley_information_outer_sector_required") is not True:
            failures.append(f"{name}: corrected metric sector not required")

    u = d.get("ungauged_timeout_route", {})
    if u.get("yaw_only_quotient_disproved") is not True:
        failures.append("ungauged route did not consume yaw-only quotient disproof")
    if "AXIAL_GYRO_BIAS" not in str(u.get("required_route", "")):
        failures.append("ungauged route does not name axial gyro-bias neutral direction")

    if d.get("gauged_full_heading_first_failure") != "SOURCE_SHAPED_CAYLEY_INFORMATION_OUTER_SECTOR_NOT_CERTIFIED":
        failures.append("wrong gauged next obstruction")
    if d.get("first_failure") != "OBSERVABLE_GRAVITY_ONLY_QUOTIENT_WORD_NOT_CERTIFIED":
        failures.append("wrong complete-family next obstruction")
    if d.get("P5_OUTER_H_BRIDGE_CERTIFICATE") != "NOT_ESTABLISHED":
        failures.append("P5 bridge promoted before corrected outer routes close")
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
        "raw_V_R_sector": out.get("raw_V_R_large_angle_sector_audit"),
        "yaw_only_quotient": out.get("yaw_only_quotient_audit"),
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
