#!/usr/bin/env python3
"""Compose the staged OU-III P5 outer-H capture bridge.

The P5 bridge now consumes four source-faithful pieces rather than the retired
normal-Live global covariance/radius surrogates:

* ``ou3_p5_go_live_covariance_stage``: P_theta,S=0 and K_theta,S=0 exactly at
  goLive, with the pseudo phase reset;
* ``ou3_p5_first_s_gain_certificate``: a small source-staged first-due
  K_theta,S bound from the persistent constructor covariance;
* ``ou3_p5_first_s_state_prefix_certificate``: a deterministic first-S state
  bound, conditional only on the declared outer-node bootstrap, retaining every
  source-possible accepted physical correction;
* ``ou3_p5_heading_handoff_contract``: P1 gravity cosines are tilt bounds, not
  full SO(3) bounds.  Gauged quality/timeout subbranches receive composed
  full-attitude Cayley nodes, while the ungauged timeout branch is routed to the
  gravity/yaw quotient until a magnetic gauge exists.

For the gauged full-heading branches the early S->attitude obstruction is now
closed conditionally: the staged K_theta,S coefficient multiplied by the
source-staged S prefix remains inside the validated finite-correction helper.
The next full-heading obligation is therefore the paper's exact large-angle
source correction sector and the associated outer-node prefix invariance.  The
old comparison of O(1) vector geometry against the translation-limited ~1e-35
P3 gap is diagnostic only.

The complete startup family also contains an ungauged timeout branch.  This file
does not silently assign that branch a full-heading Cayley radius: it remains a
separate yaw-quotient capture obligation.  P5 is not promoted and N_H is not set
until both routes are composed correctly.
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
import ou3_p5_outer_h_word_certificate as OUTER
import ou3_startup_stability_certificate as P1

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 3


def _gauged_node(name: str, heading_row: dict, outer_name: str,
                 outer: dict, p4: dict, first_s: dict, s_prefix: dict) -> dict:
    o = outer["node_word_tests"][outer_name]
    cfull = float(heading_row["full_attitude_cosine_lower"])
    qfull = float(heading_row["full_attitude_cayley_norm_upper"])
    V_R_upper = math.nextafter(1.0 - cfull, math.inf)
    if not (0.0 <= V_R_upper < 2.0 and 0.0 < qfull < 1.0):
        raise RuntimeError(f"{name}: invalid gauged full-attitude handoff node")

    perturb = float(o["outer_vector_nonlinear_information_ratio_upper"])
    p3sqrt = float(o["P3_homogeneous_sqrt_decrease_lower"])
    if not (math.isfinite(perturb) and perturb >= 0.0 and p3sqrt > 0.0):
        raise RuntimeError(f"{name}: invalid retired outer diagnostic ratios")

    kS = float(first_s["K_thetaS_operator_norm_upper_first_due"])
    Sbound = float(s_prefix["first_due_S_error_norm_upper_m_s"])
    Sinjection = float(s_prefix["first_due_S_induced_attitude_correction_norm_upper_rad"])
    helper = float(s_prefix["deployed_group_helper_correction_limit_rad"])
    if not (Sinjection < helper):
        raise RuntimeError(f"{name}: staged first-S correction left validated group helper")

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
        "finite_angle_perturbation_diagnostic_ratio_upper": perturb,
        "P3_translation_limited_sqrt_gap_lower": p3sqrt,
        "perturbation_over_P3_sqrt_gap": perturb / p3sqrt,
        "perturbation_vs_P3_gap_is_P5_promotion_route": False,
        "exact_large_angle_sector_required": True,
        "large_angle_sector_form": (
            "D_R,g = V_R(R_e)-V_R(exp([d_g]_x)R_e) "
            ">= alpha_R,i V_R(R_e)-beta_R,i ||xi||^2"
        ),
        "large_angle_sector_status": "NOT_ESTABLISHED",
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

    prereq = []
    prereq += [f"P1: {x}" for x in P1.validate(p1)]
    prereq += [f"P4: {x}" for x in P4.validate(p4)]
    prereq += [f"outer-diagnostic: {x}" for x in OUTER.validate(outer)]
    prereq += [f"goLive-stage: {x}" for x in GOLIVE.validate(stage)]
    prereq += [f"first-S-gain: {x}" for x in FIRSTS.validate(first_s)]
    prereq += [f"first-S-state-prefix: {x}" for x in SPREFIX.validate(s_prefix)]
    prereq += [f"heading-handoff: {x}" for x in HEADING.validate(heading)]
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
            outer, p4, first_s, s_prefix,
        ),
        "timeout_gauged": _gauged_node(
            "timeout_gauged", heading["gauged_timeout_subbranch"], "timeout",
            outer, p4, first_s, s_prefix,
        ),
    }
    quotient = heading["ungauged_timeout_subbranch"]

    # The conditional first-S state/gain pair is closed.  What remains on each
    # gauged node is the exact finite-angle dissipation/prefix bootstrap.  The
    # complete startup family has a second route: an ungauged timeout cannot be
    # promoted into the full-heading node at all and must be captured on the
    # yaw quotient until magnetic gauge acquisition.
    gauged_first_failure = "EXACT_LARGE_ANGLE_VECTOR_DISSIPATION_SECTOR_NOT_CERTIFIED"
    complete_first_failure = "UNGAUGED_TIMEOUT_YAW_QUOTIENT_CAPTURE_NOT_CERTIFIED"

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_STAGED_OUTER_H_BRIDGE",
        "claim": "SOURCE_STAGED_OUTER_H_CAPTURE_WITH_GAUGED_AND_YAW_QUOTIENT_BRANCHES",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "local_P4_recurrence_used_as_outer_bridge": False,
        "goLive_covariance_stage": {
            "status": stage["P5_GOLIVE_COVARIANCE_STAGE_CERTIFICATE"],
            "theta_S_cross_covariance_operator_norm_upper": 0.0,
            "S_to_attitude_gain_exact_zero": True,
            "pseudo_elapsed_s": seed["pseudo_update_elapsed_s_at_goLive"],
            "P_SS_variance_per_axis": seed["P_SS_variance_per_axis"],
            "P_awaw_source_std_outward_mps2": seed["P_awaw_source_std_outward_mps2"],
        },
        "first_due_S_gain_stage": {
            "status": first_s["P5_FIRST_DUE_S_GAIN_CERTIFICATE"],
            "theta_S_canonical_correlation_upper": first_s["theta_S_canonical_correlation_upper"],
            "K_thetaS_operator_norm_upper": first_s["K_thetaS_operator_norm_upper_first_due"],
            "gain_widening_factor_vs_global_P3_bound_lower": first_s["gain_widening_factor_vs_global_P3_bound_lower"],
            "first_due_time_upper_s": first_s["timing"]["first_due_time_upper_s"],
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
            "ungauged_timeout_required_route": quotient["required_route"],
        },
        "global_normal_live_P3_covariance_bound_retired_as_outer_S_prefix_gate": True,
        "old_global_P3_S_induced_attitude_bound_is_diagnostic_only": True,
        "finite_angle_perturbation_vs_tiny_P3_gap_retired_as_P5_promotion_route": True,
        "exact_large_angle_vector_dissipation_is_primary_full_heading_outer_route": True,
        "gauged_full_heading_nodes": nodes,
        "ungauged_timeout_route": {
            "status": "NOT_ESTABLISHED",
            "full_heading_cayley_bound_available": False,
            "required_route": quotient["required_route"],
            "promotion_to_full_heading_requires": "magnetic_lock_or_regauge_hybrid_gauge_event",
        },
        "gauged_full_heading_first_failure": gauged_first_failure,
        "first_failure": complete_first_failure,
        "next_full_heading_numerical_certificate": (
            "validated source-correlated exact large-angle vector dissipation sector and prefix-invariance "
            "bootstrap on the normal-gauged and timeout-gauged H nodes"
        ),
        "next_complete_startup_family_certificate": (
            "gravity-only/yaw-quotient H capture for the ungauged timeout branch through the magnetic gauge hybrid event"
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
        failures.append("P5 still compares outer physical dissipation to tiny P3 perturbation gap")
    if d.get("exact_large_angle_vector_dissipation_is_primary_full_heading_outer_route") is not True:
        failures.append("exact large-angle sector is not primary full-heading outer route")

    g = d.get("goLive_covariance_stage", {})
    if g.get("status") != "PASS" or g.get("S_to_attitude_gain_exact_zero") is not True:
        failures.append("goLive covariance seed is not closed")
    fs = d.get("first_due_S_gain_stage", {})
    if fs.get("status") != "PASS":
        failures.append("first-due S gain stage did not pass")
    sp = d.get("first_due_S_state_prefix_stage", {})
    if sp.get("status") != "PASS_CONDITIONAL" or sp.get("conditional_on_outer_node_bootstrap") is not True:
        failures.append("conditional first-S state prefix not consumed correctly")
    if not (float(sp.get("first_due_S_induced_attitude_correction_norm_upper_rad", math.inf))
            < float(sp.get("group_helper_limit_rad", -math.inf))):
        failures.append("staged first-S correction exceeds helper")

    h = d.get("heading_handoff_contract", {})
    if h.get("P1_gravity_cosines_are_tilt_only") is not True:
        failures.append("bridge still treats gravity cosine as full attitude")
    if h.get("ungauged_timeout_full_heading_node_available") is not False:
        failures.append("bridge invented a full-heading node for ungauged timeout")
    if "YAW_QUOTIENT" not in str(h.get("ungauged_timeout_required_route", "")):
        failures.append("bridge did not route ungauged timeout to yaw quotient")

    nodes = d.get("gauged_full_heading_nodes", {})
    if set(nodes) != {"normal_gauged", "timeout_gauged"}:
        failures.append("gauged normal/timeout nodes missing")
    for name, n in nodes.items():
        if n.get("P1_gravity_tilt_cosine_not_used_as_full_attitude_cosine") is not True:
            failures.append(f"{name}: gravity cosine reused as full attitude")
        if n.get("inside_antipodal_exclusion") is not True:
            failures.append(f"{name}: gauged handoff reaches antipodal set")
        if n.get("inside_candidate_outer_cayley_bootstrap") is not True:
            failures.append(f"{name}: gauged handoff leaves candidate q<=1 node")
        if n.get("first_due_S_gain_certificate") != "PASS":
            failures.append(f"{name}: first-due S gain missing")
        if n.get("first_due_S_state_prefix_certificate") != "PASS_CONDITIONAL":
            failures.append(f"{name}: first-S state prefix missing")
        if n.get("S_induced_correction_inside_group_helper") is not True:
            failures.append(f"{name}: S-induced correction outside helper")
        if n.get("large_angle_sector_status") != "NOT_ESTABLISHED":
            failures.append(f"{name}: large-angle sector promoted without proof")
        if n.get("outer_prefix_domain_bootstrap_status") != "NOT_ESTABLISHED":
            failures.append(f"{name}: outer prefix bootstrap promoted without proof")
        if n.get("perturbation_vs_P3_gap_is_P5_promotion_route") is not False:
            failures.append(f"{name}: retired P3 perturbation route active")

    if d.get("gauged_full_heading_first_failure") != "EXACT_LARGE_ANGLE_VECTOR_DISSIPATION_SECTOR_NOT_CERTIFIED":
        failures.append("wrong next gauged full-heading obstruction")
    if d.get("first_failure") != "UNGAUGED_TIMEOUT_YAW_QUOTIENT_CAPTURE_NOT_CERTIFIED":
        failures.append("wrong complete-startup-family obstruction")
    if d.get("P5_OUTER_H_BRIDGE_CERTIFICATE") != "NOT_ESTABLISHED":
        failures.append("P5 bridge promoted before both routes close")
    if d.get("N_H_words") is not None:
        failures.append("finite H word count set before outer bridge closes")
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
        "first_due_S_state_prefix": out.get("first_due_S_state_prefix_stage"),
        "heading": out.get("heading_handoff_contract"),
        "nodes": out.get("gauged_full_heading_nodes"),
        "ungauged_timeout": out.get("ungauged_timeout_route"),
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
