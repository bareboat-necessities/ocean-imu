#!/usr/bin/env python3
"""Compose the staged OU-III P5 outer-H capture bridge.

The earlier outer-H diagnostic intentionally started from the source-uniform
normal-Live P3 covariance envelope.  CI showed why that is not the right object
at startup: it produced an 8.5e7-rad S->attitude prefix bound.  The shipping
handoff has much more structure.  :mod:`ou3_p5_go_live_covariance_stage`
proves from source that P_theta,S is exactly zero at goLive, the a_w cross
covariances are reset to zero, and the pseudo timer starts from zero.

This composer therefore retires the global-P3 S-prefix number as a P5 gate and
replaces it with the actual finite source-stage obligations:

1. propagate the exact goLive covariance seed to the first due S pseudo through
   every admissible accepted/rejected physical-measurement prefix;
2. certify the exact large-angle correction sector from the paper on the same
   source-correlated covariance/gain tuples;
3. prove the complete outer prefix remains in the chosen finite Cayley node and
   overlaps the already-certified P4 inner level.

Nothing in this file promotes P5 before those obligations are numerical.  Its
purpose is to make the implementation theorem chain consume the correct next
objects rather than the deliberately crude diagnostic bound.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_nonlinear_word_certificate as P4
import ou3_p5_go_live_covariance_stage as GOLIVE
import ou3_p5_outer_h_word_certificate as OUTER
import ou3_startup_stability_certificate as P1

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def _node(name: str, p1: dict, outer: dict, p4: dict) -> dict:
    o = outer["node_word_tests"][name]
    handoff = outer["handoff_nodes"][name]
    cos_lower = (
        float(p1["normal_handoff"]["true_gravity_cosine_lower"])
        if name == "normal"
        else float(p1["timeout_handoff"]["combined_true_gravity_cosine_lower"])
    )
    V_R_upper = math.nextafter(1.0 - cos_lower, math.inf)
    if not (0.0 <= V_R_upper < 2.0):
        raise RuntimeError(f"{name}: invalid finite-angle handoff energy")
    perturb = float(o["outer_vector_nonlinear_information_ratio_upper"])
    p3sqrt = float(o["P3_homogeneous_sqrt_decrease_lower"])
    if not (math.isfinite(perturb) and perturb >= 0.0 and p3sqrt > 0.0):
        raise RuntimeError(f"{name}: invalid outer diagnostic ratios")
    return {
        "node": name,
        "handoff_cayley_norm_upper": float(handoff["rotation"]["cayley_norm_upper"]),
        "handoff_group_energy_V_R_upper": V_R_upper,
        "inside_antipodal_exclusion": cos_lower > -1.0,
        "goLive_S_to_attitude_gain_exact_zero": True,
        "global_P3_S_prefix_bound_used_for_P5_gate": False,
        "pre_first_due_S_cross_covariance_status": "NOT_ESTABLISHED",
        "finite_angle_perturbation_diagnostic_ratio_upper": perturb,
        "P3_translation_limited_sqrt_gap_lower": p3sqrt,
        "perturbation_over_P3_sqrt_gap": perturb / p3sqrt,
        "perturbation_vs_P3_gap_is_P5_promotion_route": False,
        "exact_large_angle_sector_required": True,
        "large_angle_sector_form": (
            "D_R,g = V_R(R_e)-V_R(exp([d_g]_x)R_e) "
            ">= alpha_R,i V_R(R_e)-beta_R,i ||xi||^2"
        ),
        "large_angle_sector_status": "PENDING_SOURCE_CORRELATED_GAIN_TUPLES",
        "P4_inner_overlap_target_W": float(p4["modes"]["H"]["certified_level_W"]),
        "outer_prefix_domain_bootstrap_status": "NOT_ESTABLISHED",
        "pass": False,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    p1 = P1.build(domain_path)
    p4 = P4.build(domain_path)
    outer = OUTER.build(domain_path)
    stage = GOLIVE.build(domain_path)

    prereq = []
    prereq += [f"P1: {x}" for x in P1.validate(p1)]
    prereq += [f"P4: {x}" for x in P4.validate(p4)]
    prereq += [f"outer: {x}" for x in OUTER.validate(outer)]
    prereq += [f"goLive-stage: {x}" for x in GOLIVE.validate(stage)]
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
    first_stage = stage["pre_first_S_stage"]
    if seed["theta_S_cross_covariance_operator_norm_upper"] != 0.0:
        raise RuntimeError("goLive theta-S covariance is not exact zero")
    if seed["S_to_attitude_gain_at_goLive_exact_zero"] is not True:
        raise RuntimeError("goLive S-to-attitude gain is not exact zero")

    nodes = {name: _node(name, p1, outer, p4) for name in ("normal", "timeout")}
    # Both the old finite-angle perturbation comparison and the global P3
    # S-prefix bound are now diagnostics only.  The first *actual* unclosed P5
    # obligation is the finite covariance propagation to the first due pseudo.
    first_failure = "PRE_FIRST_DUE_THETA_S_CROSS_COVARIANCE_NOT_ENCLOSED"

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_STAGED_OUTER_H_BRIDGE",
        "claim": "SOURCE_STAGED_OUTER_H_CAPTURE_BRIDGE_WITH_EXACT_LARGE_ANGLE_ROUTE",
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
        "pre_first_due_stage": first_stage,
        "global_normal_live_P3_covariance_bound_retired_as_outer_S_prefix_gate": True,
        "old_global_P3_S_induced_attitude_bound_is_diagnostic_only": True,
        "finite_angle_perturbation_vs_tiny_P3_gap_retired_as_P5_promotion_route": True,
        "exact_large_angle_vector_dissipation_is_primary_outer_route": True,
        "nodes": nodes,
        "first_failure": first_failure,
        "first_required_numerical_certificate": (
            "source-correlated covariance/gain propagation from the exact goLive H seed through the finite "
            "pre-first-S stage, including accepted accelerometer branches, followed by direct K_thetaS enclosure"
        ),
        "second_required_numerical_certificate": (
            "exact source-correlated large-angle vector dissipation sector on normal and timeout H nodes"
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
    if d.get("exact_large_angle_vector_dissipation_is_primary_outer_route") is not True:
        failures.append("exact large-angle sector is not primary outer route")
    g = d.get("goLive_covariance_stage", {})
    if g.get("status") != "PASS" or g.get("S_to_attitude_gain_exact_zero") is not True:
        failures.append("goLive covariance seed is not closed")
    if g.get("theta_S_cross_covariance_operator_norm_upper") != 0.0:
        failures.append("goLive theta-S covariance lost exact zero")
    if d.get("first_failure") != "PRE_FIRST_DUE_THETA_S_CROSS_COVARIANCE_NOT_ENCLOSED":
        failures.append("wrong next P5 obstruction")
    if d.get("P5_OUTER_H_BRIDGE_CERTIFICATE") != "NOT_ESTABLISHED":
        failures.append("P5 bridge promoted before first-due/large-angle closure")
    if d.get("N_H_words") is not None:
        failures.append("finite H word count set before outer bridge closes")
    nodes = d.get("nodes", {})
    if set(nodes) != {"normal", "timeout"}:
        failures.append("normal/timeout outer nodes missing")
    for name, n in nodes.items():
        if n.get("inside_antipodal_exclusion") is not True:
            failures.append(f"{name}: handoff reaches antipodal set")
        if n.get("goLive_S_to_attitude_gain_exact_zero") is not True:
            failures.append(f"{name}: exact goLive S gain not consumed")
        if n.get("exact_large_angle_sector_required") is not True:
            failures.append(f"{name}: large-angle sector not required")
        ratio = n.get("perturbation_over_P3_sqrt_gap")
        if not (isinstance(ratio, (int, float)) and math.isfinite(float(ratio)) and float(ratio) > 1.0):
            failures.append(f"{name}: old perturbation obstruction not quantified")
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
        "goLive": out.get("goLive_covariance_stage"),
        "nodes": out.get("nodes"),
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
