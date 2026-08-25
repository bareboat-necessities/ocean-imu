#!/usr/bin/env python3
"""Exact fixed-mode nonlinear source-word map contract for OU-III P4.

This module does not simulate selected trajectories.  It defines the operations
that every source-complete normal-Live H/A word must compose, in the shipping
order extracted by :mod:`ou3_implementation_proof_manifest`.

A correction is represented in the MEKF's actual form:

    r = y - yhat,
    dx = K r,
    x <- x + dx,
    R_e <- R_inj(dx_theta) R_e,
    delta_theta <- 0,
    P <- G_reset P G_reset',

where ``R_inj`` is the exact deployed normalized quaternion map supplied by
:mod:`ou3_p4_group_algebra`.  The full K is retained, including S->attitude and
all attitude/non-attitude cross gains.  Rejected/not-due branches are identity
corrections.  In A mode the source's accelerometer-bias ball projection follows
every state injection and is retained as an exact convex projection operation.

The numerical P4 enclosure backend consumes this semantic map together with
validated covariance/gain/source-node boxes.  Keeping semantics in a separate
source-bound producer prevents a later interval optimization from changing the
filter it is proving.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_implementation_proof_manifest as MANIFEST
import ou3_implementation_word_language as WORDS

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1

EXPECTED_ORDER = [
    "commit_previous_tune",
    "prediction",
    "apply_pending_aw_covariance_psd_increment",
    "periodic_S_zero_when_due_then_immediate_quaternion_injection_and_left_error_reset",
    "accelerometer_correction_or_rejection_then_immediate_quaternion_injection_and_left_error_reset_if_accepted",
    "source_tuner_evolution_and_stage_next_tune",
    "periodic_aw_covariance_sync_tick_stages_future_psd_increment",
    "asynchronous_magnetometer_correction_or_rejection_then_immediate_quaternion_injection_and_left_error_reset_if_accepted",
]


def _mode(mode: str, dimension: int, coordinates: list[str]) -> dict:
    active_bias = mode == "A"
    return {
        "mode": mode,
        "dimension": dimension,
        "coordinates": coordinates,
        "fixed_dimension_inside_word": True,
        "operators": [
            {
                "name": "commit_previous_tune",
                "kind": "source_state_commit",
                "state_map": "identity_on_error_state",
                "covariance_map": "parameter_commit_only",
            },
            {
                "name": "prediction",
                "kind": "exact_ME KF_prediction".replace(" ", ""),
                "attitude_map": "shipping quaternion gyro propagation with bias error",
                "linear_chain": "exact IntegratedOUChain [v,p,S,a_w] transition",
                "bias_dynamics": "held" if mode == "H" else "active_first_order_Gauss_Markov",
                "covariance_map": "F P F^T + Q",
            },
            {
                "name": "pending_aw_covariance_sync",
                "kind": "PSD_covariance_increment_or_identity",
                "state_map": "identity",
                "covariance_map": "P_awaw <- P_awaw + Delta, Delta >= 0",
            },
            {
                "name": "S_zero",
                "branch_family": ["not_due", "due"],
                "innovation": "r_S=-delta_S",
                "gain": "full implemented K_S=P(:,S)(P_SS+R_S)^-1",
                "accepted_map": "x+=K_S r_S; immediate deployed quaternion injection; immediate left-error covariance reset",
                "S_to_attitude_cross_gain_retained": True,
            },
            {
                "name": "accelerometer",
                "branch_family": ["accepted", "rejected"],
                "innovation": "exact physical specific-force residual, including R_e, a_w, b_a and configured lever-arm terms",
                "gain": "full implemented K_a",
                "accepted_map": "x+=K_a r_a; immediate deployed quaternion injection; immediate left-error covariance reset",
                "rejected_map": "identity correction",
                "bias_projection_after_injection": active_bias,
            },
            {
                "name": "source_evolution",
                "kind": "joint_source_recurrence",
                "map": "measurement-only tuner/EMA evolves candidates; current sample cannot retroactively choose its own gain schedule",
            },
            {
                "name": "aw_sync_tick",
                "kind": "source_scheduler",
                "map": "may stage a future PSD a_w covariance floor; no immediate error-state correction",
            },
            {
                "name": "magnetometer",
                "branch_family": ["not_due", "accepted", "rejected"],
                "innovation": "exact vector residual r_m=y_m-R_wb m_world",
                "gain": "full implemented K_m",
                "accepted_map": "x+=K_m r_m; immediate deployed quaternion injection; immediate left-error covariance reset",
                "rejected_or_not_due_map": "identity correction",
                "bias_projection_after_injection": active_bias,
            },
        ],
        "correction_policy": {
            "deployed_normalized_quaternion_map": True,
            "series_branch_threshold_norm": 1.0e-2,
            "linearized_attitude_injection_allowed": False,
            "one_shared_reset_after_all_measurements": False,
            "left_error_reset_after_each_accepted_correction": True,
            "full_S_to_attitude_cross_gain": True,
        },
        "active_bias_projection": ({
            "kind": "exact_Euclidean_projection_onto_closed_ball",
            "radius_source": "acc_bias_limit_",
            "nonexpansive": True,
        } if active_bias else None),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    manifest = MANIFEST.build()
    mf = MANIFEST.validate(manifest)
    words = WORDS.build(domain_path)
    wf = WORDS.validate(words)
    failures = [f"manifest: {x}" for x in mf] + [f"word-language: {x}" for x in wf]
    order = manifest.get("normal_live_update_order", [])
    if order != EXPECTED_ORDER:
        failures.append("source manifest normal-Live operation order differs from P4 word map")
    reset = manifest.get("same_sample_reset_policy", {})
    if reset.get("single_shared_end_of_sample_reset") is not False:
        failures.append("same-sample reset policy is not source-faithful")
    dims = manifest["state_coordinates"]
    return {
        "schema": SCHEMA,
        "qualification": "SOURCE_BOUND_EXACT_NONLINEAR_OU3_FIXED_MODE_WORD_MAP",
        "source_generated_not_trajectory_fit": True,
        "joint_source_reachability_required": True,
        "shipping_operation_order": order,
        "H": _mode("H", dims["H_dimension"], dims["H"]),
        "A": _mode("A", dims["A_dimension"], dims["A"]),
        "source_word_horizon_s": words["word_contract"]["conditional_word_language"]["word_horizon_lower_s"],
        "word_samples_upper": words["word_contract"]["conditional_word_language"]["word_samples_upper_at_configured_dt"],
        "arbitrary_admissible_rejections_between_required_PE_events": True,
        "dimension_changing_events_inside_word": False,
        "hybrid_events_separate": manifest["hybrid_events"],
        "continuous_gain_covariance_enclosure_supplied_here": False,
        "theorem_promotion": "SEMANTIC_MAP_ONLY",
        "failures": failures,
        "pass": not failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("word map is not source generated")
    if d.get("joint_source_reachability_required") is not True:
        failures.append("word map permits Cartesian edge mixing")
    if d.get("shipping_operation_order") != EXPECTED_ORDER:
        failures.append("shipping operation order mismatch")
    for mode in ("H", "A"):
        m = d.get(mode, {})
        cp = m.get("correction_policy", {})
        if cp.get("deployed_normalized_quaternion_map") is not True:
            failures.append(f"{mode}: deployed quaternion injection missing")
        if cp.get("linearized_attitude_injection_allowed") is not False:
            failures.append(f"{mode}: linearized attitude fallback enabled")
        if cp.get("left_error_reset_after_each_accepted_correction") is not True:
            failures.append(f"{mode}: per-correction reset missing")
        if cp.get("full_S_to_attitude_cross_gain") is not True:
            failures.append(f"{mode}: S-to-attitude cross gain discarded")
    if d.get("dimension_changing_events_inside_word") is not False:
        failures.append("dimension-changing hybrid event placed inside fixed-mode word")
    if not failures and d.get("pass") is not True:
        failures.append("word map did not pass")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain.resolve())
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"pass": not failures, "order": d["shipping_operation_order"], "failures": failures}, indent=2))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
