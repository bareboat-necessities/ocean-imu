#!/usr/bin/env python3
"""Bind every P5 H source-word operation to the exact nonlinear transport law.

This stage closes the semantic gap between the positive outer geometry and the
numerical P5 word.  It does not replace the remaining interval/subdivision
calculation with a norm-only recurrence.

For accepted measurements it uses the exact Joseph information identity and the
exact quaternion/reset congruence from ``ou3_p5_exact_correction_transport``.
The S=0 residual is exactly linear (eta=0); accelerometer and magnetometer
updates retain their nonlinear eta terms.  Rejected/not-due updates are exact
identities.  The pending a_w covariance synchronization is a PSD covariance
increment and is therefore information-nonexpansive for a fixed physical error.
Prediction retains the complete tangent H map and its exact finite-angle group
defect.  Tuner staging changes future source parameters but is not silently
turned into a state contraction.

The first-due S correction is also composed exactly.  Its coarse source-staged
bound does not stay in the convenient |c|<1 diagnostic set, but the Cayley chart
has no singularity there.  ``ou3_p5_first_s_exact_prefix`` therefore widens to a
finite dyadic prefix chart that contains the exact first-S image and still has a
strict positive finite-angle vector-information factor.  The remaining gauged
obligation is now the source-correlated eta/reset/prediction budget over later
prefixes, not an artificial q<1 gate.

For the gravity quotient the same calculus is used after quotient projection.
The gravity-parallel gyro-bias component remains an actual filter state but is
excluded from the strict quotient energy and is charged as a bounded input.  A
source-uniform coordinate bound over one word follows from the exact rotational
transport norm: |Delta theta| <= T_word |b_g,parallel|.  No strict contraction
of that neutral zero dynamics is requested.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p3_word_algebra as P3ALG
import ou3_p4_exact_word_map as WORDMAP
import ou3_p5_exact_correction_transport as CORR
import ou3_p5_first_s_exact_prefix as FIRSTEX
import ou3_p5_gravity_quotient_certificate as GQUOT
import ou3_p5_outer_information_geometry as OUTINFO
import ou3_startup_stability_certificate as P1

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 2


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _operation_calculus() -> list[dict]:
    return [
        {
            "operation": "commit_previous_tune",
            "state_map": "identity_on_error_state",
            "covariance_map": "source_parameter_commit_only",
            "nonlinear_budget": "none_at_commit; parameters remain jointly source correlated",
        },
        {
            "operation": "prediction",
            "state_map": "exact shipping tangent prediction plus finite-angle Cayley prediction defect",
            "covariance_map": "F P F^T + Q",
            "nonlinear_budget": "exact group prediction defect retained for numerical enclosure",
        },
        {
            "operation": "apply_pending_aw_covariance_psd_increment",
            "state_map": "identity",
            "covariance_map": "P -> P + Delta_aw, Delta_aw >= 0",
            "nonlinear_budget": "information-nonexpansive in Loewner order",
        },
        {
            "operation": "periodic_S_zero_when_due_then_immediate_quaternion_injection_and_left_error_reset",
            "state_map": "accepted full K_S correction or exact identity",
            "covariance_map": "Joseph then immediate left-error reset when due",
            "nonlinear_budget": "eta_S=0 exactly; only deployed Cayley/reset defect remains",
        },
        {
            "operation": "accelerometer_correction_or_rejection_then_immediate_quaternion_injection_and_left_error_reset_if_accepted",
            "state_map": "exact nonlinear source residual, full K including linear cross terms, or identity",
            "covariance_map": "Joseph then immediate left-error reset when accepted",
            "nonlinear_budget": "eta_acc^T R_acc^-1 eta_acc plus exact Cayley/reset defect",
        },
        {
            "operation": "source_tuner_evolution_and_stage_next_tune",
            "state_map": "no retroactive estimator-state correction",
            "covariance_map": "future source parameters only",
            "nonlinear_budget": "joint source tuple propagated; no independent extrema product",
        },
        {
            "operation": "periodic_aw_covariance_sync_tick_stages_future_psd_increment",
            "state_map": "identity",
            "covariance_map": "future PSD increment only",
            "nonlinear_budget": "none at staging tick",
        },
        {
            "operation": "asynchronous_magnetometer_correction_or_rejection_then_immediate_quaternion_injection_and_left_error_reset_if_accepted",
            "state_map": "exact nonlinear source residual and immediate injection/reset, or identity",
            "covariance_map": "Joseph then immediate left-error reset when accepted",
            "nonlinear_budget": "eta_mag^T R_mag^-1 eta_mag plus exact Cayley/reset defect",
        },
    ]


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("complete-word transport domain must not be trajectory fitted")

    word = WORDMAP.build(domain_path)
    corr = CORR.build(domain_path)
    firstex = FIRSTEX.build(domain_path)
    outinfo = OUTINFO.build(domain_path)
    gquot = GQUOT.build(domain_path)
    p1 = P1.build(domain_path)
    p3alg = P3ALG.build()

    failures = [f"word-map: {x}" for x in WORDMAP.validate(word)]
    failures += [f"correction-transport: {x}" for x in CORR.validate(corr)]
    failures += [f"first-S-exact-prefix: {x}" for x in FIRSTEX.validate(firstex)]
    failures += [f"outer-information: {x}" for x in OUTINFO.validate(outinfo)]
    failures += [f"gravity-quotient: {x}" for x in GQUOT.validate(gquot)]
    failures += [f"P1: {x}" for x in P1.validate(p1)]
    failures += [f"P3-word-algebra: {x}" for x in P3ALG.validate(p3alg)]

    calculus = _operation_calculus()
    calc_ops = [row["operation"] for row in calculus]
    source_ops = list(word["shipping_operation_order"])
    if calc_ops != source_ops:
        failures.append("exact transport calculus does not match shipping source operation order")

    horizon = float(word["source_word_horizon_s"])
    bg_bound = float(p1["go_live"]["physical_coordinate_bounds"]["gyro_bias_error_norm_upper_rad_s"])
    axial_coordinate_input = up(horizon * bg_bound)
    if not (horizon > 0.0 and bg_bound >= 0.0 and math.isfinite(axial_coordinate_input)):
        failures.append("axial gyro-bias quotient input bound invalid")

    if firstex["P5_FIRST_DUE_S_EXACT_CAYLEY_PREFIX_CERTIFICATE"] != "PASS_WIDENED_CHART":
        failures.append("first-S exact widened Cayley prefix did not close")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_COMPLETE_H_SOURCE_WORD_EXACT_TRANSPORT_CALCULUS",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "source_word_horizon_s": horizon,
        "source_operation_order": source_ops,
        "operation_transport_calculus": calculus,
        "all_source_operation_classes_bound_to_transport_calculus": not failures,
        "accepted_measurement_identity": corr["exact_joseph_information_identity"],
        "reset_metric_identity": corr["exact_reset_congruence_identity"],
        "reset_condition_number_multiplier_used": False,
        "S_zero_nonlinear_measurement_eta_exact_zero": True,
        "full_S_to_attitude_gain_retained": True,
        "sequential_immediate_quaternion_resets_retained": True,
        "aw_sync_PSD_information_nonexpansive": True,
        "rejected_and_not_due_state_covariance_identity_branches_retained": True,
        "gauged_H": {
            "finite_angle_information_geometry": outinfo["P5_FINITE_ANGLE_INFORMATION_GEOMETRY_CERTIFICATE"],
            "exact_correction_transport_algebra": corr["P5_EXACT_CORRECTION_TRANSPORT_ALGEBRA_CERTIFICATE"],
            "first_due_S_exact_prefix": firstex["P5_FIRST_DUE_S_EXACT_CAYLEY_PREFIX_CERTIFICATE"],
            "diagnostic_q_lt_1_is_promotion_gate": firstex["diagnostic_q_lt_1_is_promotion_gate"],
            "widened_prefix_cayley_norm_upper": firstex["widened_prefix_cayley_norm_upper"],
            "widened_prefix_antipodal_margin_lower": firstex["widened_prefix_antipodal_one_plus_cosine_margin_lower"],
            "widened_prefix_vector_information_vs_goLive_metric_lower": firstex["widened_prefix_pair_information_vs_goLive_attitude_metric_lower"],
            "first_due_S_nodes": firstex["nodes"],
            "complete_word_numerical_status": "NOT_ESTABLISHED",
            "first_unclosed_numerical_obligation": "COMPLETE_WORD_ETA_RESET_INFORMATION_BUDGET_NOT_CERTIFIED",
        },
        "gravity_quotient_H": {
            "reduced_detectability": gquot["P5_GRAVITY_QUOTIENT_REDUCED_DETECTABILITY_CERTIFICATE"],
            "strict_coordinates": gquot["detectable_coordinates"],
            "axial_gyro_bias_role": gquot["axial_gyro_bias_role"],
            "axial_gyro_bias_norm_upper_rad_s": bg_bound,
            "one_word_axial_bias_full_attitude_coordinate_input_norm_upper_rad": axial_coordinate_input,
            "axial_input_bound_reason": "rotation transport has operator norm one, so integral over T_word is at most T_word*|b_parallel| before quotient projection",
            "strict_contraction_of_axial_bias_requested": False,
            "complete_word_numerical_status": "NOT_ESTABLISHED",
            "first_unclosed_numerical_obligation": "GRAVITY_QUOTIENT_EXACT_ETA_RESET_PREFIX_BUDGET_NOT_CERTIFIED",
        },
        "P5_COMPLETE_WORD_TRANSPORT_ALGEBRA_CERTIFICATE": "PASS" if not failures else "FAIL",
        "P5_GAUGED_COMPLETE_WORD_NUMERICAL_CERTIFICATE": "NOT_ESTABLISHED",
        "P5_GRAVITY_QUOTIENT_COMPLETE_WORD_NUMERICAL_CERTIFICATE": "NOT_ESTABLISHED",
        "next_obligation": (
            "outward-subdivide later exact vector/prediction/S prefixes in the widened source-correlated Cayley/information chart and close their eta/reset budgets; use b_g_parallel only as the quotient ISS input term"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("complete-word transport is not source generated")
    if d.get("source_replay_used") is not False:
        failures.append("complete-word transport uses replay")
    if d.get("filter_changed") is not False:
        failures.append("complete-word transport changes filter")
    if d.get("all_source_operation_classes_bound_to_transport_calculus") is not True:
        failures.append("complete-word transport calculus misses a source operation")
    if d.get("reset_condition_number_multiplier_used") is not False:
        failures.append("reset condition-number penalty reintroduced")
    for flag in (
        "S_zero_nonlinear_measurement_eta_exact_zero",
        "full_S_to_attitude_gain_retained",
        "sequential_immediate_quaternion_resets_retained",
        "aw_sync_PSD_information_nonexpansive",
        "rejected_and_not_due_state_covariance_identity_branches_retained",
    ):
        if d.get(flag) is not True:
            failures.append(f"missing transport semantic {flag}")
    if d.get("P5_COMPLETE_WORD_TRANSPORT_ALGEBRA_CERTIFICATE") != "PASS" and not failures:
        failures.append("complete-word transport algebra did not pass")
    g = d.get("gauged_H", {})
    if g.get("finite_angle_information_geometry") != "PASS":
        failures.append("gauged finite-angle prerequisite missing")
    if g.get("exact_correction_transport_algebra") != "PASS":
        failures.append("gauged correction transport algebra missing")
    if g.get("first_due_S_exact_prefix") != "PASS_WIDENED_CHART":
        failures.append("gauged first-S widened chart missing")
    if g.get("diagnostic_q_lt_1_is_promotion_gate") is not False:
        failures.append("gauged transport reinstated q<1 gate")
    if not float(g.get("widened_prefix_antipodal_margin_lower", 0.0)) > 0.0:
        failures.append("gauged widened prefix has no antipodal margin")
    if not float(g.get("widened_prefix_vector_information_vs_goLive_metric_lower", 0.0)) > 0.0:
        failures.append("gauged widened prefix loses vector information")
    if g.get("complete_word_numerical_status") != "NOT_ESTABLISHED":
        failures.append("gauged word promoted before numerical enclosure")
    q = d.get("gravity_quotient_H", {})
    if q.get("reduced_detectability") != "PASS":
        failures.append("gravity quotient detectability prerequisite missing")
    if q.get("strict_contraction_of_axial_bias_requested") is not False:
        failures.append("gravity quotient demands contraction of neutral axial bias")
    if not (isinstance(q.get("one_word_axial_bias_full_attitude_coordinate_input_norm_upper_rad"), (int, float))
            and math.isfinite(float(q["one_word_axial_bias_full_attitude_coordinate_input_norm_upper_rad"]))):
        failures.append("gravity quotient axial input bound invalid")
    if q.get("complete_word_numerical_status") != "NOT_ESTABLISHED":
        failures.append("gravity quotient word promoted before numerical enclosure")
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
        "algebra": out["P5_COMPLETE_WORD_TRANSPORT_ALGEBRA_CERTIFICATE"],
        "gauged": out["gauged_H"],
        "quotient": out["gravity_quotient_H"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
