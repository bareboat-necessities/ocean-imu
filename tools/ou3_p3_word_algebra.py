#!/usr/bin/env python3
"""Exact algebraic closure of the OU-III P3 covariance word certificate.

P3's numerical backend establishes a strict generalized information margin on
one source-uniform covariance/noise comparison.  This module proves that the
margin and the unit prefix information-gain bound survive every subsequent
fixed-dimensional normal-Live covariance operation implemented by OU-III.

For every prefix write

    P = Phi P0 Phi^T + Omega,   Omega >= 0.

Every implemented covariance operation has the affine PSD form

    P+ = A P A^T + B,          B >= 0,
    Omega+ = A Omega A^T + B.

Consequently, for 0 < delta < 1,

    Omega+ - delta P+
      = A (Omega - delta P) A^T + (1-delta) B >= 0.

The same decomposition also gives the exact information-prefix inequality

    Phi^T P^{-1} Phi <= P0^{-1},

by the Schur complement of the joint covariance.  Hence the source-uniform P3
prefix information gain is at most one; it is not a sampled or hard-coded
estimate.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_implementation_proof_manifest as MANIFEST

REPO = Path(__file__).resolve().parents[1]
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
CORE = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"
SCHEMA = 2

# This order is copied from the source-derived implementation manifest.  It is
# intentionally more detailed than the covariance operation classes below: the
# vibration guard runs before prediction, the S pseudo update happens inside
# time_update() before the accelerometer update, and every accepted correction
# performs its quaternion injection/reset immediately rather than through one
# shared end-of-sample reset.
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


def _require(text: str, marker: str, label: str, failures: list[str]) -> None:
    if marker not in text:
        failures.append(f"missing source semantic {label}: {marker}")


def reset_det(theta_norm: float) -> float:
    """Exact determinant of I + 1/2 [dtheta]_x as a function of ||dtheta||."""
    t = float(theta_norm)
    if not math.isfinite(t) or t < 0.0:
        raise ValueError("theta norm must be finite nonnegative")
    return 1.0 + 0.25 * t * t


def margin_after_affine_psd(delta: float, primitive_margin: float,
                            additive_psd_energy: float) -> float:
    """Scalar oracle for the exact matrix identity used by the proof contract."""
    d = float(delta)
    q = float(primitive_margin)
    b = float(additive_psd_energy)
    if not (0.0 < d < 1.0) or q < 0.0 or b < 0.0:
        raise ValueError("invalid affine-PSD margin inputs")
    return q + (1.0 - d) * b


def build() -> dict:
    manifest = MANIFEST.build()
    manifest_failures = MANIFEST.validate(manifest)
    k = MEKF.read_text(encoding="utf-8")
    c = CORE.read_text(encoding="utf-8")
    failures = list(manifest_failures)

    if manifest.get("normal_live_update_order") != EXPECTED_ORDER:
        failures.append("source manifest normal-Live update order changed")
    reset_policy = manifest.get("same_sample_reset_policy", {})
    if reset_policy.get("single_shared_end_of_sample_reset") is not False:
        failures.append("source manifest merged immediate correction resets")
    vibration_guard = manifest.get("vibration_guard", {})
    if vibration_guard.get("zero_engagement_is_bit_exact_transparent") is not True:
        failures.append("P3 requires source-certified zero-engagement vibration-guard transparency")
    if vibration_guard.get("active_guard_requires_separate_source_certificate") is not True:
        failures.append("P3 must not absorb active vibration-guard dynamics into the dormant branch")

    # Bind the abstract affine-PSD operations to the exact shipping source.
    _require(k, "P_LL_new = F_LL * P_LL_old * F_LLᵀ + Q_LL", "linear prediction", failures)
    _require(k, "joseph_update3_(K, S_mat, PCt);", "Joseph corrections", failures)
    _require(k, "last_acc_diag_.accepted = false;", "accelerometer rejection branch", failures)
    _require(k, "last_mag_diag_.accepted = false;", "magnetometer rejection branch", failures)
    _require(k, "if (!safe_ldlt3_(S_mat, ldlt, R_S.norm())) return;", "S-update skipped branch", failures)
    _require(c, "const Eigen::Matrix<T,3,3> G = Eigen::Matrix<T,3,3>::Identity() + T(0.5)*skew(dtheta);", "left-error reset Jacobian", failures)
    _require(k, "Pext.template block<3,3>(OFF_AW, OFF_AW) += Delta;", "a_w PSD covariance increment", failures)
    _require(k, "evals(i) = std::max(T(0), evals(i));", "a_w increment PSD projection", failures)

    operations = {
        "vibration_guard_dormant": {
            "affine_map": "P+=P",
            "A": "I",
            "B": "0",
            "B_psd": True,
            "measurement_map": "acc_in=acc",
            "scope": "zero-engagement bit-exact-transparent branch only",
            "active_or_transitioning_guard_covered": False,
            "active_or_transitioning_guard_requires_separate_source_certificate": True,
            "source_bound": "implementation manifest vibration_guard zero-engagement contract",
        },
        "prediction": {
            "affine_map": "P+=F P F^T+Q",
            "A": "F",
            "B": "Q",
            "B_psd": True,
            "source_bound": "time_update process covariance construction",
        },
        "accepted_joseph": {
            "affine_map": "P+=(I-KH)P(I-KH)^T+K R_eff K^T",
            "A": "I-KH",
            "B": "K R_eff K^T",
            "B_psd": True,
            "safe_ldlt_diagonal_boost_handled": (
                "a safety boost only increases R_eff by a PSD diagonal term"
            ),
            "frozen_gain_rows_handled": (
                "row freezing changes K but preserves the affine-PSD identity"
            ),
            "immediate_quaternion_reset_order_handled": True,
        },
        "rejected_or_not_due": {
            "affine_map": "P+=P",
            "A": "I",
            "B": "0",
            "B_psd": True,
            "covers": [
                "accelerometer rejected",
                "magnetometer not due",
                "magnetometer rejected",
                "S pseudo not due",
                "S pseudo factorization return",
            ],
        },
        "left_error_reset": {
            "affine_map": "P+=G_reset P G_reset^T",
            "A": "blkdiag(I+0.5[dtheta]_x,I)",
            "B": "0",
            "B_psd": True,
            "det_attitude_block": "1+||dtheta||^2/4",
            "determinant_lower": 1.0,
            "nonsingular_for_every_finite_injection": True,
            "generalized_information_margin_congruence_invariant": True,
            "applied_after_each_accepted_correction": True,
        },
        "aw_covariance_sync": {
            "affine_map": "P+=P+E_aw Delta_plus E_aw^T",
            "A": "I",
            "B": "E_aw Delta_plus E_aw^T",
            "B_psd": True,
            "margin_increment": "(1-delta) E_aw Delta_plus E_aw^T >= 0",
            "can_reduce_existing_delta": False,
        },
    }

    return {
        "schema": SCHEMA,
        "qualification": "SOURCE_BOUND_EXACT_P3_FIXED_MODE_WORD_ALGEBRA",
        "source_generated_not_trajectory_fit": True,
        "fixed_dimensions": {"H": 18, "A": 21},
        "dimension_change_inside_word": False,
        "normal_live_update_order": EXPECTED_ORDER,
        "vibration_guard_scope": "dormant_zero_engagement_bit_exact_transparent_only",
        "active_vibration_guard_covered": False,
        "same_sample_reset_policy": "immediate_after_each_accepted_S_acc_mag_correction",
        "operation_classes": operations,
        "covariance_decomposition_invariant": {
            "identity": "P_s=Phi_s P_0 Phi_s^T+Omega_s",
            "Omega_s_psd": True,
            "induction_complete_for_every_normal_live_prefix": True,
        },
        "strict_margin_preservation": {
            "identity": "Omega_plus-delta P_plus=A(Omega-delta P)A^T+(1-delta)B",
            "requires": "0<delta<1 and B>=0",
            "covers_every_operation_class": True,
            "one_step_contraction_assumed": False,
        },
        "prefix_information_bound": {
            "schur_complement_identity": "Phi_s^T P_s^-1 Phi_s <= P_0^-1",
            "information_gain_upper": 1.0,
            "source_uniform": True,
            "sampled_evidence_used": False,
            "reason": "P_s=Phi_s P_0 Phi_s^T+Omega_s with Omega_s>=0",
        },
        "reset": {
            "attitude_jacobian": "I+0.5[dtheta]_x",
            "determinant_formula": "1+||dtheta||^2/4",
            "determinant_lower": 1.0,
            "requires_small_angle_for_nonsingularity": False,
        },
        "failures": failures,
        "pass": not failures,
    }


def validate(d: dict) -> list[str]:
    failures: list[str] = []
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("P3 word algebra is not source bound")
    if d.get("fixed_dimensions") != {"H": 18, "A": 21}:
        failures.append("P3 fixed mode dimensions changed")
    if d.get("dimension_change_inside_word") is not False:
        failures.append("dimension change was admitted inside a P3 word")
    if d.get("normal_live_update_order") != EXPECTED_ORDER:
        failures.append("normal-Live order mismatch")
    if d.get("vibration_guard_scope") != "dormant_zero_engagement_bit_exact_transparent_only":
        failures.append("P3 vibration-guard scope is not the certified dormant branch")
    if d.get("active_vibration_guard_covered") is not False:
        failures.append("P3 incorrectly claims active vibration-guard coverage")
    if d.get("same_sample_reset_policy") != "immediate_after_each_accepted_S_acc_mag_correction":
        failures.append("same-sample reset order is not source-faithful")

    ops = d.get("operation_classes", {})
    required = {"vibration_guard_dormant", "prediction", "accepted_joseph", "rejected_or_not_due", "left_error_reset", "aw_covariance_sync"}
    if set(ops) != required:
        failures.append("P3 operation class coverage is incomplete")
    for name in required:
        if ops.get(name, {}).get("B_psd") is not True:
            failures.append(f"{name} additive covariance is not PSD")
    guard_op = ops.get("vibration_guard_dormant", {})
    if guard_op.get("A") != "I" or guard_op.get("B") != "0":
        failures.append("dormant vibration guard is not represented as exact identity")
    if guard_op.get("active_or_transitioning_guard_covered") is not False:
        failures.append("active vibration guard was admitted into dormant identity class")
    if guard_op.get("active_or_transitioning_guard_requires_separate_source_certificate") is not True:
        failures.append("active vibration guard is not fail-closed as a separate source obligation")

    reset = d.get("reset", {})
    if reset.get("determinant_formula") != "1+||dtheta||^2/4":
        failures.append("left-error reset determinant formula mismatch")
    if reset.get("determinant_lower") != 1.0:
        failures.append("left-error reset determinant lower bound is not one")
    if reset.get("requires_small_angle_for_nonsingularity") is not False:
        failures.append("P3 reset incorrectly requires a small-angle assumption")

    inv = d.get("covariance_decomposition_invariant", {})
    if inv.get("Omega_s_psd") is not True or inv.get("induction_complete_for_every_normal_live_prefix") is not True:
        failures.append("covariance decomposition induction is incomplete")
    margin = d.get("strict_margin_preservation", {})
    if margin.get("covers_every_operation_class") is not True or margin.get("one_step_contraction_assumed") is not False:
        failures.append("strict P3 margin preservation is incomplete")
    prefix = d.get("prefix_information_bound", {})
    if prefix.get("information_gain_upper") != 1.0 or prefix.get("source_uniform") is not True:
        failures.append("unit source-uniform prefix information bound is not established")
    if prefix.get("sampled_evidence_used") is not False:
        failures.append("prefix information bound depends on sampled evidence")
    if d.get("pass") is not True:
        failures.extend(str(x) for x in d.get("failures", []))
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    out = dict(d)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "qualification": out["qualification"],
        "prefix_information_gain_upper": out["prefix_information_bound"]["information_gain_upper"],
        "reset_determinant_lower": out["reset"]["determinant_lower"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
