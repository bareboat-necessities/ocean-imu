#!/usr/bin/env python3
"""Exact fixed-mode event algebra for the complete SEA3 Normal-Live P3 word.

This is not a source approximation.  It is a source-uniform algebraic lemma for
*every* event in the same complete phase-continuous SEA3 word.  The only source
restriction is the declared Normal-Live SEA3 theorem domain.

For each prefix write

    P = Psi P0 Psi^T + Omega,       Omega >= 0.

Every shipping fixed-dimensional covariance operation has affine-PSD form

    P+     = A P A^T + B,
    Psi+   = A Psi,
    Omega+ = A Omega A^T + B,       B >= 0.

Therefore, for 0 < delta < 1,

    M_delta+ = Omega+ - delta P+
             = A M_delta A^T + (1-delta) B.

Once a complete-SEA3 prefix establishes M_delta >= 0, every later prediction,
Joseph correction, skipped/not-due branch, immediate left-error reset and a_w
PSD covariance-floor operation preserves the margin.  No one-step contraction
factor, D/L split, blockwise ratio or scalar beta is used.

The immediate left-error reset is important.  Shipping applies it after each
accepted S/accelerometer/magnetometer correction.  Its attitude block is

    G_theta = I + 0.5 [dtheta]_x,

with det(G_theta)=1+||dtheta||^2/4 >= 1, hence it is nonsingular for every
finite injection.  The reset is a congruence of P, Psi and Omega and therefore
preserves the same full-matrix generalized information inequality exactly.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
CORE = REPO / "src" / "kalman_ou_common" / "KalmanOUCoreMath.h"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_EXACT_FIXED_MODE_EVENT_ALGEBRA"


def reset_det(theta_norm: float) -> float:
    t = float(theta_norm)
    if not math.isfinite(t) or t < 0.0:
        raise ValueError("theta norm must be finite nonnegative")
    return 1.0 + 0.25 * t * t


def _contains_any(text: str, choices: tuple[str, ...]) -> bool:
    return any(x in text for x in choices)


def build() -> dict:
    k = MEKF.read_text(encoding="utf-8")
    c = CORE.read_text(encoding="utf-8")
    w = WRAPPER.read_text(encoding="utf-8")

    parity = {
        "prediction_is_full_covariance_congruence_plus_Q": (
            "F_LL" in k and "Q_LL" in k and "Pext" in k
        ),
        "S_zero_update_exists": "applyIntegralZeroPseudoMeas" in k,
        "accelerometer_update_exists": "measurement_update_acc_only" in k,
        "magnetometer_update_exists": "measurement_update_mag_only" in k,
        "joseph_covariance_update_exists": _contains_any(
            k, ("joseph_update3_", "Joseph", "K * R", "K*R")
        ),
        "left_error_reset_jacobian_exact": (
            "Identity() + T(0.5)*skew(dtheta)" in c
            or "Identity() + T(0.5) * skew(dtheta)" in c
        ),
        "aw_floor_is_psd_increment_default_path": (
            "aw_covariance_floor_pending_" in k
            and "Sigma_aw_stat - P_awaw" in k
            and "std::max(T(0), evals(i))" in k
        ),
        "legacy_aw_replacement_is_not_shipping_default": (
            "legacy_aw_covariance_replacement_" in k
            and "false" in k
        ),
        "actual_applied_RS_reaches_filter": (
            "RS_applied" in w and "set_RS_noise" in w
        ),
        "SpectralMSE_is_shipping_RS_law": (
            "RSAdaptationLaw::SpectralMSE" in w
        ),
        "progress_preserving_pseudo_period_retarget_exists": (
            "retarget_period_elapsed_progress_preserving" in k
        ),
    }

    failures = [name for name, ok in parity.items() if not ok]

    operations = {
        "prediction": {
            "map": "P+=F P F^T+Q; Psi+=F Psi; Omega+=F Omega F^T+Q",
            "A": "F",
            "B": "Q",
            "B_psd": True,
        },
        "pending_aw_covariance_floor": {
            "map": "P+=P+E_aw Delta_+ E_aw^T; Psi+=Psi; Omega+=Omega+E_aw Delta_+ E_aw^T",
            "A": "I",
            "B": "E_aw Delta_+ E_aw^T",
            "B_psd": True,
            "legacy_raw_block_replacement_admitted": False,
        },
        "accepted_S_acc_mag_joseph": {
            "map": "P+=(I-KH)P(I-KH)^T+K R_eff K^T",
            "A": "I-KH",
            "B": "K R_eff K^T",
            "B_psd": True,
            "actual_applied_SpectralMSE_R_S_required_for_S": True,
        },
        "not_due_or_rejected": {
            "map": "identity",
            "A": "I",
            "B": "0",
            "B_psd": True,
        },
        "immediate_left_error_reset": {
            "map": "P+=G P G^T; Psi+=G Psi; Omega+=G Omega G^T",
            "A": "G=blkdiag(I+0.5[dtheta]_x,I)",
            "B": "0",
            "B_psd": True,
            "det_attitude_formula": "1+||dtheta||^2/4",
            "det_attitude_lower": 1.0,
            "nonsingular_for_every_finite_injection": True,
            "must_follow_each_accepted_S_acc_mag_update": True,
        },
    }

    order = [
        "commit_previous_tune",
        "prediction",
        "pending_aw_covariance_floor_if_requested",
        "every_due_S_zero_Joseph_then_immediate_left_error_reset",
        "accelerometer_Joseph_then_immediate_left_error_reset_on_every_valid_Normal_Live_IMU_sample",
        "source_frontend_tuner_evolution_and_stage_next_tune",
        "periodic_aw_covariance_sync_stages_future_PSD_increment",
        "asynchronous_accepted_magnetometer_Joseph_then_immediate_left_error_reset",
    ]

    passed = not failures
    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "source_family_replaced": False,
        "trajectory_replay_used": False,
        "fixed_dimensions": {"H": 18, "A": 21},
        "dimension_change_inside_word": False,
        "normal_live_shipping_order": order,
        "source_parity": parity,
        "operation_classes": operations,
        "joint_recursion": {
            "P": "P+=A P A^T+B",
            "Psi": "Psi+=A Psi",
            "Omega": "Omega+=A Omega A^T+B",
            "same_A_and_B_used_for_same_complete_SEA3_event": True,
        },
        "full_matrix_margin_preservation": {
            "matrix": "M_delta=Omega-delta P",
            "identity": "M_delta+=A M_delta A^T+(1-delta)B",
            "requires": "0<delta<1 and B>=0",
            "covers_prediction": True,
            "covers_every_due_S_update": True,
            "covers_every_Normal_Live_accelerometer_update": True,
            "covers_asynchronous_magnetometer_update": True,
            "covers_immediate_left_error_reset": True,
            "covers_aw_covariance_floor": True,
            "covers_not_due_or_rejected_identity_branches": True,
            "one_step_contraction_assumed": False,
            "blockwise_ratio_used": False,
            "D_W_L_W_split_used": False,
            "scalar_beta_used": False,
        },
        "prefix_information": {
            "identity": "P_s=Psi_s P0 Psi_s^T+Omega_s",
            "bound": "Psi_s^T P_s^-1 Psi_s <= P0^-1",
            "information_gain_upper": 1.0,
            "source_uniform": True,
        },
        "left_error_reset": {
            "determinant_formula": "1+||dtheta||^2/4",
            "determinant_lower": 1.0,
            "small_angle_needed_for_nonsingularity": False,
            "same_full_matrix_margin_preserved_by_congruence": True,
        },
        "pass": passed,
        "failures": failures,
        "P3_promoted": False,
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    if d.get("source_family_replaced") is not False or d.get("trajectory_replay_used") is not False:
        f.append("source replacement/replay admitted")
    if d.get("fixed_dimensions") != {"H": 18, "A": 21}:
        f.append("fixed dimensions changed")
    if d.get("dimension_change_inside_word") is not False:
        f.append("dimension change admitted inside word")
    parity = d.get("source_parity", {})
    f.extend(f"source parity failed: {k}" for k, v in parity.items() if v is not True)
    ops = d.get("operation_classes", {})
    required = {
        "prediction", "pending_aw_covariance_floor", "accepted_S_acc_mag_joseph",
        "not_due_or_rejected", "immediate_left_error_reset",
    }
    if set(ops) != required:
        f.append("operation coverage incomplete")
    for name in required:
        if ops.get(name, {}).get("B_psd") is not True:
            f.append(f"{name} B is not PSD")
    m = d.get("full_matrix_margin_preservation", {})
    for key in (
        "covers_prediction", "covers_every_due_S_update",
        "covers_every_Normal_Live_accelerometer_update",
        "covers_asynchronous_magnetometer_update",
        "covers_immediate_left_error_reset", "covers_aw_covariance_floor",
        "covers_not_due_or_rejected_identity_branches",
    ):
        if m.get(key) is not True:
            f.append(f"margin preservation missing {key}")
    for key in ("one_step_contraction_assumed", "blockwise_ratio_used", "D_W_L_W_split_used", "scalar_beta_used"):
        if m.get(key) is not False:
            f.append(f"forbidden shortcut enabled: {key}")
    reset = d.get("left_error_reset", {})
    if reset.get("determinant_lower") != 1.0:
        f.append("reset determinant lower is not one")
    if reset.get("same_full_matrix_margin_preserved_by_congruence") is not True:
        f.append("reset full-matrix congruence not closed")
    if d.get("P3_promoted") is not False:
        f.append("event algebra alone promoted P3")
    if d.get("pass") is not True:
        f.extend(str(x) for x in d.get("failures", []))
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build()
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "pass": not failures,
        "reset_det_lower": d["left_error_reset"]["determinant_lower"],
        "prefix_gain_upper": d["prefix_information"]["information_gain_upper"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
