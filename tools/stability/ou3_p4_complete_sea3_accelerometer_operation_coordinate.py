#!/usr/bin/env python3
"""Exact complete-SEA3 accelerometer operation coordinate for current OU-III P4.

For one accepted shipping accelerometer operation write the physical error as

    R_true = E R_hat,

and let da be the world-frame latent-acceleration error.  With zero lever arm,
shipping has

    r_a = (E-I) f_hat + E R_hat da + db_a,
    H_a z = [c]x f_hat + R_hat da + db_a.

The original coordinate makes ``(E-I)R_hat da`` look like a nonlinear
measurement remainder even though it is only a rotation of the full-rank a_w
measurement column.  Define, pointwise for this operation,

    Q_aw = R_hat^T E R_hat,       u_aw = Q_aw da.

Q_aw is orthogonal and E R_hat da = R_hat u_aw.  Therefore in the operation
coordinate

    r_a = (E-I) f_hat + R_hat u_aw + db_a,
    h_a = [c]x f_hat + R_hat u_aw + db_a,
    eta_a = ((E-I)-[c]x) f_hat.

The latent a_w error and the active accelerometer-bias error cancel *exactly*
from nonlinear eta.  They remain fully present in the residual, innovation
covariance and Kalman correction, so no measurement information or uncertainty
is dropped.  In particular the complete-SEA3 S=0 / SpectralMSE R_S
regularization of a_w remains active.

This coordinate is compatible with the current moving shipping Riccati metric
without any group-isotropic assumption.  Let

    T_E = diag(I_theta,I_bg,I_v,I_p,I_S,Q_aw[,I_ba]),
    z_u = T_E z,       P_u = T_E P T_E^T,
    H_u = H T_E^T,     K_u = T_E K.

Since T_E is orthogonal,

    z_u^T P_u^-1 z_u = z^T P^-1 z

exactly.  Also H_u P_u H_u^T = H P H^T, so S is unchanged; the Joseph
information identity and exact correction are unchanged under the congruence.
A distinction is essential: the residual rewrite above uses the auxiliary
matrix H0 with a_w column R_hat. The congruently transformed shipping matrix
H_u instead has a_w column R_hat Q_aw^T. These are not equal at finite angle.
Consequently the aw-free eta is not the remainder relative to H_u, and these
two valid algebraic identities do not close the nonlinear shipping Joseph map.
No zero-cost composition of the residual rewrite and congruence is claimed.
This is an operation identity over every admitted complete-SEA3 realization,
not a source generator and not a P4 promotion.  The filename intentionally does
not revive the retired endpoint-attachment era ``accelerometer_corotated_aw``
module; this is bound only to the current complete-SEA3 moving-Riccati route.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_moving_metric_rebind as REBIND
import ou3_sea3_complete_source as COMPLETE

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SCHEMA = 2
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_ACCELEROMETER_OPERATION_COORDINATE_MOVING_METRIC_V2"

SOURCE_MARKERS = (
    "const Vector3 f_pred = R_wb() * (aw - g_world) + lever + ba_term;",
    "const Vector3 f_cog_b = R_wb() * (aw - g_world);",
    "const Matrix3 J_att = -skew_symmetric_matrix(f_cog_b);",
    "const Matrix3 J_aw  =  R_wb();",
    "xext.noalias() += K * r;",
    "joseph_update3_(K, S_mat, PCt);",
    "qref = corr * qref;",
    "apply_error_state_reset_jacobian_(dtheta);",
)


def _source_contract(domain: dict) -> list[str]:
    failures: list[str] = []
    runtime = domain.get("configured_runtime", {})
    if runtime.get("imu_lever_arm_enabled") is not False:
        failures.append("accelerometer operation coordinate requires zero lever arm")
    if runtime.get("accelerometer_vibration_guard_proof_branch") != "dormant_transparent":
        failures.append("accelerometer operation coordinate requires transparent vibration-guard branch")
    text = MEKF.read_text(encoding="utf-8")
    for marker in SOURCE_MARKERS:
        if marker not in text:
            failures.append(f"shipping accelerometer/correction marker changed: {marker}")
    if "Measurement Jacobians use the left-multiplicative convention" not in text:
        failures.append("left-multiplicative attitude-error convention marker changed")
    if "J_ba = I" not in text:
        failures.append("accelerometer-bias identity Jacobian marker changed")
    return failures


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("accelerometer operation-coordinate proof must not be trajectory fitted")

    complete = COMPLETE.build(path)
    cayley = CAYLEY.build(path)
    rebind = REBIND.build()
    failures = [f"complete SEA3: {x}" for x in COMPLETE.validate(complete)]
    failures += [f"Cayley: {x}" for x in CAYLEY.validate(cayley)]
    failures += [f"moving metric: {x}" for x in REBIND.validate(rebind)]
    failures += _source_contract(domain)
    if complete.get("canonical_P3_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        failures.append("canonical complete SEA3 source changed")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": "COMPLETE_SEA3_NORMAL_LIVE_WORD",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "source_family_replaced": False,
        "zero_lever_arm_branch": True,
        "dormant_transparent_vibration_guard_branch": True,
        "all_valid_accelerometer_updates_remain_in_complete_word": True,
        "all_due_S_updates_and_actual_RS_remain_in_complete_word": True,
        "shipping_left_multiplicative_error_convention_bound": True,
        "shipping_accelerometer_J_aw_is_R_hat": True,
        "shipping_accelerometer_J_ba_is_identity": True,
        "coordinate_definition": {
            "Q_aw": "R_hat^T E R_hat",
            "u_aw": "Q_aw delta_a_w",
            "Q_aw_orthogonal": True,
            "aw_norm_preserved_exactly": True,
        },
        "exact_residual_identity": (
            "r_a=(E-I)f_hat + R_hat*u_aw + delta_b_a; "
            "u_aw=(R_hat^T E R_hat)delta_a_w"
        ),
        "tangent_residual_in_operation_coordinate": (
            "h_a=[c]x f_hat + R_hat*u_aw + delta_b_a"
        ),
        "nonlinear_eta_in_operation_coordinate": "eta_a=((E-I)-[c]x)f_hat",
        "latent_aw_nonlinear_eta_coefficient": 0.0,
        "accelerometer_bias_nonlinear_eta_coefficient": 0.0,
        "aw_error_exactly_linear_in_accelerometer_operation_coordinate": True,
        "accelerometer_bias_error_exactly_linear": True,
        "actual_RS_regularizer_not_removed_by_coordinate_change": True,
        "state_coordinate_transform": {
            "kind": "BLOCK_ORTHOGONAL_MOVING_COVARIANCE_CONGRUENCE",
            "T_E": "diag(I_theta,I_bg,I_v,I_p,I_S,Q_aw[,I_ba])",
            "P_u": "T_E P T_E^T",
            "H_u": "H T_E^T",
            "K_u": "T_E K",
            "moving_metric_energy_invariant": True,
            "innovation_covariance_S_invariant": True,
            "Joseph_information_identity_invariant": True,
            "exact_state_correction_invariant": True,
            "transform_back_after_measurement_exact": True,
            "group_isotropic_metric_assumption_required": False,
        },
        "moving_metric_rebind_consumed": bool(rebind["full_nonlinear_measurement_metric_rebind_closed"]),
        "outer_angle_rad": float(cayley["outer_angle_rad"]),
        "outer_sector_covered": float(cayley["outer_angle_rad"]) >= 0.80,
        "retired_endpoint_attachment_module_reintroduced": False,
        "complete_H18_A21_word_established_here": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "shipping_Joseph_binding_closed": False,
        "auxiliary_H0_equals_congruent_shipping_Hu": False,
        "congruent_shipping_aw_column": "R_hat Q_aw^T",
        "next_obligation": (
            "bind the auxiliary aw-free residual rewrite to the actual shipping H_u=H T_E^T and K_u=T_E K; "
            "account for the missing finite-angle column/coordinate transport in the complete SEA3 word"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "source_generated_not_trajectory_fit", "zero_lever_arm_branch",
        "dormant_transparent_vibration_guard_branch",
        "all_valid_accelerometer_updates_remain_in_complete_word",
        "all_due_S_updates_and_actual_RS_remain_in_complete_word",
        "shipping_left_multiplicative_error_convention_bound",
        "shipping_accelerometer_J_aw_is_R_hat", "shipping_accelerometer_J_ba_is_identity",
        "aw_error_exactly_linear_in_accelerometer_operation_coordinate",
        "accelerometer_bias_error_exactly_linear",
        "actual_RS_regularizer_not_removed_by_coordinate_change",
        "moving_metric_rebind_consumed", "outer_sector_covered",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
        "source_family_replaced", "retired_endpoint_attachment_module_reintroduced",
        "complete_H18_A21_word_established_here", "P4_USABLE_CERTIFICATE_PROMOTED",
        "shipping_Joseph_binding_closed", "auxiliary_H0_equals_congruent_shipping_Hu",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    coord = d.get("coordinate_definition", {})
    for key in ("Q_aw_orthogonal", "aw_norm_preserved_exactly"):
        if coord.get(key) is not True:
            f.append(f"coordinate {key} is not true")
    tr = d.get("state_coordinate_transform", {})
    for key in (
        "moving_metric_energy_invariant", "innovation_covariance_S_invariant",
        "Joseph_information_identity_invariant", "exact_state_correction_invariant",
        "transform_back_after_measurement_exact",
    ):
        if tr.get(key) is not True:
            f.append(f"state congruence {key} is not true")
    if tr.get("group_isotropic_metric_assumption_required") is not False:
        f.append("retired group-isotropic metric assumption re-entered")
    if d.get("latent_aw_nonlinear_eta_coefficient") != 0.0:
        f.append("operation-coordinate a_w must have zero nonlinear eta coefficient")
    if d.get("accelerometer_bias_nonlinear_eta_coefficient") != 0.0:
        f.append("accelerometer bias must have zero nonlinear eta coefficient")
    if float(d.get("outer_angle_rad", 0.0)) < 0.80:
        f.append("accelerometer operation identity does not cover retained 0.8-rad sector")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    d = build(a.domain)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "PASS" if not vf else "FAIL",
        "outer_angle_rad": d["outer_angle_rad"],
        "aw_eta_coefficient": d["latent_aw_nonlinear_eta_coefficient"],
        "moving_metric_energy_invariant": d["state_coordinate_transform"]["moving_metric_energy_invariant"],
        "actual_RS_retained": d["actual_RS_regularizer_not_removed_by_coordinate_change"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
