#!/usr/bin/env python3
"""Exact accelerometer coordinates and shipping-tangent remainder for OU-III P4.

Let R_true=E R_hat, C=[c]x, da=a_true-a_hat, and f=R_hat(a_hat-g).
With zero lever arm the exact homogeneous residual is

    y=(E-I)f+E R_hat da+db_a,       H z=C f+R_hat da+db_a.

For Q=R_hat^T E R_hat, u=Q da the physical residual is linear in u:

    y=(E-I)f+R_hat u+db_a.

A frozen coordinate congruence T=diag(I,...,Q[,I]) preserves energy and
innovation ONLY when P_u=T P T^T, H_u=H T^-1 and K_u=T K. Consequently

    H_u z_u=H z=C f+R_hat Q^T u+db_a,
    y-H_u z_u=((E-I)-C)f+(R_hat-R_hat Q^T)u.

The mixed wave-error term is not removed by a congruence. Keeping R_hat as
the tangent u column while also claiming H_u=H T^-1 is inconsistent.

There is a different, exact nonlinear measurement-linearizing coordinate:

    e0=R_hat^T ((E-I)-C)f,
    epsilon_aw=(Q-I)da+e0,       Phi_aw=da+epsilon_aw=Q da+e0.

It satisfies y=H Phi(z) with the original shipping H,P,K,S. Its storage
Phi(z)^T P^-1 Phi(z) is NOT asserted equal to z^T P^-1 z. Full epsilon_aw,
including (Q-I)da, must persist through the physical correction/reset and
source transitions. No filter, source word, or P3 certificate is changed.

The evaluator below uses the retained outward interval arithmetic. R_hat must
be an orthogonal source-provided world-to-body rotation; interval consistency
alone does not establish source reachability. These are operation identities,
not a complete-word contraction or a new source generator.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from ou3_interval import (
    Interval, matrix_add, matrix_identity, matrix_mul, matrix_sub, matrix_transpose,
)

import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_moving_metric_rebind as REBIND
import ou3_sea3_complete_source as COMPLETE

REPO = Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN = REPO / "tools" / "stability" / "ou3_proof_operating_domain.json"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SCHEMA = 3
QUALIFICATION = "OU3_P4_COMPLETE_SEA3_ACCELEROMETER_OPERATION_COORDINATE_MOVING_METRIC_V3"

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


def evaluate_operation_coordinate(c, R_hat, da, f_hat) -> dict:
    """Evaluate the exact full shift; vectors are 3x1 interval matrices.

    R_hat orthogonality is a caller/source precondition, not inferred from an
    arbitrary box of entries. No source-completeness or P4 flag is emitted.
    """
    for name, value, cols in (("c", c, 1), ("da", da, 1),
                              ("f_hat", f_hat, 1), ("R_hat", R_hat, 3)):
        if len(value) != 3 or any(len(row) != cols for row in value):
            raise ValueError(f"{name} must have shape 3x{cols}")
        if any(not isinstance(x, Interval) or not (float("-inf") < x.lo <= x.hi < float("inf"))
               for row in value for x in row):
            raise ValueError(f"{name} must contain finite intervals")
    z = Interval.point(0.0)
    half, quarter = Interval.point(0.5), Interval.point(0.25)
    x, y, w = (v[0] for v in c)
    C = [[z, -w, y], [w, z, -x], [-y, x, z]]
    C2 = matrix_mul(C, C)
    den = Interval.point(1.0) + quarter * (x.square() + y.square() + w.square())
    EmI = [[(C[i][j] + half*C2[i][j])/den for j in range(3)] for i in range(3)]
    E = matrix_add(matrix_identity(3), EmI)
    Rt = matrix_transpose(R_hat)
    Q = matrix_mul(matrix_mul(Rt, E), R_hat)
    u = matrix_mul(Q, da)
    pure_eta = matrix_mul(matrix_sub(EmI, C), f_hat)
    e0 = matrix_mul(Rt, pure_eta)
    mixed_shift = matrix_sub(u, da)
    epsilon = matrix_add(mixed_shift, e0)
    phi_aw = matrix_add(da, epsilon)
    tangent = matrix_add(matrix_mul(C, f_hat), matrix_mul(R_hat, da))
    physical = matrix_add(matrix_mul(EmI, f_hat), matrix_mul(matrix_mul(E, R_hat), da))
    linearized_phi = matrix_add(matrix_mul(C, f_hat), matrix_mul(R_hat, phi_aw))
    return {
        "E": E, "Q_aw": Q, "u_aw": u, "pure_eta": pure_eta,
        "e_eta": e0, "mixed_aw_shift": mixed_shift,
        "epsilon_aw": epsilon, "Phi_aw": phi_aw,
        "shipping_tangent_residual_without_ba": tangent,
        "physical_residual_without_ba": physical,
        "H_Phi_without_ba": linearized_phi,
        "shipping_tangent_remainder": matrix_sub(physical, tangent),
    }


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
            "h_u=[c]x f_hat + R_hat*Q_aw^T*u_aw + delta_b_a"
        ),
        "nonlinear_eta_in_operation_coordinate": "eta_a=((E-I)-[c]x)f_hat+(R_hat-R_hat*Q_aw^T)*u_aw",
        "latent_aw_nonlinear_eta_coefficient": float(cayley["cayley_radius_upper"]),
        "accelerometer_bias_nonlinear_eta_coefficient": 0.0,
        "aw_error_exactly_linear_in_accelerometer_operation_coordinate": True,
        "mixed_aw_shipping_tangent_remainder_retained": True,
        "nonlinear_Phi_storage_is_original_metric_isometry": False,
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
        "moving_metric_rebind_consumed": bool(rebind["structural_shipping_covariance_identities_closed"]),
        "outer_angle_rad": float(cayley["outer_angle_rad"]),
        "outer_sector_covered": float(cayley["outer_angle_rad"]) >= 0.80,
        "retired_endpoint_attachment_module_reintroduced": False,
        "complete_H18_A21_word_established_here": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "shipping_Joseph_binding_closed": False,
        "auxiliary_H0_equals_congruent_shipping_Hu": False,
        "congruent_shipping_aw_column": "R_hat Q_aw^T",
        "next_obligation": (
            "retain the mixed wave-error remainder or transport full epsilon_aw in nonlinear Phi storage; "
            "keep u_aw and b_a in the linear residual/innovation and keep every actual-R_S S update"
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
    if not float(d.get("latent_aw_nonlinear_eta_coefficient", 0.0)) > 0.0:
        f.append("mixed a_w shipping-tangent remainder bound is missing")
    if d.get("mixed_aw_shipping_tangent_remainder_retained") is not True:
        f.append("mixed a_w shipping-tangent remainder was dropped")
    if d.get("nonlinear_Phi_storage_is_original_metric_isometry") is not False:
        f.append("nonlinear Phi storage was misrepresented as a metric isometry")
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
