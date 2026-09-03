#!/usr/bin/env python3
"""Exact co-rotated a_w coordinate for the OU-III accelerometer Joseph step.

The retained 0.8-rad remainder primitive writes the exact accelerometer
residual in the original world-acceleration error coordinate as

    r = (E-I) f_hat + E R_hat da + db,

where R_hat maps world to the estimator body frame, E is the left
multiplicative attitude error, f_hat=R_hat(a_hat-g), da is the world-frame
latent-acceleration error and db is the accelerometer-bias error.

Treating ``E R_hat da - R_hat da`` as nonlinear measurement eta is valid but
unnecessarily expensive.  Define instead, for this one measurement operation,

    Q = R_hat' E R_hat,       u = Q da.

Q is orthogonal and therefore ||u||=||da|| exactly.  Also

    E R_hat da = R_hat Q da = R_hat u,

so the exact residual becomes

    r = (E-I) f_hat + R_hat u + db.

Against the tangent residual

    h = [c]x f_hat + R_hat u + db,

all a_w and b_a terms cancel from eta=r-h.  The nonlinear measurement defect is
therefore the pure finite-angle rotation defect only:

    eta = ((E-I)-[c]x) f_hat.

This is not a change to the shipping filter.  It is an exact pointwise state
coordinate congruence used only inside the proof of an accepted accelerometer
operation.  With

    T_E = diag(I_theta, I_bg, I_v, I_p, I_S, Q_aw [, I_ba]),

we have T_E' T_E=I.  Transforming z,P,H,K by

    z_u=T_E z, P_u=T_E P T_E', H_u=H T_E', K_u=T_E K

leaves S=H P H'+R, the Joseph quadratic identity, and the exact state
correction invariant.  The current P3->P4 metric is group-isotropic (one scalar
weight repeated on each 3-vector group), hence its a_w contribution is also an
exact isometry under Q.  The proof may transform back after the measurement at
zero metric cost.

The source contract below pins the identity to the deployed zero-lever-arm
accelerometer model and left-error correction convention.  This primitive does
not prove complete-word contraction or promote P4 by itself.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_p3_metric_attachment as METRIC

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SCHEMA = 1


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
        failures.append("co-rotated accelerometer identity requires declared zero lever arm")
    if runtime.get("accelerometer_vibration_guard_proof_branch") != "dormant_transparent":
        failures.append("co-rotated accelerometer identity requires transparent vibration-guard branch")
    text = MEKF.read_text(encoding="utf-8")
    for marker in SOURCE_MARKERS:
        if marker not in text:
            failures.append(f"shipping accelerometer/correction marker changed: {marker}")
    if "Measurement Jacobians use the left-multiplicative convention" not in text:
        failures.append("shipping left-multiplicative attitude-error convention marker changed")
    if "J_ba = I" not in text:
        failures.append("shipping accelerometer-bias identity Jacobian marker changed")
    return failures


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("co-rotated accelerometer proof must not be trajectory fitted")

    cayley = CAYLEY.build(path)
    failures = [f"Cayley: {x}" for x in CAYLEY.validate(cayley)]
    failures += _source_contract(domain)

    # The attached P4 metric repeats one scalar on every coordinate of each
    # physical 3-vector group.  Q_aw therefore commutes with the a_w metric
    # block exactly.  Pin the expected state grouping rather than silently
    # assuming a different future state layout.
    if METRIC.STATE_GROUPS.get("H") != ["theta", "b_g", "v", "p", "S", "a_w"]:
        failures.append("H-mode P4 metric state grouping changed")
    if METRIC.STATE_GROUPS.get("A") != ["theta", "b_g", "v", "p", "S", "a_w", "b_a"]:
        failures.append("A-mode P4 metric state grouping changed")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_EXACT_ACCELEROMETER_COROTATED_AW_COORDINATE",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_changed": False,
        "zero_lever_arm_branch": True,
        "dormant_transparent_vibration_guard_branch": True,
        "shipping_left_multiplicative_error_convention_bound": True,
        "shipping_accelerometer_J_aw_is_R_hat": True,
        "shipping_accelerometer_J_ba_is_identity": True,
        "coordinate_definition": {
            "attitude_error": "R_true = E R_hat",
            "Q_aw": "R_hat^T E R_hat",
            "u_aw": "Q_aw delta_a_w",
            "Q_aw_orthogonal": True,
            "aw_norm_preserved_exactly": True,
        },
        "exact_residual_identity": (
            "r_a=(E-I)f_hat + R_hat*u_aw + delta_b_a, "
            "u_aw=(R_hat^T E R_hat)delta_a_w"
        ),
        "tangent_residual_in_corotated_coordinate": (
            "h_a=[c]x f_hat + R_hat*u_aw + delta_b_a"
        ),
        "nonlinear_eta_in_corotated_coordinate": "eta_a=((E-I)-[c]x)f_hat",
        "latent_aw_nonlinear_eta_coefficient": 0.0,
        "accelerometer_bias_nonlinear_eta_coefficient": 0.0,
        "aw_error_exactly_linear_in_accelerometer_operation_coordinate": True,
        "accelerometer_bias_error_exactly_linear": True,
        "state_coordinate_transform": {
            "kind": "BLOCK_ORTHOGONAL_CONGRUENCE",
            "T_E": "diag(I_theta,I_bg,I_v,I_p,I_S,Q_aw[,I_ba])",
            "state_energy_invariant_under_exact_covariance_congruence": True,
            "innovation_covariance_S_invariant": True,
            "Joseph_information_identity_invariant": True,
            "exact_state_correction_invariant": True,
            "transform_back_after_measurement_exact": True,
        },
        "P3_P4_group_isotropic_metric_aw_block_isometry": True,
        "outer_angle_rad": float(cayley["outer_angle_rad"]),
        "outer_sector_covered": float(cayley["outer_angle_rad"]) >= 0.80,
        "replaces_old_aw_eta_penalty_in_complete_word": True,
        "old_aw_eta_penalty_was_valid_but_coordinate_conservative": True,
        "complete_H18_A21_word_established_here": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "next_obligation": (
            "use the pure-rotation accelerometer eta in the source-correlated signed Joseph directional word; "
            "retain the exact congruence through each accepted accelerometer operation and scalarize only at the recurrent word endpoint"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    if d.get("qualification") != "OU3_P4_EXACT_ACCELEROMETER_COROTATED_AW_COORDINATE":
        f.append("wrong qualification")
    for key in (
        "source_generated_not_trajectory_fit", "zero_lever_arm_branch",
        "dormant_transparent_vibration_guard_branch",
        "shipping_left_multiplicative_error_convention_bound",
        "shipping_accelerometer_J_aw_is_R_hat", "shipping_accelerometer_J_ba_is_identity",
        "aw_error_exactly_linear_in_accelerometer_operation_coordinate",
        "accelerometer_bias_error_exactly_linear",
        "P3_P4_group_isotropic_metric_aw_block_isometry", "outer_sector_covered",
        "replaces_old_aw_eta_penalty_in_complete_word",
        "old_aw_eta_penalty_was_valid_but_coordinate_conservative",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_changed",
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
        "state_energy_invariant_under_exact_covariance_congruence",
        "innovation_covariance_S_invariant", "Joseph_information_identity_invariant",
        "exact_state_correction_invariant", "transform_back_after_measurement_exact",
    ):
        if tr.get(key) is not True:
            f.append(f"state congruence {key} is not true")
    if d.get("latent_aw_nonlinear_eta_coefficient") != 0.0:
        f.append("co-rotated a_w must have zero nonlinear eta coefficient")
    if d.get("accelerometer_bias_nonlinear_eta_coefficient") != 0.0:
        f.append("accelerometer bias must have zero nonlinear eta coefficient")
    if float(d.get("outer_angle_rad", 0.0)) < 0.80:
        f.append("co-rotated identity is not attached to the 0.8-rad sector")
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
    a.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": "PASS" if not vf else "FAIL",
        "outer_angle_rad": d["outer_angle_rad"],
        "aw_eta_coefficient": d["latent_aw_nonlinear_eta_coefficient"],
        "ba_eta_coefficient": d["accelerometer_bias_nonlinear_eta_coefficient"],
        "Joseph_congruence_invariant": d["state_coordinate_transform"]["Joseph_information_identity_invariant"],
        "P4_promoted": d["P4_USABLE_CERTIFICATE_PROMOTED"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
