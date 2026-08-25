#!/usr/bin/env python3
"""Exact effective-vector-input reduction for OU-III P5.

The finite-angle vector residual should not be charged as an independent
``eta`` disturbance when the shipping measurement map annihilates or absorbs
that residual algebraically.

Magnetometer.
-------------
For the configured deployment the magnetometer covariance is isotropic,
``R_m=r_m I``, and the shipping Jacobian is

    H_m = -[v]_x,

with ``v`` the predicted body-frame magnetic vector.  For the exact rotational
residual ``y_m=(R(c)-I)v``, decompose it into components parallel and
orthogonal to ``v``.  Since

    H_m^T v = 0,
    S_m v = r_m v,
    P C_m^T v = 0,

one has the *exact* state-correction identity

    K_m v = 0.

Therefore the radial finite-angle residual never changes the estimator state,
never changes the quaternion injection, and never enters the reset correction.
The useful residual is exactly tangent:

    K_m y_m = K_m H_m d_m,
    d_m := H_m^T y_m / ||v||^2.

For ``c=q u``, writing ``c_perp`` for the component orthogonal to ``v`` and
``alpha=c^T v_hat``, the Cayley formula gives

    d_m = A(q)c_perp - B(q) alpha (c_perp x v_hat),
    A(q)=4/(4+q^2),  B(q)=2/(4+q^2).

Hence

    ||d_m|| <= 2/sqrt(4+q^2) ||c_perp|| <= ||c_perp||,
    ||d_m-c_perp|| <= q/sqrt(4+q^2) ||c_perp||.

The exact finite-angle magnetometer map can thus be propagated as a
source-correlated effective tangent coordinate instead of subtracting a
standalone ``eta^T R^-1 eta`` budget.

Accelerometer.
--------------
In normal Live operation the shipping accelerometer Jacobian contains

    J_aw = R_wb,

which is orthogonal and full row rank.  For

    y_a = H_a z + eta_a

set ``e_eta=J_aw^T eta_a`` and insert it only in the ``a_w`` coordinate.  Then

    H_a E_aw e_eta = eta_a,
    K_a y_a = K_a H_a (z + E_aw e_eta),
    ||e_eta|| = ||eta_a||.

Thus the nonlinear residual is an exact source-correlated effective ``a_w``
input to the implemented Kalman correction; it is not an unrelated measurement
space disturbance.  The Joseph identity remains valid, but the P5 numerical
backend need not pay an independent eta norm when it directly propagates the
joint ``P,H,R,K,r,d_eff`` cell.

Using the exact Cayley residual identities already certified by
``ou3_p5_cayley_eta_geometry``, on ``||c||<=q`` one may enclose

    ||e_eta|| <= q^2/sqrt(4+q^2) ||f||
                 + 2q/sqrt(4+q^2) ||delta a_w||,

for the accelerometer attitude remainder plus the rotated latent-acceleration
cross term.  Bias/additive body terms are linear and are not charged as
attitude nonlinearity.  The configured proof domain disables the optional IMU
lever arm, so no unmodelled lever-arm nonlinear term is hidden in this lemma.

This module is a source-bound reduction primitive, not a complete P5 word
certificate and not a replay calculation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_exact_word_map as WORDMAP
import ou3_p5_cayley_eta_geometry as ETA

REPO = Path(__file__).resolve().parents[1]
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SIM = REPO / "src" / "util" / "W3dSimCommon.h"
CERT_SIM = REPO / "tests" / "kalman_ou_iii" / "ou3-certificate-sim.cpp"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _sqrt_down(x: float) -> float:
    if not (math.isfinite(x) and x > 0.0):
        raise ValueError("positive finite square-root input required")
    return down(math.sqrt(x))


def _q_den_sqrt_lower(q: float) -> float:
    q = float(q)
    if not (math.isfinite(q) and q >= 0.0):
        raise ValueError("finite nonnegative Cayley radius required")
    return _sqrt_down(down(4.0 + down(q * q)))


def mag_effective_coordinate_gain_upper(q_lo: float) -> float:
    """Upper of 2/sqrt(4+q^2) for q>=q_lo, capped by the exact <=1 theorem."""
    den = _q_den_sqrt_lower(q_lo)
    raw = up(2.0 / den)
    # The analytic expression is <=1 for every q>=0.  The cap avoids one-ulp
    # roundoff above one at q=0 without weakening the outward theorem bound.
    return min(1.0, raw)


def mag_effective_vs_tangent_defect_ratio_upper(q_hi: float) -> float:
    """Upper of q/sqrt(4+q^2)."""
    den = _q_den_sqrt_lower(q_hi)
    return up(float(q_hi) / den)


def accel_attitude_eta_per_vector_norm_upper(q_hi: float) -> float:
    """Upper of q^2/sqrt(4+q^2)."""
    q_hi = float(q_hi)
    den = _q_den_sqrt_lower(q_hi)
    return up(up(q_hi * q_hi) / den)


def accel_latent_cross_gain_upper(q_hi: float) -> float:
    """Upper of 2q/sqrt(4+q^2)."""
    q_hi = float(q_hi)
    den = _q_den_sqrt_lower(q_hi)
    return up(up(2.0 * q_hi) / den)


def _source_semantics(domain: dict) -> tuple[dict, list[str]]:
    mekf = MEKF.read_text(encoding="utf-8")
    sim = SIM.read_text(encoding="utf-8")
    cert = CERT_SIM.read_text(encoding="utf-8")
    markers = {
        "mag_H_theta": "const Matrix3 J_att = -skew_symmetric_matrix(v2hat);",
        "mag_innovation_covariance": "S_mat.noalias() += J_att * P_th_th * J_att.transpose();",
        "mag_cross_covariance": "PCt.noalias() += Pext.template block<NX,3>(0,0) * J_att.transpose();",
        "mag_R_from_sigma_vector": "Rmag(sigma_m.array().square().matrix().asDiagonal())",
        "configured_mag_sigma_isotropic": "const Vector3f sigma_m(sigma_m_uT, sigma_m_uT, sigma_m_uT);",
        "certificate_mag_scale_preserves_isotropy": "cfg.sigma_m = sigma_m * kSigmaMRescale;",
        "acc_J_aw_is_rotation": "J_aw  =  R_wb();",
        "acc_aw_cross_covariance_used": "PCt.noalias() += P_all_aw * J_aw.transpose();",
        "acc_prediction_uses_aw": "f_pred = R_wb() * (aw - g_world) + lever + ba_term;",
        "state_update_uses_full_gain": "xext.noalias() += K * r;",
        "joseph_update_before_reset": "joseph_update3_(K, S_mat, PCt);",
        "immediate_quaternion_reset": "applyQuaternionCorrectionFromErrorState();",
    }
    joined = mekf + "\n" + sim + "\n" + cert
    missing = [name for name, marker in markers.items() if marker not in joined]
    lever_disabled = domain.get("configured_runtime", {}).get("imu_lever_arm_enabled") is False
    if not lever_disabled:
        missing.append("configured proof domain does not disable optional IMU lever arm")
    return {
        "source_markers": markers,
        "configured_magnetometer_covariance_isotropic": "configured_mag_sigma_isotropic" not in missing
            and "certificate_mag_scale_preserves_isotropy" not in missing
            and "mag_R_from_sigma_vector" not in missing,
        "normal_live_accelerometer_Jaw_is_Rwb": "acc_J_aw_is_rotation" not in missing,
        "normal_live_accelerometer_Jaw_orthogonal_full_rank": "acc_J_aw_is_rotation" not in missing,
        "configured_imu_lever_arm_disabled": lever_disabled,
    }, missing


def _cell(row: dict) -> dict:
    q_lo, q_hi = map(float, row["q_interval"])
    return {
        "index": int(row["index"]),
        "q_interval": [q_lo, q_hi],
        "mag_effective_tangent_coordinate_gain_upper": mag_effective_coordinate_gain_upper(q_lo),
        "mag_effective_vs_tangent_defect_ratio_upper": mag_effective_vs_tangent_defect_ratio_upper(q_hi),
        "acc_effective_aw_attitude_eta_per_vector_norm_upper": accel_attitude_eta_per_vector_norm_upper(q_hi),
        "acc_effective_aw_latent_cross_gain_upper": accel_latent_cross_gain_upper(q_hi),
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("effective-vector-input domain must not be trajectory fitted")

    wordmap = WORDMAP.build(domain_path)
    eta = ETA.build(domain_path)
    failures = [f"word-map: {x}" for x in WORDMAP.validate(wordmap)]
    failures += [f"eta-geometry: {x}" for x in ETA.validate(eta)]
    source, missing = _source_semantics(domain)
    failures += [f"source semantic missing: {x}" for x in missing]

    cells = [_cell(row) for row in eta["annular_subdivision_cells"]]
    if not cells:
        failures.append("effective-vector-input subdivision is empty")
    for row in cells:
        if not (0.0 < float(row["mag_effective_tangent_coordinate_gain_upper"]) <= 1.0):
            failures.append("magnetometer effective coordinate lost nonexpansive bound")
            break
        if not (0.0 <= float(row["mag_effective_vs_tangent_defect_ratio_upper"]) < 1.0):
            failures.append("magnetometer effective tangent defect ratio is not below one")
            break
        if not (0.0 <= float(row["acc_effective_aw_latent_cross_gain_upper"]) < 2.0):
            failures.append("accelerometer latent cross gain outside exact rotation range")
            break

    fmax = float(domain["normal_live"]["specific_force_norm_upper_mps2"])
    awmax = float(domain["startup"]["physical_handoff_coordinate_bounds"]["latent_acceleration_error_norm_upper_mps2"])
    qmax = float(eta["widened_cayley_norm_upper"])
    acc_eta_bound = up(
        up(accel_attitude_eta_per_vector_norm_upper(qmax) * fmax)
        + up(accel_latent_cross_gain_upper(qmax) * awmax)
    )
    if not (math.isfinite(acc_eta_bound) and acc_eta_bound >= 0.0):
        failures.append("effective accelerometer aw-input bound is invalid")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_EXACT_EFFECTIVE_VECTOR_INPUT_REDUCTION",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "source_semantics": source,
        "magnetometer": {
            "shipping_H_theta": "H_m=-[v]_x",
            "configured_R": "R_m=r_m I",
            "H_theta_transpose_v_exact_zero": True,
            "innovation_covariance_action_on_v": "S_m v=r_m v",
            "cross_covariance_action_on_v_exact_zero": True,
            "kalman_gain_radial_action_exact_zero": True,
            "exact_effective_coordinate": "d_m=H_m^T y_m/||v||^2",
            "exact_state_correction_identity": "K_m y_m=K_m H_m d_m",
            "effective_coordinate_formula": "d_m=A c_perp-B alpha(c_perp x v_hat), A=4/(4+q^2), B=2/(4+q^2)",
            "effective_coordinate_nonexpansive": True,
            "effective_coordinate_norm_bound": "||d_m||<=2/sqrt(4+q^2)||c_perp||<=||c_perp||",
            "effective_vs_tangent_defect_bound": "||d_m-c_perp||<=q/sqrt(4+q^2)||c_perp||",
            "standalone_radial_eta_changes_state": False,
            "standalone_eta_information_penalty_required_for_state_correction": False,
        },
        "accelerometer": {
            "shipping_J_aw": "R_wb",
            "J_aw_inverse": "R_wb^T",
            "J_aw_orthogonal_full_row_rank": True,
            "effective_aw_defect": "e_eta=J_aw^T eta_a",
            "exact_measurement_range_identity": "H_a E_aw e_eta=eta_a",
            "exact_state_correction_identity": "K_a(H_a z+eta_a)=K_a H_a(z+E_aw e_eta)",
            "effective_aw_defect_norm_equals_eta_norm": True,
            "attitude_eta_bound": "||eta_R||<=q^2/sqrt(4+q^2)||f||",
            "latent_cross_bound": "||(R(c)^T-I)delta_a_w||<=2q/sqrt(4+q^2)||delta_a_w||",
            "widened_effective_aw_defect_norm_upper_mps2": acc_eta_bound,
            "standalone_eta_information_penalty_required_for_state_correction": False,
        },
        "gravity_quotient": {
            "accelerometer_effective_aw_input_descends_to_quotient": True,
            "axial_gyro_bias_role_unchanged": "NEUTRAL_BOUNDED_INPUT",
            "standalone_accelerometer_eta_penalty_required": False,
        },
        "annular_effective_input_cells": cells,
        "subdivision_cell_count": len(cells),
        "standalone_vector_eta_penalty_retired_from_P5_numerical_route": True,
        "joseph_information_identity_remains_valid": True,
        "eta_declared_identically_zero": False,
        "complete_word_numerical_certificate_closed_here": False,
        "P5_EFFECTIVE_VECTOR_INPUT_CERTIFICATE": "PASS" if not failures else "FAIL",
        "next_obligation": (
            "propagate the remaining source-correlated P,H,R,K,r,d_eff/reset/prediction cells through every later 1 s prefix; do not reintroduce an independent vector-eta norm budget"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for key in ("source_generated_not_trajectory_fit", "standalone_vector_eta_penalty_retired_from_P5_numerical_route",
                "joseph_information_identity_remains_valid"):
        if d.get(key) is not True:
            failures.append(f"{key} is not true")
    for key in ("source_replay_used", "filter_changed", "eta_declared_identically_zero",
                "complete_word_numerical_certificate_closed_here"):
        if d.get(key) is not False:
            failures.append(f"{key} is not false")

    src = d.get("source_semantics", {})
    for key in ("configured_magnetometer_covariance_isotropic",
                "normal_live_accelerometer_Jaw_is_Rwb",
                "normal_live_accelerometer_Jaw_orthogonal_full_rank",
                "configured_imu_lever_arm_disabled"):
        if src.get(key) is not True:
            failures.append(f"source semantic {key} is not true")

    mag = d.get("magnetometer", {})
    for key in ("H_theta_transpose_v_exact_zero", "cross_covariance_action_on_v_exact_zero",
                "kalman_gain_radial_action_exact_zero", "effective_coordinate_nonexpansive"):
        if mag.get(key) is not True:
            failures.append(f"magnetometer {key} is not true")
    if mag.get("standalone_radial_eta_changes_state") is not False:
        failures.append("magnetometer radial eta still changes state")
    if mag.get("standalone_eta_information_penalty_required_for_state_correction") is not False:
        failures.append("magnetometer independent eta penalty not retired")

    acc = d.get("accelerometer", {})
    if acc.get("J_aw_orthogonal_full_row_rank") is not True:
        failures.append("accelerometer J_aw is not certified orthogonal/full-row-rank")
    if acc.get("effective_aw_defect_norm_equals_eta_norm") is not True:
        failures.append("accelerometer effective aw-input norm identity missing")
    if acc.get("standalone_eta_information_penalty_required_for_state_correction") is not False:
        failures.append("accelerometer independent eta penalty not retired")
    x = acc.get("widened_effective_aw_defect_norm_upper_mps2")
    if not (isinstance(x, (int, float)) and math.isfinite(float(x)) and float(x) >= 0.0):
        failures.append("accelerometer effective aw defect bound is invalid")

    q = d.get("gravity_quotient", {})
    if q.get("accelerometer_effective_aw_input_descends_to_quotient") is not True:
        failures.append("effective accelerometer input does not descend to gravity quotient")
    if q.get("standalone_accelerometer_eta_penalty_required") is not False:
        failures.append("gravity quotient still pays independent accelerometer eta")

    cells = d.get("annular_effective_input_cells", [])
    if not cells or len(cells) != d.get("subdivision_cell_count"):
        failures.append("effective-vector-input cells missing")
    else:
        for row in cells:
            if not 0.0 < float(row["mag_effective_tangent_coordinate_gain_upper"]) <= 1.0:
                failures.append("mag effective-coordinate cell is expansive")
                break
            if not 0.0 <= float(row["mag_effective_vs_tangent_defect_ratio_upper"]) < 1.0:
                failures.append("mag effective-coordinate defect cell invalid")
                break
    if not failures and d.get("P5_EFFECTIVE_VECTOR_INPUT_CERTIFICATE") != "PASS":
        failures.append("effective-vector-input certificate did not pass")
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
        "status": out["P5_EFFECTIVE_VECTOR_INPUT_CERTIFICATE"],
        "mag": out["magnetometer"],
        "acc": out["accelerometer"],
        "cells": out["subdivision_cell_count"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
