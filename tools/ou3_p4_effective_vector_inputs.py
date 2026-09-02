#!/usr/bin/env python3
"""Exact effective-vector-input identities for the broad OU-III P4 sector.

The 0.8-rad P4 route must not charge every finite-angle vector residual as an
independent adverse ``eta`` term.  Two exact range/nullspace identities from the
shipping measurement model remove that artificial penalty.

Magnetometer
------------
For the configured isotropic magnetic covariance and source Jacobian

    H_m = -[v]x,

any residual component parallel to v lies in the innovation direction on which
K_m is exactly zero.  The useful part of y=(R(c)-I)v can therefore be written

    K_m y = K_m H_m d_m,
    d_m = H_m' y / ||v||^2.

For c=2 tan(theta/2) u, q=||c||, vhat=v/||v||,

    d_m = 4/(4+q^2) c_perp
          - 2 alpha/(4+q^2) (c_perp x vhat),

and hence

    ||d_m|| <= ||c_perp||,
    ||d_m-c_perp|| <= q/sqrt(4+q^2) ||c_perp||.

Accelerometer
-------------
With the proof-scope lever arm disabled, the source uses

    J_aw = R_wb,
    J_aw' J_aw = I.

Thus any nonlinear measurement-space remainder eta_a is exactly in the range of
that state block.  Defining e_eta=J_aw' eta_a gives

    H_a E_aw e_eta = eta_a,
    K_a(H_a z + eta_a) = K_a H_a(z + E_aw e_eta),
    ||e_eta|| = ||eta_a||.

On the Cayley sector the exact rotational/latent-acceleration part obeys

    ||e_eta|| <= q^2/sqrt(4+q^2) * ||f||
                 + 2q/sqrt(4+q^2) * ||delta a_w||.

These identities do not themselves prove a complete H/A word.  They are the
source-faithful input reduction consumed by the retained source-correlated
signed-Joseph word composition.  In particular, this producer never invents a
condition-number conversion between the P3 information metric and P4.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_cayley_sector_certificate as CAYLEY
import ou3_p4_vector_remainder_sector as REMAINDER
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
SCHEMA = 1
DEFAULT_OUTER_ANGLE_RAD = 0.80


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def build(domain_path: Path = DEFAULT_DOMAIN,
          outer_angle_rad: float = DEFAULT_OUTER_ANGLE_RAD) -> dict:
    path = Path(domain_path).resolve()
    domain = json.loads(path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("effective-vector P4 domain must not be trajectory fitted")

    runtime = domain.get("configured_runtime", {})
    if runtime.get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("effective accelerometer range identity currently requires lever arm disabled")

    source = MEKF.read_text(encoding="utf-8")
    source_markers = (
        "const Matrix3 J_att = -skew_symmetric_matrix(f_cog_b);",
        "const Matrix3 J_aw  =  R_wb();",
        "const Matrix3 J_att = -skew_symmetric_matrix(v2hat);",
    )
    missing = [m for m in source_markers if m not in source]
    if missing:
        raise RuntimeError(f"shipping vector Jacobian semantics changed: {missing}")

    cayley = CAYLEY.build(path, outer_angle_rad)
    remainder = REMAINDER.build(path, outer_angle_rad)
    vector = VECTOR.build()
    prereq = [f"cayley: {x}" for x in CAYLEY.validate(cayley)]
    prereq += [f"remainder: {x}" for x in REMAINDER.validate(remainder)]
    prereq += [f"vector: {x}" for x in VECTOR.validate(vector)]

    q = float(cayley["cayley_radius_upper"])
    if not (math.isfinite(q) and q >= 0.0):
        raise RuntimeError("invalid Cayley radius")
    den_lo = down(math.sqrt(4.0 + q * q))
    den_hi = up(math.sqrt(4.0 + q * q))
    if den_lo <= 0.0:
        raise RuntimeError("effective-vector denominator lost positivity")

    # The first magnetic factor 2/sqrt(4+q^2) is maximized at q=0, so the
    # source-uniform nonexpansive upper bound is exactly one.  The defect factor
    # is monotone increasing and is evaluated at the sector boundary.
    mag_nonexpansive_factor_upper = 1.0
    mag_tangent_defect_factor_upper = up(q / den_lo)

    # Both accelerometer coefficients are monotone for q>=0 on this chart.
    acc_force_remainder_factor_upper = up((q * q) / den_lo)
    acc_aw_rotation_factor_upper = up((2.0 * q) / den_lo)

    live = domain["normal_live"]
    f_upper = float(live["specific_force_norm_upper_mps2"])
    aw_upper = float(
        domain["startup"]["physical_handoff_coordinate_bounds"]
        ["latent_acceleration_error_norm_upper_mps2"]
    )
    effective_aw_input_norm_upper = up(
        acc_force_remainder_factor_upper * f_upper
        + acc_aw_rotation_factor_upper * aw_upper
    )

    pass_ = bool(
        not prereq
        and cayley.get("declared_filter_entrance_covered") is True
        and remainder.get("accelerometer_bias_cancels_exactly_from_eta") is True
        and q < 1.0
        and 0.0 <= mag_tangent_defect_factor_upper < 1.0
        and mag_nonexpansive_factor_upper == 1.0
        and math.isfinite(effective_aw_input_norm_upper)
    )

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P4_EXACT_EFFECTIVE_VECTOR_INPUT_SECTOR",
        "source_generated_not_trajectory_fit": True,
        "trajectory_replay_used": False,
        "filter_changed": False,
        "outer_angle_rad": float(outer_angle_rad),
        "outer_angle_deg_display": float(outer_angle_rad) * 180.0 / math.pi,
        "cayley_radius_upper": q,
        "shipping_accelerometer_J_aw_is_R_wb": True,
        "shipping_magnetometer_J_att_is_minus_skew_predicted_vector": True,
        "lever_arm_disabled_in_proof_scope": True,
        "configured_vector_covariance_scope": "ISOTROPIC_THREE_AXIS",
        "mag_radial_residual_gain_null_exact": True,
        "mag_effective_coordinate_identity_exact": True,
        "mag_effective_coordinate_nonexpansive_factor_upper": mag_nonexpansive_factor_upper,
        "mag_effective_coordinate_tangent_defect_factor_upper": mag_tangent_defect_factor_upper,
        "acc_eta_in_aw_measurement_range_exact": True,
        "acc_effective_aw_input_isometry_exact": True,
        "acc_force_remainder_factor_upper": acc_force_remainder_factor_upper,
        "acc_aw_rotation_factor_upper": acc_aw_rotation_factor_upper,
        "specific_force_norm_upper_mps2": f_upper,
        "latent_acceleration_error_norm_upper_mps2": aw_upper,
        "acc_effective_aw_input_norm_upper_mps2": effective_aw_input_norm_upper,
        "accelerometer_bias_standalone_nonlinear_penalty": 0.0,
        "standalone_vector_eta_penalty_retired_from_active_word_route": True,
        "condition_number_conversion_inserted_between_P3_and_P4": False,
        "complete_HA_signed_Joseph_word_established_here": False,
        "P4_USABLE_CERTIFICATE_PROMOTED": False,
        "prerequisite_failures": prereq,
        "pass": pass_,
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("prerequisite_failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for key in (
        "source_generated_not_trajectory_fit",
        "shipping_accelerometer_J_aw_is_R_wb",
        "shipping_magnetometer_J_att_is_minus_skew_predicted_vector",
        "lever_arm_disabled_in_proof_scope",
        "mag_radial_residual_gain_null_exact",
        "mag_effective_coordinate_identity_exact",
        "acc_eta_in_aw_measurement_range_exact",
        "acc_effective_aw_input_isometry_exact",
        "standalone_vector_eta_penalty_retired_from_active_word_route",
        "pass",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "trajectory_replay_used",
        "filter_changed",
        "condition_number_conversion_inserted_between_P3_and_P4",
        "complete_HA_signed_Joseph_word_established_here",
        "P4_USABLE_CERTIFICATE_PROMOTED",
    ):
        if d.get(key) is not False:
            f.append(f"{key} is not false")
    if d.get("accelerometer_bias_standalone_nonlinear_penalty") != 0.0:
        f.append("accelerometer bias must not create a standalone nonlinear eta penalty")
    for key in (
        "mag_effective_coordinate_tangent_defect_factor_upper",
        "acc_force_remainder_factor_upper",
        "acc_aw_rotation_factor_upper",
        "acc_effective_aw_input_norm_upper_mps2",
    ):
        x = d.get(key)
        if not isinstance(x, (int, float)) or not math.isfinite(float(x)) or float(x) < 0.0:
            f.append(f"{key} is invalid")
    if float(d.get("mag_effective_coordinate_tangent_defect_factor_upper", math.inf)) >= 1.0:
        f.append("magnetic effective tangent defect is not contractive on the declared sector")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--outer-angle-rad", type=float, default=DEFAULT_OUTER_ANGLE_RAD)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain, args.outer_angle_rad)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "qualification": d["qualification"],
        "q_upper": d["cayley_radius_upper"],
        "mag_tangent_defect_upper": d["mag_effective_coordinate_tangent_defect_factor_upper"],
        "acc_effective_aw_input_norm_upper_mps2": d["acc_effective_aw_input_norm_upper_mps2"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
