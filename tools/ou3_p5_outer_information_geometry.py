#!/usr/bin/env python3
"""Finite-angle source-information geometry prerequisite for OU-III P5.

This producer closes two analytical pieces needed by the gauged outer-H word
certificate without reverting to the false isotropic V_R sector.

1. Exact vector-pair geometry in Cayley coordinates.  For c=q u and any source
   vector v,

       ||(R(c)^T-I)v||^2
         = 4 q^2/(4+q^2) * (||v||^2-(u^T v)^2).

   Therefore a non-collinear accepted accelerometer/magnetometer pair satisfying
   the declared source packet bounds has a strictly positive finite-angle
   residual-information lower bound over each q<=q_* handoff node.

2. The exact Joseph/Kalman tangent information identity.  If y=H z+eta is the
   signed nonlinear measurement mismatch used by the local correction,

       z^T P^-1 z - (z-Ky)^T (P+)^-1 (z-Ky)
          = y^T S^-1 y - eta^T R^-1 eta.

   This is an identity, not a condition-number bound.  It explains why the
   outer proof must keep residual, covariance, gain and nonlinear defect from
   one source tuple together.  The complete word still has to transport the
   exact deployed Cayley/quaternion/reset defect and the non-attitude state
   through all source prefixes; this module deliberately does not promote that
   final step.

The output supplies separate positive constants for the q~=0.27 and q~=0.60
H nodes and pins them to the actual goLive attitude covariance.  It uses no
replay and changes no filter/tuning value.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p5_go_live_covariance_stage as GOLIVE
import ou3_p5_heading_handoff_contract as HEADING
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _angular_factor(s: float) -> float:
    # Stable evaluation of 1-sqrt(1-s^2).
    root = up(math.sqrt(max(0.0, down(1.0 - up(s*s)))))
    return down(down(s*s) / up(1.0 + root))


def _node(name: str, qmax: float, mu_pair: float, qinfo_max: float) -> dict:
    finite_factor = down(4.0 / up(4.0 + up(qmax*qmax)))
    residual_per_q2 = down(finite_factor * mu_pair)
    alpha_vs_initial_info = down(residual_per_q2 / qinfo_max)
    return {
        "node": name,
        "cayley_norm_upper": qmax,
        "exact_cayley_residual_factor_lower": finite_factor,
        "exact_pair_residual_information_per_cayley_norm_sq_lower": residual_per_q2,
        "goLive_attitude_information_lambda_max_upper": qinfo_max,
        "exact_pair_residual_information_vs_goLive_attitude_metric_lower": alpha_vs_initial_info,
        "strict": alpha_vs_initial_info > 0.0,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("outer information geometry domain must not be trajectory fitted")

    heading = HEADING.build(domain_path)
    hf = HEADING.validate(heading)
    stage = GOLIVE.build(domain_path)
    sf = GOLIVE.validate(stage)
    vector = VECTOR.build()
    vf = VECTOR.validate(vector)
    failures = [f"heading: {x}" for x in hf] + [f"goLive: {x}" for x in sf] + [f"vector: {x}" for x in vf]

    live = domain["normal_live"]
    fmin = float(live["specific_force_norm_lower_mps2"])
    mmin = float(live["magnetic_vector_norm_lower_uT"])
    sine = float(live["vector_sine_separation_lower"])
    if not (fmin > 0.0 and mmin > 0.0 and 0.0 < sine < 1.0):
        failures.append("declared full-heading vector packet geometry is invalid")

    cm = vector["configured_measurement_bounds"]
    ra = float(cm["acc_measurement_variance_upper"])
    rm = float(cm["mag_measurement_variance_upper"])
    if not (ra > 0.0 and rm > 0.0):
        failures.append("configured vector measurement variance upper is invalid")

    af = down(down(fmin*fmin) / ra)
    am = down(down(mmin*mmin) / rm)
    angular = _angular_factor(sine)
    mu_pair = down(min(af, am) * angular)
    if not (math.isfinite(mu_pair) and mu_pair > 0.0):
        failures.append("finite-angle vector-pair information lower lost positivity")

    seed = stage["goLive_H_covariance_seed"]["attitude_covariance_seed"]
    tilt_var = float(seed["tilt_variance"])
    yaw_var = float(seed["gauged_yaw_variance"])
    pmin = min(tilt_var, yaw_var)
    pmax = max(tilt_var, yaw_var)
    qinfo_max = up(1.0 / down(pmin))
    qinfo_min = down(1.0 / up(pmax))
    if not (0.0 < qinfo_min <= qinfo_max < math.inf):
        failures.append("goLive attitude information eigenvalue enclosure invalid")

    q_normal = float(heading["gauged_quality_handoff"]["full_attitude_cayley_norm_upper"])
    q_timeout = float(heading["gauged_timeout_subbranch"]["full_attitude_cayley_norm_upper"])
    nodes = {
        "normal_gauged": _node("normal_gauged", q_normal, mu_pair, qinfo_max),
        "timeout_gauged": _node("timeout_gauged", q_timeout, mu_pair, qinfo_max),
    }
    for name, row in nodes.items():
        if row["strict"] is not True:
            failures.append(f"{name}: exact finite-angle residual information is not strict")

    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_EXACT_FINITE_ANGLE_SOURCE_INFORMATION_GEOMETRY",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "actual_goLive_covariance_tuple_used": True,
        "joint_source_tuple_required": True,
        "independent_gain_extrema_product_used": False,
        "exact_cayley_vector_residual_identity": (
            "||(R(c)^T-I)v||^2 = 4||c||^2/(4+||c||^2) * (||v||^2-(u^T v)^2)"
        ),
        "exact_joseph_tangent_information_identity": (
            "z^T P^-1 z-(z-Ky)^T(Pplus)^-1(z-Ky)=y^T S^-1 y-eta^T R^-1 eta"
        ),
        "identity_requires_source_correlated_H_P_R_K_S": True,
        "packet_geometry": {
            "specific_force_norm_lower_mps2": fmin,
            "magnetic_vector_norm_lower_uT": mmin,
            "vector_sine_separation_lower": sine,
            "acc_measurement_variance_upper": ra,
            "mag_measurement_variance_upper": rm,
            "a_f_lower": af,
            "a_m_lower": am,
            "angular_factor_lower": angular,
            "linear_pair_information_mu_lower": mu_pair,
        },
        "goLive_attitude_information": {
            "tilt_variance": tilt_var,
            "gauged_yaw_variance": yaw_var,
            "lambda_min_lower": qinfo_min,
            "lambda_max_upper": qinfo_max,
        },
        "nodes": nodes,
        "P5_FINITE_ANGLE_INFORMATION_GEOMETRY_CERTIFICATE": "PASS" if not failures else "FAIL",
        "P5_GAUGED_OUTER_CAYLEY_INFORMATION_WORD_SECTOR": "NOT_ESTABLISHED",
        "remaining_word_term": (
            "validated source-correlated transport of the exact nonlinear eta term, deployed Cayley/quaternion injection, "
            "Joseph covariance/reset metric change, S-to-attitude prefix, and all non-attitude coordinates over the complete 1 s source word"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for flag in ("source_generated_not_trajectory_fit", "actual_goLive_covariance_tuple_used",
                 "joint_source_tuple_required", "identity_requires_source_correlated_H_P_R_K_S"):
        if d.get(flag) is not True:
            failures.append(f"{flag} is not true")
    if d.get("source_replay_used") is not False:
        failures.append("outer information geometry uses replay")
    if d.get("filter_changed") is not False:
        failures.append("outer information geometry changes filter")
    if d.get("independent_gain_extrema_product_used") is not False:
        failures.append("independent gain extrema are used")
    for name in ("normal_gauged", "timeout_gauged"):
        row = d.get("nodes", {}).get(name, {})
        if row.get("strict") is not True:
            failures.append(f"{name}: finite-angle information lower is not strict")
        a = row.get("exact_pair_residual_information_vs_goLive_attitude_metric_lower")
        if not (isinstance(a, (int, float)) and math.isfinite(float(a)) and float(a) > 0.0):
            failures.append(f"{name}: finite-angle metric coefficient is not finite positive")
    if not failures and d.get("P5_FINITE_ANGLE_INFORMATION_GEOMETRY_CERTIFICATE") != "PASS":
        failures.append("finite-angle information geometry did not pass")
    if d.get("P5_GAUGED_OUTER_CAYLEY_INFORMATION_WORD_SECTOR") != "NOT_ESTABLISHED":
        failures.append("complete gauged outer word sector promoted before exact word transport")
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
        "status": out["P5_FINITE_ANGLE_INFORMATION_GEOMETRY_CERTIFICATE"],
        "packet": out["packet_geometry"],
        "nodes": out["nodes"],
        "remaining_word_term": out["remaining_word_term"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
