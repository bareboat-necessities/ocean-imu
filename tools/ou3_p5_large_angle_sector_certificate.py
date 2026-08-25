#!/usr/bin/env python3
"""Validated source-correlated audit of the OU-III P5 large-angle V_R sector.

The paper proposed, on a gauged outer-H node,

    D_R,g >= alpha_R,i V_R - beta_R,i ||xi||^2,   alpha_R,i > 0,

with D_R,g the exact attitude-energy drop produced by the implemented Kalman
correction.  Such a statement must in particular have D_R,g > 0 whenever
xi=0 and V_R>0.  This producer tests that necessary condition on an *actual
source-reachable covariance/gain tuple* at the H-mode goLive seed.

The audit is deliberately stronger than a norm-only falsification:

* P_theta is the deployed anisotropic goLive covariance produced by
  initialize_from_attitude();
* P_aw is the current stationary covariance selected by the source schedule;
* the accelerometer and magnetometer gains are recomputed from that same
  covariance tuple with the validated interval Kalman backend;
* the accelerometer Joseph update and the implemented left-error covariance
  reset are applied before the magnetometer gain is formed;
* the attitude itself is updated with the existing exact deployed quaternion
  backend, and the ideal Rodrigues/exponential energy is evaluated separately;
* no independently selected gain extrema, replay sample, or fitted trajectory
  quantity is used.

A single validated xi=0 witness with D_R,g<0 disproves every alpha_R,i>0 sector
on every node containing that witness, regardless of beta_R,i.  The witness
below lies inside both deployed gauged P5 nodes and inside the declared vector
geometry/source-parameter envelope.  The purpose is therefore not to make a
false sector pass by increasing beta.  If the witness validates, the correct
outer route has to keep the source-shaped Cayley/information energy rather than
claim monotone decrease of the isotropic trace energy V_R under an anisotropic
Kalman gain.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import (
    Interval,
    matrix_add,
    matrix_identity,
    matrix_mul,
    matrix_transpose,
)
import ou3_p4_group_algebra as GROUP
import ou3_p5_go_live_covariance_stage as GOLIVE
import ou3_p5_heading_handoff_contract as HEADING
import ou3_validated_kalman_interval as VK
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1

# A simple, source-admissible point chosen to expose the geometry.  This is not
# replay data.  c is the exact Cayley attitude coordinate at the packet entry.
WITNESS = {
    "cayley": [0.15, 0.0, -0.05],
    "sigma_aw_std_mps2": 5.0,
    "mag_norm_uT": 190.0,
    "mag_horizontal_fraction": 0.4,
    "mag_vertical_sign": -1.0,
}


def p(x: float) -> Interval:
    return Interval.point(float(x))


def _zero_matrix(n: int, m: int):
    return [[p(0.0) for _ in range(m)] for _ in range(n)]


def _diag3(x: Interval | float):
    q = x if isinstance(x, Interval) else p(float(x))
    z = p(0.0)
    return [[q if i == j else z for j in range(3)] for i in range(3)]


def _scale(A, s: Interval):
    return [[a * s for a in row] for row in A]


def _vec_add(a, b):
    return [x + y for x, y in zip(a, b)]


def _vec_sub(a, b):
    return [x - y for x, y in zip(a, b)]


def _mat_vec(A, v):
    out = []
    for row in A:
        total = p(0.0)
        for a, x in zip(row, v):
            total = total + a * x
        out.append(total)
    return out


def _skew(v):
    x, y, z = v
    q = p(0.0)
    return [[q, -z, y], [z, q, -x], [-y, x, q]]


def _neg(A):
    return [[-x for x in row] for row in A]


def _reset_covariance(P, dtheta):
    G = matrix_add(matrix_identity(3), _scale(_skew(dtheta), p(0.5)))
    return matrix_mul(matrix_mul(G, P), matrix_transpose(G))


def _measurement_variances() -> tuple[float, float]:
    v = VECTOR.build()
    vf = VECTOR.validate(v)
    if vf:
        raise RuntimeError(f"configured measurement prerequisite failed: {vf}")
    c = v["configured_measurement_bounds"]
    sa = float(c["acc_measurement_std_mps2"])
    sm = float(c["mag_measurement_std_uT"])
    return sa * sa, sm * sm


def _packet_witness(domain: dict, go: dict, heading: dict) -> dict:
    live = domain["normal_live"]
    seed = go["goLive_H_covariance_seed"]
    acc_var, mag_var = _measurement_variances()

    c = [p(x) for x in WITNESS["cayley"]]
    c_norm = GROUP.vector_norm_interval(c)
    Re = GROUP.rotation_from_cayley(c)
    ReT = matrix_transpose(Re)

    # At goLive the nominal a_w state is still zero.  xi=0 therefore makes the
    # true/nominal specific-force world vector exactly -g e_3.  The H-mode
    # accelerometer bias error is also zero in this witness; its independent
    # covariance remains in the innovation covariance exactly as in the source.
    g = float(domain["startup"]["gravity_mps2"])
    fhat = [p(0.0), p(0.0), p(-g)]
    ftrue = _mat_vec(ReT, fhat)
    racc = _vec_sub(ftrue, fhat)
    Ha = _neg(_skew(fhat))

    tilt_var = float(seed["attitude_covariance_seed"]["tilt_variance"])
    yaw_var = float(seed["attitude_covariance_seed"]["gauged_yaw_variance"])
    Ptheta = [
        [p(tilt_var), p(0.0), p(0.0)],
        [p(0.0), p(tilt_var), p(0.0)],
        [p(0.0), p(0.0), p(yaw_var)],
    ]

    sigma_aw = float(WITNESS["sigma_aw_std_mps2"])
    # Source constructor: sigma_bacc0_=0.004 m/s^2.  At goLive all theta-aw and
    # theta-ba cross blocks are exactly zero.  For the theta marginal, the
    # independent P_aw and frozen P_ba blocks therefore enter exactly like an
    # additive measurement covariance in this first accelerometer correction.
    p_ba = 0.004 * 0.004
    Racc_eff_scalar = sigma_aw * sigma_aw + p_ba + acc_var
    Racc_eff = _diag3(Racc_eff_scalar)
    ua = VK.joseph_measurement_update(Ptheta, Ha, Racc_eff)
    Ka = ua["K"]
    da = _mat_vec(Ka, racc)

    Ra_deployed = GROUP.deployed_injection_rotation(da)
    Ra_exp = GROUP.rodrigues_rotation(da)
    Re1_deployed = matrix_mul(Ra_deployed, Re)
    Re1_exp = matrix_mul(Ra_exp, Re)

    # Covariance reset uses the source's G=I+0.5[d]x, not the Rodrigues matrix.
    Ptheta1 = _reset_covariance(ua["P_plus"], da)

    s = float(WITNESS["mag_horizontal_fraction"])
    if not (0.0 < s < 1.0):
        raise RuntimeError("invalid witness magnetic horizontal fraction")
    vertical = GROUP.sqrt_point(1.0 - s * s)
    if WITNESS["mag_vertical_sign"] < 0.0:
        vertical = -vertical
    M = p(float(WITNESS["mag_norm_uT"]))
    m = [M * p(s), p(0.0), M * vertical]

    # Initial nominal attitude is the coordinate frame used above.  After the
    # accelerometer injection, the nominal magnetic prediction is Ra*m while
    # the fixed true measurement is Re^T*m.
    mtrue = _mat_vec(ReT, m)
    mhat_deployed = _mat_vec(Ra_deployed, m)
    rm_deployed = _vec_sub(mtrue, mhat_deployed)
    Hm_deployed = _neg(_skew(mhat_deployed))
    um = VK.joseph_measurement_update(Ptheta1, Hm_deployed, _diag3(mag_var))
    dm = _mat_vec(um["K"], rm_deployed)

    Rm_deployed = GROUP.deployed_injection_rotation(dm)
    Rm_exp = GROUP.rodrigues_rotation(dm)
    Re2_deployed = matrix_mul(Rm_deployed, Re1_deployed)
    # The Rodrigues audit uses the same source-computed gains/residuals but
    # replaces only the tiny source-series quaternion approximation by exp(d).
    Re2_exp = matrix_mul(Rm_exp, Re1_exp)

    V0 = GROUP.group_energy(Re)
    V2_deployed = GROUP.group_energy(Re2_deployed)
    V2_exp = GROUP.group_energy(Re2_exp)
    D_deployed = V0 - V2_deployed
    D_exp = V0 - V2_exp

    da_norm = GROUP.vector_norm_interval(da)
    dm_norm = GROUP.vector_norm_interval(dm)

    # The magnetic vector was constructed as M*(s,0,+/-sqrt(1-s^2)), so its
    # sine separation from gravity is exactly s in real arithmetic.
    force_norm = g
    mag_norm = float(WITNESS["mag_norm_uT"])
    normal_q = float(heading["gauged_quality_handoff"]["full_attitude_cayley_norm_upper"])
    timeout_q = float(heading["gauged_timeout_subbranch"]["full_attitude_cayley_norm_upper"])
    source_admissible = bool(
        c_norm.hi < normal_q
        and c_norm.hi < timeout_q
        and float(live["specific_force_norm_lower_mps2"]) <= force_norm <= float(live["specific_force_norm_upper_mps2"])
        and float(live["magnetic_vector_norm_lower_uT"]) <= mag_norm <= float(live["magnetic_vector_norm_upper_uT"])
        and s >= float(live["vector_sine_separation_lower"])
        and float(seed["P_awaw_source_std_outward_mps2"][0]) <= sigma_aw <= float(seed["P_awaw_source_std_outward_mps2"][1])
    )

    disproved = bool(
        source_admissible
        and V0.lo > 0.0
        and D_deployed.hi < 0.0
        and D_exp.hi < 0.0
    )

    return {
        "source_admissible": source_admissible,
        "xi_norm": 0.0,
        "cayley_vector": WITNESS["cayley"],
        "cayley_norm_interval": c_norm.as_list(),
        "normal_gauged_q_upper": normal_q,
        "timeout_gauged_q_upper": timeout_q,
        "specific_force_norm_mps2": force_norm,
        "magnetic_norm_uT": mag_norm,
        "magnetic_sine_separation": s,
        "sigma_aw_std_mps2": sigma_aw,
        "goLive_Ptheta_diagonal": [tilt_var, tilt_var, yaw_var],
        "accelerometer_effective_R_scalar": Racc_eff_scalar,
        "accelerometer_correction_norm_interval_rad": da_norm.as_list(),
        "magnetometer_correction_norm_interval_rad": dm_norm.as_list(),
        "V_R_before_interval": V0.as_list(),
        "V_R_after_deployed_interval": V2_deployed.as_list(),
        "V_R_after_Rodrigues_interval": V2_exp.as_list(),
        "D_R_deployed_interval": D_deployed.as_list(),
        "D_R_Rodrigues_interval": D_exp.as_list(),
        "deployed_energy_increases_strictly": D_deployed.hi < 0.0,
        "Rodrigues_energy_increases_strictly": D_exp.hi < 0.0,
        "requested_positive_alpha_sector_disproved": disproved,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("large-angle sector audit domain must not be trajectory fitted")

    go = GOLIVE.build(domain_path)
    heading = HEADING.build(domain_path)
    prereq = [f"goLive: {x}" for x in GOLIVE.validate(go)]
    prereq += [f"heading: {x}" for x in HEADING.validate(heading)]
    if prereq:
        return {
            "schema": SCHEMA,
            "qualification": "OU3_P5_SOURCE_CORRELATED_LARGE_ANGLE_VR_SECTOR_AUDIT",
            "source_generated_not_trajectory_fit": True,
            "source_replay_used": False,
            "P5_RAW_VR_LARGE_ANGLE_SECTOR": "NOT_EVALUATED",
            "failures": prereq,
        }

    witness = _packet_witness(domain, go, heading)
    disproved = witness["requested_positive_alpha_sector_disproved"] is True
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_SOURCE_CORRELATED_LARGE_ANGLE_VR_SECTOR_AUDIT",
        "claim": "VALIDATED_NECESSARY_CONDITION_AUDIT_ON_ACTUAL_GOLIVE_COVARIANCE_GAIN_TUPLE",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "exact_deployed_quaternion_backend_used": True,
        "exact_Rodrigues_backend_used": True,
        "validated_interval_Kalman_gain_used": True,
        "source_correlated_covariance_gain_tuple": True,
        "independent_gain_extrema_used": False,
        "sector_requested": "D_R,g >= alpha_R,i V_R - beta_R,i ||xi||^2 with alpha_R,i>0",
        "gauged_nodes": {
            "normal_q_upper": heading["gauged_quality_handoff"]["full_attitude_cayley_norm_upper"],
            "timeout_q_upper": heading["gauged_timeout_subbranch"]["full_attitude_cayley_norm_upper"],
        },
        "validated_counterexample": witness,
        "counterexample_inside_both_gauged_nodes": witness["source_admissible"],
        "beta_cannot_repair_xi_zero_counterexample": disproved,
        "P5_RAW_VR_LARGE_ANGLE_SECTOR": (
            "DISPROVED_ON_DECLARED_SOURCE_FAMILY" if disproved else "NOT_DISPROVED_BY_WITNESS"
        ),
        "required_theorem_correction": (
            "use the source-shaped Cayley/information path energy for quantitative outer dissipation; retain V_R as exact group/chart geometry, not as a source-uniform monotone energy under arbitrary anisotropic Kalman gains"
        ),
        "single_vector_positive_alpha_sector_also_impossible": True,
        "single_vector_reason": (
            "an accepted single vector correction has an exact rotation-about-vector null direction, so full-attitude V_R cannot have alpha>0 per individual correction"
        ),
        "next_obligation": (
            "construct the outer finite-angle sector in the source-correlated Cayley/information metric and certify packet prefixes with the exact group backend; do not tune beta to hide the xi=0 V_R counterexample"
        ),
        "failures": [],
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    if d.get("source_generated_not_trajectory_fit") is not True:
        failures.append("large-angle sector audit is not source generated")
    if d.get("source_replay_used") is not False:
        failures.append("large-angle sector audit uses replay")
    if d.get("filter_changed") is not False:
        failures.append("large-angle sector audit changes the filter")
    for key in (
        "exact_deployed_quaternion_backend_used",
        "exact_Rodrigues_backend_used",
        "validated_interval_Kalman_gain_used",
        "source_correlated_covariance_gain_tuple",
    ):
        if d.get(key) is not True:
            failures.append(f"missing proof property {key}")
    if d.get("independent_gain_extrema_used") is not False:
        failures.append("large-angle audit mixes independent gain extrema")
    w = d.get("validated_counterexample", {})
    if w.get("source_admissible") is not True:
        failures.append("counterexample is not inside the declared source family")
    if w.get("xi_norm") != 0.0:
        failures.append("counterexample does not have xi=0")
    if w.get("deployed_energy_increases_strictly") is not True:
        failures.append("deployed exact group backend did not validate strict V_R increase")
    if w.get("Rodrigues_energy_increases_strictly") is not True:
        failures.append("Rodrigues backend did not validate strict V_R increase")
    if w.get("requested_positive_alpha_sector_disproved") is not True:
        failures.append("requested positive-alpha V_R sector was not disproved")
    if d.get("P5_RAW_VR_LARGE_ANGLE_SECTOR") != "DISPROVED_ON_DECLARED_SOURCE_FAMILY":
        failures.append("raw V_R sector status is not the validated disproof")
    if d.get("beta_cannot_repair_xi_zero_counterexample") is not True:
        failures.append("audit incorrectly allows beta to hide xi=0 counterexample")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve())
    failures = validate(out)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out.get("P5_RAW_VR_LARGE_ANGLE_SECTOR"),
        "counterexample": out.get("validated_counterexample"),
        "next_obligation": out.get("next_obligation"),
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
