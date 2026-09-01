#!/usr/bin/env python3
"""Joint attitude/a_w first-accelerometer innovation contraction for P5.

A large first accelerometer correction is not an isolated attitude kick.  The
same Kalman update simultaneously changes attitude and latent acceleration.
Bounding only ||dtheta|| and then forgetting the correlated a_w correction loses
that mechanism and makes the generic 45 deg P5 route appear much more violent
than the implemented filter.

At first Live in H mode the source-audited geometry gives

* J_aw=R_wb, hence J_aw J_aw^T=I;
* P_aw=p_aw I after prediction / optional first S=0 update;
* P_theta,aw=0 exactly;
* R_acc=r I.

For the modeled/effective measurement state z, with innovation r_e=H z, the
accepted state correction z+=z-K r_e has the exact measurement-space map

    H z+ = (I-HK) r_e = R_acc S^-1 r_e,
    S = H P H^T + R_acc.

Because H P H^T >= p_aw I and R_acc=r I,

    ||I-HK||_2 <= r/(r+p_aw).

This is a source-correlated *joint* attitude/a_w contraction; it neither says
attitude alone contracts nor ignores the finite-angle effective-a_w embedding.
The exact accelerometer finite-angle residual is already representable as an
additional a_w tangent input in the certified P5 effective-vector lemma, so it
belongs to the state-driven residual above.  The held physical accelerometer
bias is different: H mode has no active b_a state, so its declared deterministic
error is carried additively and is not multiplied by the contraction factor.

The producer evaluates the factor over every first source/phase cell and reports
a conservative post-correction effective-residual bound for the generic 45 deg
P5 entrance.  It is an algebraic first-update certificate only.  Quaternion
injection/reset and one-step re-evaluation of the physical residual remain the
next obligation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_candidate_first_accel_exact_source as FIRST
import ou3_p5_45deg_first_accel_q8_bridge as Q8
import ou3_p5_effective_vector_input as VEFF
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = FIRST.DEFAULT_DOMAIN
SCHEMA = 1


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2) -> dict:
    path = Path(domain_path).resolve()
    dom = json.loads(path.read_text(encoding="utf-8"))
    first = FIRST.build(path, source_pieces=source_pieces)
    q8 = Q8.build(path, source_pieces=source_pieces)
    veff = VEFF.build(path)
    vector = VECTOR.build()
    failures = [f"first: {x}" for x in FIRST.validate(first)]
    failures += [f"q8: {x}" for x in Q8.validate(q8)]
    failures += [f"effective-vector: {x}" for x in VEFF.validate(veff)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    if veff.get("accelerometer", {}).get("effective_aw_defect_norm_equals_eta_norm") is not True:
        failures.append("finite-angle accelerometer defect is not certified as effective a_w input")
    if first.get("first_accel_yaw_covariance_axis_aligned_with_force_axis") is not True:
        failures.append("first accelerometer source alignment lost")

    q = float(q8["P5_45deg_entrance_first_accel"]["post_prediction_q_upper"])
    gravity = float(dom["startup"]["gravity_mps2"])
    ba = float(dom["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    acc_std = float(vector["configured_measurement_bounds"]["acc_measurement_std_mps2"])
    racc = up(acc_std * acc_std)
    if not racc > 0.0:
        failures.append("accelerometer covariance is not positive")

    # The source/phase family and p_aw are independent of the candidate angle;
    # use the first candidate row solely as a convenient complete inventory.
    crows = first.get("candidate_rows", [])
    source_rows = crows[0].get("source_rows", []) if crows else []
    if not source_rows:
        failures.append("first source/phase rows missing")

    rotational_state = up(gravity * q)
    rows = []
    worst_gamma = 0.0
    worst_pre = 0.0
    worst_post = 0.0
    first_bad = None
    for row in source_rows:
        try:
            paw_lo = float(row["P_aw_variance_interval"][0])
            aw = float(row["predicted_aw_error_norm_upper_mps2"])
            if not (math.isfinite(paw_lo) and paw_lo >= 0.0):
                raise RuntimeError("invalid p_aw lower bound")
            den = down(racc + paw_lo)
            if not den > 0.0:
                raise RuntimeError("joint innovation contraction denominator lost positivity")
            gamma = min(1.0, up(racc / den))
            state_pre = up(rotational_state + aw)
            total_pre = up(state_pre + ba)
            state_post = up(gamma * state_pre)
            total_post = up(state_post + ba)
            strict = gamma < 1.0 and state_post < state_pre
            r = {
                "source_phase_cell": int(row["source_phase_cell"]),
                "pseudo_phase": row["pseudo_phase"],
                "tau_s": row["tau_s"],
                "sigma_aw_mps2": row["sigma_aw_mps2"],
                "R_S_filter_std": row["R_S_filter_std"],
                "P_aw_variance_lower": paw_lo,
                "Racc_variance_upper": racc,
                "joint_state_residual_contraction_factor_upper": gamma,
                "generic_45deg_rotational_effective_residual_upper_mps2": rotational_state,
                "pre_update_state_driven_effective_residual_upper_mps2": state_pre,
                "held_accel_bias_additive_upper_mps2": ba,
                "pre_update_total_effective_residual_upper_mps2": total_pre,
                "post_state_correction_state_driven_effective_residual_upper_mps2": state_post,
                "post_state_correction_total_effective_residual_upper_mps2": total_post,
                "state_driven_residual_strictly_contracts": strict,
            }
            rows.append(r)
            worst_gamma = max(worst_gamma, gamma)
            worst_pre = max(worst_pre, total_pre)
            worst_post = max(worst_post, total_post)
            if not strict and first_bad is None:
                first_bad = r
        except Exception as exc:
            first_bad = {
                "source_phase_cell": row.get("source_phase_cell"),
                "pseudo_phase": row.get("pseudo_phase"),
                "reason": f"{type(exc).__name__}: {exc}",
            }
            break

    complete = len(rows) == len(source_rows) and first_bad is None and bool(rows)
    strict_all = complete and all(r["state_driven_residual_strictly_contracts"] for r in rows)
    if not complete:
        failures.append("joint first-accelerometer source family incomplete")
    if not strict_all:
        failures.append("joint state-driven first-accelerometer residual did not strictly contract on every cell")

    passed = not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_FIRST_ACCEL_JOINT_ATTITUDE_AW_INNOVATION_CONTRACTION",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "generic_P5_45deg_entrance_used": True,
        "additional_P1_tilt_bound_used": False,
        "H_first_live_mode": True,
        "J_aw_orthogonality_used": True,
        "first_theta_aw_cross_covariance_exact_zero_used": True,
        "isotropic_aw_marginal_used": True,
        "isotropic_Racc_used": True,
        "exact_measurement_space_identity": "I-HK=Racc*S^-1",
        "finite_angle_accel_defect_embedded_as_effective_aw_input": True,
        "held_accel_bias_falsely_contracted": False,
        "attitude_alone_claimed_contracting": False,
        "quaternion_reset_composed_here": False,
        "source_phase_cell_count": len(source_rows),
        "evaluated_source_phase_cells": len(rows),
        "worst_joint_state_residual_contraction_factor_upper": worst_gamma,
        "worst_pre_update_total_effective_residual_upper_mps2": worst_pre,
        "worst_post_state_correction_total_effective_residual_upper_mps2": worst_post,
        "total_effective_residual_improvement_factor_lower": (
            worst_pre / worst_post if worst_post > 0.0 else math.inf),
        "source_rows": rows,
        "first_unclosed_child": first_bad,
        "P5_FIRST_ACCEL_JOINT_INNOVATION_CONTRACTION": "PASS" if passed else "NOT_ESTABLISHED",
        "P5_FINITE_CAPTURE_TO_P4_ESTABLISHED_HERE": False,
        "next_obligation": (
            "compose the exact quaternion injection/reset with this joint residual contraction and the physical H group-norm child; re-evaluate the source-correlated effective accelerometer residual at sample1 before applying magnetometer yaw information"
        ),
        "failures": list(dict.fromkeys(failures)),
    }


def validate(d: dict) -> list[str]:
    f = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        f.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit", "generic_P5_45deg_entrance_used",
        "H_first_live_mode", "J_aw_orthogonality_used",
        "first_theta_aw_cross_covariance_exact_zero_used", "isotropic_aw_marginal_used",
        "isotropic_Racc_used", "finite_angle_accel_defect_embedded_as_effective_aw_input",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "additional_P1_tilt_bound_used",
        "held_accel_bias_falsely_contracted", "attitude_alone_claimed_contracting",
        "quaternion_reset_composed_here", "P5_FINITE_CAPTURE_TO_P4_ESTABLISHED_HERE",
    ):
        if d.get(k) is not False:
            f.append(f"{k} is not false")
    if d.get("exact_measurement_space_identity") != "I-HK=Racc*S^-1":
        f.append("measurement-space identity mismatch")
    if int(d.get("evaluated_source_phase_cells", 0)) != int(d.get("source_phase_cell_count", -1)):
        f.append("source family incomplete")
    if d.get("first_unclosed_child") is not None:
        f.append("joint source family retains unclosed child")
    g = float(d.get("worst_joint_state_residual_contraction_factor_upper", math.inf))
    if not (0.0 <= g < 1.0):
        f.append("worst joint residual contraction factor is not strict")
    pre = float(d.get("worst_pre_update_total_effective_residual_upper_mps2", -math.inf))
    post = float(d.get("worst_post_state_correction_total_effective_residual_upper_mps2", math.inf))
    if not (0.0 < post < pre < math.inf):
        f.append("worst total effective residual did not improve")
    if d.get("P5_FIRST_ACCEL_JOINT_INNOVATION_CONTRACTION") == "PASS" and f:
        f.append("PASS carries validation failures")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--output", type=Path, required=True)
    x = ap.parse_args()
    d = build(x.domain.resolve(), source_pieces=x.source_pieces)
    vf = validate(d)
    d["validation_pass"] = not vf
    d["validation_failures"] = vf
    x.output.parent.mkdir(parents=True, exist_ok=True)
    x.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": d["P5_FIRST_ACCEL_JOINT_INNOVATION_CONTRACTION"],
        "cells": d["evaluated_source_phase_cells"],
        "gamma_max": d["worst_joint_state_residual_contraction_factor_upper"],
        "pre_total": d["worst_pre_update_total_effective_residual_upper_mps2"],
        "post_total": d["worst_post_state_correction_total_effective_residual_upper_mps2"],
        "improvement_factor": d["total_effective_residual_improvement_factor_lower"],
        "validation_failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
