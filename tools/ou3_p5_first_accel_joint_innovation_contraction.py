#!/usr/bin/env python3
"""Joint attitude/a_w first-accelerometer innovation contraction for P5.

A large first accelerometer correction is not an isolated attitude kick. The
same Kalman update simultaneously changes attitude and latent acceleration.
Bounding only ||dtheta|| and forgetting the correlated a_w correction loses the
implemented filter's strongest first-packet mechanism.

At first Live in H mode the source-audited geometry gives J_aw=R_wb orthogonal,
P_aw=p_aw I, P_theta,aw=0, and R_acc=r I. For the *current effective input* u
seen by the linear Kalman correction,

    (I-HK) u = R_acc S^-1 u,     S=H P H^T+R_acc,

hence, because HPH^T >= p_aw I,

    ||I-HK||_2 <= r/(r+p_aw).

The finite-angle accelerometer defect is permitted here because the established
P5 effective-vector lemma represents its action on the current correction as an
effective a_w tangent input. This producer does NOT call the expression above
the next physical accelerometer residual: quaternion injection/reset changes the
attitude and therefore changes that effective embedding. It is a same-
linearization innovation-cancellation bound which must be composed with the
exact reset and then re-embedded at sample 1.

A due first S=0 update needs care. Direct interval subtraction for
p_aw - p_Saw^2/(p_SS+R_S) may touch zero even though the exact conditional
covariance is positive. The exact 2x2 covariance block is PSD, so
p_aw p_SS-p_Saw^2>=0 and therefore

 p_aw|S >= p_aw^- R_S^- /(p_SS^+ + R_S^+) > 0.

That Schur lower is used only for the contraction factor. The ordinary interval
posterior is retained for all other covariance propagation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import ou3_p4_candidate_first_accel_exact_source as FIRST
import ou3_p5_45deg_first_accel_q8_bridge as Q8
import ou3_p5_effective_vector_input as VEFF
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_first_accel_rotation_gauge_v3 as RG3
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = FIRST.DEFAULT_DOMAIN
SCHEMA = 2


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _due_paw_positive_lower(paw_pred, pss, rs2) -> float:
    """Positive lower from exact PSD Schur geometry, avoiding interval cancellation."""
    num = down(float(paw_pred.lo) * float(rs2.lo))
    den = up(float(pss.hi) + float(rs2.hi))
    if not (num > 0.0 and den > 0.0):
        raise RuntimeError("due first-S Schur lower lost positivity")
    return down(num / den)


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
    racc_lo = down(acc_std * acc_std)
    racc_hi = up(acc_std * acc_std)
    if not racc_lo > 0.0:
        failures.append("accelerometer covariance is not positive")

    RG3._install_backend(path, source_pieces)
    FULL3._install_backend()
    src_phases = RG._source_phase_children(source_pieces)
    if not src_phases:
        failures.append("first source/phase family missing")

    rotational_state = up(gravity * q)
    rows = []
    worst_gamma = 0.0
    worst_pre = 0.0
    worst_post = 0.0
    min_cell_improvement = math.inf
    first_bad = None
    for si, (src, phase) in enumerate(src_phases):
        try:
            P0 = FULL._initial_covariance(src, path)
            F, Q, _ = FULL._transition_and_Q(src, dom)
            Pp = FULL._psd_tighten(FULL.matrix_add(
                FULL.matrix_mul(FULL.matrix_mul(F, P0), FULL.matrix_transpose(F)), Q))
            pss, _psa, paw_pred = RG._scalar_axis_structure(Pp)
            aw_pred, eS_pred = RG._prediction_norms(src, dom)
            if phase == "due":
                rs2 = src["R_S_filter_std"].square()
                paw_interval, aw_norm = RG._due_paw_and_error_norm(Pp, src, aw_pred, eS_pred)
                paw_lo = _due_paw_positive_lower(paw_pred, pss, rs2)
                paw_route = "PSD_SCHUR_LOWER_pawR_over_pSSplusR"
            else:
                paw_interval = paw_pred
                aw_norm = aw_pred
                paw_lo = float(paw_pred.lo)
                paw_route = "PREDICTED_ISOTROPIC_PAW_LOWER"

            if not (math.isfinite(paw_lo) and paw_lo > 0.0):
                raise RuntimeError("invalid positive p_aw lower bound")
            den = down(racc_lo + paw_lo)
            if not den > 0.0:
                raise RuntimeError("joint innovation contraction denominator lost positivity")
            gamma = min(1.0, up(racc_hi / den))

            state_pre = up(rotational_state + aw_norm)
            total_pre = up(state_pre + ba)
            # This is the same-linearization current-effective-input remainder,
            # not yet the next physical measurement residual after reset.
            state_remainder = up(gamma * state_pre)
            total_remainder = up(state_remainder + ba)
            strict = gamma < 1.0 and state_remainder < state_pre
            improvement = total_pre / total_remainder if total_remainder > 0.0 else math.inf
            r = {
                "source_phase_cell": si,
                "pseudo_phase": phase,
                "tau_s": src["tau_s"].as_list(),
                "sigma_aw_mps2": src["sigma_aw_mps2"].as_list(),
                "R_S_filter_std": src["R_S_filter_std"].as_list(),
                "P_aw_interval_after_optional_S": paw_interval.as_list(),
                "P_aw_strict_lower_for_joint_factor": paw_lo,
                "P_aw_lower_route": paw_route,
                "Racc_variance_interval": [racc_lo, racc_hi],
                "joint_current_effective_input_remainder_factor_upper": gamma,
                "generic_45deg_rotational_effective_input_upper_mps2": rotational_state,
                "pre_update_state_driven_effective_input_upper_mps2": state_pre,
                "held_accel_bias_additive_upper_mps2": ba,
                "pre_update_total_effective_input_upper_mps2": total_pre,
                "same_linearization_state_driven_remainder_upper_mps2": state_remainder,
                "same_linearization_total_remainder_upper_mps2": total_remainder,
                "cell_total_remainder_improvement_factor_lower": improvement,
                "state_driven_current_effective_input_strictly_reduced": strict,
            }
            rows.append(r)
            worst_gamma = max(worst_gamma, gamma)
            worst_pre = max(worst_pre, total_pre)
            worst_post = max(worst_post, total_remainder)
            min_cell_improvement = min(min_cell_improvement, improvement)
            if not strict and first_bad is None:
                first_bad = r
        except Exception as exc:
            first_bad = {
                "source_phase_cell": si,
                "pseudo_phase": phase,
                "reason": f"{type(exc).__name__}: {exc}",
            }
            break

    complete = len(rows) == len(src_phases) and first_bad is None and bool(rows)
    strict_all = complete and all(r["state_driven_current_effective_input_strictly_reduced"] for r in rows)
    if not complete:
        failures.append("joint first-accelerometer source family incomplete")
    if not strict_all:
        failures.append("joint current effective input did not strictly reduce on every cell")

    passed = not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_FIRST_ACCEL_JOINT_ATTITUDE_AW_CURRENT_INNOVATION_CANCELLATION",
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
        "due_S_positive_conditional_aw_Schur_lower_used": True,
        "finite_angle_accel_defect_embedded_as_current_effective_aw_input": True,
        "held_accel_bias_falsely_contracted": False,
        "attitude_alone_claimed_contracting": False,
        "same_linearization_remainder_claimed_as_next_physical_residual": False,
        "quaternion_reset_composed_here": False,
        "source_phase_cell_count": len(src_phases),
        "evaluated_source_phase_cells": len(rows),
        "worst_joint_current_effective_input_remainder_factor_upper": worst_gamma,
        "worst_pre_update_total_effective_input_upper_mps2": worst_pre,
        "worst_same_linearization_total_remainder_upper_mps2": worst_post,
        "minimum_cell_total_remainder_improvement_factor_lower": min_cell_improvement,
        "source_rows": rows,
        "first_unclosed_child": first_bad,
        "P5_FIRST_ACCEL_JOINT_INNOVATION_CONTRACTION": "PASS" if passed else "NOT_ESTABLISHED",
        "P5_FINITE_CAPTURE_TO_P4_ESTABLISHED_HERE": False,
        "next_obligation": (
            "compose exact quaternion injection/reset with the joint current-input cancellation and physical H norm child, then re-embed the finite-angle accelerometer input at sample1; only that re-embedded quantity may be called the next physical residual"
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
        "isotropic_Racc_used", "due_S_positive_conditional_aw_Schur_lower_used",
        "finite_angle_accel_defect_embedded_as_current_effective_aw_input",
    ):
        if d.get(k) is not True:
            f.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "additional_P1_tilt_bound_used",
        "held_accel_bias_falsely_contracted", "attitude_alone_claimed_contracting",
        "same_linearization_remainder_claimed_as_next_physical_residual",
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
    g = float(d.get("worst_joint_current_effective_input_remainder_factor_upper", math.inf))
    if not (0.0 <= g < 1.0):
        f.append("worst joint current-input remainder factor is not strict")
    imp = float(d.get("minimum_cell_total_remainder_improvement_factor_lower", 0.0))
    if not (math.isfinite(imp) and imp > 1.0):
        f.append("no strict total-remainder improvement on every source cell")
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
        "gamma_max": d["worst_joint_current_effective_input_remainder_factor_upper"],
        "pre_total": d["worst_pre_update_total_effective_input_upper_mps2"],
        "same_linearization_remainder": d["worst_same_linearization_total_remainder_upper_mps2"],
        "min_cell_improvement": d["minimum_cell_total_remainder_improvement_factor_lower"],
        "validation_failures": vf,
        "next": d["next_obligation"],
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
