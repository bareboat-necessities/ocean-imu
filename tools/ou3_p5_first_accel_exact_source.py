#!/usr/bin/env python3
"""Exact source-reachable first accelerometer gate for OU-III P5.

The preceding structured-gain stage intentionally still allowed an arbitrary
normal-Live force vector and an independently oriented yaw covariance axis at
sample zero.  Those combinations are not source reachable at the deployed
Mahony-proxy handoff.

At the first Live accelerometer correction:

* the MEKF extended state was constructed with xext=0;
* MahonyProxy warmup advances the measurement-only front end and leaves the
  MEKF state untouched;
* initialize_from_attitude changes attitude/error covariance, not v,p,S,a_w;
* enabling the linear block resets a_w covariance and pseudo phase, not its
  mean;
* the first homogeneous OU prediction preserves a_w=v=p=S=0;
* an optional first S=0 pseudo update has residual exactly zero and therefore
  also preserves those means.

Consequently the first shipping accelerometer Jacobian is built from

    f_cog_b = R_wb (0-g_world),

so ||f_cog_b|| is exactly g and its axis is the body image of world down.  The
anisotropic goLive attitude covariance uses that same body-down axis for its yaw
rank-one term.  Prediction transports the rank-one axis with the body frame;
the gyro-bias/process contribution is already carried by the small PSD remainder
E of the structured-gain lemma.  Hence the source-reachable rank-one alignment
coordinate is x=0, not an arbitrary x in [0,1].

The physical attitude error is also not charged by its full Cayley norm in the
accelerometer residual.  If c_t is the Cayley component tangent to gravity,
then exactly

    cos(gamma)=1-2||c_t||^2/(4+||c||^2),

where gamma is the true-gravity misalignment.  P1 supplies a lower cosine at
handoff; one 5 ms bounded transport step gives a post-prediction lower cosine,
and the already certified full q bound then gives a finite ||c_t|| bound.
For a gravity vector the exact rotational measurement residual satisfies

    ||y_R|| <= g ||c_t||.

Finally the latent linear term and its finite-angle cross term must not be
triangled independently:

    e_aw + (R(c)^T-I)e_aw = R(c)^T e_aw,

so their combined norm is exactly ||e_aw||.  This removes the artificial
(1+||R-I||)||e_aw|| penalty that dominated the previous 11.51-rad result.

The output certifies only that every source-reachable first accelerometer
correction remains inside the already validated [0,6] deployed-quaternion
range.  It does not promote the complete q<=8 word or set N_H_words.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

from ou3_interval import Interval
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_first_accel_rotation_gauge_v3 as RG3
import ou3_p5_first_accel_structured_gain as SG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_p5_heading_handoff_contract as HEADING
import ou3_startup_stability_certificate as P1
import ou3_vector_uco_certificate as VECTOR

REPO = Path(__file__).resolve().parents[1]
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = RG.DEFAULT_DOMAIN
SCHEMA = 1
DEPLOYED_CORRECTION_LIMIT_RAD = 6.0


def _source_semantics() -> tuple[dict, list[str]]:
    k = MEKF.read_text(encoding="utf-8")
    w = WRAPPER.read_text(encoding="utf-8")
    markers = {
        "constructor_zeroes_extended_mean": "xext.setZero();",
        "front_end_does_not_drive_mekf": "updateCore_(dt, gyro, acc, /*tempC=*/35.0f, /*drive_mekf=*/false);",
        "attitude_init_only_zeroes_attitude_error_mean": "set_quaternion_boat(q_bw);",
        "linear_enable_resets_aw_covariance": "reset_aw_covariance_to_stationary();",
        "linear_enable_resets_pseudo_phase": "pseudo_update_elapsed_s_ = T(0);",
        "first_accel_mean_uses_aw_state": "const Vector3 aw = xext.template segment<3>(OFF_AW);",
        "first_accel_force_is_aw_minus_gravity": "f_cog_b = R_wb() * (aw - g_world);",
        "yaw_covariance_axis_is_body_world_down": "const Vector3 u_down_body = R_wb() * world_down;",
        "S_pseudo_residual_is_minus_S": "const Vector3 r = -Sstate;",
    }
    joined = k + "\n" + w
    missing = [name for name, marker in markers.items() if marker not in joined]
    return {
        "source_markers": markers,
        "constructor_extended_mean_zero": "constructor_zeroes_extended_mean" not in missing,
        "mahony_proxy_front_end_leaves_mekf_state_untouched": "front_end_does_not_drive_mekf" not in missing,
        "linear_enable_changes_aw_covariance_not_mean": "linear_enable_resets_aw_covariance" not in missing,
        "first_accel_force_reads_aw_minus_gravity": "first_accel_force_is_aw_minus_gravity" not in missing,
        "goLive_yaw_covariance_axis_is_body_gravity_axis": "yaw_covariance_axis_is_body_world_down" not in missing,
        "zero_S_mean_makes_first_due_pseudo_mean_correction_zero": "S_pseudo_residual_is_minus_S" not in missing,
    }, missing


def _post_prediction_gravity_cosine_lower(domain: dict, handoff_cos: float, h: float) -> tuple[float, float]:
    b = domain["startup"]["physical_handoff_coordinate_bounds"]
    bg = float(b["gyro_bias_error_norm_upper_rad_s"])
    wd = float(domain["startup"]["effective_deterministic_gyro_transport_disturbance_upper_rad_s"])
    a = FULL.up(h * (bg + wd))
    if not (0.0 <= a < 1.0 and -1.0 < handoff_cos <= 1.0):
        raise RuntimeError("invalid first-step gravity transport bound")
    cos_a_lo = FULL.down(math.sqrt(max(0.0, FULL.down(1.0 - FULL.up(a * a)))))
    sin_gamma_hi = FULL.up(math.sqrt(max(0.0, FULL.up(1.0 - FULL.down(handoff_cos * handoff_cos)))))
    c = FULL.down(FULL.down(handoff_cos * cos_a_lo) - FULL.up(sin_gamma_hi * a))
    return max(-1.0, c), a


def _cayley_tangent_upper(q: float, gravity_cos_lower: float) -> float:
    if not (q >= 0.0 and -1.0 < gravity_cos_lower <= 1.0):
        raise RuntimeError("invalid Cayley/gravity bound")
    d = FULL.up(4.0 + FULL.up(q * q))
    r2 = FULL.up(FULL.up(0.5 * d) * FULL.up(1.0 - gravity_cos_lower))
    return FULL.up(math.sqrt(max(0.0, r2)))


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2) -> dict:
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("first-accelerometer exact-source domain must not be trajectory fitted")
    if domain.get("configured_runtime", {}).get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("first-accelerometer exact-source stage requires lever arm disabled")

    source, missing = _source_semantics()
    p1 = P1.build(domain_path)
    heading = HEADING.build(domain_path)
    vector = VECTOR.build()
    failures = [f"source semantic missing: {x}" for x in missing]
    failures += [f"P1: {x}" for x in P1.validate(p1)]
    failures += [f"heading: {x}" for x in HEADING.validate(heading)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]

    RG3._install_backend(domain_path, source_pieces)
    FULL3._install_backend()
    src_phases = RG._source_phase_children(source_pieces)
    if not src_phases:
        failures.append("no first-prefix source phase cells")

    h = float(FULL._source_cell()["dt_s"])
    gravity = float(domain["startup"]["gravity_mps2"])
    m = Interval.outward_bounds(gravity, gravity)
    q0 = float(heading["gauged_timeout_subbranch"]["full_attitude_cayley_norm_upper"])
    qpred = RG._q_after_first_prediction(q0, domain, h)
    handoff_tilt_cos = float(heading["gauged_timeout_subbranch"]["tilt_cosine_lower"])
    pred_gravity_cos, transport_angle = _post_prediction_gravity_cosine_lower(domain, handoff_tilt_cos, h)
    c_tan = _cayley_tangent_upper(qpred, pred_gravity_cos)

    tilt, yaw, eps = RG._attitude_covariance_epsilon(domain_path, h)
    vc = vector["configured_measurement_bounds"]
    Racc = FULL._R_diag(float(vc["acc_measurement_std_mps2"]))
    racc_var = Racc[0][0]
    ba = float(domain["startup"]["physical_handoff_coordinate_bounds"]["accelerometer_bias_error_norm_upper_mps2"])
    x_aligned = Interval.outward_bounds(0.0, 0.0)

    max_d = 0.0
    min_margin = math.inf
    max_k = 0.0
    first_over = None
    over = 0
    rows = []

    for si, (src, phase) in enumerate(src_phases):
        P0 = FULL._initial_covariance(src, domain_path)
        F, Q, _ = FULL._transition_and_Q(src, domain)
        Pp = FULL._psd_tighten(FULL.matrix_add(FULL.matrix_mul(FULL.matrix_mul(F, P0), FULL.matrix_transpose(F)), Q))
        _pss, _psa, paw_pred = RG._scalar_axis_structure(Pp)
        aw_pred, eS_pred = RG._prediction_norms(src, domain)
        if phase == "due":
            paw, aw_norm = RG._due_paw_and_error_norm(Pp, src, aw_pred, eS_pred)
        else:
            paw, aw_norm = paw_pred, aw_pred

        k, _kh, detail = SG._structured_gain_bounds(
            tilt=tilt,
            yaw=yaw,
            eps=eps,
            x=x_aligned,
            m=m,
            paw=paw,
            racc_var=racc_var,
        )
        # Exact source combination:
        # 1) rotational vector residual <= g*|c_tangent|;
        # 2) e_aw + (R^T-I)e_aw = R^T e_aw, norm exactly |e_aw|;
        # 3) accelerometer-bias error remains additive.
        rotational_residual = FULL.up(gravity * c_tan)
        useful_residual = FULL.up(rotational_residual + FULL.up(aw_norm + ba))
        d = FULL.up(k * useful_residual)
        margin = DEPLOYED_CORRECTION_LIMIT_RAD - d
        max_k = max(max_k, k)
        max_d = max(max_d, d)
        min_margin = min(min_margin, margin)
        row = {
            "source_phase_cell": si,
            "pseudo_phase": phase,
            "tau_s": src["tau_s"].as_list(),
            "sigma_aw_mps2": src["sigma_aw_mps2"].as_list(),
            "R_S_filter_std": src["R_S_filter_std"].as_list(),
            "P_aw_variance_interval": paw.as_list(),
            "predicted_aw_error_norm_upper_mps2": aw_norm,
            "Ktheta_norm_upper": k,
            "rotational_residual_norm_upper_mps2": rotational_residual,
            "combined_useful_residual_norm_upper_mps2": useful_residual,
            "correction_norm_upper_rad": d,
            "gain_detail": detail,
        }
        rows.append(row)
        if not math.isfinite(d) or d > DEPLOYED_CORRECTION_LIMIT_RAD:
            over += 1
            if first_over is None:
                first_over = row

    closed = bool(rows) and over == 0 and not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_FIRST_ACCEL_EXACT_STARTUP_SOURCE_GEOMETRY",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "source_semantics": source,
        "first_accel_aw_mean_exact_zero_before_measurement": True,
        "first_due_S_mean_correction_exact_zero": True,
        "first_accel_specific_force_magnitude_exact_gravity": True,
        "first_accel_specific_force_magnitude_mps2": gravity,
        "first_accel_yaw_covariance_axis_aligned_with_force_axis": True,
        "yaw_alignment_x_equals_zero": True,
        "post_prediction_full_cayley_norm_upper": qpred,
        "handoff_true_gravity_cosine_lower": handoff_tilt_cos,
        "first_prediction_transport_angle_upper_rad": transport_angle,
        "post_prediction_true_gravity_cosine_lower": pred_gravity_cos,
        "post_prediction_cayley_tangent_norm_upper": c_tan,
        "exact_gravity_cayley_tangent_identity_used": True,
        "exact_rotational_residual_tangent_bound_used": True,
        "latent_linear_plus_rotation_cross_combined_before_norm": True,
        "latent_combined_norm_identity": "||e+(R^T-I)e||=||R^T e||=||e||",
        "independent_latent_rotation_cross_norm_added": False,
        "analytic_structured_gain_no_matrix_inverse": True,
        "attitude_PSD_remainder_retained": True,
        "deployed_correction_limit_rad": DEPLOYED_CORRECTION_LIMIT_RAD,
        "deployed_correction_limit_increased": False,
        "evaluated_source_phase_cells": len(rows),
        "children_above_validated_correction_limit": over,
        "max_Ktheta_norm_upper": max_k,
        "max_first_accelerometer_correction_norm_upper_rad": max_d,
        "minimum_correction_range_margin_rad": min_margin,
        "all_first_accelerometer_source_cells_inside_validated_correction_range": closed,
        "first_unclosed_child": first_over,
        "source_cells": rows,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_FIRST_ACCEL_EXACT_SOURCE_CERTIFICATE": "PASS" if closed else "NOT_ESTABLISHED",
        "next_obligation": (
            "PROPAGATE_FIRST_ACCEL_EXACT_SOURCE_CELLS_THROUGH_JOSEPH_RESET_AND_NEXT_PREFIX"
            if closed else
            "DERIVE_SMALLER_SOURCE_REACHABLE_FIRST_ACCEL_STATE_CHART"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "first_accel_aw_mean_exact_zero_before_measurement",
        "first_due_S_mean_correction_exact_zero",
        "first_accel_specific_force_magnitude_exact_gravity",
        "first_accel_yaw_covariance_axis_aligned_with_force_axis",
        "yaw_alignment_x_equals_zero",
        "exact_gravity_cayley_tangent_identity_used",
        "exact_rotational_residual_tangent_bound_used",
        "latent_linear_plus_rotation_cross_combined_before_norm",
        "analytic_structured_gain_no_matrix_inverse",
        "attitude_PSD_remainder_retained",
    ):
        if d.get(k) is not True:
            failures.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "independent_latent_rotation_cross_norm_added",
        "deployed_correction_limit_increased", "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            failures.append(f"{k} is not false")
    s = d.get("source_semantics", {})
    for k in (
        "constructor_extended_mean_zero",
        "mahony_proxy_front_end_leaves_mekf_state_untouched",
        "linear_enable_changes_aw_covariance_not_mean",
        "first_accel_force_reads_aw_minus_gravity",
        "goLive_yaw_covariance_axis_is_body_gravity_axis",
        "zero_S_mean_makes_first_due_pseudo_mean_correction_zero",
    ):
        if s.get(k) is not True:
            failures.append(f"source semantic {k} is not true")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        failures.append("deployed correction range changed")
    if int(d.get("evaluated_source_phase_cells", 0)) <= 0:
        failures.append("no first-accel exact source cells evaluated")
    status = d.get("P5_FIRST_ACCEL_EXACT_SOURCE_CERTIFICATE")
    if status == "PASS":
        if d.get("all_first_accelerometer_source_cells_inside_validated_correction_range") is not True:
            failures.append("PASS without complete first-accelerometer source closure")
        if d.get("first_unclosed_child") is not None:
            failures.append("PASS retains an unclosed first-accelerometer child")
    elif status == "NOT_ESTABLISHED":
        if d.get("first_unclosed_child") is None:
            failures.append("nonclosure missing first source witness")
    else:
        failures.append("invalid first-accelerometer exact-source status")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--source-pieces", type=int, default=2)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    out = build(args.domain.resolve(), source_pieces=args.source_pieces)
    vf = validate(out)
    out["validation_pass"] = not vf
    out["validation_failures"] = vf
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": out["P5_FIRST_ACCEL_EXACT_SOURCE_CERTIFICATE"],
        "cells": out["evaluated_source_phase_cells"],
        "q_full": out["post_prediction_full_cayley_norm_upper"],
        "q_tangent": out["post_prediction_cayley_tangent_norm_upper"],
        "force": out["first_accel_specific_force_magnitude_mps2"],
        "max_K": out["max_Ktheta_norm_upper"],
        "max_d": out["max_first_accelerometer_correction_norm_upper_rad"],
        "margin": out["minimum_correction_range_margin_rad"],
        "over_limit": out["children_above_validated_correction_limit"],
        "first_unclosed": out["first_unclosed_child"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
