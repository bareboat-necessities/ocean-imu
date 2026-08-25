#!/usr/bin/env python3
"""Propagate the certified first accelerometer family through Joseph/reset to sample 1.

This stage starts from the source-reachable sample-zero family certified by
``ou3_p5_first_accel_exact_source_v2``.  It does not return to the broad V3
normal-Live initial hull.  For every first-prefix tau/sigma_aw/R_S/pseudo-phase
cell it:

* forms the first predicted 18x18 H covariance;
* installs the exact first-prefix gravity/yaw-axis covariance alignment already
  certified by the exact-source stage;
* applies the due S=0 covariance update with its exact zero estimator residual,
  or the not-due identity branch;
* forms the canonical first accelerometer Jacobian J_att=-[g e3]_x, J_aw=I;
* applies the shipping Joseph covariance update and immediate left-error reset;
* hulls the accepted/reset result with the conservative identity branch;
* propagates that state/covariance family through the next 5 ms prediction,
  producing the entry enclosure for sample 1.

The exact-source scalar correction certificate is used only to cap each
attitude-correction component by the already certified norm bound.  It is not
used to replace Joseph or reset covariance propagation.  The post-reset Cayley
chart is advanced one prediction with the same bounded gyro-bias/disturbance
transport law.  This producer certifies entry to sample 1 only; it does not yet
evaluate the sample-1 S/accelerometer/magnetometer packet and cannot set
N_H_words.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, hull, matrix_add, matrix_mul, matrix_transpose
import ou3_p5_first_accel_exact_source_v2 as FIRST
import ou3_p5_first_accel_post_reset as POST
import ou3_p5_first_accel_rotation_gauge as RG
import ou3_p5_full_h_prefix_cells as FULL
import ou3_p5_full_h_prefix_cells_v3 as FULL3
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = FIRST.DEFAULT_DOMAIN
SCHEMA = 1
N = FULL.N


def _canonical_first_H(gravity: float):
    H = FULL._zero(3, N)
    g = FULL.I(gravity)
    H[0][1] = g
    H[1][0] = -g
    for i in range(3):
        H[i][15 + i] = FULL.I(1.0)
    return H


def _canonicalize_first_attitude_covariance(Pm, domain_path: Path, h: float):
    tilt, yaw, eps = RG._attitude_covariance_epsilon(domain_path, h)
    v = [FULL.I(0.0), FULL.I(0.0), FULL.I(1.0)]
    Pt = RG._ptheta_cell(v, tilt, yaw, eps)
    out = [[Pm[i][j] for j in range(N)] for i in range(N)]
    for i in range(3):
        for j in range(3):
            out[i][j] = Pt[i][j]
    return FULL._psd_tighten(out)


def _zero_residual_S_covariance(Pm, src: dict):
    H = FULL._H_S()
    R = FULL._R_S(src)
    PHt, S = FULL._innovation(Pm, H, R)
    Sinv, backend = FULL._spd_inverse_enclosure(S, R)
    K = matrix_mul(PHt, Sinv)
    # The filter mean S is exactly zero on this first prefix, so the shipping
    # estimator correction is exactly zero even though the physical S error is
    # not.  Joseph still contracts covariance.
    Pj = FULL._shipping_joseph(Pm, K, S, PHt)
    return Pj, backend


def _first_accel_covariance_and_state(Pm, e, H, R, residual_norm: float, d_norm_cap: float):
    r = FULL._vec_box(residual_norm)
    PHt, S = FULL._innovation(Pm, H, R)
    Sinv, backend = FULL._spd_inverse_enclosure(S, R)
    K = matrix_mul(PHt, Sinv)
    dx = FULL._mat_vec(K, r)

    cap = Interval(-FULL.up(d_norm_cap), FULL.up(d_norm_cap))
    dx_theta = [FULL._intersect(dx[i], cap) for i in range(3)]
    dx_capped = list(dx)
    dx_capped[0:3] = dx_theta

    Pj = FULL._shipping_joseph(Pm, K, S, PHt)
    Pr = FULL._reset_covariance(Pj, dx_theta)

    e_acc = list(e)
    for i in range(3, N):
        e_acc[i] = e[i] - dx_capped[i]

    # Retain the identity branch as the same conservative source-complete hull
    # used by the established V3 backend.
    Pout = FULL._psd_tighten(FULL._mat_hull(Pm, Pr))
    eout = FULL._vec_hull(e, e_acc)
    return Pout, eout, {
        "inverse_backend": backend,
        "K": K,
        "S": S,
        "r": r,
        "dx": dx_capped,
        "P_accepted_reset": Pr,
    }


def _group_norm(e, idxs) -> float:
    return FULL._norm_upper([e[i] for i in idxs])


def build(domain_path: Path = DEFAULT_DOMAIN, *, source_pieces: int = 2) -> dict:
    FULL3._install_backend()
    domain_path = Path(domain_path).resolve()
    domain = json.loads(domain_path.read_text(encoding="utf-8"))
    if domain.get("trajectory_fit") is not False:
        raise RuntimeError("sample-1 entry domain must not be trajectory fitted")
    if domain.get("configured_runtime", {}).get("imu_lever_arm_enabled") is not False:
        raise RuntimeError("sample-1 entry stage requires lever arm disabled")

    first = FIRST.build(domain_path, source_pieces=source_pieces)
    post = POST.build(domain_path, source_pieces=source_pieces)
    vector = VECTOR.build()
    failures = [f"first-accel: {x}" for x in FIRST.validate(first)]
    failures += [f"post-reset: {x}" for x in POST.validate(post)]
    failures += [f"vector: {x}" for x in VECTOR.validate(vector)]
    if first.get("P5_FIRST_ACCEL_EXACT_SOURCE_CERTIFICATE") != "PASS":
        failures.append("first accelerometer prerequisite did not pass")
    if post.get("P5_FIRST_ACCEL_POST_RESET_PREFIX_CERTIFICATE") != "PASS":
        failures.append("first post-reset prerequisite did not pass")

    src_phases = RG._source_phase_children(source_pieces)
    first_rows = first.get("source_cells", [])
    if len(src_phases) != len(first_rows):
        failures.append("first source cell ordering/count mismatch")

    h = float(FULL._source_cell()["dt_s"])
    gravity = float(domain["startup"]["gravity_mps2"])
    Hacc = _canonical_first_H(gravity)
    vc = vector["configured_measurement_bounds"]
    Racc = FULL._R_diag(float(vc["acc_measurement_std_mps2"]))
    e0 = FULL._initial_error(domain)
    dcap = float(first["max_first_accelerometer_correction_norm_upper_rad"])
    qpost = float(post["post_accel_cayley_norm_upper"])
    q1 = RG._q_after_first_prediction(qpost, domain, h)
    q1_safe = math.isfinite(q1) and q1 < float(post["q_chart_target"])
    if not q1_safe:
        failures.append("sample-1 prediction leaves declared q<8 chart")

    backend_counts = {
        "FIXED_PIVOT_INTERVAL_GAUSS_JORDAN": 0,
        "SPD_S_GE_R_SPECTRAL_ENTRY_ENCLOSURE": 0,
    }
    rows = []
    first_failure = None
    max_diag = 0.0
    max_cross = 0.0
    max_state = {"gyro_bias": 0.0, "velocity": 0.0, "position": 0.0, "S": 0.0, "aw": 0.0}

    for si, item in enumerate(src_phases):
        src, phase = item
        try:
            P0 = FULL._initial_covariance(src, domain_path)
            F, Q, _Rstep = FULL._transition_and_Q(src, domain)
            Pp = FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, P0), matrix_transpose(F)), Q))
            Pp = _canonicalize_first_attitude_covariance(Pp, domain_path, h)
            ep = FULL._predict_error(e0, F)

            s_backend = "NOT_DUE_IDENTITY"
            if phase == "due":
                Ppre, s_backend = _zero_residual_S_covariance(Pp, src)
                backend_counts[s_backend] += 1
            else:
                Ppre = Pp

            source_row = first_rows[si]
            rho = float(source_row["combined_useful_residual_norm_upper_mps2"])
            Pa, ea, acell = _first_accel_covariance_and_state(Ppre, ep, Hacc, Racc, rho, dcap)
            backend_counts[acell["inverse_backend"]] += 1

            # This is the actual next IMU prediction: sample-1 entry before its
            # optional S/vector corrections.
            P1 = FULL._psd_tighten(matrix_add(matrix_mul(matrix_mul(F, Pa), matrix_transpose(F)), Q))
            e1 = FULL._predict_error(ea, F)
            summary = FULL._matrix_summary(P1)
            max_diag = max(max_diag, max(float(x[1]) for x in summary["diagonal_intervals"]))
            max_cross = max(max_cross, float(summary["max_offdiagonal_abs_upper"]))
            group = {
                "gyro_bias": _group_norm(e1, FULL.BG),
                "velocity": _group_norm(e1, FULL.V),
                "position": _group_norm(e1, FULL.P),
                "S": _group_norm(e1, FULL.SS),
                "aw": _group_norm(e1, FULL.AW),
            }
            for k, v in group.items():
                max_state[k] = max(max_state[k], v)

            rows.append({
                "source_phase_cell": si,
                "pseudo_phase_at_sample0": phase,
                "tau_s": src["tau_s"].as_list(),
                "sigma_aw_mps2": src["sigma_aw_mps2"].as_list(),
                "R_S_filter_std": src["R_S_filter_std"].as_list(),
                "sample0_S_inverse_backend": s_backend,
                "sample0_acc_inverse_backend": acell["inverse_backend"],
                "sample0_residual_norm_upper_mps2": rho,
                "sample0_attitude_dx_component_intervals": [x.as_list() for x in acell["dx"][0:3]],
                "sample1_entry_covariance": summary,
                "sample1_entry_state_group_norm_uppers": group,
            })
        except Exception as exc:
            first_failure = {
                "source_phase_cell": si,
                "pseudo_phase_at_sample0": phase,
                "reason": f"{type(exc).__name__}: {exc}",
            }
            break

    closed = bool(rows) and len(rows) == len(src_phases) and first_failure is None and q1_safe and not failures
    return {
        "schema": SCHEMA,
        "qualification": "OU3_P5_EXACT_SOURCE_SAMPLE1_ENTRY_AFTER_FIRST_ACCEL_JOSEPH_RESET",
        "source_generated_not_trajectory_fit": True,
        "source_replay_used": False,
        "filter_changed": False,
        "starts_from_exact_source_sample0_family_not_broad_v3_initial_hull": True,
        "full_18x18_covariance_propagated": True,
        "shipping_Joseph_update_used_for_first_accelerometer": True,
        "immediate_left_error_reset_congruence_used": True,
        "first_due_S_zero_mean_but_covariance_update_retained": True,
        "first_accel_canonical_Jatt_gravity_and_Jaw_identity_used": True,
        "first_accel_exact_norm_cap_used_only_to_intersect_attitude_dx_components": True,
        "source_complete_identity_branch_hull_retained": True,
        "sample1_entry_is_before_sample1_measurements": True,
        "deployed_correction_limit_rad": float(first["deployed_correction_limit_rad"]),
        "deployed_correction_limit_increased": False,
        "sample0_post_reset_cayley_norm_upper": qpost,
        "sample1_pre_measurement_cayley_norm_upper": q1,
        "q_chart_target": float(post["q_chart_target"]),
        "sample1_entry_inside_q8": q1_safe,
        "evaluated_source_phase_cells": len(rows),
        "expected_source_phase_cells": len(src_phases),
        "inverse_backend_counts": backend_counts,
        "max_sample1_covariance_diagonal_upper": max_diag,
        "max_sample1_covariance_offdiagonal_abs_upper": max_cross,
        "max_sample1_state_group_norm_uppers": max_state,
        "source_cells": rows,
        "first_failure": first_failure,
        "sample1_S_accel_mag_prefix_evaluated_here": False,
        "whole_word_promoted_here": False,
        "N_H_words_set_here": False,
        "P5_SAMPLE1_ENTRY_CERTIFICATE": "PASS" if closed else "NOT_ESTABLISHED",
        "next_obligation": (
            "CLASSIFY_SAMPLE1_PSEUDO_PHASE_FROM_SAMPLE0_PHASE_AND_EVALUATE_SAMPLE1_S_ACCEL_PREFIX"
            if closed else
            "REFINE_SAMPLE0_JOSEPH_RESET_DIRECTIONAL_CELL_BEFORE_SAMPLE1_PREDICTION"
        ),
        "failures": failures,
    }


def validate(d: dict) -> list[str]:
    failures = list(d.get("failures", []))
    if d.get("schema") != SCHEMA:
        failures.append("schema mismatch")
    for k in (
        "source_generated_not_trajectory_fit",
        "starts_from_exact_source_sample0_family_not_broad_v3_initial_hull",
        "full_18x18_covariance_propagated",
        "shipping_Joseph_update_used_for_first_accelerometer",
        "immediate_left_error_reset_congruence_used",
        "first_due_S_zero_mean_but_covariance_update_retained",
        "first_accel_canonical_Jatt_gravity_and_Jaw_identity_used",
        "first_accel_exact_norm_cap_used_only_to_intersect_attitude_dx_components",
        "source_complete_identity_branch_hull_retained",
        "sample1_entry_is_before_sample1_measurements",
        "sample1_entry_inside_q8",
    ):
        if d.get(k) is not True:
            failures.append(f"{k} is not true")
    for k in (
        "source_replay_used", "filter_changed", "deployed_correction_limit_increased",
        "sample1_S_accel_mag_prefix_evaluated_here", "whole_word_promoted_here", "N_H_words_set_here",
    ):
        if d.get(k) is not False:
            failures.append(f"{k} is not false")
    if float(d.get("deployed_correction_limit_rad", 0.0)) != 6.0:
        failures.append("deployed correction limit changed")
    if int(d.get("evaluated_source_phase_cells", 0)) != int(d.get("expected_source_phase_cells", -1)):
        failures.append("not all source/phase cells reached sample 1")
    if d.get("first_failure") is not None:
        failures.append("sample-1 entry retains a source-cell failure")
    q1 = d.get("sample1_pre_measurement_cayley_norm_upper")
    if not isinstance(q1, (int, float)) or not math.isfinite(float(q1)) or float(q1) >= 8.0:
        failures.append("sample-1 pre-measurement Cayley bound is not finite below 8")
    if d.get("P5_SAMPLE1_ENTRY_CERTIFICATE") != "PASS" and not failures:
        failures.append("sample-1 entry stage did not pass")
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
        "status": out["P5_SAMPLE1_ENTRY_CERTIFICATE"],
        "cells": out["evaluated_source_phase_cells"],
        "q0_post_reset": out["sample0_post_reset_cayley_norm_upper"],
        "q1_pre_measurement": out["sample1_pre_measurement_cayley_norm_upper"],
        "inverse_backends": out["inverse_backend_counts"],
        "max_cov_diag": out["max_sample1_covariance_diagonal_upper"],
        "max_cov_cross": out["max_sample1_covariance_offdiagonal_abs_upper"],
        "max_state": out["max_sample1_state_group_norm_uppers"],
        "first_failure": out["first_failure"],
        "next": out["next_obligation"],
        "validation_failures": vf,
    }, indent=2, sort_keys=True))
    return 0 if not vf else 2


if __name__ == "__main__":
    raise SystemExit(main())
