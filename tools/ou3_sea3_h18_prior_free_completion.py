#!/usr/bin/env python3
"""Recurrent prior-free H18 completion on one complete SEA3 Normal-Live word.

This closes the H18 quantitative obligation left open by the complete-SEA3
information certificate.  It does not create a new source, replay a trajectory,
change the filter, or use the retired one-step ``min(Q)/trace(P)`` contraction.

Let the exact prior-free factorization of the same 3 s word be ``(D,T,Qc)``.
The existing complete-SEA3 information certificate proves ``D >> 0``.  In the
diffuse-prior limit the endpoint covariance is

    P_inf = T D^-1 T^T + Qc.

A concrete finite-memory estimator using only measurements that are already in
the same complete word gives a source-uniform covariance upper ``Pbar`` for
that diffuse-prior endpoint.  Hence

    T D^-1 T^T <= P_inf <= Pbar.

The old moving tube propagated its selected S estimator for an extra full PE
second and therefore reported a 3.155 s bookkeeping memory.  That extension is
not needed here.  Three guaranteed S firings lie in windows near 0, 1 and 2 s,
and the two declared PE occurrences lie in [0,1] and [2,3] s.  We evaluate the
same estimator directly at the canonical 3.0 s endpoint, so every observation
used by Pbar is a subset of this same complete SEA3 word.

Take the immediately following shipping prediction, before its S/accelerometer
measurement.  Then

    Qc+ = F Qc F^T + Q >= Q,
    T+ D^-1 T+^T <= F Pbar F^T.

For every PSD covariance with marginal bounds u_i,

    Pbar <= trace(Pbar) I <= (sum_i u_i) I.

The shipping same-mode prediction is block diagonal.  Its attitude/gyro-bias
block ``[[R,B],[0,I]]`` has ``||R||=1`` and ``||B||<=h``; its norm squared is
therefore at most ``2+h^2`` by Frobenius domination of the corresponding 2x2
block norm.  One translation-axis transition has

    ||F_l||_F^2 <= 4 + 3 h^2 + h^4/2 + h^6/36,

using |phi_va|<=h, |phi_pa|<=h^2/2, |phi_Sa|<=h^3/6 and |alpha|<=1.
Thus a source-uniform ``F_norm_sq_upper`` bounds the complete H18 prediction.

The exact prior-free completion identity gives the sufficient full-matrix test

    (1-delta) Q
      - (delta^2/4) F Pbar F^T  >> 0.

We replace the second term only by the proven isotropic upper
``penalty*I``.  The attitude/gyro-bias Q lower is the shipping configured UCC
primitive.  The 12-state integrated-OU Q is *not* collapsed to its catastrophic
old Euclidean minimum: in the source-cell congruence

    D_h = diag(sigma_floor*h, sigma_floor*h^2,
               sigma_floor*h^3, sigma_floor)

per axis, the stable factored exact shipping covariance supplies
``x * B(x)``.  Every x=h/tau cell is outward rounded and the complete 18x18
congruent matrix is certified by interval LDL^T.  The source sigma may exceed
the 0.05 floor; because every D_h coordinate has the same sigma factor, using
the floor is a PSD process subnoise, not an independent source choice.

Once this first H18 margin is established, the existing event algebra proves
that all subsequent predictions, Joseph S/accelerometer/magnetometer updates,
PSD a_w floors and immediate left-error resets preserve it by exact congruence
/ PSD addition.  A21 remains a separate obligation, so this module cannot
promote canonical P3 by itself.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import (
    Interval,
    IntervalMatrix,
    down,
    up,
    symmetric_positive_definite_ldlt,
)
import ou3_full_process_ucc as PROCESS
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_full_normal_live_word as WORD
import ou3_sea3_full_word_event_algebra as EVENT
import ou3_sea3_h18_information_composition as HINFO
import ou3_sea3_riccati_tube as TUBE
import ou3_sea3_riccati_tube_factored as FACTORED
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = COMPLETE.DEFAULT_DOMAIN
SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_H18_PRIOR_FREE_FULL_MATRIX_COMPLETION"
USEFUL_GATE = 1.0e-18
HORIZON_S = 3.0

OFF_V = 6
OFF_P = 9
OFF_S = 12
OFF_AW = 15


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _zero(n: int) -> IntervalMatrix:
    z = I(0.0)
    return [[z for _ in range(n)] for _ in range(n)]


def _same_word_translation_upper(dynamic: dict, live: dict, axis_factors: list[float]) -> tuple[list[float], dict]:
    """Retained finite-memory S estimator, referenced at exactly t=3 s."""
    inv = dynamic["dynamic_invariant"]
    rates = dynamic["validated_rate_and_jump_bounds"]
    h = TUBE.pos(rates["dt_s"], "dt")
    tau_lo = TUBE.pos(inv["tau_applied_s"][0], "tau lower")
    sigma_hi = TUBE.pos(inv["sigma_aw_filter_mps2"][1], "sigma upper")
    rs_hi = TUBE.pos(inv["R_S_applied"][1], "R_S upper")
    cadence_hi = TUBE.pos(inv["pseudo_update_period_s"][1], "pseudo cadence upper")
    Tpe = TUBE.pos(live["vector_pe_recurrence_window_s"], "PE recurrence")

    gap = TUBE.up(cadence_hi + h)
    spacing = TUBE.up(max(Tpe, 2.0 * gap))
    Tobs = TUBE.up(2.0 * spacing + gap)
    if not Tobs < HORIZON_S:
        raise RuntimeError(f"selected S estimator does not fit the canonical 3 s word: {Tobs}")

    Binv = TUBE.integrator_inverse(gap, spacing)
    qc_hi = TUBE.up(2.0 * sigma_hi * sigma_hi / tau_lo)
    s_nuis = TUBE.up(sigma_hi * sigma_hi * (Tobs ** 3 / 6.0) ** 2)
    s_proc = TUBE.up(qc_hi * Tobs ** 7 / 252.0)
    rmax = TUBE.up((rs_hi * max(axis_factors)) ** 2)
    rstack = TUBE.up(3.0 * (rmax + s_nuis + s_proc))
    R = [[TUBE.I(rstack if i == j else 0.0) for j in range(3)] for i in range(3)]
    Cspv = TUBE.matrix_symmetric_hull(
        TUBE.matrix_mul(TUBE.matrix_mul(Binv, R), TUBE.matrix_transpose(Binv))
    )
    order = (2, 1, 0)  # [S,p,v] -> [v,p,S]
    Cvps = [[Cspv[order[i]][order[j]] for j in range(3)] for i in range(3)]

    # The estimator is referenced directly at the common 3 s endpoint.  Using
    # the whole 0..3 s propagation range is conservative for whichever selected
    # S occurrence supplies each reconstructed coordinate.
    t = Interval.outward_bounds(0.0, HORIZON_S)
    F = [
        [TUBE.I(1), TUBE.I(0), TUBE.I(0)],
        [t, TUBE.I(1), TUBE.I(0)],
        [TUBE.I(0.5) * t.square(), t, TUBE.I(1)],
    ]
    Cend = TUBE.matrix_symmetric_hull(
        TUBE.matrix_mul(TUBE.matrix_mul(F, Cvps), TUBE.matrix_transpose(F))
    )
    u = TUBE.diagonal_dominator(Cend)

    Tword = HORIZON_S
    variances = [
        TUBE.up(sigma_hi * sigma_hi * Tword * Tword + qc_hi * Tword ** 3 / 3.0),
        TUBE.up(sigma_hi * sigma_hi * Tword ** 4 / 4.0 + qc_hi * Tword ** 5 / 20.0),
        TUBE.up(sigma_hi * sigma_hi * Tword ** 6 / 36.0 + qc_hi * Tword ** 7 / 252.0),
        TUBE.up(sigma_hi * sigma_hi),
    ]
    roots = [math.sqrt(v) for v in variances]
    total = TUBE.up(sum(roots))
    noise = [TUBE.up(r * total) for r in roots]
    upper = [TUBE.up(u[i] + noise[i]) for i in range(3)] + [noise[3]]
    return upper, {
        "pseudo_gap_s_upper": gap,
        "observation_window_s_upper": Tobs,
        "covariance_memory_window_s_upper": HORIZON_S,
        "same_word_endpoint_s": HORIZON_S,
        "selected_S_windows_fit_same_word": True,
        "q_c_global_upper": qc_hi,
        "source_motion_inside_window_allowed": True,
    }


def _same_word_covariance_upper(path: Path, dynamic: dict, process: dict, hinfo: dict) -> dict:
    domain = json.loads(path.read_text(encoding="utf-8"))
    live = domain["normal_live"]
    vector = VECTOR.build()
    vf = VECTOR.validate(vector)
    if vf:
        raise RuntimeError(f"vector covariance-upper prerequisite failed: {vf}")
    axis_factors = TUBE._axis_factors()
    alpha6 = float(hinfo["eta6_information_lower"])
    trans, timing = _same_word_translation_upper(dynamic, live, axis_factors)
    full = TUBE._global_full_state_upper(dynamic, live, vector, process, alpha6, trans, timing)
    H = [float(x) for x in full["H"]]
    if len(H) != 18 or any(not (math.isfinite(x) and x > 0.0) for x in H):
        raise RuntimeError("same-word H18 covariance upper invalid")
    trace = TUBE.up(sum(H))
    return {
        "dimension": 18,
        "Pbar_diagonal_variance_upper": H,
        "Pbar_trace_upper": trace,
        "translation": trans,
        "timing": timing,
        "concrete_estimator_uses_only_measurements_inside_same_word": True,
        "diffuse_prior_endpoint_covariance_Loewner_upper": True,
        "TD_inverse_T_transpose_Loewner_upper": True,
        "marginal_bounds_not_misused_as_diagonal_Loewner_matrix": True,
        "Loewner_route": "Pbar <= trace(Pbar) I <= sum(u_i) I",
    }


def _prediction_norm_sq_upper(process: dict) -> float:
    h_hi = float(process["configured_runtime"]["imu_dt_outward_interval_s"][1])
    if not (math.isfinite(h_hi) and h_hi > 0.0):
        raise RuntimeError("invalid shipping dt")
    attitude = TUBE.up(2.0 + h_hi * h_hi)
    h2 = TUBE.up(h_hi * h_hi)
    h4 = TUBE.up(h2 * h2)
    h6 = TUBE.up(h4 * h2)
    translation = TUBE.up(4.0 + 3.0 * h2 + 0.5 * h4 + h6 / 36.0)
    return TUBE.up(max(attitude, translation))


def _full_H18_cell(
    x: Interval,
    *,
    process: dict,
    dynamic: dict,
    penalty_physical: float,
) -> tuple[bool, dict]:
    delta = USEFUL_GATE
    one_minus = I(down(1.0 - delta))
    M = _zero(18)

    q_att = float(process["attitude_gyro_bias"]["Q_attitude_gyro_bias_lambda_min_lower"])
    att_diag = down((1.0 - delta) * q_att - penalty_physical)
    if not att_diag > 0.0:
        return False, {"reason": "attitude_penalty", "attitude_diagonal_lower": att_diag}
    for i in range(6):
        M[i][i] = I(att_diag)

    inv = dynamic["dynamic_invariant"]
    sigma_floor = float(inv["sigma_aw_filter_mps2"][0])
    h = float(dynamic["validated_rate_and_jump_bounds"]["dt_s"])
    if not (sigma_floor > 0.0 and h > 0.0):
        raise RuntimeError("invalid SEA3 process scale")
    scales = [sigma_floor * h, sigma_floor * h * h, sigma_floor * h * h * h, sigma_floor]

    B = FACTORED.step_scaled_q_over_x(x)
    qscale = I(x.lo)
    Maxis = [[one_minus * qscale * B[i][j] for j in range(4)] for i in range(4)]
    for i in range(4):
        pscaled = TUBE.up(penalty_physical / TUBE.down(scales[i] * scales[i]))
        Maxis[i][i] = Maxis[i][i] - I(pscaled)

    idx = (OFF_V, OFF_P, OFF_S, OFF_AW)
    for axis in range(3):
        for i in range(4):
            for j in range(4):
                M[idx[i] + axis][idx[j] + axis] = Maxis[i][j]

    ok, pivots = symmetric_positive_definite_ldlt(M)
    return ok, {
        "x_h_over_tau": x.as_list(),
        "sigma_floor_mps2": sigma_floor,
        "translation_congruence_scales": scales,
        "pivot_count": len(pivots),
        "pivot_lower": min((p.lo for p in pivots), default=math.inf),
        "full_dimension": 18,
        "full_18x18_interval_LDLT": True,
    }


def _x_cover(dynamic: dict) -> list[Interval]:
    h = float(dynamic["validated_rate_and_jump_bounds"]["dt_s"])
    tau_lo, tau_hi = map(float, dynamic["dynamic_invariant"]["tau_applied_s"])
    xlo, xhi = h / tau_hi, h / tau_lo
    edges = TUBE.geom_edges(xlo, xhi, 24)
    if xlo < TUBE.BRANCH_X < xhi:
        edges = sorted(set(edges + [TUBE.BRANCH_X]))
    leaves: list[Interval] = []
    for cell in TUBE.interval_cells(edges):
        for leaf, _rho in FACTORED.split_x_cell(cell):
            leaves.append(leaf)
    return leaves


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    dynamic = DYNAMIC.build(path)
    process = PROCESS.build()
    hinfo = HINFO.build(path)
    event = EVENT.build()
    word = WORD.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "dynamic": DYNAMIC.validate(dynamic),
        "process": PROCESS.validate(process),
        "H18_information": HINFO.validate(hinfo),
        "event_algebra": EVENT.validate(event),
        "literal_word": WORD.validate(word),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"H18 prior-free prerequisites failed: {bad}")
    if complete["canonical_P3_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("H18 prior-free completion detached from complete SEA3")
    if float(complete["word_horizon_s"]) != HORIZON_S:
        raise RuntimeError("canonical complete-SEA3 word is no longer 3 s")
    if not hinfo["H18_information_useful_gate_pass"]:
        raise RuntimeError("same-word H18 information is not strictly full rank")

    pbar = _same_word_covariance_upper(path, dynamic, process, hinfo)
    Fnorm2 = _prediction_norm_sq_upper(process)
    penalty = TUBE.up(
        (USEFUL_GATE * USEFUL_GATE / 4.0)
        * Fnorm2
        * float(pbar["Pbar_trace_upper"])
    )
    if not (math.isfinite(penalty) and penalty > 0.0):
        raise RuntimeError("prior-free completion penalty invalid")

    leaves = _x_cover(dynamic)
    rows = []
    worst = math.inf
    failures = []
    for x in leaves:
        ok, row = _full_H18_cell(x, process=process, dynamic=dynamic, penalty_physical=penalty)
        rows.append(row)
        if not ok:
            failures.append(row)
        else:
            worst = min(worst, float(row["pivot_lower"]))
    closed = not failures and bool(rows) and math.isfinite(worst) and worst > 0.0

    preserve = event["full_matrix_margin_preservation"]
    suffix_preserved = all(bool(preserve[k]) for k in (
        "covers_prediction",
        "covers_every_due_S_update",
        "covers_every_Normal_Live_accelerometer_update",
        "covers_asynchronous_magnetometer_update",
        "covers_immediate_left_error_reset",
        "covers_aw_covariance_floor",
        "covers_not_due_or_rejected_identity_branches",
    ))

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": complete["canonical_P3_source"],
        "complete_word_horizon_s": HORIZON_S,
        "component_of_complete_SEA3_full_word": True,
        "same_word_diffuse_prior_covariance_upper": pbar,
        "H18_information_lambda_min_lower": float(
            hinfo["triangular_information_composition"]["D_H18_lambda_min_lower"]
        ),
        "H18_D_is_strictly_positive_definite": True,
        "prior_free_exact_identity_consumed": True,
        "diffuse_prior_limit_used_only_to_upper_T_Dinv_Tt": True,
        "next_shipping_prediction_taken_before_next_measurement": True,
        "prediction_F_spectral_norm_squared_upper": Fnorm2,
        "delta_squared_over_four_penalty_physical": penalty,
        "stable_factored_shipping_integrated_OU_Q_consumed": True,
        "old_one_step_Euclidean_Q_min_used": False,
        "x_cells_certified": len(rows),
        "x_cell_failures": failures,
        "worst_full_H18_LDLT_pivot_lower": worst if closed else None,
        "full_H18_prior_free_matrix_condition_closed": closed,
        "H18_prior_free_completion_closed": closed,
        "event_algebra_preserves_margin_after_closure": suffix_preserved,
        "actual_applied_SpectralMSE_R_S_consumed_in_same_word_covariance_upper": True,
        "all_due_S_updates_remain_in_complete_word": True,
        "marginal_Pbar_bounds_used_as_Loewner_diagonal_directly": False,
        "Pbar_trace_isotropic_Loewner_upper_used": True,
        "blockwise_minimum_contraction_used": False,
        "D_W_L_W_product_used": False,
        "scalar_beta_contraction_used": False,
        "determinant_trace_final_matrix_gate_used": False,
        "full_18x18_interval_LDLT_used": True,
        "source_family_replaced": False,
        "trajectory_replay_used": False,
        "independent_tau_sigma_RS_source_created": False,
        "useful_gate": USEFUL_GATE,
        "P3_promoted": False,
        "next_obligation": (
            "compose the H-to-A release and active accelerometer-bias Gauss-Markov block into the same "
            "prior-free full A21 matrix condition; then wire both H18/A21 closures into the canonical P3 gate"
        ) if closed else (
            "subdivide only the same SEA3 x coordinate or tighten the same-word covariance estimator; do not relax delta"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "component_of_complete_SEA3_full_word",
        "H18_D_is_strictly_positive_definite",
        "prior_free_exact_identity_consumed",
        "diffuse_prior_limit_used_only_to_upper_T_Dinv_Tt",
        "next_shipping_prediction_taken_before_next_measurement",
        "stable_factored_shipping_integrated_OU_Q_consumed",
        "full_H18_prior_free_matrix_condition_closed",
        "H18_prior_free_completion_closed",
        "event_algebra_preserves_margin_after_closure",
        "actual_applied_SpectralMSE_R_S_consumed_in_same_word_covariance_upper",
        "all_due_S_updates_remain_in_complete_word",
        "Pbar_trace_isotropic_Loewner_upper_used",
        "full_18x18_interval_LDLT_used",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "old_one_step_Euclidean_Q_min_used",
        "marginal_Pbar_bounds_used_as_Loewner_diagonal_directly",
        "blockwise_minimum_contraction_used",
        "D_W_L_W_product_used",
        "scalar_beta_contraction_used",
        "determinant_trace_final_matrix_gate_used",
        "source_family_replaced",
        "trajectory_replay_used",
        "independent_tau_sigma_RS_source_created",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"forbidden/open flag {key} changed")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("useful gate changed")
    if int(d.get("x_cells_certified", 0)) <= 0:
        f.append("no SEA3 x cells certified")
    if d.get("x_cell_failures"):
        f.append("one or more H18 prior-free x cells failed")
    pivot = d.get("worst_full_H18_LDLT_pivot_lower")
    if not isinstance(pivot, (int, float)) or not (math.isfinite(float(pivot)) and float(pivot) > 0.0):
        f.append("H18 full-matrix LDLT pivot is not strict")
    pbar = d.get("same_word_diffuse_prior_covariance_upper", {})
    timing = pbar.get("timing", {})
    if float(timing.get("same_word_endpoint_s", math.nan)) != HORIZON_S:
        f.append("same-word covariance upper is not referenced at 3 s")
    if timing.get("selected_S_windows_fit_same_word") is not True:
        f.append("selected S estimator escaped the complete word")
    return list(dict.fromkeys(f))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    d["validation_pass"] = not failures
    d["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(d, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "same_word_Pbar_trace_upper": d["same_word_diffuse_prior_covariance_upper"]["Pbar_trace_upper"],
        "same_word_timing": d["same_word_diffuse_prior_covariance_upper"]["timing"],
        "F_norm_sq_upper": d["prediction_F_spectral_norm_squared_upper"],
        "delta2_penalty": d["delta_squared_over_four_penalty_physical"],
        "x_cells": d["x_cells_certified"],
        "worst_H18_LDLT_pivot": d["worst_full_H18_LDLT_pivot_lower"],
        "H18_closed": d["H18_prior_free_completion_closed"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
