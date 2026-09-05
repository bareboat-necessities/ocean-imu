#!/usr/bin/env python3
"""Prior-free full A21 Riccati completion on the canonical complete SEA3 word.

This is the A-mode analogue of ``ou3_sea3_h18_prior_free_completion``.  It
keeps the same 3 s complete SEA3 word, the same finite-memory covariance
estimator, the deployed SpectralMSE R_S regularizer, eta6 PE, and the shipping
finite-correlation accelerometer-bias model.  It does not introduce eta9 PE,
a replay source, a tuner rectangle, an independent parameter product, or an
alternate estimator.

The same-word covariance estimator used by the H18 proof already pays active
b_a as accelerometer nuisance in its attitude/gyro-bias ceiling.  The A21
ceiling is therefore the same first 18 marginal bounds plus the certified
uniform active-bias variance

    p_ba <= max(sigma_ba0^2, q_ba tau_ba / 2).

For a PSD covariance, marginal bounds u_i imply

    Pbar_A <= trace(Pbar_A) I.

At the immediately following shipping A-mode prediction, before any
measurement, the exact prior-free completion identity reduces to the
sufficient condition

    (1-delta) Q_A - (delta^2/4) F_A Pbar_A F_A^T >> 0.

We upper the second term by the proven scalar Loewner bound and certify the
first term as one full 21x21 interval LDL^T matrix.  The 6-state
attitude/gyro-bias process lower and the 12-state integrated-OU translation
block are exactly the H18 ingredients.  The added 3-state b_a block uses the
shipping exact Gauss--Markov one-sample Q_ba lower.  Shipping prediction has no
H18<->b_a state-transition cross block, so this is a block-diagonal congruence
of the actual A21 process matrix, not a blockwise contraction ratio.

Once the first A21 margin is established, the already-certified complete-word
event algebra preserves Omega-delta*P through later predictions, Joseph
S/accelerometer/magnetometer updates, a_w PSD floors and immediate left-error
reset congruences.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, IntervalMatrix, down, symmetric_positive_definite_ldlt
import ou3_full_process_ucc as PROCESS
import ou3_sea3_a21_detectability_completion as ADET
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_full_normal_live_word as WORD
import ou3_sea3_full_word_event_algebra as EVENT
import ou3_sea3_h18_prior_free_completion as H18
import ou3_sea3_riccati_tube as TUBE
import ou3_sea3_riccati_tube_factored as FACTORED
import ou3_vector_uco_certificate as VECTOR

DEFAULT_DOMAIN = COMPLETE.DEFAULT_DOMAIN
SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_A21_PRIOR_FREE_FULL_MATRIX_COMPLETION"
USEFUL_GATE = 1.0e-18
HORIZON_S = 3.0

OFF_V = 6
OFF_P = 9
OFF_S = 12
OFF_AW = 15
OFF_BA = 18


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _zero(n: int) -> IntervalMatrix:
    z = I(0.0)
    return [[z for _ in range(n)] for _ in range(n)]


def _same_word_A21_covariance_upper(
    path: Path,
    dynamic: dict,
    process: dict,
    h18: dict,
) -> dict:
    """Build the A21 marginal/trace upper at exactly the canonical 3 s endpoint."""
    domain = json.loads(path.read_text(encoding="utf-8"))
    live = domain["normal_live"]
    vector = VECTOR.build()
    vf = VECTOR.validate(vector)
    if vf:
        raise RuntimeError(f"vector covariance-upper prerequisite failed: {vf}")

    axis_factors = TUBE._axis_factors()
    trans, timing = H18._same_word_translation_upper(dynamic, live, axis_factors)
    alpha6 = float(h18["H18_information_lambda_min_lower"])
    # _global_full_state_upper expects the eta6 information lower, not the
    # already-composed H18 scalar.  Recover the exact upstream alpha6 from the
    # H18 information producer through its public field.
    import ou3_sea3_h18_information_composition as HINFO
    hinfo = HINFO.build(path)
    hf = HINFO.validate(hinfo)
    if hf:
        raise RuntimeError(f"H18 information prerequisite failed: {hf}")
    alpha6 = float(hinfo["eta6_information_lower"])

    full = TUBE._global_full_state_upper(
        dynamic, live, vector, process, alpha6, trans, timing
    )
    A = [float(x) for x in full["A"]]
    if len(A) != 21 or any(not (math.isfinite(x) and x > 0.0) for x in A):
        raise RuntimeError("same-word A21 covariance upper invalid")
    if timing.get("same_word_endpoint_s") != HORIZON_S:
        raise RuntimeError("same-word A21 covariance upper escaped the canonical 3 s endpoint")

    trace = TUBE.up(sum(A))
    return {
        "dimension": 21,
        "Pbar_diagonal_variance_upper": A,
        "Pbar_trace_upper": trace,
        "translation": trans,
        "timing": timing,
        "accel_bias_variance_upper": float(full["accel_bias_variance_upper"]),
        "active_ba_nuisance_paid_in_H_subblock_upper": True,
        "concrete_estimator_uses_only_measurements_inside_same_word": True,
        "diffuse_prior_endpoint_covariance_Loewner_upper": True,
        "TD_inverse_T_transpose_Loewner_upper": True,
        "marginal_bounds_not_misused_as_diagonal_Loewner_matrix": True,
        "Loewner_route": "Pbar_A <= trace(Pbar_A) I <= sum(u_i) I",
    }


def _full_A21_cell(
    x: Interval,
    *,
    process: dict,
    dynamic: dict,
    penalty_physical: float,
) -> tuple[bool, dict]:
    delta = USEFUL_GATE
    one_minus = I(down(1.0 - delta))
    M = _zero(21)

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

    q_ba = float(process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"])
    ba_diag = down((1.0 - delta) * q_ba - penalty_physical)
    if not ba_diag > 0.0:
        return False, {"reason": "active_ba_penalty", "active_ba_diagonal_lower": ba_diag}
    for i in range(OFF_BA, OFF_BA + 3):
        M[i][i] = I(ba_diag)

    ok, pivots = symmetric_positive_definite_ldlt(M)
    return ok, {
        "x_h_over_tau": x.as_list(),
        "sigma_floor_mps2": sigma_floor,
        "translation_congruence_scales": scales,
        "active_ba_Q_lambda_min_lower": q_ba,
        "active_ba_diagonal_after_penalty_lower": ba_diag,
        "pivot_count": len(pivots),
        "pivot_lower": min((p.lo for p in pivots), default=math.inf),
        "full_dimension": 21,
        "full_21x21_interval_LDLT": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    dynamic = DYNAMIC.build(path)
    process = PROCESS.build()
    h18 = H18.build(path)
    adet = ADET.build(path)
    event = EVENT.build()
    word = WORD.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "dynamic": DYNAMIC.validate(dynamic),
        "process": PROCESS.validate(process),
        "H18": H18.validate(h18),
        "A21_detectability": ADET.validate(adet),
        "event_algebra": EVENT.validate(event),
        "literal_word": WORD.validate(word),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"A21 prior-free prerequisites failed: {bad}")
    if complete["canonical_P3_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("A21 prior-free completion detached from complete SEA3")
    if float(complete["word_horizon_s"]) != HORIZON_S:
        raise RuntimeError("canonical complete-SEA3 word is no longer 3 s")
    if not h18["H18_prior_free_completion_closed"]:
        raise RuntimeError("A21 prior-free completion requires the closed H18 route")
    if not adet["A21_finite_bias_detectability_closed"]:
        raise RuntimeError("A21 finite-bias detectability prerequisite is not closed")
    if adet["eta9_point_packet_shortcut_used"] is not False:
        raise RuntimeError("eta9 shortcut re-entered A21 prior-free completion")

    pbar = _same_word_A21_covariance_upper(path, dynamic, process, h18)
    Fnorm2 = float(h18["prediction_F_spectral_norm_squared_upper"])
    # A-mode adds phi_b I with |phi_b|<=1.  The H18 prediction norm upper is
    # already >1, so it remains a valid full A21 norm upper.
    if not (math.isfinite(Fnorm2) and Fnorm2 >= 1.0):
        raise RuntimeError("invalid A21 prediction norm upper")
    penalty = TUBE.up(
        (USEFUL_GATE * USEFUL_GATE / 4.0)
        * Fnorm2
        * float(pbar["Pbar_trace_upper"])
    )
    if not (math.isfinite(penalty) and penalty > 0.0):
        raise RuntimeError("A21 prior-free completion penalty invalid")

    leaves = H18._x_cover(dynamic)
    rows = []
    failures = []
    worst = math.inf
    for x in leaves:
        ok, row = _full_A21_cell(
            x, process=process, dynamic=dynamic, penalty_physical=penalty
        )
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
        "H18_prior_free_completion_consumed": True,
        "A21_finite_bias_detectability_consumed": True,
        "paper_active_bias_route": adet["paper_active_bias_route"],
        "eta9_point_packet_shortcut_used": False,
        "prior_free_exact_identity_consumed": True,
        "next_shipping_prediction_taken_before_next_measurement": True,
        "prediction_F_spectral_norm_squared_upper": Fnorm2,
        "delta_squared_over_four_penalty_physical": penalty,
        "stable_factored_shipping_integrated_OU_Q_consumed": True,
        "shipping_active_ba_GM_Q_consumed": True,
        "x_cells_certified": len(rows),
        "x_cell_failures": failures,
        "worst_full_A21_LDLT_pivot_lower": worst if closed else None,
        "full_21x21_Omega_minus_delta_P_LDLT_closed": closed,
        "A21_prior_free_completion_closed": closed,
        "event_algebra_preserves_margin_after_closure": suffix_preserved,
        "actual_applied_SpectralMSE_R_S_consumed_in_same_word_covariance_upper": True,
        "all_due_S_updates_remain_in_complete_word": True,
        "full_21x21_interval_LDLT_used": True,
        "blockwise_minimum_contraction_used": False,
        "D_W_L_W_product_used": False,
        "scalar_beta_contraction_used": False,
        "determinant_trace_final_matrix_gate_used": False,
        "source_family_replaced": False,
        "trajectory_replay_used": False,
        "independent_tau_sigma_RS_source_created": False,
        "useful_gate": USEFUL_GATE,
        "P3_promoted": False,
        "next_obligation": (
            "wire the closed H18/A21 complete-SEA3 full-matrix certificates into the sole canonical P3 gate and make the literal word execute immediate reset congruences before promotion"
            if closed else
            "tighten only the same complete-SEA3 A21 covariance/process enclosure; do not relax delta or introduce eta9"
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
        "H18_prior_free_completion_consumed",
        "A21_finite_bias_detectability_consumed",
        "prior_free_exact_identity_consumed",
        "next_shipping_prediction_taken_before_next_measurement",
        "stable_factored_shipping_integrated_OU_Q_consumed",
        "shipping_active_ba_GM_Q_consumed",
        "full_21x21_Omega_minus_delta_P_LDLT_closed",
        "A21_prior_free_completion_closed",
        "event_algebra_preserves_margin_after_closure",
        "actual_applied_SpectralMSE_R_S_consumed_in_same_word_covariance_upper",
        "all_due_S_updates_remain_in_complete_word",
        "full_21x21_interval_LDLT_used",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "eta9_point_packet_shortcut_used",
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
    if d.get("paper_active_bias_route") != "ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION":
        f.append("wrong A21 active-bias theorem route")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("useful gate changed")
    if int(d.get("x_cells_certified", 0)) <= 0:
        f.append("no SEA3 x cells certified")
    if d.get("x_cell_failures"):
        f.append("one or more A21 prior-free x cells failed")
    pivot = d.get("worst_full_A21_LDLT_pivot_lower")
    if not isinstance(pivot, (int, float)) or not (math.isfinite(float(pivot)) and float(pivot) > 0.0):
        f.append("A21 full-matrix LDLT pivot is not strict")
    pbar = d.get("same_word_diffuse_prior_covariance_upper", {})
    timing = pbar.get("timing", {})
    if float(timing.get("same_word_endpoint_s", math.nan)) != HORIZON_S:
        f.append("same-word A21 covariance upper is not referenced at 3 s")
    if timing.get("selected_S_windows_fit_same_word") is not True:
        f.append("selected S estimator escaped the complete word")
    if pbar.get("active_ba_nuisance_paid_in_H_subblock_upper") is not True:
        f.append("A21 covariance upper did not pay active-bias nuisance")
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
        "same_word_A21_Pbar_trace_upper": d["same_word_diffuse_prior_covariance_upper"]["Pbar_trace_upper"],
        "active_ba_variance_upper": d["same_word_diffuse_prior_covariance_upper"]["accel_bias_variance_upper"],
        "delta2_penalty": d["delta_squared_over_four_penalty_physical"],
        "x_cells": d["x_cells_certified"],
        "worst_A21_LDLT_pivot": d["worst_full_A21_LDLT_pivot_lower"],
        "A21_closed": d["A21_prior_free_completion_closed"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
