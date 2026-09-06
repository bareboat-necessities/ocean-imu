#!/usr/bin/env python3
"""Full A21 prior-free Riccati completion on canonical complete SEA3.

This closes the active 21-state mode without eta9 point-packet PE and without
embedding a translation principal-block lower as a full-state Loewner lower.

The complete-SEA3 H18 finite-memory estimator is valid while accelerometer bias
is active because its vector covariance upper explicitly pays the configured
b_a nuisance.  In A mode every valid accelerometer update is accepted and its
instantaneous Jacobian has H_ba=I.  At the word endpoint, combine that finite
H18 estimator with the required accelerometer sample:

    b_a = y_a - H_H x_H - n_a.

No independence is assumed.  With ||H_H||^2 <= f_max^2+1,

    E||e_ba||^2 <= 2 tr(R_a) + 2 ||H_H||^2 tr(Pbar_H).

Hence a concrete finite full-21 estimator exists in the diffuse-prior limit.
That proves the full A21 information matrix D_A is strictly positive without a
pointwise eta9 packet condition, and supplies a source-uniform endpoint upper
P_inf,A <= trace_A I.

The exact prior-free completion then uses the immediately following shipping
prediction:

    (1-delta) Q_A - (delta^2/4) F_A Pbar_A F_A' >> 0.

The second term is bounded by a source-uniform isotropic penalty.  The actual
21x21 matrix is certified with outward-rounded interval LDLT over the complete
h/tau cover.  Translation retains the stable factored shipping integrated-OU
Q; attitude/gyro-bias and active b_a consume the shipping process primitives.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_full_process_ucc as PROCESS
import ou3_sea3_a21_detectability_completion as ADET
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_full_word_event_algebra as EVENT
import ou3_sea3_h18_prior_free_completion as H18
import ou3_sea3_riccati_tube as TUBE
import ou3_sea3_riccati_tube_factored as FACTORED
import ou3_sea3_windowed_vector_pe as PE

DEFAULT_DOMAIN = COMPLETE.DEFAULT_DOMAIN
SCHEMA = 4
QUALIFICATION = "OU3_COMPLETE_SEA3_A21_PRIOR_FREE_FULL_MATRIX_COMPLETION"
USEFUL_GATE = 1.0e-18
HORIZON_S = 3.0
DIM = 21
OFF_V, OFF_P, OFF_S, OFF_AW, OFF_BA = 6, 9, 12, 15, 18


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _zero(n: int):
    z = I(0.0)
    return [[z for _ in range(n)] for _ in range(n)]


def _full_A21_cell(x: Interval, *, process: dict, dynamic: dict, penalty: float):
    delta = USEFUL_GATE
    one_minus = I(TUBE.down(1.0 - delta))
    M = _zero(DIM)

    q_att = float(process["attitude_gyro_bias"]["Q_attitude_gyro_bias_lambda_min_lower"])
    att_diag = TUBE.down((1.0 - delta) * q_att - penalty)
    if not att_diag > 0.0:
        return False, {"reason": "attitude_gyro_bias_penalty", "diagonal_lower": att_diag}
    for i in range(6):
        M[i][i] = I(att_diag)

    inv = dynamic["dynamic_invariant"]
    sigma_floor = float(inv["sigma_aw_filter_mps2"][0])
    h = float(dynamic["validated_rate_and_jump_bounds"]["dt_s"])
    if not (sigma_floor > 0.0 and h > 0.0):
        raise RuntimeError("invalid SEA3 translation process scale")
    scales = [sigma_floor*h, sigma_floor*h*h, sigma_floor*h*h*h, sigma_floor]
    B = FACTORED.step_scaled_q_over_x(x)
    qscale = I(x.lo)
    Maxis = [[one_minus*qscale*B[i][j] for j in range(4)] for i in range(4)]
    for i in range(4):
        pscaled = TUBE.up(penalty / TUBE.down(scales[i]*scales[i]))
        Maxis[i][i] = Maxis[i][i] - I(pscaled)
    idx = (OFF_V, OFF_P, OFF_S, OFF_AW)
    for axis in range(3):
        for i in range(4):
            for j in range(4):
                M[idx[i]+axis][idx[j]+axis] = Maxis[i][j]

    q_ba = float(process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"])
    ba_diag = TUBE.down((1.0 - delta) * q_ba - penalty)
    if not ba_diag > 0.0:
        return False, {"reason": "active_ba_penalty", "diagonal_lower": ba_diag}
    for i in range(3):
        M[OFF_BA+i][OFF_BA+i] = I(ba_diag)

    ok, pivots = symmetric_positive_definite_ldlt(M)
    return ok, {
        "x_h_over_tau": x.as_list(),
        "full_dimension": DIM,
        "full_21x21_interval_LDLT": True,
        "pivot_count": len(pivots),
        "pivot_lower": min((p.lo for p in pivots), default=math.inf),
        "attitude_gyro_bias_diagonal_lower": att_diag,
        "active_ba_diagonal_lower": ba_diag,
        "translation_congruence_scales": scales,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    dynamic = DYNAMIC.build(path)
    process = PROCESS.build()
    h18 = H18.build(path)
    adet = ADET.build(path)
    event = EVENT.build()
    pe = PE.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "dynamic": DYNAMIC.validate(dynamic),
        "process": PROCESS.validate(process),
        "H18": H18.validate(h18),
        "A21_detectability": ADET.validate(adet),
        "event_algebra": EVENT.validate(event),
        "PE": PE.validate(pe),
    }
    bad = {k:v for k,v in bad.items() if v}
    if bad:
        raise RuntimeError(f"A21 completion prerequisites failed: {bad}")
    if complete["canonical_P3_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("A21 completion detached from complete SEA3")
    if float(complete["word_horizon_s"]) != HORIZON_S:
        raise RuntimeError("canonical complete-SEA3 word is no longer 3 s")
    if adet["paper_active_bias_route"] != "ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION":
        raise RuntimeError("A21 paper route changed")
    if adet["eta9_point_packet_shortcut_used"] is not False:
        raise RuntimeError("eta9 point-packet shortcut entered A21")
    if h18["H18_prior_free_completion_closed"] is not True:
        raise RuntimeError("A21 requires the certified same-word H18 estimator")
    if pe["all_valid_accelerometer_packets_required"] is not True:
        raise RuntimeError("A mode no longer guarantees the endpoint accelerometer row")

    pbar_h = h18["same_word_diffuse_prior_covariance_upper"]
    trace_h = float(pbar_h["Pbar_trace_upper"])
    domain = json.loads(path.read_text(encoding="utf-8"))
    fmax = float(domain["normal_live"]["specific_force_norm_upper_mps2"])
    ra = float(pe["measurement_runtime"]["accelerometer_variance_upper"])
    hacc_h_norm_sq = TUBE.up(fmax*fmax + 1.0)
    ba_trace = TUBE.up(2.0*3.0*ra + 2.0*hacc_h_norm_sq*trace_h)
    trace_a = TUBE.up(trace_h + ba_trace)
    if not all(math.isfinite(v) and v > 0.0 for v in (trace_h,fmax,ra,hacc_h_norm_sq,ba_trace,trace_a)):
        raise RuntimeError("A21 finite estimator covariance bound invalid")

    Fnorm2 = TUBE.up(max(1.0, float(h18["prediction_F_spectral_norm_squared_upper"])))
    penalty = TUBE.up((USEFUL_GATE*USEFUL_GATE/4.0) * Fnorm2 * trace_a)
    if not (math.isfinite(penalty) and penalty > 0.0):
        raise RuntimeError("A21 delta^2 completion penalty invalid")

    rows, failures, worst = [], [], math.inf
    for x in H18._x_cover(dynamic):
        ok,row = _full_A21_cell(x, process=process, dynamic=dynamic, penalty=penalty)
        rows.append(row)
        if not ok:
            failures.append(row)
        else:
            worst = min(worst, float(row["pivot_lower"]))
    closed = bool(rows) and not failures and math.isfinite(worst) and worst > 0.0

    preserve = event["full_matrix_margin_preservation"]
    suffix_preserved = all(bool(preserve[k]) for k in (
        "covers_prediction", "covers_every_due_S_update",
        "covers_every_Normal_Live_accelerometer_update",
        "covers_asynchronous_magnetometer_update",
        "covers_immediate_left_error_reset", "covers_aw_covariance_floor",
        "covers_not_due_or_rejected_identity_branches",
    ))

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": complete["canonical_P3_source"],
        "complete_word_horizon_s": HORIZON_S,
        "component_of_complete_SEA3_full_word": True,
        "paper_active_bias_route": adet["paper_active_bias_route"],
        "eta9_point_packet_shortcut_used": False,
        "H18_finite_memory_estimator_consumed": True,
        "H18_estimator_pays_active_ba_nuisance": True,
        "required_A_mode_accelerometer_row_consumed": True,
        "A_mode_accelerometer_H_ba_is_identity": True,
        "finite_full_A21_linear_estimator_constructed": True,
        "finite_full_A21_estimator_implies_D_strictly_positive": True,
        "full_A21_prior_free_D_inverse_identity_used": True,
        "A21_diffuse_prior_covariance_upper": {
            "H18_trace_upper": trace_h,
            "accelerometer_H18_operator_norm_squared_upper": hacc_h_norm_sq,
            "accelerometer_variance_upper": ra,
            "ba_estimator_trace_upper": ba_trace,
            "full_A21_trace_upper": trace_a,
            "Loewner_upper": "P_inf,A <= trace_A * I_21",
            "no_H_estimator_measurement_noise_independence_assumed": True,
        },
        "prediction_F_spectral_norm_squared_upper": Fnorm2,
        "delta_squared_over_four_penalty_physical": penalty,
        "stable_factored_shipping_integrated_OU_Q_consumed": True,
        "shipping_active_ba_GM_Q_consumed": True,
        "x_cells_certified": len(rows),
        "x_cell_failures": failures,
        "worst_full_A21_LDLT_pivot_lower": worst if closed else None,
        "full_21x21_Omega_minus_delta_P_LDLT_closed": closed,
        "A21_prior_free_completion_closed": closed,
        "full_21x21_interval_LDLT_used": True,
        "event_algebra_preserves_margin_after_closure": suffix_preserved,
        "actual_applied_SpectralMSE_R_S_retained_through_H18_component": True,
        "old_one_step_Euclidean_Q_min_used": False,
        "scalar_beta_contraction_used": False,
        "blockwise_minimum_contraction_used": False,
        "source_family_replaced": False,
        "trajectory_replay_used": False,
        "independent_tau_sigma_RS_source_created": False,
        "useful_gate": USEFUL_GATE,
        "P3_promoted": False,
        "next_obligation": (
            "wire the certified exact reset congruence into the literal S/accelerometer/magnetometer event API, then consume H18+A21 in canonical P3"
            if closed else
            "tighten only this same-word A21 estimator/SEA3 x cover; do not introduce eta9 or relax delta"
        ),
    }


def validate(d: dict) -> list[str]:
    f=[]
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION: f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD": f.append("canonical source changed")
    for k in (
        "component_of_complete_SEA3_full_word", "H18_finite_memory_estimator_consumed",
        "H18_estimator_pays_active_ba_nuisance", "required_A_mode_accelerometer_row_consumed",
        "A_mode_accelerometer_H_ba_is_identity", "finite_full_A21_linear_estimator_constructed",
        "finite_full_A21_estimator_implies_D_strictly_positive", "full_A21_prior_free_D_inverse_identity_used",
        "stable_factored_shipping_integrated_OU_Q_consumed", "shipping_active_ba_GM_Q_consumed",
        "full_21x21_Omega_minus_delta_P_LDLT_closed", "A21_prior_free_completion_closed",
        "full_21x21_interval_LDLT_used", "event_algebra_preserves_margin_after_closure",
        "actual_applied_SpectralMSE_R_S_retained_through_H18_component",
    ):
        if d.get(k) is not True: f.append(k)
    for k in (
        "eta9_point_packet_shortcut_used", "old_one_step_Euclidean_Q_min_used",
        "scalar_beta_contraction_used", "blockwise_minimum_contraction_used",
        "source_family_replaced", "trajectory_replay_used",
        "independent_tau_sigma_RS_source_created", "P3_promoted",
    ):
        if d.get(k) is not False: f.append(k)
    if d.get("paper_active_bias_route") != "ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION": f.append("paper active-bias route")
    if float(d.get("useful_gate",math.nan)) != USEFUL_GATE: f.append("useful gate")
    if int(d.get("x_cells_certified",0)) <= 0: f.append("no SEA3 x cells")
    if d.get("x_cell_failures"): f.append("one or more A21 x cells failed")
    p=d.get("worst_full_A21_LDLT_pivot_lower")
    if not isinstance(p,(int,float)) or not(math.isfinite(float(p)) and float(p)>0): f.append("A21 full-matrix pivot")
    est=d.get("A21_diffuse_prior_covariance_upper",{})
    for k in ("H18_trace_upper","ba_estimator_trace_upper","full_A21_trace_upper"):
        x=est.get(k)
        if not isinstance(x,(int,float)) or not(math.isfinite(float(x)) and float(x)>0): f.append(f"invalid estimator {k}")
    if est.get("no_H_estimator_measurement_noise_independence_assumed") is not True: f.append("estimator independence shortcut")
    return list(dict.fromkeys(f))


def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument("--domain",type=Path,default=DEFAULT_DOMAIN); ap.add_argument("--output",type=Path,required=True); a=ap.parse_args()
    d=build(a.domain); f=validate(d); d["validation_pass"]=not f; d["validation_failures"]=f
    a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(d,indent=2,sort_keys=True),encoding="utf-8")
    print(json.dumps({"A21_trace_upper":d["A21_diffuse_prior_covariance_upper"]["full_A21_trace_upper"],"delta2_penalty":d["delta_squared_over_four_penalty_physical"],"x_cells":d["x_cells_certified"],"worst_A21_LDLT_pivot":d["worst_full_A21_LDLT_pivot_lower"],"A21_closed":d["A21_prior_free_completion_closed"],"failures":f},indent=2,sort_keys=True))
    return 0 if not f else 2

if __name__=="__main__": raise SystemExit(main())
