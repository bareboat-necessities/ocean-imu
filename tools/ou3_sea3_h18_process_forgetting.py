#!/usr/bin/env python3
"""Full-matrix H18 process-forgetting prefix on the complete SEA3 source.

This certificate closes a different obligation from the four-S/vector-PE UCO
lemma.  It proves that, after Live entry, the shipping process noise has made
the covariance forget at least ``delta=1e-18`` of the source-generated Live
covariance.  The UCO/R_S lemma is still required separately for boundedness and
is not replaced here.

No source word is generated.  The proof uses facts true of every admitted
``COMPLETE_SEA3_NORMAL_LIVE_WORD``:

* Live starts from the shipping structured covariance seed.
* Every valid IMU sample advances the tuner EMAs, candidates are staged only
  after ``ADAPT_EVERY_SECS`` and committed on the next IMU boundary.
* After any active commit, two further active commits cannot occur in the next
  0.2 s because successive staging instants are strictly more than 0.1 s apart.
  The first following commit is enclosed by the existing source-derived commit
  gap upper.  Hence the 0.2 s prefix contains at most two constant-active OU
  pieces.  We allow arbitrary legal tau on both pieces; this is a cover of the
  same SEA3 path, not an independent source language.
* sigma_aw is never below the deployed 0.05 m/s^2 process floor.  For a fixed
  active tau the shipping integrated-OU Q is linear in sigma_aw^2, so evaluating
  the exact shipping Q at 0.05 is a PSD process subnoise of every admitted word.

Condition on the word-start error and on every unselected process/noise
coordinate.  The selected process subnoise and configured measurement noise
then form a smaller linear-Gaussian estimation problem.  Its optimal posterior
covariance is a lower bound on the actual noise contribution Omega, because the
shipping gains are merely one admissible estimator after that extra
conditioning.  If beta bounds the complete selected-mode measurement
information, then

    Omega >= (1/(1+beta)) Q_selected.

The beta bound admits every valid accelerometer packet and an S packet at every
IMU slot.  Actual applied SpectralMSE R_S is retained through its deployed
positive minimum and XY factors.  No hardware magnetometer rate is required:
closure is taken immediately after the final accelerometer/reset and before
any asynchronous magnetometer event.

For P, conditioning can only reduce covariance.  Before measurement/reset
conditioning, the no-measurement H18 covariance remains block diagonal between
[theta,b_g] and each translation axis.  Translation diagonal variances are
bounded with the undamped integrator kernel and the complete-SEA3 q_c upper.
For any PSD 4x4 translation covariance with diagonal <= U,

    P_axis <= 4 diag(U),

by diagonal normalization and trace(C)<=4.  The attitude/gyro-bias block is
bounded by its no-measurement trace.  Left-error resets are pure nonsingular
congruences; the generalized inequality Omega-delta P>=0 is evaluated in the
reset-free transported chart and therefore needs no reset-norm bound.

For every tau1 x tau2 x commit-phase cell we assemble the complete 18x18
interval matrix

    M = Omega_lower - delta P_upper

and certify it directly by interval LDL^T.  There is no scalar block ratio,
D/L product, determinant/trace eigenvalue conversion, replay, predecessor
history graph, or alternate estimator.
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
    matrix_add,
    matrix_mul,
    matrix_sub,
    matrix_transpose,
    symmetric_positive_definite_ldlt,
)
import ou3_full_process_ucc as PROCESS
import ou3_sea3_complete_source as COMPLETE
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_live_covariance_seed as LIVE
import ou3_sea3_shipping_prediction_primitives as PRED
import ou3_sea3_windowed_vector_pe as PE

REPO = Path(__file__).resolve().parents[1]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_H18_FULL_MATRIX_PROCESS_FORGETTING"
USEFUL_GATE = 1.0e-18
PREFIX_HORIZON_S = 0.2
TAU_BASE_CELLS = 32
PHASE_BASE_CELLS = 2
MAX_REFINEMENT_DEPTH = 5
MAX_ACCEPTED_CELLS = 50000

OFF_V = 6
OFF_P = 9
OFF_S = 12
OFF_AW = 15


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _zero(n: int) -> IntervalMatrix:
    return [[I(0.0) for _ in range(n)] for _ in range(n)]


def _scale(A: IntervalMatrix, c: float) -> IntervalMatrix:
    ci = I(c)
    return [[ci * x for x in row] for row in A]


def _diag(values: list[float]) -> IntervalMatrix:
    n = len(values)
    z = I(0.0)
    return [[I(values[i]) if i == j else z for j in range(n)] for i in range(n)]


def _translation_diag_upper(H: float, live: dict, dynamic: dict) -> list[float]:
    """No-measurement [v,p,S,a] diagonal upper through physical time H."""
    tr = live["translation_seed"]
    Pv = float(tr["P_v"])
    Pp = float(tr["P_p"])
    PS = float(tr["P_S"])
    sigma_hi = float(dynamic["dynamic_invariant"]["sigma_aw_filter_mps2"][1])
    Pa = up(sigma_hi * sigma_hi)
    tau_lo = float(dynamic["dynamic_invariant"]["tau_applied_s"][0])
    qc_hi = up(2.0 * Pa / tau_lo)
    h2 = up(H * H)
    h3 = up(h2 * H)
    h4 = up(h2 * h2)
    h5 = up(h4 * H)
    h6 = up(h3 * h3)
    h7 = up(h6 * H)
    Uv = up(Pv + up(h2 * Pa) + up(qc_hi * h3 / 3.0))
    Up = up(Pp + up(h2 * Pv) + up(0.25 * h4 * Pa) + up(qc_hi * h5 / 20.0))
    US = up(
        PS
        + up(h2 * Pp)
        + up(0.25 * h4 * Pv)
        + up((h6 / 36.0) * Pa)
        + up(qc_hi * h7 / 252.0)
    )
    Ua = up(Pa + up(qc_hi * H))
    return [Uv, Up, US, Ua]


def _eta_trace_upper(H: float, live: dict, process: dict) -> float:
    """Crude but source-uniform no-measurement trace upper for [theta,b_g]."""
    att = live["full_heading_gauged_live_attitude_seed"]
    tilt = float(att["tilt_std_rad"])
    yaw = float(att["yaw_std_rad"])
    Pb = float(live["constructor"]["P_bg_variance"])
    trace = up(2.0 * tilt * tilt + yaw * yaw + 3.0 * Pb)

    rt = process["configured_runtime"]["imu_dt_outward_interval_s"]
    h_lo, h_hi = map(float, rt)
    if not (0.0 < h_lo <= h_hi):
        raise RuntimeError("invalid IMU dt for eta trace upper")
    steps = int(math.ceil(H / h_lo)) + 2
    qg = max(float(x) ** 2 for x in process["source_constants"]["gyro_noise_density_rad_sqrt_s_per_axis"])
    qb = float(process["source_constants"]["gyro_bias_rw_variance_density"])
    gain2 = up((1.0 + h_hi) ** 2)
    qtrace = up(3.0 * qg * h_hi + qb * h_hi ** 3 + 3.0 * qb * h_hi)
    for _ in range(steps):
        trace = up(gain2 * trace + qtrace)
    return trace


def _measurement_attenuation_lower(
    H: float, dynamic: dict, complete: dict, pe: dict, process: dict
) -> dict:
    """One-shot lifted information bound for selected process modes."""
    dt = float(dynamic["validated_rate_and_jump_bounds"]["dt_s"])
    packets = int(math.ceil(H / dt)) + 2
    tau_lo = float(dynamic["dynamic_invariant"]["tau_applied_s"][0])
    sigma_floor = 0.05
    # Selected subnoise uses sigma_floor at the actual tau.  This is its largest
    # possible driving intensity and therefore bounds all selected-mode
    # measurement sensitivities from above.
    qc_selected_hi = up(2.0 * sigma_floor * sigma_floor / tau_lo)

    acc_std = float(pe["measurement_runtime"]["accelerometer_std_mps2"][0])
    racc_lo = down(acc_std * acc_std)
    if not racc_lo > 0.0:
        raise RuntimeError("accelerometer variance lower lost positivity")

    # Sum t_k <= dt*N(N+1)/2 for the aw marginal of the selected base process.
    sum_t = up(dt * packets * (packets + 1) / 2.0)
    beta_acc_translation = up(3.0 * qc_selected_hi * sum_t / racc_lo)

    rs_lo = float(dynamic["dynamic_invariant"]["R_S_applied"][0])
    factors = [float(x) for x in complete["R_S_regularizer"]["axis_std_factors"]]
    inv_rs_var_sum = 0.0
    for f in factors:
        s = down(rs_lo * f)
        inv_rs_var_sum = up(inv_rs_var_sum + up(1.0 / down(s * s)))
    # At most one S packet per IMU slot.  Undamped S response gives
    # Var(S(t)) <= q_c t^7/252 for one axis.
    beta_S_translation = up(
        packets * qc_selected_hi * up(H ** 7) * inv_rs_var_sum / 252.0
    )

    q_eta = float(process["attitude_gyro_bias"]["Q_attitude_gyro_bias_lambda_min_lower"])
    domain = json.loads(DEFAULT_DOMAIN.read_text(encoding="utf-8"))
    f_hi = float(domain["normal_live"]["specific_force_norm_upper_mps2"])
    # eta selected modes are injected by the final prediction and closure is
    # taken immediately after that sample's accelerometer/reset.  Only that
    # accelerometer can observe them.  ||skew(f)||_F^2=2||f||^2.
    beta_eta = up(2.0 * q_eta * f_hi * f_hi / racc_lo)

    beta = up(beta_acc_translation + beta_S_translation + beta_eta)
    gamma = down(1.0 / up(1.0 + beta))
    if not (math.isfinite(gamma) and gamma > 0.0):
        raise RuntimeError("selected-process posterior attenuation lost positivity")
    return {
        "measurement_slots_upper": packets,
        "selected_sigma_floor_mps2": sigma_floor,
        "selected_qc_upper_for_information": qc_selected_hi,
        "accelerometer_variance_lower": racc_lo,
        "actual_applied_R_S_base_lower": rs_lo,
        "actual_applied_R_S_axis_factors": factors,
        "beta_acc_translation_upper": beta_acc_translation,
        "beta_S_translation_upper": beta_S_translation,
        "beta_final_eta_acc_upper": beta_eta,
        "beta_total_upper": beta,
        "posterior_attenuation_lower": gamma,
        "all_accelerometer_packets_admitted": True,
        "S_packet_admitted_at_every_IMU_slot": True,
        "asynchronous_magnetometer_needed_before_closure": False,
    }


def _tau_bounds(lo: float, hi: float, count: int) -> list[tuple[float, float]]:
    if not (0.0 < lo < hi and count >= 1):
        raise ValueError("invalid tau cover")
    ratio = hi / lo
    edges = [lo * (ratio ** (k / count)) for k in range(count + 1)]
    edges[0], edges[-1] = lo, hi
    return [(down(edges[k]), up(edges[k + 1])) for k in range(count)]


def _phase_bounds(lo: float, hi: float, count: int) -> list[tuple[float, float]]:
    w = (hi - lo) / count
    return [
        (down(lo + k * w), up(hi if k + 1 == count else lo + (k + 1) * w))
        for k in range(count)
    ]


def _selected_translation_Q(
    tau1: Interval, tau2: Interval, h1: Interval, H: float
) -> IntervalMatrix:
    h2 = Interval.outward_bounds(down(H - h1.hi), up(H - h1.lo))
    s2 = I(0.05 * 0.05)
    Q1 = PRED.translation_axis_process(tau1, h1, s2)
    F2 = PRED.translation_axis_transition(tau2, h2)
    Q2 = PRED.translation_axis_process(tau2, h2, s2)
    return matrix_add(matrix_mul(matrix_mul(F2, Q1), matrix_transpose(F2)), Q2)


def _full_matrices(
    Qaxis: IntervalMatrix,
    gamma: float,
    q_eta: float,
    eta_trace_upper: float,
    Utrans: list[float],
    delta: float,
) -> tuple[IntervalMatrix, IntervalMatrix, IntervalMatrix]:
    n = 18
    Om = _zero(n)
    Pu = _zero(n)
    for i in range(6):
        Om[i][i] = I(down(gamma * q_eta))
        Pu[i][i] = I(eta_trace_upper)
    idx = (OFF_V, OFF_P, OFF_S, OFF_AW)
    for axis in range(3):
        for i in range(4):
            Pu[idx[i] + axis][idx[i] + axis] = I(up(4.0 * Utrans[i]))
            for j in range(4):
                Om[idx[i] + axis][idx[j] + axis] = I(gamma) * Qaxis[i][j]
    M = matrix_sub(Om, _scale(Pu, delta))
    return Om, Pu, M


def _cell_certificate(
    tau1: tuple[float, float],
    tau2: tuple[float, float],
    h1: tuple[float, float],
    *,
    H: float,
    gamma: float,
    q_eta: float,
    eta_trace: float,
    Utrans: list[float],
    delta: float,
) -> tuple[bool, dict]:
    t1 = Interval.outward_bounds(*tau1)
    t2 = Interval.outward_bounds(*tau2)
    hp = Interval.outward_bounds(*h1)
    Qaxis = _selected_translation_Q(t1, t2, hp, H)
    _, _, M = _full_matrices(Qaxis, gamma, q_eta, eta_trace, Utrans, delta)
    ok, pivots = symmetric_positive_definite_ldlt(M)
    return ok, {
        "tau1": [t1.lo, t1.hi],
        "tau2": [t2.lo, t2.hi],
        "first_piece_s": [hp.lo, hp.hi],
        "second_piece_s": [down(H - hp.hi), up(H - hp.lo)],
        "pivot_lower": min((p.lo for p in pivots), default=math.inf),
        "pivot_count": len(pivots),
    }


def _split_pair(x: tuple[float, float]) -> tuple[tuple[float, float], tuple[float, float]]:
    lo, hi = x
    mid = 0.5 * (lo + hi)
    return (lo, up(mid)), (down(mid), hi)


def _certify_cover(
    tau_cells: list[tuple[float, float]],
    phase_cells: list[tuple[float, float]],
    **kwargs,
) -> dict:
    accepted = 0
    worst_pivot = math.inf
    worst = None
    refinements = 0
    stack: list[tuple[tuple[float,float], tuple[float,float], tuple[float,float], int]] = [
        (a, b, p, 0) for a in tau_cells for b in tau_cells for p in phase_cells
    ]
    while stack:
        a, b, p, depth = stack.pop()
        ok, row = _cell_certificate(a, b, p, **kwargs)
        if ok:
            accepted += 1
            if row["pivot_lower"] < worst_pivot:
                worst_pivot = row["pivot_lower"]
                worst = row
            if accepted > MAX_ACCEPTED_CELLS:
                raise RuntimeError("H18 process-prefix cell cap exceeded")
            continue
        if depth >= MAX_REFINEMENT_DEPTH:
            return {
                "pass": False,
                "accepted_cells": accepted,
                "refinements": refinements,
                "failed_cell": row,
            }
        # Split the widest relative parameter.  Phase width is normalized by H.
        scores = [
            (a[1] - a[0]) / max(a[0], 1e-12),
            (b[1] - b[0]) / max(b[0], 1e-12),
            (p[1] - p[0]) / PREFIX_HORIZON_S,
        ]
        which = max(range(3), key=lambda i: scores[i])
        parts = _split_pair((a, b, p)[which])
        for part in parts:
            vals = [a, b, p]
            vals[which] = part
            stack.append((vals[0], vals[1], vals[2], depth + 1))
        refinements += 1
    return {
        "pass": True,
        "accepted_cells": accepted,
        "refinements": refinements,
        "worst_full_H18_LDLT_pivot_lower": worst_pivot,
        "worst_cell": worst,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    complete = COMPLETE.build(path)
    dynamic = DYNAMIC.build(path)
    live = LIVE.build(path)
    process = PROCESS.build()
    pe = PE.build(path)
    bad = {
        "complete": COMPLETE.validate(complete),
        "dynamic": DYNAMIC.validate(dynamic),
        "live": LIVE.validate(live),
        "process": PROCESS.validate(process),
        "PE": PE.validate(pe),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"H18 process-forgetting prerequisites failed: {bad}")
    if complete["canonical_P3_source"] != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        raise RuntimeError("process forgetting detached from complete SEA3")

    wrapper = WRAPPER.read_text(encoding="utf-8")
    parity = {
        "staging_uses_strict_adaptation_cadence": (
            "if (time_ - last_adapt_time_sec_ > adapt_every_secs_)" in wrapper
            and "last_adapt_time_sec_ = time_;" in wrapper
        ),
        "commit_is_next_sample_boundary": (
            "void apply_pending_online_tune_()" in wrapper
            and "online_tune_apply_pending_ = false;" in wrapper
        ),
        "same_candidate_commits_tau_and_RS": (
            "apply_ou_tune_(false);" in wrapper and "apply_RS_tune_();" in wrapper
        ),
    }
    if not all(parity.values()):
        raise RuntimeError(f"shipping commit parity failed: {parity}")

    adapt = float(dynamic["validated_rate_and_jump_bounds"]["active_commit_gap_s_upper"])
    cadence = 0.1
    if not (cadence < adapt <= 0.111):
        raise RuntimeError(f"unexpected active commit gap upper {adapt}")
    # Start at the first post-Live active commit (which occurs within adapt).
    # The next commit lies strictly after cadence and no later than adapt.
    # A third commit cannot fit before 2*cadence, so the 0.2 s prefix has at
    # most one interior active-schedule change.
    phase_lo = down(cadence)
    phase_hi = up(min(adapt, PREFIX_HORIZON_S - 0.005))
    if not phase_lo < phase_hi < PREFIX_HORIZON_S:
        raise RuntimeError("invalid one-change phase cover")

    total_from_live = up(adapt + PREFIX_HORIZON_S)
    Utrans = _translation_diag_upper(total_from_live, live, dynamic)
    eta_trace = _eta_trace_upper(total_from_live, live, process)
    attenuation = _measurement_attenuation_lower(
        PREFIX_HORIZON_S, dynamic, complete, pe, process
    )
    gamma = float(attenuation["posterior_attenuation_lower"])
    q_eta = float(process["attitude_gyro_bias"]["Q_attitude_gyro_bias_lambda_min_lower"])
    if not (q_eta > 0.0 and gamma > 0.0):
        raise RuntimeError("H18 selected-process lower lost positivity")

    tau_lo, tau_hi = map(float, dynamic["dynamic_invariant"]["tau_applied_s"])
    tau_cells = _tau_bounds(tau_lo, tau_hi, TAU_BASE_CELLS)
    phase_cells = _phase_bounds(phase_lo, phase_hi, PHASE_BASE_CELLS)
    cover = _certify_cover(
        tau_cells,
        phase_cells,
        H=PREFIX_HORIZON_S,
        gamma=gamma,
        q_eta=q_eta,
        eta_trace=eta_trace,
        Utrans=Utrans,
        delta=USEFUL_GATE,
    )

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_source": complete["canonical_P3_source"],
        "theorem_prefix_starts_from_shipping_Live_seed": True,
        "first_post_Live_commit_wait_s_upper": adapt,
        "process_prefix_horizon_s": PREFIX_HORIZON_S,
        "closure_time_from_Live_s_upper": total_from_live,
        "at_most_one_active_schedule_change_inside_prefix": True,
        "interior_commit_phase_s": [phase_lo, phase_hi],
        "two_piece_tau_values_are_same_SEA3_path_coordinates": True,
        "independent_tau_sigma_RS_source_created": False,
        "trajectory_replay_used": False,
        "source_history_graph_used": False,
        "sigma_floor_selected_as_PSD_subnoise": 0.05,
        "shipping_integrated_OU_process_primitive_consumed": True,
        "measurement_conditioning": attenuation,
        "actual_applied_SpectralMSE_R_S_consumed": True,
        "left_error_resets_quotiented_as_nonsingular_congruences": True,
        "no_measurement_covariance_upper": {
            "translation_axis_order": ["v", "p", "S", "a_w"],
            "translation_diagonal_upper": Utrans,
            "translation_Loewner_upper": "P_axis <= 4*diag(U)",
            "attitude_gyro_bias_trace_upper": eta_trace,
            "measurement_updates_can_only_reduce_conditional_covariance": True,
            "aw_floor_is_covered_by_sigma_safety_upper": True,
        },
        "selected_eta_last_prediction_Q_lambda_min_lower": q_eta,
        "useful_gate": USEFUL_GATE,
        "parameter_cover": {
            "tau_interval_s": [tau_lo, tau_hi],
            "tau_base_cells_per_piece": TAU_BASE_CELLS,
            "phase_base_cells": PHASE_BASE_CELLS,
            "max_refinement_depth": MAX_REFINEMENT_DEPTH,
            **cover,
        },
        "full_H18_Omega_minus_delta_P_LDLT_closed": bool(cover.get("pass")),
        "blockwise_minimum_ratio_used": False,
        "scalar_information_beta_used_for_contraction": False,
        "determinant_trace_eigenvalue_conversion_used": False,
        "D_W_L_W_product_used": False,
        "one_step_process_Q_used_as_contraction_ratio": False,
        "H18_PROCESS_FORGETTING_PASS": bool(cover.get("pass")),
        "P3_promoted": False,
        "next_obligation": (
            "combine this full-matrix H18 forgetting prefix with the corrected actual-R_S H18 UCO gate; "
            "then prove the H->A release and A21 full-matrix forgetting/detectability transition before promoting P3"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        f.append("canonical source changed")
    for key in (
        "theorem_prefix_starts_from_shipping_Live_seed",
        "at_most_one_active_schedule_change_inside_prefix",
        "two_piece_tau_values_are_same_SEA3_path_coordinates",
        "shipping_integrated_OU_process_primitive_consumed",
        "actual_applied_SpectralMSE_R_S_consumed",
        "left_error_resets_quotiented_as_nonsingular_congruences",
        "full_H18_Omega_minus_delta_P_LDLT_closed",
        "H18_PROCESS_FORGETTING_PASS",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    for key in (
        "independent_tau_sigma_RS_source_created",
        "trajectory_replay_used",
        "source_history_graph_used",
        "blockwise_minimum_ratio_used",
        "scalar_information_beta_used_for_contraction",
        "determinant_trace_eigenvalue_conversion_used",
        "D_W_L_W_product_used",
        "one_step_process_Q_used_as_contraction_ratio",
        "P3_promoted",
    ):
        if d.get(key) is not False:
            f.append(f"forbidden/open flag {key} changed")
    cover = d.get("parameter_cover", {})
    if cover.get("pass") is not True:
        f.append(f"full H18 LDLT cover failed: {cover.get('failed_cell')}")
    pivot = cover.get("worst_full_H18_LDLT_pivot_lower")
    if not isinstance(pivot, (int, float)) or not (math.isfinite(float(pivot)) and float(pivot) > 0.0):
        f.append("H18 LDLT pivot lower is not strict")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        f.append("useful gate changed")
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
        "closure_time_from_Live_s_upper": d["closure_time_from_Live_s_upper"],
        "posterior_attenuation_lower": d["measurement_conditioning"]["posterior_attenuation_lower"],
        "beta_total_upper": d["measurement_conditioning"]["beta_total_upper"],
        "translation_diagonal_upper": d["no_measurement_covariance_upper"]["translation_diagonal_upper"],
        "eta_trace_upper": d["no_measurement_covariance_upper"]["attitude_gyro_bias_trace_upper"],
        "accepted_cells": d["parameter_cover"].get("accepted_cells"),
        "refinements": d["parameter_cover"].get("refinements"),
        "worst_H18_LDLT_pivot_lower": d["parameter_cover"].get("worst_full_H18_LDLT_pivot_lower"),
        "pass": d["H18_PROCESS_FORGETTING_PASS"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
