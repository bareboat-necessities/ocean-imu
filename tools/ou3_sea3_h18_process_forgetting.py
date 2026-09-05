#!/usr/bin/env python3
"""Full-matrix H18 process-forgetting prefix on the complete SEA3 source.

This certificate closes a different obligation from the four-S/vector-PE UCO
lemma. It proves that, after Live entry, shipping process noise makes the
covariance forget at least ``delta=1e-18`` of the source-generated Live
covariance. The UCO/R_S lemma remains a separate required same-source result.

No source word is generated. The proof uses facts true of every admitted
``COMPLETE_SEA3_NORMAL_LIVE_WORD``:

* Live starts from the shipping structured covariance seed.
* Tuner candidates are staged only after ``ADAPT_EVERY_SECS`` and committed at
  the next IMU boundary.
* Starting at the first post-Live active commit, a 0.2 s prefix has at most one
  interior active-schedule change because successive staging instants are
  strictly more than 0.1 s apart. The possible change phase is bounded by the
  existing binary32 source-derived active-commit gap certificate.
* sigma_aw is never below the deployed 0.05 m/s^2 process floor. For fixed
  active tau, shipping integrated-OU Q is linear in sigma_aw^2, so the exact
  shipping Q at 0.05 is a PSD process subnoise of every admitted word.

Condition on the word-start error and every unselected process/noise coordinate.
For the selected process subnoise, the complete sequence of configured
accelerometer and S measurements is one linear-Gaussian observation operator.
Bounding its Hilbert-Schmidt information by beta gives the operator inequality

    (I + A^T R^-1 A)^-1 >= I/(1+beta),

and therefore the endpoint selected process covariance obeys

    Omega >= gamma Q_selected, gamma=1/(1+beta).

The beta bound admits every valid accelerometer packet and an S packet at every
IMU slot. Actual applied SpectralMSE R_S is retained through its deployed
positive minimum and XY factors. Closure is taken immediately after the final
accelerometer/reset and before any asynchronous magnetometer event.

For P, conditioning can only reduce covariance. Before measurement/reset
conditioning, the no-measurement H18 covariance remains block diagonal between
[theta,b_g] and each translation axis. Translation diagonal variances are
bounded with the undamped integrator kernel and the complete-SEA3 q_c upper.
For any PSD 4x4 translation covariance with diagonal <= U,

    P_axis <= 4 diag(U),

by diagonal normalization and trace(C)<=4. The attitude/gyro-bias block is
bounded by its no-measurement trace. Left-error resets are nonsingular
congruences; the generalized inequality is certified in the reset-free
transported chart and is invariant under the shipping reset congruences.

For every tau1 x tau2 x commit-phase cell this module assembles the complete
18x18 interval matrix

    M = Omega_lower - delta P_upper

and certifies it directly by interval LDL^T. There is no scalar block ratio,
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
    """Source-uniform no-measurement trace upper for [theta,b_g]."""
    att = live["full_heading_gauged_live_attitude_seed"]
    tilt = float(att["tilt_std_rad"])
    yaw = float(att["yaw_std_rad"])
    Pb = float(live["constructor"]["P_bg_variance"])
    trace = up(2.0 * tilt * tilt + yaw * yaw + 3.0 * Pb)
    h_lo, h_hi = map(float, process["configured_runtime"]["imu_dt_outward_interval_s"])
    if not (0.0 < h_lo <= h_hi):
        raise RuntimeError("invalid IMU dt for eta trace upper")
    steps = int(math.ceil(H / h_lo)) + 2
    qg = max(
        float(x) ** 2
        for x in process["source_constants"]["gyro_noise_density_rad_sqrt_s_per_axis"]
    )
    qb = float(process["source_constants"]["gyro_bias_rw_variance_density"])
    gain2 = up((1.0 + h_hi) ** 2)
    qtrace = up(3.0 * qg * h_hi + qb * h_hi ** 3 + 3.0 * qb * h_hi)
    for _ in range(steps):
        trace = up(gain2 * trace + qtrace)
    return trace


def _measurement_attenuation_lower(
    H: float,
    dynamic: dict,
    complete: dict,
    pe: dict,
    process: dict,
    domain_path: Path,
) -> dict:
    """One-shot lifted information bound for selected process modes."""
    dt = float(dynamic["validated_rate_and_jump_bounds"]["dt_s"])
    packets = int(math.ceil(H / dt)) + 2
    tau_lo = float(dynamic["dynamic_invariant"]["tau_applied_s"][0])
    sigma_floor = 0.05
    qc_selected_hi = up(2.0 * sigma_floor * sigma_floor / tau_lo)

    acc_std = float(pe["measurement_runtime"]["accelerometer_std_mps2"][0])
    racc_lo = down(acc_std * acc_std)
    if not racc_lo > 0.0:
        raise RuntimeError("accelerometer variance lower lost positivity")

    # The selected OU driving noise has a_w variance <= q_c t when damping is
    # dropped. Summing trace(H Cov H^T R^-1) over all possible accelerometer
    # slots bounds the Hilbert-Schmidt information operator.
    sum_t = up(dt * packets * (packets + 1) / 2.0)
    beta_acc_translation = up(3.0 * qc_selected_hi * sum_t / racc_lo)

    rs_lo = float(dynamic["dynamic_invariant"]["R_S_applied"][0])
    factors = [float(x) for x in complete["R_S_regularizer"]["axis_std_factors"]]
    inv_rs_var_sum = 0.0
    for factor in factors:
        std = down(rs_lo * factor)
        inv_rs_var_sum = up(inv_rs_var_sum + up(1.0 / down(std * std)))
    # For S, dropping damping gives Var(S(t)) <= q_c t^7/252 per axis.
    beta_S_translation = up(
        packets * qc_selected_hi * up(H ** 7) * inv_rs_var_sum / 252.0
    )

    q_eta = float(process["attitude_gyro_bias"]["Q_attitude_gyro_bias_lambda_min_lower"])
    domain = json.loads(Path(domain_path).read_text(encoding="utf-8"))
    f_hi = float(domain["normal_live"]["specific_force_norm_upper_mps2"])
    # eta selected noise is the final prediction's q_eta I6 subnoise. Only the
    # final accelerometer exists after that injection before closure.
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
    width = (hi - lo) / count
    return [
        (
            down(lo + k * width),
            up(hi if k + 1 == count else lo + (k + 1) * width),
        )
        for k in range(count)
    ]


def _selected_translation_Q(
    tau1: Interval,
    tau2: Interval,
    h1: Interval,
    H: float,
) -> IntervalMatrix:
    h2 = Interval.outward_bounds(down(H - h1.hi), up(H - h1.lo))
    sigma2 = I(0.05 * 0.05)
    Q1 = PRED.translation_axis_process(tau1, h1, sigma2)
    F2 = PRED.translation_axis_transition(tau2, h2)
    Q2 = PRED.translation_axis_process(tau2, h2, sigma2)
    return matrix_add(matrix_mul(matrix_mul(F2, Q1), matrix_transpose(F2)), Q2)


def _full_matrix(
    Qaxis: IntervalMatrix,
    gamma: float,
    q_eta: float,
    eta_trace_upper: float,
    Utrans: list[float],
    delta: float,
) -> IntervalMatrix:
    n = 18
    omega_lower = _zero(n)
    p_upper = _zero(n)
    for i in range(6):
        omega_lower[i][i] = I(down(gamma * q_eta))
        p_upper[i][i] = I(eta_trace_upper)
    idx = (OFF_V, OFF_P, OFF_S, OFF_AW)
    for axis in range(3):
        for i in range(4):
            p_upper[idx[i] + axis][idx[i] + axis] = I(up(4.0 * Utrans[i]))
            for j in range(4):
                omega_lower[idx[i] + axis][idx[j] + axis] = I(gamma) * Qaxis[i][j]
    return matrix_sub(omega_lower, _scale(p_upper, delta))


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
    M = _full_matrix(Qaxis, gamma, q_eta, eta_trace, Utrans, delta)
    ok, pivots = symmetric_positive_definite_ldlt(M)
    return ok, {
        "tau1": [t1.lo, t1.hi],
        "tau2": [t2.lo, t2.hi],
        "first_piece_s": [hp.lo, hp.hi],
        "second_piece_s": [down(H - hp.hi), up(H - hp.lo)],
        "pivot_lower": min((p.lo for p in pivots), default=math.inf),
        "pivot_count": len(pivots),
    }


def _split_pair(
    bounds: tuple[float, float],
) -> tuple[tuple[float, float], tuple[float, float]]:
    lo, hi = bounds
    mid = 0.5 * (lo + hi)
    return (lo, up(mid)), (down(mid), hi)


def _certify_cover(
    tau_cells: list[tuple[float, float]],
    phase_cells: list[tuple[float, float]],
    **kwargs,
) -> dict:
    accepted = 0
    refinements = 0
    worst_pivot = math.inf
    worst = None
    stack: list[
        tuple[tuple[float, float], tuple[float, float], tuple[float, float], int]
    ] = [(a, b, p, 0) for a in tau_cells for b in tau_cells for p in phase_cells]
    while stack:
        a, b, phase, depth = stack.pop()
        ok, row = _cell_certificate(a, b, phase, **kwargs)
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
        scores = [
            (a[1] - a[0]) / max(a[0], 1e-12),
            (b[1] - b[0]) / max(b[0], 1e-12),
            (phase[1] - phase[0]) / PREFIX_HORIZON_S,
        ]
        which = max(range(3), key=lambda i: scores[i])
        pieces = _split_pair((a, b, phase)[which])
        for piece in pieces:
            values = [a, b, phase]
            values[which] = piece
            stack.append((values[0], values[1], values[2], depth + 1))
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
    bad = {key: value for key, value in bad.items() if value}
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

    rate = dynamic["validated_rate_and_jump_bounds"]
    active_gap_upper = float(rate["active_commit_gap_s_upper"])
    dt = float(rate["dt_s"])
    cadence = 0.1
    # Binary32 dt is slightly below 5 ms, so ceil(0.1/dt) can be 21 rather
    # than 20. The dynamic source certificate deliberately adds two further
    # samples of source-independent padding. Derive the admissible upper from
    # that exact construction instead of guessing a decimal threshold.
    derived_gap_guard = up(cadence + 3.0 * dt)
    if not (cadence < active_gap_upper <= derived_gap_guard):
        raise RuntimeError(
            f"active commit gap upper {active_gap_upper} exceeds binary32 cadence guard "
            f"{derived_gap_guard}"
        )

    # Start at the first post-Live active commit, reached no later than the
    # source-certified active gap upper. A second commit may occur strictly
    # after 0.1 s; a third cannot occur before the 0.2 s closure boundary.
    phase_lo = down(cadence)
    phase_hi = up(min(active_gap_upper, PREFIX_HORIZON_S - dt))
    if not phase_lo < phase_hi < PREFIX_HORIZON_S:
        raise RuntimeError("invalid one-change phase cover")

    total_from_live = up(active_gap_upper + PREFIX_HORIZON_S)
    Utrans = _translation_diag_upper(total_from_live, live, dynamic)
    eta_trace = _eta_trace_upper(total_from_live, live, process)
    attenuation = _measurement_attenuation_lower(
        PREFIX_HORIZON_S, dynamic, complete, pe, process, path
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
        "first_post_Live_commit_wait_s_upper": active_gap_upper,
        "binary32_active_commit_gap_guard_s": derived_gap_guard,
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
    failures: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        failures.append("schema/qualification mismatch")
    if d.get("canonical_source") != "COMPLETE_SEA3_NORMAL_LIVE_WORD":
        failures.append("canonical source changed")
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
            failures.append(f"{key} is not true")
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
            failures.append(f"forbidden/open flag {key} changed")
    cover = d.get("parameter_cover", {})
    if cover.get("pass") is not True:
        failures.append(f"full H18 LDLT cover failed: {cover.get('failed_cell')}")
    pivot = cover.get("worst_full_H18_LDLT_pivot_lower")
    if not isinstance(pivot, (int, float)) or not (
        math.isfinite(float(pivot)) and float(pivot) > 0.0
    ):
        failures.append("H18 LDLT pivot lower is not strict")
    if float(d.get("useful_gate", math.nan)) != USEFUL_GATE:
        failures.append("useful gate changed")
    return list(dict.fromkeys(failures))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
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
