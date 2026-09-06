#!/usr/bin/env python3
"""Literal shipping Normal-Live H18/A21 word assembler for canonical OU-III P3.

This module is deliberately *not* another reduced certificate.  It exposes the
actual full-state operations used by shipping Normal Live and applies them to the
same :class:`JointWordState` from ``ou3_sea3_full_word_riccati_backend``.

Per valid IMU sample the source order is

  1. commit a previously staged tuner candidate (if armed),
  2. prediction with the committed tau/sigma schedule,
  3. apply any queued PSD a_w covariance-floor increment,
  4. service every due S=0 pseudo update,
  5. apply the required accelerometer Joseph update,
  6. update the measurement-only front end / tuner candidate,
  7. stage the next-sample commit and periodic a_w-floor request as due.

Magnetometer Joseph updates are asynchronous external calls and are inserted at
their actual source-reachable positions between IMU samples.  The theorem's
finite-window PE premise constrains their accepted source language; hardware ODR
is not substituted for that premise.

The routines below pack the *complete* H18/A21 F,Q,H,R matrices.  In particular
accelerometer H retains theta<->a_w and, in A mode, theta<->a_w<->b_a coupling;
S=0 and magnetometer updates act through the complete covariance and hence every
cross covariance via the backend Joseph recursion.

No scalar information beta, determinant/trace gate, blockwise contraction,
independent tau/sigma/R_S/T_S extrema product, arbitrary P0 rectangle, or old
P2 predecessor/history graph is available here.  P3 promotion remains false
until a source-reachable 3 s family has actually executed this API and the full
Omega-delta*P LDLT closes for both H18 and A21.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
import json
import math
from pathlib import Path
from typing import Sequence

import ou3_p3_pseudo_scheduler_progress_certificate as SCHED
import ou3_sea3_dynamic_source_certificate as DYNAMIC
import ou3_sea3_full_word_riccati_backend as BACKEND
import ou3_sea3_live_covariance_seed as LIVE
import ou3_sea3_p3_full_preconditions as FULL
import ou3_sea3_windowed_vector_pe as PE
from ou3_interval import (
    Interval,
    IntervalMatrix,
    matrix_identity,
    matrix_point,
    symmetric_positive_definite_ldlt,
)

REPO = Path(__file__).resolve().parents[2]
WRAPPER = REPO / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
MEKF = REPO / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
DEFAULT_DOMAIN = REPO / "tools" / "ou3_proof_operating_domain.json"
SCHEMA = 1
QUALIFICATION = "OU3_SEA3_LITERAL_FULL_NORMAL_LIVE_WORD_ASSEMBLER"
USEFUL_GATE = 1.0e-18

H_DIM = 18
A_DIM = 21
OFF_TH = 0
OFF_BG = 3
OFF_V = 6
OFF_P = 9
OFF_S = 12
OFF_AW = 15
OFF_BA = 18


def _I(x: float) -> Interval:
    return Interval.point(float(x))


def _zero(rows: int, cols: int) -> IntervalMatrix:
    return [[_I(0.0) for _ in range(cols)] for _ in range(rows)]


def _shape(A: Sequence[Sequence[Interval]]) -> tuple[int, int]:
    r = len(A)
    c = len(A[0]) if r else 0
    if any(len(row) != c for row in A):
        raise ValueError("ragged interval matrix")
    return r, c


def _body(text: str, signature: str, next_marker: str | None = None) -> str:
    start = text.find(signature)
    if start < 0:
        return ""
    if next_marker is not None:
        end = text.find(next_marker, start + len(signature))
        if end > start:
            return text[start:end]
    # Brace-counting fallback: begin at first opening brace after signature.
    b = text.find("{", start)
    if b < 0:
        return ""
    depth = 0
    for i in range(b, len(text)):
        ch = text[i]
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                return text[start:i + 1]
    return ""


def shipping_event_order_parity() -> dict[str, bool]:
    """Bind assembler ordering to executable method bodies, not file comments."""
    wrapper = WRAPPER.read_text(encoding="utf-8")
    mekf = MEKF.read_text(encoding="utf-8")
    core = _body(wrapper, "void updateCore_(")
    tu = _body(mekf, "::time_update(")
    mag = _body(wrapper, "void updateMag(")
    commit = _body(wrapper, "void apply_pending_online_tune_()")
    floor_tick = _body(wrapper, "void periodic_aw_cov_sync_tick_()")

    p_commit = core.find("apply_pending_online_tune_();")
    p_predict = core.find("mekf_->time_update(gyro, dt);")
    p_acc = core.find("mekf_->measurement_update_acc_only(acc_in, tempC);")
    p_tuner = core.find("update_tuner(dt, a_vert_measurement, tuner_frequency_hz_());")
    p_period = core.find("wave_period_.update(dt, wave_period_input_ms2_(direction_accel));")
    p_floor_tick = core.find("periodic_aw_cov_sync_tick_();")

    p_floor_apply = tu.find("apply_pending_aw_covariance_inflation_();")
    p_s_due = tu.find("periodic_update_due(")
    p_s_apply = tu.find("applyIntegralZeroPseudoMeas();")

    return {
        "updateCore_body_scoped": bool(core),
        "time_update_body_scoped": bool(tu),
        "mag_update_body_scoped": bool(mag),
        "pending_commit_body_scoped": bool(commit),
        "periodic_floor_tick_body_scoped": bool(floor_tick),
        "commit_precedes_prediction": 0 <= p_commit < p_predict,
        "prediction_precedes_accelerometer": 0 <= p_predict < p_acc,
        "accelerometer_precedes_tuner_candidate_update": 0 <= p_acc < p_tuner,
        "tuner_candidate_update_precedes_wave_period_current_sample_update": (
            0 <= p_tuner < p_period
        ),
        "periodic_floor_request_occurs_after_current_accelerometer": (
            0 <= p_acc < p_floor_tick
        ),
        "pending_aw_floor_applied_inside_prediction_before_S_service": (
            0 <= p_floor_apply < p_s_due < p_s_apply
        ),
        "pending_commit_applies_ou_then_RS": (
            commit.find("apply_ou_tune_(false);") >= 0
            and commit.find("apply_RS_tune_();") > commit.find("apply_ou_tune_(false);")
        ),
        "magnetometer_is_separate_external_call": (
            "measurement_update_mag_only" in mag
            and "measurement_update_mag_only" not in core
        ),
        "normal_live_accelerometer_is_full_Joseph": (
            "joseph_update3_(K, S_mat, PCt);" in _body(mekf, "::measurement_update_acc_only(")
        ),
        "S_zero_is_full_Joseph": (
            "joseph_update3_(K, S_mat, PCt);" in _body(mekf, "::applyIntegralZeroPseudoMeas()")
        ),
    }


def state_dimension(mode: str) -> int:
    if mode == "H":
        return H_DIM
    if mode == "A":
        return A_DIM
    raise ValueError("mode must be H or A")


def _embed3(n: int, offset: int, block: Sequence[Sequence[Interval]]) -> IntervalMatrix:
    if _shape(block) != (3, 3):
        raise ValueError("3x3 block required")
    M = _zero(n, n)
    for i in range(3):
        for j in range(3):
            M[offset + i][offset + j] = block[i][j]
    return M


def _skew(v: Sequence[Interval]) -> IntervalMatrix:
    if len(v) != 3:
        raise ValueError("three-vector required")
    z = _I(0.0)
    return [
        [z, -v[2], v[1]],
        [v[2], z, -v[0]],
        [-v[1], v[0], z],
    ]


def pack_prediction(
    mode: str,
    F_att_bg: Sequence[Sequence[Interval]],
    Q_att_bg: Sequence[Sequence[Interval]],
    F_translation: Sequence[Sequence[Interval]],
    Q_translation: Sequence[Sequence[Interval]],
    *,
    phi_ba: Interval | None = None,
    Q_ba: Sequence[Sequence[Interval]] | None = None,
) -> tuple[IntervalMatrix, IntervalMatrix]:
    """Pack shipping block-diagonal same-mode prediction without dropping blocks."""
    n = state_dimension(mode)
    if _shape(F_att_bg) != (6, 6) or _shape(Q_att_bg) != (6, 6):
        raise ValueError("attitude/gyro-bias F/Q must be 6x6")
    if _shape(F_translation) != (12, 12) or _shape(Q_translation) != (12, 12):
        raise ValueError("translation F/Q must be 12x12")
    F = matrix_identity(n)
    Q = _zero(n, n)
    for i in range(6):
        for j in range(6):
            F[i][j] = F_att_bg[i][j]
            Q[i][j] = Q_att_bg[i][j]
    for i in range(12):
        for j in range(12):
            F[OFF_V + i][OFF_V + j] = F_translation[i][j]
            Q[OFF_V + i][OFF_V + j] = Q_translation[i][j]
    if mode == "A":
        if phi_ba is None or Q_ba is None or _shape(Q_ba) != (3, 3):
            raise ValueError("active A mode requires exact b_a GM F/Q")
        for i in range(3):
            F[OFF_BA + i][OFF_BA + i] = phi_ba
            for j in range(3):
                Q[OFF_BA + i][OFF_BA + j] = Q_ba[i][j]
    return F, Q


def H_S_zero(mode: str) -> IntervalMatrix:
    n = state_dimension(mode)
    H = _zero(3, n)
    for i in range(3):
        H[i][OFF_S + i] = _I(1.0)
    return H


def R_S_zero(rs_std_xyz: Sequence[Interval]) -> IntervalMatrix:
    if len(rs_std_xyz) != 3:
        raise ValueError("three applied R_S standard deviations required")
    R = _zero(3, 3)
    for i, s in enumerate(rs_std_xyz):
        if not s.lo > 0.0:
            raise ValueError("R_S std must stay positive")
        R[i][i] = s.square()
    return R


def H_accelerometer(
    mode: str,
    f_cog_body: Sequence[Interval],
    R_wb: Sequence[Sequence[Interval]],
) -> IntervalMatrix:
    """Exact default no-lever-arm accelerometer Jacobian packing.

    H_theta=-[f]_x, H_aw=R_wb, and H_ba=I in A mode.  H_bg is exactly zero
    because the canonical configured runtime has the lever arm disabled.
    """
    n = state_dimension(mode)
    if len(f_cog_body) != 3 or _shape(R_wb) != (3, 3):
        raise ValueError("accelerometer geometry must be 3-vector plus 3x3 rotation")
    H = _zero(3, n)
    Jth = _skew(f_cog_body)
    for i in range(3):
        for j in range(3):
            H[i][OFF_TH + j] = -Jth[i][j]
            H[i][OFF_AW + j] = R_wb[i][j]
        if mode == "A":
            H[i][OFF_BA + i] = _I(1.0)
    return H


def H_magnetometer(mode: str, m_body: Sequence[Interval]) -> IntervalMatrix:
    n = state_dimension(mode)
    if len(m_body) != 3:
        raise ValueError("magnetic body vector must have length three")
    H = _zero(3, n)
    Jth = _skew(m_body)
    for i in range(3):
        for j in range(3):
            H[i][OFF_TH + j] = -Jth[i][j]
    return H


def diagonal_R(std_xyz: Sequence[float]) -> IntervalMatrix:
    if len(std_xyz) != 3:
        raise ValueError("three measurement standard deviations required")
    R = _zero(3, 3)
    for i, s in enumerate(std_xyz):
        x = float(s)
        if not (math.isfinite(x) and x > 0.0):
            raise ValueError("measurement std must be finite positive")
        R[i][i] = Interval.outward_bounds(x, x).square()
    return R


def aw_floor_increment(mode: str, Delta_aw: Sequence[Sequence[Interval]]) -> IntervalMatrix:
    """Embed the source-computed PSD Pi_+(Sigma_aw-P_awaw^-) increment exactly."""
    n = state_dimension(mode)
    if _shape(Delta_aw) != (3, 3):
        raise ValueError("a_w floor increment must be 3x3")
    M = _zero(n, n)
    for i in range(3):
        for j in range(3):
            M[OFF_AW + i][OFF_AW + j] = Delta_aw[i][j]
    return M


@dataclass
class LiteralWordState:
    mode: str
    riccati: BACKEND.JointWordState
    imu_samples: int = 0
    S_updates: int = 0
    accel_updates: int = 0
    mag_updates: int = 0
    aw_floor_applications: int = 0
    event_log: list[str] = field(default_factory=list)

    @property
    def dimension(self) -> int:
        return state_dimension(self.mode)


def initialize_word(mode: str, P0: Sequence[Sequence[Interval]]) -> LiteralWordState:
    n = state_dimension(mode)
    if _shape(P0) != (n, n):
        raise ValueError(f"{mode} P0 must be {n}x{n}")
    return LiteralWordState(mode=mode, riccati=BACKEND.initialize(P0))


def apply_imu_sample(
    word: LiteralWordState,
    *,
    F: Sequence[Sequence[Interval]],
    Q: Sequence[Sequence[Interval]],
    f_cog_body: Sequence[Interval],
    R_wb: Sequence[Sequence[Interval]],
    Racc: Sequence[Sequence[Interval]],
    due_S: bool,
    rs_std_xyz: Sequence[Interval] | None,
    Delta_aw: Sequence[Sequence[Interval]] | None,
) -> None:
    """Execute one shipping IMU Kalman sequence after the source commit decision.

    Candidate/front-end/staging state is advanced by the source-family assembler
    around this call; it never changes F,Q,R_S used by this same sample.
    """
    BACKEND.predict(word.riccati, F, Q)
    word.event_log.append("prediction")
    if Delta_aw is not None:
        BACKEND.add_psd_floor(word.riccati, aw_floor_increment(word.mode, Delta_aw))
        word.aw_floor_applications += 1
        word.event_log.append("aw_floor")
    if due_S:
        if rs_std_xyz is None:
            raise ValueError("due S update requires the actual applied per-axis R_S")
        BACKEND.joseph_measurement(word.riccati, H_S_zero(word.mode), R_S_zero(rs_std_xyz))
        word.S_updates += 1
        word.event_log.append("S_zero")
    BACKEND.joseph_measurement(
        word.riccati,
        H_accelerometer(word.mode, f_cog_body, R_wb),
        Racc,
    )
    word.accel_updates += 1
    word.imu_samples += 1
    word.event_log.append("accelerometer")


def apply_magnetometer(
    word: LiteralWordState,
    *,
    m_body: Sequence[Interval],
    Rmag: Sequence[Sequence[Interval]],
) -> None:
    BACKEND.joseph_measurement(word.riccati, H_magnetometer(word.mode, m_body), Rmag)
    word.mag_updates += 1
    word.event_log.append("magnetometer")


def certify_literal_endpoint(word: LiteralWordState, delta: float = USEFUL_GATE) -> dict:
    cert = BACKEND.certify_contraction(word.riccati, delta)
    return {
        **cert,
        "mode": word.mode,
        "imu_samples": word.imu_samples,
        "S_updates": word.S_updates,
        "accelerometer_updates": word.accel_updates,
        "magnetometer_updates": word.mag_updates,
        "aw_floor_applications": word.aw_floor_applications,
    }


def _point_kernel_self_test(mode: str) -> dict:
    """Exercise every literal operation type; never used as a theorem result."""
    n = state_dimension(mode)
    P0 = matrix_point([[2.0 if i == j else 0.0 for j in range(n)] for i in range(n)])
    w = initialize_word(mode, P0)
    Faa = matrix_identity(6)
    Qaa = matrix_point([[0.02 if i == j else 0.0 for j in range(6)] for i in range(6)])
    Fll = matrix_identity(12)
    # Preserve the actual integrated-chain structure in the smoke operation.
    h = 0.005
    for a in range(3):
        Fll[a][9 + a] = _I(h)
        Fll[3 + a][a] = _I(h)
        Fll[6 + a][3 + a] = _I(h)
    Qll = matrix_point([[0.01 if i == j else 0.0 for j in range(12)] for i in range(12)])
    if mode == "A":
        F, Q = pack_prediction(
            mode, Faa, Qaa, Fll, Qll,
            phi_ba=_I(0.999999),
            Q_ba=matrix_point([[1e-8 if i == j else 0.0 for j in range(3)] for i in range(3)]),
        )
    else:
        F, Q = pack_prediction(mode, Faa, Qaa, Fll, Qll)
    f = [_I(0.0), _I(0.0), _I(-9.80665)]
    Rwb = matrix_identity(3)
    Racc = diagonal_R([0.2, 0.2, 0.2])
    Rmag = diagonal_R([0.3, 0.3, 0.3])
    Delta = matrix_point([[0.001 if i == j else 0.0 for j in range(3)] for i in range(3)])
    apply_imu_sample(
        w, F=F, Q=Q, f_cog_body=f, R_wb=Rwb, Racc=Racc,
        due_S=True, rs_std_xyz=[_I(0.72), _I(0.72), _I(1.0)], Delta_aw=Delta,
    )
    apply_magnetometer(w, m_body=[_I(20.0), _I(0.0), _I(40.0)], Rmag=Rmag)
    return {
        "mode": mode,
        "dimension": n,
        "decomposition_identity_enclosed": BACKEND.decomposition_identity_enclosed(w.riccati),
        "event_log": w.event_log,
        "every_literal_event_type_executed": (
            w.imu_samples == 1 and w.S_updates == 1 and w.accel_updates == 1
            and w.mag_updates == 1 and w.aw_floor_applications == 1
        ),
        "kernel_self_test_only_not_P3": True,
    }


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    path = Path(domain_path).resolve()
    full = FULL.build(path)
    dynamic = DYNAMIC.build(path)
    live = LIVE.build(path)
    sched = SCHED.build(path)
    pe = PE.build(path)
    bad = {
        "full": FULL.validate(full),
        "dynamic": DYNAMIC.validate(dynamic),
        "live": LIVE.validate(live),
        "scheduler": SCHED.validate(sched),
        "PE": PE.validate(pe),
        "backend": BACKEND.validate_backend(),
    }
    bad = {k: v for k, v in bad.items() if v}
    if bad:
        raise RuntimeError(f"literal full-word prerequisites failed: {bad}")

    parity = shipping_event_order_parity()
    horizon = float(full["final_numeric_contract"]["common_word_horizon_s"])
    samples = int(full["word"]["samples_upper"])
    max_s_gap = int(sched["certified_uniform_max_gap_samples"])
    guaranteed_s_updates = max(0, (samples - 1) // max_s_gap)
    # enterLive/apply_ou_tune_(true) seats the covariance and last sync time.
    # The periodic tick requests another floor only after >0.1 s, and the
    # request is applied by the following prediction.  At 200 Hz that is one
    # application per 21-sample service interval after the first request.
    floor_service_samples = int(math.floor(0.1 / float(sched["dt_binary32_s"]))) + 1
    guaranteed_floor_applications = max(0, (samples - 1) // floor_service_samples)

    measurement_runtime = full["measurement_runtime"]
    htest = _point_kernel_self_test("H")
    atest = _point_kernel_self_test("A")

    return {
        "schema": SCHEMA,
        "qualification": QUALIFICATION,
        "canonical_architecture": "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD",
        "trajectory_replay_used": False,
        "filter_changed": False,
        "declared_domain_shrunk": False,
        "shipping_event_order_parity": parity,
        "shipping_event_order_parity_pass": all(parity.values()),
        "state_orders": {
            "H18": ["theta", "b_g", "v", "p", "S", "a_w"],
            "A21": ["theta", "b_g", "v", "p", "S", "a_w", "b_a"],
        },
        "word_horizon_s": horizon,
        "imu_samples_upper": samples,
        "every_valid_imu_sample_requires_prediction": True,
        "every_valid_imu_sample_requires_accelerometer_Joseph": True,
        "S_scheduler_is_executed_not_replaced_by_selected_four": True,
        "certified_S_gap_samples_upper": max_s_gap,
        "guaranteed_S_updates_lower_over_word": guaranteed_s_updates,
        "periodic_aw_floor_request_and_next_prediction_application_retained": True,
        "aw_floor_service_samples": floor_service_samples,
        "guaranteed_aw_floor_applications_lower_over_word": guaranteed_floor_applications,
        "magnetometer_is_asynchronous_external_event_family": True,
        "hardware_magnetometer_ODR_used_as_PE_recurrence": False,
        "PE_recurrence_window_s": pe["declared_normal_live_PE"]["recurrence_window_s"],
        "PE_accumulated_information_lower": pe["eta6_information"]["alpha_6_information_lower"],
        "same_source_state_must_feed_F_Q_pseudo_period_and_applied_RS": True,
        "candidate_update_occurs_after_current_Kalman_innovation": True,
        "candidate_commit_only_affects_following_sample": True,
        "full_front_end_state_manifest": full["front_end_state_manifest"],
        "measurement_runtime": measurement_runtime,
        "live_covariance_seed": live,
        "numeric_execution_api": {
            "prediction": "pack_prediction + BACKEND.predict",
            "S_zero": "H_S_zero/R_S_zero + BACKEND.joseph_measurement",
            "accelerometer": "H_accelerometer + BACKEND.joseph_measurement",
            "magnetometer": "H_magnetometer + BACKEND.joseph_measurement",
            "aw_floor": "aw_floor_increment + BACKEND.add_psd_floor",
            "final_gate": "BACKEND.certify_contraction on full H18/A21 Omega-delta*P",
        },
        "H_kernel_self_test": htest,
        "A_kernel_self_test": atest,
        "no_reduced_promotion_routes": {
            "D_W_L_W_split": False,
            "prior_free_batch_as_canonical_gate": False,
            "one_sample_process_ratio": False,
            "scalar_information_beta": False,
            "determinant_trace_scalarization": False,
            "blockwise_minimum_ratio": False,
            "independent_tau_sigma_RS_TS_extrema_product": False,
            "arbitrary_P0_rectangle": False,
            "old_P2_800_state_graph": False,
            "predecessor_history_enumeration": False,
        },
        # Fail closed: the API is literal, but theorem promotion waits for the
        # forward source-coordinate enclosure to execute all samples/events.
        "SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED": False,
        "FULL_H18_WORD_EXECUTED": False,
        "FULL_A21_WORD_EXECUTED": False,
        "FULL_H18_A21_LDLT_CLOSED": False,
        "P3_CANONICAL_PASS": False,
        "first_open_obligation": (
            "materialize the forward source-reachable 3 s family carrying the complete "
            "front-end/tuner/commit/scheduler/vector geometry state, then execute every "
            "shipping event above through the literal H18/A21 API"
        ),
    }


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("schema") != SCHEMA or d.get("qualification") != QUALIFICATION:
        f.append("schema/qualification mismatch")
    if d.get("canonical_architecture") != "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD":
        f.append("literal assembler is not bound to canonical P3")
    if d.get("shipping_event_order_parity_pass") is not True:
        f.append("shipping event order parity failed")
    if not all(d.get("shipping_event_order_parity", {}).values()):
        f.append("one or more method-body source-order checks failed")
    for key in (
        "every_valid_imu_sample_requires_prediction",
        "every_valid_imu_sample_requires_accelerometer_Joseph",
        "S_scheduler_is_executed_not_replaced_by_selected_four",
        "periodic_aw_floor_request_and_next_prediction_application_retained",
        "magnetometer_is_asynchronous_external_event_family",
        "same_source_state_must_feed_F_Q_pseudo_period_and_applied_RS",
        "candidate_update_occurs_after_current_Kalman_innovation",
        "candidate_commit_only_affects_following_sample",
    ):
        if d.get(key) is not True:
            f.append(f"{key} is not true")
    if d.get("hardware_magnetometer_ODR_used_as_PE_recurrence") is not False:
        f.append("magnetometer ODR shortcut re-entered")
    if int(d.get("imu_samples_upper", 0)) < 600:
        f.append("literal word no longer spans the 3 s sample family")
    if int(d.get("guaranteed_S_updates_lower_over_word", 0)) < 4:
        f.append("literal word lost recurrent S correction")
    if int(d.get("guaranteed_aw_floor_applications_lower_over_word", 0)) <= 0:
        f.append("literal word lost recurrent a_w floor events")
    for mode in ("H", "A"):
        row = d.get(f"{mode}_kernel_self_test", {})
        if row.get("decomposition_identity_enclosed") is not True:
            f.append(f"{mode} literal kernel lost P/Psi/Omega identity")
        if row.get("every_literal_event_type_executed") is not True:
            f.append(f"{mode} literal kernel did not execute every event type")
        if row.get("kernel_self_test_only_not_P3") is not True:
            f.append(f"{mode} kernel self-test was misclassified")
    reduced = d.get("no_reduced_promotion_routes", {})
    if not reduced or any(v is not False for v in reduced.values()):
        f.append("a reduced/dead-end promotion route re-entered")
    for key in (
        "trajectory_replay_used", "filter_changed", "declared_domain_shrunk",
        "SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED", "FULL_H18_WORD_EXECUTED",
        "FULL_A21_WORD_EXECUTED", "FULL_H18_A21_LDLT_CLOSED", "P3_CANONICAL_PASS",
    ):
        if d.get(key) is not False:
            f.append(f"open literal assembler flag {key} must remain false until execution closes")
    if not d.get("first_open_obligation"):
        f.append("literal assembler does not name its first unresolved source enclosure")
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
        "shipping_event_order": d["shipping_event_order_parity_pass"],
        "word_horizon_s": d["word_horizon_s"],
        "samples": d["imu_samples_upper"],
        "S_updates_lower": d["guaranteed_S_updates_lower_over_word"],
        "aw_floor_applications_lower": d["guaranteed_aw_floor_applications_lower_over_word"],
        "H_kernel": d["H_kernel_self_test"],
        "A_kernel": d["A_kernel_self_test"],
        "P3_CANONICAL_PASS": d["P3_CANONICAL_PASS"],
        "open": d["first_open_obligation"],
        "validation_failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
