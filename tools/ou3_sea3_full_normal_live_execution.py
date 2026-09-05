#!/usr/bin/env python3
"""Literal 3 s shipping-order H18/A21 execution for canonical OU-III P3.

This module is an execution engine for ``ou3_sea3_full_normal_live_word``.  It
does not define a second proof architecture.  It drives the existing exact
``P/Psi/Omega`` backend through every event in one complete Normal-Live word.

The point word used here is deliberately source/domain consistent and is only an
execution witness for the assembler machinery:

* 200 Hz effective binary32 sample time;
* a stationary SpectralMSE operating point whose target equals the committed
  tau/sigma/R_S point, so the actual candidate EMA/commit code path executes
  without inventing source jumps;
* the shipping progress-preserving pseudo scheduler;
* every due S=0 update with the committed per-axis R_S;
* every accelerometer update;
* periodic state-dependent a_w covariance floors;
* two asynchronous accepted magnetic events in the paper's [0,W] and [2W,3W]
  recurrence cells, with physical vectors satisfying the declared PE domain.

The word is *not* a source-family certificate.  In particular, executing this
point word cannot set ``SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED`` or promote
P3.  Canonical promotion still requires the same literal execution over a
validated enclosure of every admitted source/front-end/event path.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import struct
from typing import Any

import ou3_p3_pseudo_scheduler_progress_certificate as SCHED
import ou3_p3_scaled_process as SCALED_Q
import ou3_sea3_full_word_riccati_backend as BACKEND
from ou3_interval import Interval, IntervalMatrix, matrix_identity, matrix_point

USEFUL_GATE = 1.0e-18
G = 9.80665
TAU_BA = 5000.0
Q_BA_DENSITY = 2.5e-7
Q_BG_DENSITY = 1.0e-11
SIGMA_G = 0.01
TAU0 = 1.1
RS0 = 0.5
SIGMA_COEFF = 0.9
CJ = 0.0538
R_A_REDUCED = 0.0148 * 0.0148 * 0.005
PSEUDO_RATIO = 0.015 / 1.1
PSEUDO_MIN = 0.005
PSEUDO_MAX = 0.15
RS_XY = 0.72


def _f32(x: float) -> float:
    return struct.unpack("!f", struct.pack("!f", float(x)))[0]


def _I(x: float) -> Interval:
    return Interval.point(float(x))


def _zero(r: int, c: int) -> IntervalMatrix:
    return [[_I(0.0) for _ in range(c)] for _ in range(r)]


def _mid(x: Interval) -> float:
    return 0.5 * (x.lo + x.hi)


def _matrix_mid(A: IntervalMatrix) -> list[list[float]]:
    return [[_mid(x) for x in row] for row in A]


def _spectral_mse_sigma_for_rs(tau: float, rs: float) -> float:
    """Solve the deployed SpectralMSE target exactly for sigma at a point."""
    ts = min(max(PSEUDO_RATIO * tau, PSEUDO_MIN), PSEUDO_MAX)
    qpow = (2.0 * R_A_REDUCED) ** (1.0 / 14.0)
    coeff = CJ * qpow * ((tau ** 4) / SIGMA_COEFF) ** (6.0 / 7.0) / math.sqrt(ts)
    sigma = (rs / coeff) ** (7.0 / 6.0)
    if not (0.05 <= sigma <= 4.0):
        raise RuntimeError("source-consistent SpectralMSE equilibrium left deployed sigma box")
    return sigma


SIGMA0 = _spectral_mse_sigma_for_rs(TAU0, RS0)


@dataclass
class PointSource:
    dt: float
    candidate_tau: float
    candidate_sigma: float
    candidate_rs: float
    active_tau: float
    active_sigma: float
    active_rs: float
    pseudo_period: float
    pseudo_elapsed: float = 0.0
    time: float = 0.0
    last_adapt_time: float = 0.0
    last_floor_time: float = 0.0
    pending_commit: bool = False
    floor_pending: bool = False
    commits: int = 0
    floor_requests: int = 0

    @classmethod
    def equilibrium(cls, dt: float) -> "PointSource":
        sched = SCHED.BASE.source_schedule()
        p = SCHED._period_from_tau(TAU0, sched)
        return cls(
            dt=_f32(dt), candidate_tau=TAU0, candidate_sigma=SIGMA0,
            candidate_rs=RS0, active_tau=TAU0, active_sigma=SIGMA0,
            active_rs=RS0, pseudo_period=p,
        )

    def commit_if_pending(self) -> None:
        if not self.pending_commit:
            return
        sched = SCHED.BASE.source_schedule()
        self.active_tau = self.candidate_tau
        self.active_sigma = self.candidate_sigma
        self.active_rs = self.candidate_rs
        new_period = SCHED._period_from_tau(self.active_tau, sched)
        self.pseudo_elapsed = SCHED._retarget(self.pseudo_elapsed, new_period)
        self.pseudo_period = new_period
        self.pending_commit = False
        self.commits += 1

    def advance_scheduler(self) -> bool:
        due, elapsed = SCHED._due(self.dt, self.pseudo_period, self.pseudo_elapsed)
        self.pseudo_elapsed = elapsed
        return bool(due)

    def finish_sample(self) -> None:
        """Execute the source-only part after the current Kalman innovation."""
        self.time += float(self.dt)
        sea_time = self.candidate_tau
        adapt_sec = max(0.05, min(35.0, 0.40 * sea_time))
        alpha = 1.0 - math.exp(-float(self.dt) / adapt_sec)
        rs_sec = max(0.05, min(35.0, 1.5 * self.candidate_tau))
        alpha_rs = 1.0 - math.exp(-float(self.dt) / rs_sec)
        self.candidate_tau += alpha * (TAU0 - self.candidate_tau)
        self.candidate_sigma += alpha * (SIGMA0 - self.candidate_sigma)
        self.candidate_rs += alpha_rs * (RS0 - self.candidate_rs)
        if self.time - self.last_adapt_time > 0.1:
            self.pending_commit = True
            self.last_adapt_time = self.time
        if self.time - self.last_floor_time > 0.1:
            self.floor_pending = True
            self.last_floor_time = self.time
            self.floor_requests += 1


def _attitude_prediction(dt: float) -> tuple[IntervalMatrix, IntervalMatrix]:
    """Exact shipping attitude/bg prediction at the admitted omega=0 point."""
    F = matrix_identity(6)
    Q = _zero(6, 6)
    h = float(dt)
    qg = SIGMA_G * SIGMA_G
    qb = Q_BG_DENSITY
    for i in range(3):
        F[i][3 + i] = _I(h)
        Q[i][i] = _I(qg * h + qb * h ** 3 / 3.0)
        Q[i][3 + i] = _I(qb * h * h / 2.0)
        Q[3 + i][i] = _I(qb * h * h / 2.0)
        Q[3 + i][3 + i] = _I(qb * h)
    return F, Q


def _scaled_process_matrix(x: Interval) -> IntervalMatrix:
    if x.hi < SCALED_Q.BRANCH_X:
        B = SCALED_Q.small_normalized_matrix(x)
        return [[x * B[i][j] for j in range(4)] for i in range(4)]
    if x.lo >= SCALED_Q.BRANCH_X and x.hi <= SCALED_Q.NEAR_EXACT_SERIES_MAX_X:
        B = SCALED_Q.near_exact_normalized_matrix(x)
        return [[x * B[i][j] for j in range(4)] for i in range(4)]
    return SCALED_Q._large_scaled_matrix(x)


def _translation_prediction(tau: float, sigma: float, dt: float) -> tuple[IntervalMatrix, IntervalMatrix]:
    """Shipping integrated-OU F/Q at one committed point, all three axes."""
    h = float(dt)
    x = h / float(tau)
    a = math.exp(-x)
    em1 = math.expm1(-x)
    phi_va = -tau * em1
    phi_pa = tau * tau * (x + em1)
    phi_Sa = tau ** 3 * (0.5 * x * x - x - em1)
    F4 = [
        [1.0, 0.0, 0.0, phi_va],
        [h, 1.0, 0.0, phi_pa],
        [0.5 * h * h, h, 1.0, phi_Sa],
        [0.0, 0.0, 0.0, a],
    ]
    xi = Interval.outward_bounds(x, x)
    qscaled = _scaled_process_matrix(xi)
    d = [
        Interval.outward_bounds(sigma * h, sigma * h),
        Interval.outward_bounds(sigma * h * h, sigma * h * h),
        Interval.outward_bounds(sigma * h ** 3, sigma * h ** 3),
        Interval.outward_bounds(sigma, sigma),
    ]
    Q4 = [[d[i] * qscaled[i][j] * d[j] for j in range(4)] for i in range(4)]
    F = matrix_identity(12)
    Q = _zero(12, 12)
    for axis in range(3):
        idx = [axis, 3 + axis, 6 + axis, 9 + axis]
        for i in range(4):
            for j in range(4):
                F[idx[i]][idx[j]] = _I(F4[i][j])
                Q[idx[i]][idx[j]] = Q4[i][j]
    return F, Q


def _active_ba_prediction(dt: float) -> tuple[Interval, IntervalMatrix]:
    h = float(dt)
    phi = math.exp(-h / TAU_BA)
    qd = -0.5 * TAU_BA * math.expm1(-2.0 * h / TAU_BA)
    return _I(phi), matrix_point([
        [Q_BA_DENSITY * qd, 0.0, 0.0],
        [0.0, Q_BA_DENSITY * qd, 0.0],
        [0.0, 0.0, Q_BA_DENSITY * qd],
    ])


def _point_live_P0(word: Any, mode: str, sigma: float) -> IntervalMatrix:
    n = word.state_dimension(mode)
    P = [[0.0 for _ in range(n)] for _ in range(n)]
    tilt2, yaw2 = 0.035 ** 2, 0.087 ** 2
    P[word.OFF_TH][word.OFF_TH] = tilt2
    P[word.OFF_TH + 1][word.OFF_TH + 1] = tilt2
    P[word.OFF_TH + 2][word.OFF_TH + 2] = yaw2
    for i in range(3):
        P[word.OFF_BG + i][word.OFF_BG + i] = 1e-6
        P[word.OFF_V + i][word.OFF_V + i] = 1.0
        P[word.OFF_P + i][word.OFF_P + i] = 400.0
        P[word.OFF_S + i][word.OFF_S + i] = 2500.0
        P[word.OFF_AW + i][word.OFF_AW + i] = sigma * sigma
        if mode == "A":
            P[word.OFF_BA + i][word.OFF_BA + i] = 1.6e-5
    return matrix_point(P)


def _jacobi_eigh3(A: list[list[float]]) -> tuple[list[float], list[list[float]]]:
    """Small deterministic point-only symmetric eigensolver for Pi_+ audit."""
    V = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
    M = [[0.5 * (A[i][j] + A[j][i]) for j in range(3)] for i in range(3)]
    for _ in range(32):
        p, q = 0, 1
        best = abs(M[p][q])
        for i, j in ((0, 2), (1, 2)):
            if abs(M[i][j]) > best:
                p, q, best = i, j, abs(M[i][j])
        if best <= 1e-14 * max(1.0, max(abs(M[i][i]) for i in range(3))):
            break
        app, aqq, apq = M[p][p], M[q][q], M[p][q]
        phi = 0.5 * math.atan2(2.0 * apq, aqq - app)
        c, s = math.cos(phi), math.sin(phi)
        for k in range(3):
            mkp, mkq = M[k][p], M[k][q]
            M[k][p] = c * mkp - s * mkq
            M[k][q] = s * mkp + c * mkq
        for k in range(3):
            mpk, mqk = M[p][k], M[q][k]
            M[p][k] = c * mpk - s * mqk
            M[q][k] = s * mpk + c * mqk
        M[p][q] = M[q][p] = 0.0
        for k in range(3):
            vkp, vkq = V[k][p], V[k][q]
            V[k][p] = c * vkp - s * vkq
            V[k][q] = s * vkp + c * vkq
    return [M[i][i] for i in range(3)], V


def _point_aw_floor(word_state: Any, word: Any, target_var: float) -> IntervalMatrix:
    """Point transcription of shipping Delta=Pi_+(Sigma-P_awaw^-)."""
    P = _matrix_mid(word_state.riccati.P)
    A = [[0.0 for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            A[i][j] = (target_var if i == j else 0.0) - P[word.OFF_AW + i][word.OFF_AW + j]
    vals, V = _jacobi_eigh3(A)
    pos = [max(0.0, x) for x in vals]
    D = [[sum(V[i][k] * pos[k] * V[j][k] for k in range(3)) for j in range(3)] for i in range(3)]
    for i in range(3):
        for j in range(i + 1, 3):
            x = 0.5 * (D[i][j] + D[j][i])
            D[i][j] = D[j][i] = x
    return matrix_point(D)


def _prediction(word: Any, mode: str, src: PointSource) -> tuple[IntervalMatrix, IntervalMatrix]:
    Faa, Qaa = _attitude_prediction(src.dt)
    Fll, Qll = _translation_prediction(src.active_tau, src.active_sigma, src.dt)
    if mode == "A":
        phi, Qba = _active_ba_prediction(src.dt)
        return word.pack_prediction(mode, Faa, Qaa, Fll, Qll, phi_ba=phi, Q_ba=Qba)
    return word.pack_prediction(mode, Faa, Qaa, Fll, Qll)


def _mag_samples(samples: int) -> set[int]:
    return {min(samples, 100), min(samples, 500)}


def execute_point_word(word: Any, mode: str, horizon_s: float) -> dict:
    """Execute one complete source/domain-consistent point word through backend."""
    sched_cert = SCHED.build()
    dt = float(sched_cert["dt_binary32_s"])
    samples = int(math.ceil(float(horizon_s) / dt))
    src = PointSource.equilibrium(dt)
    P0 = _point_live_P0(word, mode, src.active_sigma)
    w = word.initialize_word(mode, P0)
    Racc = word.diagonal_R([0.2, 0.2, 0.2])
    Rmag = word.diagonal_R([0.3, 0.3, 0.3])
    Rwb = matrix_identity(3)
    f = [_I(0.0), _I(0.0), _I(-G)]
    mag = [_I(30.0), _I(0.0), _I(40.0)]
    mags = _mag_samples(samples)

    decomposition_fail_sample = None
    for k in range(1, samples + 1):
        src.commit_if_pending()
        F, Q = _prediction(word, mode, src)
        Delta = None
        if src.floor_pending:
            Delta = _point_aw_floor(w, word, src.active_sigma * src.active_sigma)
            src.floor_pending = False
        due = src.advance_scheduler()
        rs = [_I(src.active_rs * RS_XY), _I(src.active_rs * RS_XY), _I(src.active_rs)]
        word.apply_imu_sample(
            w, F=F, Q=Q, f_cog_body=f, R_wb=Rwb, Racc=Racc,
            due_S=due, rs_std_xyz=rs, Delta_aw=Delta,
        )
        if k in mags:
            word.apply_magnetometer(w, m_body=mag, Rmag=Rmag)
        if not BACKEND.decomposition_identity_enclosed(w.riccati):
            decomposition_fail_sample = k
            break
        src.finish_sample()

    completed = decomposition_fail_sample is None and w.imu_samples == samples
    cert = word.certify_literal_endpoint(w, USEFUL_GATE) if completed else {
        "pass": False, "Omega_minus_delta_P_ldlt_closed": False,
        "dimension": word.state_dimension(mode), "delta": USEFUL_GATE,
    }
    return {
        "mode": mode,
        "source_word_role": "FULL_LITERAL_POINT_EXECUTION_NOT_SOURCE_FAMILY_CERTIFICATE",
        "source_consistent_equilibrium": {
            "tau_s": src.active_tau, "sigma_aw_mps2": src.active_sigma,
            "R_S_applied": src.active_rs, "pseudo_period_s": src.pseudo_period,
            "specific_force_body_mps2": [0.0, 0.0, -G],
            "magnetic_body_uT": [30.0, 0.0, 40.0],
            "body_rate_rad_s": [0.0, 0.0, 0.0],
            "SpectralMSE_target_equals_applied_RS_by_construction": True,
        },
        "horizon_requested_s": horizon_s,
        "dt_binary32_s": dt,
        "samples_executed": w.imu_samples,
        "samples_required": samples,
        "word_time_s": w.imu_samples * dt,
        "predictions": w.riccati.predictions,
        "accelerometer_updates": w.accel_updates,
        "S_zero_updates": w.S_updates,
        "magnetometer_updates": w.mag_updates,
        "aw_floor_applications": w.aw_floor_applications,
        "candidate_commits": src.commits,
        "aw_floor_requests": src.floor_requests,
        "decomposition_identity_preserved_every_sample": completed,
        "first_decomposition_failure_sample": decomposition_fail_sample,
        "complete_3s_word_executed": completed and w.imu_samples * dt >= 3.0,
        "point_endpoint_contraction": cert,
        "point_endpoint_may_promote_P3": False,
        "event_log_prefix": w.event_log[:16],
        "event_log_suffix": w.event_log[-16:],
    }


def build_execution(word: Any, horizon_s: float) -> dict:
    H = execute_point_word(word, "H", horizon_s)
    A = execute_point_word(word, "A", horizon_s)
    return {
        "qualification": "OU3_SEA3_LITERAL_FULL_WORD_EXECUTION_WITNESS",
        "canonical_architecture": "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD",
        "H18": H, "A21": A,
        "FULL_H18_WORD_EXECUTED": bool(H["complete_3s_word_executed"]),
        "FULL_A21_WORD_EXECUTED": bool(A["complete_3s_word_executed"]),
        "both_modes_exact_joint_backend_executed": bool(
            H["complete_3s_word_executed"] and A["complete_3s_word_executed"]
        ),
        "SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED": False,
        "FULL_H18_A21_LDLT_CLOSED": False,
        "P3_CANONICAL_PASS": False,
        "no_point_or_reduced_result_may_promote": True,
        "next_obligation": (
            "replace the point source in this same literal event loop by a validated "
            "cover of the complete admitted front-end/tuner/scheduler/PE source family; "
            "only then may the full H18/A21 endpoint LDLT promote P3"
        ),
    }
