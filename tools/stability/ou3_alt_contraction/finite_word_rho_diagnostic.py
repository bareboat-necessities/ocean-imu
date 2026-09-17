#!/usr/bin/env python3
"""Non-promoting complete-word rho feasibility diagnostic for the ALT master.

``AGENTS.md`` requires this measurement before any rigorous enclosure work on a
new complete-word formulation, and it is deliberately the cheapest object that
can tell a false theorem apart from a bad enclosure.  It preserves the exact
matrix/group structure of the shipping word and introduces essentially no
interval pessimism: the shipping covariance recursion runs in binary64 point
arithmetic, and each literal event Jacobian is the same
``physical_lineage``/``joint24_events`` composition the rigorous route uses,
evaluated on point cells.

The controlling object is the complete-word ratio

    rho_w = V(F_w(x)) / V(x)

for the coercive joint24 storage V with the admitted bounded neutral supply on
coordinates ``[e_ba, beta_true] = 18:24``.  ``common_storage_master`` reduces the
existence of a finite supply completion to the projected Finsler restriction

    Z' (A_w' M A_w - rho M) Z < 0,      Z = ker(C) injection.

The floor that holds in both modes comes from an invariant subspace.  If V lies
inside ker(C) and A_w maps V into itself, then for x in V the restriction reads
``(A_w x)' M (A_w x) < rho x' M x`` with A_w x again in V, so

    rho > spectral_radius(A_w restricted to V)^2

for every storage M.  That bound is *metric independent*: no storage, no
multiplier and no sharper enclosure can go below it, and the infimum is not even
attained when the restriction is not diagonalizable.  A floor at or above one
falsifies the formulation itself rather than an enclosure of it.

H18 additionally has an invariant kernel.  Its accelerometer bias is held, so no
literal event writes a motion coordinate from a neutral one, the word map is
block lower triangular,

    A_w = [[A_mm, A_mn], [0, A_nn]],

``ker(C)`` is itself invariant, and ``spectral_radius(A_mm)^2`` is a floor.  This
module checks that entrywise rather than assuming it.  A21 releases the bias and
corrects it from motion coordinates, so its kernel is not invariant and its
motion-block spectrum is reported as an indicative ratio, not a floor.

Nothing here promotes a gate.  The module never searches storage, never fits a
metric to a trajectory, and reports ``storage_search_allowed=False``
unconditionally.  A source phase named below is an analytic source-reachable
member of the admitted language, not a captured replay.
"""
from __future__ import annotations

import copy
import math
from dataclasses import dataclass
from types import SimpleNamespace
from typing import Sequence

import numpy as np

from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_brmm_live_covariance_seed as SEED
import ou3_brmm_shipping_prediction_primitives as SHIP
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_complete_brmm_source_cover_contract as COVER
from tools.stability.ou3_alt_contraction import bias_families as BIAS
from tools.stability.ou3_alt_contraction import common_storage_master as MASTER
from tools.stability.ou3_alt_contraction import physical_lineage as LINEAGE
from tools.stability.ou3_alt_contraction import proof_plan as PLAN

QUALIFICATION = "OU3_ALT_COMPLETE_WORD_RHO_FEASIBILITY_DIAGNOSTIC_V1"
EVIDENCE_KIND = "analytic_source_reachable_word_family"
PHASE = "feasibility_diagnostic"

CANONICAL_WORD_SAMPLES = 600
# The unipotent obstruction is established per sample and composes by induction,
# so the readiness guard probes a short word that still contains both the
# with-S and without-S literal sample shapes.
GUARD_PROBE_SAMPLES = 8
MOTION_DIM = MASTER.MOTION_DIM
JOINT_DIM = MASTER.JOINT_DIM

# Proof-state layout, from ou3_p4_complete_brmm_differential_prediction.  The
# first 18 coordinates are the unsupplied motion subspace in both modes; A21
# adds the active accelerometer bias at 18:21, which the joint24 lift then
# follows with physical beta_true at 21:24.
STATE_NAMES = (
    "theta_x", "theta_y", "theta_z",
    "bg_x", "bg_y", "bg_z",
    "v_x", "v_y", "v_z",
    "p_x", "p_y", "p_z",
    "S_x", "S_y", "S_z",
    "aw_x", "aw_y", "aw_z",
    "ba_x", "ba_y", "ba_z",
    "beta_x", "beta_y", "beta_z",
)
YAW_INDEX = 2
YAW_BIAS_INDEX = 5
OFF_BA = 18

MODES = ("H", "A")


def mode_dimension(mode: str) -> int:
    if mode == "H":
        return 18
    if mode == "A":
        return 21
    raise ValueError("word mode must be H or A")


# Tuner schedule committed at the retained Live entrance, matching the
# ActiveSchedule carried by ou3_brmm_frontend_state_step's point state.
ACTIVE_TAU_S = 1.1
ACTIVE_SIGMA_AW = 0.5
ACTIVE_RS_STD = 2.0
# 25 Hz magnetometer service against the 200 Hz IMU cadence.
GAUGED_MAG_STRIDE = 8
# Committed pseudo-period 0.015 s at the 0.005 s IMU sample.
S_STRIDE = 3
WORLD_MAGNETIC_BODY = (20.0, 0.0, 40.0)
GRAVITY = 9.80665
# Shipping active accelerometer-bias projection radius, from the operating
# domain; the true-bias envelope is the BIAS family component bound.
BIAS_PROJECTION_LIMIT = 0.4
TRUE_BIAS_COMPONENT = 0.35
# Declared physical envelope every composed sample must stay inside.  A word
# outside it is not an admitted history and would prove nothing.
MAX_NON_GRAVITATIONAL_ACCEL = 8.8
MAX_BODY_RATE_RAD_S = 35.0 * math.pi / 180.0


def _I(x: float) -> Interval:
    return Interval.point(float(x))


def _mid(A) -> np.ndarray:
    return np.array([[(x.lo + x.hi) / 2 for x in row] for row in A], dtype=float)


def _point(A) -> list[list[Interval]]:
    return [[_I(float(x)) for x in row] for row in np.asarray(A, dtype=float)]


def _skew(v) -> np.ndarray:
    return np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])


@dataclass(frozen=True)
class SourcePhase:
    """One analytic source-reachable member of the admitted Live language."""

    name: str
    quiet: bool
    gauged: bool
    description: str

    def sample(self, k: int, h: float):
        if self.quiet:
            # Quiet-zero COMPLETE-BRMM member: level boat at rest, zero wave
            # displacement and every centered primitive zero.  Admitted by the
            # corrected physical envelope and by BIAS0.
            omega, force = (0.0, 0.0, 0.0), (0.0, 0.0, -GRAVITY)
        else:
            t = k * h
            omega = (0.05 * math.sin(2 * math.pi * 0.15 * t),
                     -0.04 * math.cos(2 * math.pi * 0.12 * t),
                     0.02 * math.sin(2 * math.pi * 0.08 * t))
            force = (0.5 * math.sin(2 * math.pi * 0.12 * t),
                     0.4 * math.cos(2 * math.pi * 0.10 * t),
                     -GRAVITY + 0.8 * math.sin(2 * math.pi * 0.11 * t))
        assert_admitted_sample(omega, force)
        return omega, force


def assert_admitted_sample(omega, force) -> None:
    """Keep every composed sample inside the declared physical envelope.

    A measurement taken outside the admitted language would prove nothing, so
    this is checked on each sample rather than asserted in a comment.
    """
    rate = math.sqrt(sum(x * x for x in omega))
    if rate > MAX_BODY_RATE_RAD_S:
        raise ValueError(f"body rate {rate} rad/s leaves the declared envelope")
    non_gravitational = math.sqrt(
        force[0] ** 2 + force[1] ** 2 + (force[2] + GRAVITY) ** 2)
    if non_gravitational > MAX_NON_GRAVITATIONAL_ACCEL:
        raise ValueError(
            f"non-gravitational acceleration {non_gravitational} m/s^2 "
            "leaves the declared envelope")


SOURCE_PHASES = (
    SourcePhase(
        "quiet_ungauged", True, False,
        "quiet-zero COMPLETE-BRMM member with no magnetic event in the word; "
        "admitted because the ungauged timeout path leaves mag_ref_set_ false "
        "and MAG-CALL-SCHEDULE-v1 imposes no pre-gauge acquisition deadline"),
    SourcePhase(
        "wave_ungauged", False, False,
        "bounded oscillatory wave response inside the declared envelope with no "
        "magnetic event in the word"),
    SourcePhase(
        "quiet_gauged", True, True,
        "quiet-zero member with 25 Hz magnetic service after north lock"),
    SourcePhase(
        "wave_gauged", False, True,
        "bounded oscillatory wave response with 25 Hz magnetic service"),
)


def live_entry_covariance(mode: str = "H") -> np.ndarray:
    """Live-entry seed, read from the shipping covariance-seed producer."""
    seed = SEED.build()
    failures = SEED.validate(seed)
    if failures:
        raise RuntimeError("live covariance seed prerequisite failed: " + repr(failures))
    attitude = seed["full_heading_gauged_live_attitude_seed"]
    translation = seed["translation_seed"]
    constructor = seed["constructor"]
    tilt = float(attitude["tilt_std_rad"])
    yaw = float(attitude["yaw_std_rad"])
    down = np.array([0.0, 0.0, 1.0])
    n = mode_dimension(mode)
    P = np.zeros((n, n))
    P[0:3, 0:3] = tilt * tilt * (np.eye(3) - np.outer(down, down)) + yaw * yaw * np.outer(down, down)
    P[3:6, 3:6] = float(constructor["P_bg_variance"]) * np.eye(3)
    P[6:9, 6:9] = float(translation["P_v"]) * np.eye(3)
    P[9:12, 9:12] = float(translation["P_p"]) * np.eye(3)
    P[12:15, 12:15] = float(translation["P_S"]) * np.eye(3)
    P[15:18, 15:18] = ACTIVE_SIGMA_AW * ACTIVE_SIGMA_AW * np.eye(3)
    if mode == "A":
        # Released b_a starts at the shipping diagonal seed floor, with the
        # held cross-covariances zeroed exactly as the release does.
        P[OFF_BA:OFF_BA + 3, OFF_BA:OFF_BA + 3] = float(
            seed["held_ba"]["seed_variance"]) * np.eye(3)
    return P


def _prediction_F_Q(omega, constants, mode: str) -> tuple[np.ndarray, np.ndarray]:
    Faa, Qaa = SHIP.attitude_gyro_bias_F_Q(
        tuple(_I(x) for x in omega), constants.h,
        constants.gyro_variance_density_xyz, constants.gyro_bias_variance_density)
    Fll, Qll = SHIP.translation_F_Q(_I(ACTIVE_TAU_S), constants.h, (_I(ACTIVE_SIGMA_AW),) * 3)
    n = mode_dimension(mode)
    F = np.eye(n)
    Q = np.zeros((n, n))
    F[:6, :6] = _mid(Faa)
    Q[:6, :6] = _mid(Qaa)
    F[6:MOTION_DIM, 6:MOTION_DIM] = _mid(Fll)
    Q[6:MOTION_DIM, 6:MOTION_DIM] = _mid(Qll)
    if mode == "A":
        phi_ba, Qba = SHIP.active_accel_bias_F_Q(
            constants.h, constants.accel_bias_tau_s,
            constants.accel_bias_process_variance_density)
        phi = float((phi_ba.lo + phi_ba.hi) / 2)
        F[OFF_BA:OFF_BA + 3, OFF_BA:OFF_BA + 3] = phi * np.eye(3)
        Q[OFF_BA:OFF_BA + 3, OFF_BA:OFF_BA + 3] = _mid(Qba)
    return F, Q


def _H_accelerometer(force, R_wb, mode: str) -> np.ndarray:
    H = np.zeros((3, mode_dimension(mode)))
    H[:, 0:3] = -_skew(force)
    H[:, 15:18] = R_wb
    if mode == "A":
        H[:, OFF_BA:OFF_BA + 3] = np.eye(3)
    return H


def _H_magnetometer(m_body, mode: str) -> np.ndarray:
    H = np.zeros((3, mode_dimension(mode)))
    H[:, 0:3] = -_skew(m_body)
    return H


def _H_S_zero(mode: str) -> np.ndarray:
    H = np.zeros((3, mode_dimension(mode)))
    H[:, 12:15] = np.eye(3)
    return H


def _joseph(P: np.ndarray, H: np.ndarray, R: np.ndarray) -> np.ndarray:
    S = H @ P @ H.T + R
    K = P @ H.T @ np.linalg.inv(S)
    A = np.eye(P.shape[0]) - K @ H
    return A @ P @ A.T + K @ R @ K.T


def _wave_primitive(k: int):
    zero = (_I(0.0), _I(0.0), _I(0.0))
    return KERNEL.WavePrimitivePayload(
        "alt-rho-diagnostic-generator", f"prim{k}", f"prim{k + 1}", "live-origin",
        zero, zero, zero, zero, zero, zero)


def _cell(k: int, ordinal: int, kind: str, P: np.ndarray, force, R_wb,
          m_body, primitive, Racc, Rmag, R_S, h: float, mode: str) -> COVER.SourceCoverCell:
    fields = dict(
        source_token=f"c{k + 1}:e{ordinal}", predecessor_token=f"c{k + 1}", mode=mode,
        sample_index=k, event_ordinal=ordinal, kind=kind,
        state=[_I(0.0)] * mode_dimension(mode), P=_point(P), dt_s=_I(h),
        tau_applied_s=_I(ACTIVE_TAU_S), sigma_aw_mps2=_I(ACTIVE_SIGMA_AW),
        pseudo_elapsed_s=_I(S_STRIDE * h), radial_scale=_I(1.0),
        estimator_source_token=f"c{k + 1}", estimator_predecessor_token="parent",
        estimator_generated_coefficients=True, wave_primitive=primitive)
    if kind == "accelerometer":
        fields.update(f_hat=[_I(x) for x in force], R_hat=_point(R_wb), R=_point(Racc))
    elif kind == "magnetometer":
        fields.update(m_body=[_I(x) for x in m_body], R=_point(Rmag))
    elif kind == "S_zero":
        fields.update(R=_point(R_S), R_provenance=EVENTS.ACTUAL_RS_PROVENANCE)
    else:
        fields.update(R=_point(np.eye(3)))
    if mode == "A":
        fields.update(true_bias=[_I(0.0)] * 3, bias_projection_limit=BIAS_PROJECTION_LIMIT)
    return COVER.SourceCoverCell(**fields)


def compose_word(phase: SourcePhase, samples: int = CANONICAL_WORD_SAMPLES,
                 mode: str = "H") -> dict:
    """Compose one legal single-mode word and return its joint24 map.

    The covariance is the shipping recursion on the same history; each literal
    event Jacobian consumes the covariance that the shipping filter actually
    holds when that event executes.
    """
    if samples <= 0:
        raise ValueError("a word must contain at least one sample")
    mode_dimension(mode)
    constants = KERNEL._process_constants()
    h = float((constants.h.lo + constants.h.hi) / 2)
    Racc = _mid(constants.Racc)
    Rmag = _mid(constants.Rmag)
    R_S = np.diag([ACTIVE_RS_STD * ACTIVE_RS_STD] * 3)
    R_wb = np.eye(3)

    contract = next(c for c in BIAS.contracts() if c.name == "BIAS0")
    held = [Interval(-0.75, 0.75)] * 3
    true_bias = [Interval(-TRUE_BIAS_COMPONENT, TRUE_BIAS_COMPONENT)] * 3

    P = live_entry_covariance(mode)
    Psi = np.eye(JOINT_DIM)
    event_kinds: list[str] = []
    per_event_growth: list[tuple[str, float]] = []

    for k in range(samples):
        omega, force = phase.sample(k, h)
        F, Q = _prediction_F_Q(omega, constants, mode)
        P = F @ P @ F.T + Q

        kinds = ["prediction"]
        if k % S_STRIDE == 0:
            kinds.append("S_zero")
        kinds.append("accelerometer")
        if phase.gauged and k % GAUGED_MAG_STRIDE == 0:
            kinds.append("magnetometer")

        primitive = _wave_primitive(k)
        cells = []
        P_event = P.copy()
        for ordinal, kind in enumerate(kinds):
            cells.append(_cell(k, ordinal, kind, P_event, force, R_wb,
                               WORLD_MAGNETIC_BODY, primitive, Racc, Rmag, R_S, h, mode))
            if kind == "S_zero":
                P_event = _joseph(P_event, _H_S_zero(mode), R_S)
            elif kind == "accelerometer":
                P_event = _joseph(P_event, _H_accelerometer(force, R_wb, mode), Racc)
            elif kind == "magnetometer":
                P_event = _joseph(P_event, _H_magnetometer(WORLD_MAGNETIC_BODY, mode), Rmag)

        selector = SimpleNamespace(
            prefix_length=k + 1,
            parent_source_cell_id="root" if k == 0 else f"c{k}",
            source_cell_id=f"c{k + 1}",
            sample_coordinates=SimpleNamespace(
                omega_body_corrected=tuple(_I(x) for x in omega)))
        lineage = SimpleNamespace(mode=mode, cells=tuple(cells), selector=selector)
        composed = LINEAGE.compose_attached_sample(
            lineage, bias_contract=contract,
            held_bias_error=held if mode == "H" else None,
            true_bias=true_bias if mode == "H" else None,
            tau_ba=constants.accel_bias_tau_s if mode == "A" else None)

        before = float(np.linalg.norm(Psi[:MOTION_DIM, :MOTION_DIM], 2))
        Psi = _mid(composed["J_joint24"]) @ Psi
        after = float(np.linalg.norm(Psi[:MOTION_DIM, :MOTION_DIM], 2))
        event_kinds.extend(kinds)
        per_event_growth.append((",".join(kinds), after / before if before > 0 else math.inf))
        P = P_event

    return {
        "phase": phase.name,
        "mode": mode,
        "samples": samples,
        "word_horizon_s": samples * h,
        "A_joint24": Psi,
        "literal_event_kinds": tuple(sorted(set(event_kinds))),
        "literal_event_count": len(event_kinds),
        "per_sample_growth": per_event_growth,
        "replay_used": False,
    }


def block_triangular_certificate(A: np.ndarray) -> dict:
    """Test whether the neutral rows carry no motion column, exactly.

    When they do not, ``ker(C)`` is itself invariant, the projected Finsler
    restriction equals the motion block, and the floor is
    ``spectral_radius(A_mm)^2``.  H18 satisfies this because its accelerometer
    bias is held.  A21 does not: the released bias is corrected from motion
    coordinates, so the whole motion subspace is not invariant and the
    motion-block spectrum is indicative rather than a floor.  The floor that
    holds in both modes comes from the invariant pair below.
    """
    lower_left = A[MOTION_DIM:, :MOTION_DIM]
    worst = float(np.max(np.abs(lower_left))) if lower_left.size else 0.0
    return {
        "neutral_rows_have_no_motion_column": worst == 0.0,
        "worst_neutral_to_motion_entry": worst,
        "projected_restriction_equals_motion_block": worst == 0.0,
    }


def ungauged_yaw_certificate(A: np.ndarray, horizon_s: float, tolerance: float = 1e-6) -> dict:
    """Exact unipotent structure on the ungauged (yaw, yaw gyro-bias) pair.

    When no magnetic event occurs in the word, the accelerometer residual
    sensitivity ``-skew(f)`` annihilates the specific-force direction, the S=0
    rows touch only S, and the attitude/gyro-bias prediction couples the pair to
    nothing else.  The pair is therefore an invariant subspace carrying
    ``[[1, T], [0, 1]]`` for the elapsed word time ``T``: a unipotent Jordan
    block whose spectral radius is exactly one for every word length.

    Invariance is checked against the whole joint24 map, not only the motion
    block, because that is what the floor argument needs: both coordinates lie
    in ``ker(C)``, and an invariant ``V`` inside ``ker(C)`` forces
    ``rho > spectral_radius(A|_V)^2`` for every storage.  Being checked at full
    width, the certificate holds in H18 and A21 alike.
    """
    idx = [YAW_INDEX, YAW_BIAS_INDEX]
    block = A[np.ix_(idx, idx)]
    rows = A[idx, :].copy()
    rows[:, idx] = 0.0
    cols = A[:, idx].copy()
    cols[idx, :] = 0.0
    off_row = float(np.max(np.abs(rows)))
    off_col = float(np.max(np.abs(cols)))
    invariant = off_row == 0.0 and off_col == 0.0
    unit_diagonal = bool(abs(block[0, 0] - 1.0) <= tolerance and abs(block[1, 1] - 1.0) <= tolerance)
    strictly_lower_zero = bool(abs(block[1, 0]) <= tolerance)
    drift = float(block[0, 1])
    return {
        "coordinates": (STATE_NAMES[YAW_INDEX], STATE_NAMES[YAW_BIAS_INDEX]),
        "subspace_is_exactly_invariant": invariant,
        "worst_row_leak_outside_block": off_row,
        "worst_column_leak_outside_block": off_col,
        "block": [[float(x) for x in row] for row in block],
        "unit_diagonal_within_tolerance": unit_diagonal,
        "strictly_lower_entry_zero_within_tolerance": strictly_lower_zero,
        "observed_bias_to_yaw_drift": drift,
        "elapsed_word_horizon_s": horizon_s,
        "drift_matches_elapsed_horizon": bool(abs(drift - horizon_s) <= 1e-3 * max(1.0, horizon_s)),
        "unipotent_jordan_block_certified": bool(invariant and unit_diagonal and strictly_lower_zero),
        "certified_spectral_radius": 1.0 if (invariant and unit_diagonal and strictly_lower_zero) else None,
        "note": (
            "The observed diagonal deviation is accumulated binary64 rounding of the "
            "composed Cayley products; the linearized recurrence has diagonal exactly 1."),
    }


def per_sample_unipotent_induction(samples: int = 4, mode: str = "H") -> dict:
    """Show the obstruction is per sample, hence independent of word length.

    Each ungauged literal sample leaves the (yaw, yaw gyro-bias) pair invariant
    and contributes ``[[1, h], [0, 1]]``.  Products of such blocks are again
    unipotent, so composing any number of samples gives ``[[1, N h], [0, 1]]``
    and never leaves spectral radius one.  Both literal sample shapes -- with
    and without the scheduled S=0 event -- are exercised.
    """
    phase = SOURCE_PHASES[0]
    rows = []
    shapes: set[str] = set()
    for k in range(samples):
        word = compose_word(phase, k + 1, mode)
        certificate = ungauged_yaw_certificate(word["A_joint24"], word["word_horizon_s"])
        sample_shapes = {kinds for kinds, _ in word["per_sample_growth"]}
        shapes |= sample_shapes
        rows.append({
            "samples": k + 1,
            "literal_sample_shapes": sorted(sample_shapes),
            "certified": certificate["unipotent_jordan_block_certified"],
            "bias_to_yaw_drift": certificate["observed_bias_to_yaw_drift"],
            "elapsed_word_horizon_s": certificate["elapsed_word_horizon_s"],
            "drift_matches_elapsed_horizon": certificate["drift_matches_elapsed_horizon"],
        })
    with_S = any("S_zero" in shape for shape in shapes)
    without_S = any("S_zero" not in shape for shape in shapes)
    return {
        "prefixes": rows,
        "literal_sample_shapes_exercised": sorted(shapes),
        "both_literal_sample_shapes_exercised": bool(with_S and without_S),
        "every_prefix_unipotent": all(row["certified"] for row in rows),
        "drift_grows_with_elapsed_horizon": all(
            row["drift_matches_elapsed_horizon"] for row in rows),
        "length_independent": bool(
            all(row["certified"] for row in rows) and with_S and without_S),
    }


def word_rho(A: np.ndarray, kernel_is_invariant: bool) -> dict:
    """Motion-block spectrum, and the rho floor it supports when it is one.

    ``kernel_is_invariant`` is the block-triangular certificate.  When it holds,
    ``ker(C)`` is invariant and ``spectral_radius(A_mm)^2`` is a metric
    independent floor.  When it does not, the same number is only indicative:
    the projected restriction then also carries the neutral rows of ``A Z``.
    """
    motion = A[:MOTION_DIM, :MOTION_DIM]
    values, vectors = np.linalg.eig(motion)
    order = np.argsort(-np.abs(values))
    dominant = int(order[0])
    radius = float(abs(values[dominant]))
    direction = vectors[:, dominant]
    direction = direction / np.linalg.norm(direction)
    weights = np.abs(direction)
    ranked = np.argsort(-weights)
    return {
        "spectral_radius_motion_block": radius,
        "motion_block_ratio": radius * radius,
        "motion_block_ratio_is_a_floor": bool(kernel_is_invariant),
        "rho_floor": radius * radius if kernel_is_invariant else None,
        "distance_to_rho_one": 1.0 - radius * radius,
        "dominant_eigenvalue": [float(values[dominant].real), float(values[dominant].imag)],
        "maximizing_state_direction": [
            {"coordinate": STATE_NAMES[i], "weight": float(weights[i])}
            for i in ranked[:6]],
        "leading_motion_spectrum": [float(abs(values[i])) for i in order[:6]],
        "motion_block_two_norm": float(np.linalg.norm(motion, 2)),
    }


def margin_consumption(word: dict, top: int = 8) -> list[dict]:
    """Operation-by-operation margin consumption over the composed word."""
    rows = [
        {"sample": i, "literal_events": kinds, "motion_two_norm_ratio": float(ratio)}
        for i, (kinds, ratio) in enumerate(word["per_sample_growth"])]
    rows.sort(key=lambda r: -r["motion_two_norm_ratio"])
    return rows[:top]


def diagnose(samples: int = CANONICAL_WORD_SAMPLES,
             phases: Sequence[SourcePhase] = SOURCE_PHASES,
             modes: Sequence[str] = MODES) -> dict:
    PLAN.require_theorem_task(
        obligation="ALT complete-word rho feasibility before rigorous enclosure",
        evidence_kind=EVIDENCE_KIND, complete_physical_word=False,
        requested_phase=PHASE)

    results = {}
    for mode in modes:
        for phase in phases:
            word = compose_word(phase, samples, mode)
            A = word["A_joint24"]
            triangular = block_triangular_certificate(A)
            rho = word_rho(A, triangular["projected_restriction_equals_motion_block"])
            entry = {
                "mode": mode,
                "source_phase": phase.name,
                "description": phase.description,
                "gauged": phase.gauged,
                "word_horizon_s": word["word_horizon_s"],
                "literal_event_kinds": list(word["literal_event_kinds"]),
                "literal_event_count": word["literal_event_count"],
                "block_triangular_certificate": triangular,
                **rho,
                "worst_margin_consumption": margin_consumption(word),
            }
            if not phase.gauged:
                entry["ungauged_yaw_certificate"] = ungauged_yaw_certificate(
                    A, word["word_horizon_s"])
            results[f"{mode}:{phase.name}"] = entry

    ratios = {name: row["motion_block_ratio"] for name, row in results.items()}
    certified = [name for name, row in results.items()
                 if row.get("ungauged_yaw_certificate", {}).get("unipotent_jordan_block_certified")]
    infeasible = bool(certified)

    # The composed word map is a product of hundreds of binary64 Cayley/Joseph
    # operations, so a measured eigenvalue carries accumulated rounding.  The
    # deviation of the certified unipotent diagonal from its exact value 1 is a
    # direct proxy for that drift, and every measured ratio is only meaningful
    # outside it.  The falsification below therefore rests on the exact
    # invariant-subspace algebra, never on a measured radius near 1.  Only a
    # certified-invariant block has a known exact diagonal, so only those words
    # measure drift; a word whose block leaks has a genuine deviation.
    drift = max((abs(results[name]["ungauged_yaw_certificate"]["block"][0][0] - 1.0)
                 for name in certified), default=0.0)
    limiting = certified[0] if certified else max(ratios, key=ratios.get)
    worst_by_mode = {
        mode: max((ratios[name] for name in ratios if name.startswith(mode + ":")),
                  default=None)
        for mode in modes}
    report = {
        "qualification": QUALIFICATION,
        "evidence_kind": EVIDENCE_KIND,
        "phase": PHASE,
        "samples_per_word": samples,
        "modes": list(modes),
        "words": results,
        "measured_motion_block_ratio_by_word": ratios,
        "worst_motion_block_ratio_by_mode": worst_by_mode,
        "numerical_drift_proxy": drift,
        "measured_ratios_are_meaningful_only_outside_drift": True,
        "rho_floor_over_legal_words": 1.0 if infeasible else max(
            (row["rho_floor"] for row in results.values() if row["rho_floor"] is not None),
            default=None),
        "rho_floor_source": ("exact_unipotent_subspace_algebra" if infeasible
                             else "measured_spectral_radius"),
        "limiting_word": limiting,
        "limiting_state_direction": (
            [{"coordinate": STATE_NAMES[YAW_INDEX], "weight": 1.0},
             {"coordinate": STATE_NAMES[YAW_BIAS_INDEX], "weight": 1.0}]
            if infeasible else results[limiting]["maximizing_state_direction"]),
        "distance_to_rho_one": 0.0 if infeasible else results[limiting]["distance_to_rho_one"],
        "ungauged_unipotent_words": certified,
        "per_sample_unipotent_induction": {
            mode: per_sample_unipotent_induction(mode=mode) for mode in modes
        } if infeasible else None,
        "declared_joint24_contraction_falsified": infeasible,
        "failure_classification": "theorem_failure" if infeasible else "open",
        "interpretation": (
            "An ungauged legal word carries an exactly invariant unipotent "
            "(yaw, yaw gyro-bias) pair inside ker(C), whose linearized diagonal is "
            "exactly one, so no common joint24 storage with the declared neutral "
            "supply attains rho<1 in either mode. The obstruction is the "
            "formulation, not an enclosure."
            if infeasible else
            "No ungauged unipotent pair was certified in the words examined; any "
            "reported floor is a lower bound only and does not authorize storage work."),
        "metric_independent": True,
        "replay_or_captured_trace_used": False,
        "storage_search_allowed": False,
        "ALT_STARTUP_PASS": False,
        "ALT_LIVE_PASS": False,
        "ALT_END_TO_END_PASS": False,
    }
    PLAN.assert_non_promoting_report(report)
    return report


_GUARD_PROBE_CACHE: list[dict] = []


def guard_probe() -> dict:
    """Short deterministic probe for readiness guards that run it repeatedly.

    The obstruction is per sample and composes by induction, so a short word
    settles it; the canonical 3 s measurement is the standalone run.  The result
    is memoised because composing it is not cheap and it takes no input.
    """
    if not _GUARD_PROBE_CACHE:
        _GUARD_PROBE_CACHE.append(diagnose(samples=GUARD_PROBE_SAMPLES,
                                           phases=(SOURCE_PHASES[0],), modes=("H",)))
    return copy.deepcopy(_GUARD_PROBE_CACHE[0])


def build(samples: int = CANONICAL_WORD_SAMPLES) -> dict:
    return diagnose(samples)


def validate(d: dict) -> list[str]:
    f: list[str] = []
    if d.get("qualification") != QUALIFICATION:
        f.append("qualification mismatch")
    if d.get("metric_independent") is not True:
        f.append("rho floor was not reported as metric independent")
    if d.get("replay_or_captured_trace_used") is not False:
        f.append("diagnostic admitted replay or captured-trace evidence")
    words = d.get("words", {})
    if not words:
        f.append("no legal word was composed")
    for name, row in words.items():
        if not row.get("literal_event_kinds"):
            f.append(f"{name}: word composed no literal event")
        if row.get("gauged") is False and "magnetometer" in row.get("literal_event_kinds", []):
            f.append(f"{name}: ungauged word contains a magnetic event")
        if row.get("gauged") is True and "magnetometer" not in row.get("literal_event_kinds", []):
            f.append(f"{name}: gauged word contains no magnetic event")
        # H18 holds its accelerometer bias, so ker(C) must be invariant there;
        # A21 releases it, so the motion block is indicative, not a floor.
        invariant = row.get("block_triangular_certificate", {}).get(
            "projected_restriction_equals_motion_block")
        if row.get("mode") == "H" and invariant is not True:
            f.append(f"{name}: H18 kernel is no longer invariant")
        if row.get("motion_block_ratio_is_a_floor") != bool(invariant):
            f.append(f"{name}: motion-block ratio claimed as a floor without an invariant kernel")
        if row.get("rho_floor") is not None and invariant is not True:
            f.append(f"{name}: reported a floor without an invariant kernel")
    if d.get("declared_joint24_contraction_falsified"):
        induction = d.get("per_sample_unipotent_induction") or {}
        if not induction:
            f.append("falsification was not shown independent of word length")
        for mode, row in induction.items():
            if row.get("length_independent") is not True:
                f.append(f"{mode}: falsification was not shown independent of word length")
        if d.get("rho_floor_source") != "exact_unipotent_subspace_algebra":
            f.append("falsification rests on a measured radius rather than exact algebra")
        if d.get("rho_floor_over_legal_words") != 1.0:
            f.append("falsified formulation did not report a floor of one")
    for k in ("storage_search_allowed", "ALT_STARTUP_PASS", "ALT_LIVE_PASS", "ALT_END_TO_END_PASS"):
        if d.get(k) is not False:
            f.append(k + " not false")
    return f


def main() -> int:
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=CANONICAL_WORD_SAMPLES,
                        help="IMU samples per composed word (600 = the canonical 3 s word)")
    args = parser.parse_args()
    report = build(args.samples)
    failures = validate(report)
    report["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    print(json.dumps({
        "rho_floor_over_legal_words": report["rho_floor_over_legal_words"],
        "rho_floor_source": report["rho_floor_source"],
        "limiting_word": report["limiting_word"],
        "distance_to_rho_one": report["distance_to_rho_one"],
        "measured_motion_block_ratio_by_word": report["measured_motion_block_ratio_by_word"],
        "worst_motion_block_ratio_by_mode": report["worst_motion_block_ratio_by_mode"],
        "numerical_drift_proxy": report["numerical_drift_proxy"],
        "ungauged_unipotent_words": report["ungauged_unipotent_words"],
        "declared_joint24_contraction_falsified": report["declared_joint24_contraction_falsified"],
        "failure_classification": report["failure_classification"],
        "storage_search_allowed": report["storage_search_allowed"],
        "validation_failures": failures,
    }, sort_keys=True))
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
