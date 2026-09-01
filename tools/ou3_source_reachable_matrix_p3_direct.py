#!/usr/bin/env python3
"""Direct generalized-matrix P3 comparison for deployed OU-III.

This layer removes the last scalarization in the source-reachable P3 proof.
The conditioned backend already certifies the integrated-OU process family in
an exact rational C=R L^{-1} congruence basis.  Here the same congruence is
applied to both the word-noise lower matrix and covariance upper matrix and we
certify directly

    Omega_word - delta * Sigma_upper  >>  0

with outward-rounded interval LDL^T.

No theorem domain, filter gain, source schedule, or usefulness gate is changed.
The translation process matrix is retained as a 4x4 interval family for each
axis.  A scalar measurement-information norm is used only to obtain a uniform
posterior *matrix* lower bound; the process matrix is not collapsed to its
smallest eigenvalue.  The already-established directional covariance upper
bounds are retained as a diagonal Loewner upper matrix rather than collapsed to
one global maximum.  H and A are finally checked as assembled 18x18 and 21x21
matrices in the same conditioned coordinates.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import functools

from ou3_interval import (
    Interval,
    symmetric_gershgorin_upper,
    symmetric_positive_definite_ldlt,
)
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan
import ou3_p3_word_noise_accumulation as WORDNOISE
import ou3_source_reachable_matrix_p3 as BASE
import ou3_source_reachable_matrix_p3_factored as FACTORED

DEFAULT_DOMAIN = FACTORED.DEFAULT_DOMAIN
SCHEMA = FACTORED.SCHEMA
MIN_USEFUL_DELTA = FACTORED.MIN_USEFUL_DELTA
# The factored backend now carries the true exp(lambda*x.hi) tail factor, so the
# exact series stays valid past the per-step domain.  Its order-22 truncation
# still loses the enclosure eventually; the cap keeps the search inside the
# region where the remainder is small, and any horizon whose enclosure is too
# wide fails its own LDLT and is skipped anyway.
WORD_SERIES_X_LIMIT = 2.5
# X sub-cells per source cell.  The congruence cancels three decades, so the
# enclosure certifies only on a narrow cell; 128 is the smallest power of two
# that carries the full word horizon on the binding source cell.  Shorter
# horizons tolerate a proportionally wider cell and use proportionally fewer.
WORD_X_SUBCELLS = 128
WORD_SUBCELL_REFERENCE_X = 0.25

_ORIGINAL_MODE_CELL = BASE.mode_cell


def _zero_matrix(n: int):
    z = Interval.point(0.0)
    return [[z for _ in range(n)] for _ in range(n)]


def _diag_matrix(values: list[float]):
    A = _zero_matrix(len(values))
    for i, v in enumerate(values):
        A[i][i] = Interval.outward_bounds(v, v)
    return A


def _scale_matrix(A, s: float):
    q = Interval.outward_bounds(s, s)
    return [[q * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def _subtract_delta(Omega, Sigma, delta: float):
    q = Interval.outward_bounds(delta, delta)
    return BASE.matrix_symmetric_hull([
        [Omega[i][j] - q * Sigma[i][j] for j in range(len(Omega))]
        for i in range(len(Omega))
    ])


def _spd_at_delta(Omega, Sigma, delta: float) -> bool:
    ok, _ = symmetric_positive_definite_ldlt(_subtract_delta(Omega, Sigma, delta))
    return ok


@functools.lru_cache(maxsize=200000)
def _word_translation_information(block):
    """``C' (C Q_scaled(X) C')^-1 C`` -- the word-noise information in the
    word-scaled coordinates, read through the basis where the inverse is tight.

    It depends only on the X sub-cell, so every ``sigma``, ``R_S`` and mode cell
    that reaches the same horizon shares it."""
    inv = _whitened_word_shape_inverse(block)
    if inv is None:
        return None
    C = FACTORED._C_INTERVAL
    A = BASE.matrix_symmetric_hull(
        BASE.matrix_mul(BASE.matrix_mul(BASE.matrix_transpose(C), [list(r) for r in inv]), C)
    )
    return tuple(tuple(A[i][j] for j in range(4)) for i in range(4))


def _word_translation_block_margin(information, sigma_root, information_diag) -> float:
    """``delta`` for ``(Omega_word^-1 + D)^-1 >= delta Sigma`` on one X sub-cell.

    Forming the posterior noise floor explicitly breaks down once the word
    assimilates real measurement information: ``Omega (I + D Omega)^-1`` cancels
    the measured directions almost exactly and no interval enclosure of it stays
    definite past a sixty-four step horizon.  The information form
    ``Omega^-1 + D <= (1/delta) Sigma^-1`` avoids that, but carrying it in the
    rational congruence basis fails for the opposite reason: ``Sigma^-1`` there
    is assembled from a diagonal spanning thirteen decades and is not
    certifiably definite at all.

    Normalising ``Sigma`` to the identity removes both problems.  With
    ``G = Sigma^{1/2}`` diagonal, the statement is exactly

        lambda_max( G (Omega^-1 + D) G )  <=  1/delta,

    an upper bound on one symmetric matrix -- no inverse of an ill-conditioned
    matrix, no definiteness test on a near-singular one.  ``Omega^-1`` is read in
    the congruence basis, where it is well conditioned, and transported back by
    the exact rational ``C``.
    """
    g = [Interval.outward_bounds(v, v) for v in sigma_root]
    M = [[g[i] * information[i][j] * g[j] for j in range(4)] for i in range(4)]
    for i in range(4):
        M[i][i] = M[i][i] + Interval.outward_bounds(information_diag[i], information_diag[i])
    top = symmetric_gershgorin_upper(BASE.matrix_symmetric_hull(M))
    if not (math.isfinite(top) and top > 0.0):
        return 0.0
    return BASE.down(1.0 / top)


def _certified_generalized_delta(Omega, Sigma, gate: float) -> float:
    """Largest located certified delta lower bound, using log-domain bracketing."""
    gate = float(gate)
    if gate <= 0.0:
        raise ValueError("positive usefulness gate required")

    if _spd_at_delta(Omega, Sigma, gate):
        lo = gate
        hi = gate
        while hi < 1.0:
            trial = min(1.0, hi * 10.0)
            if not _spd_at_delta(Omega, Sigma, trial):
                hi = trial
                break
            lo = trial
            hi = trial
            if trial == 1.0:
                return BASE.down(lo)
        if hi == lo:
            return BASE.down(lo)
    else:
        hi = gate
        lo = gate
        for _ in range(80):
            lo /= 10.0
            if _spd_at_delta(Omega, Sigma, lo):
                break
        else:
            return 0.0

    for _ in range(44):
        mid = math.sqrt(lo * hi)
        if _spd_at_delta(Omega, Sigma, mid):
            lo = mid
        else:
            hi = mid
    return BASE.down(lo)


def _matrix_abs_row_sum_upper(A) -> float:
    best = 0.0
    for row in A:
        total = 0.0
        for a in row:
            total = BASE.up(total + a.abs_upper())
        best = max(best, total)
    return BASE.up(best)


def _measurement_beta_upper(mode: str, sigma: Interval, rs: Interval,
                            live: dict, vector: dict, process: dict,
                            sched: dict) -> float:
    h = sched["dt_s"]
    vc = vector["configured_measurement_bounds"]
    ra = BASE.down(BASE.pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rm = BASE.down(BASE.pos(vc["mag_measurement_std_uT"], "mag std") ** 2)
    fhi = BASE.pos(live["specific_force_norm_upper_mps2"], "force upper")
    mhi = BASE.pos(live["magnetic_vector_norm_upper_uT"], "mag upper")
    qtheta = BASE.pos(process["attitude_gyro_bias"]["theta_diagonal_lower"], "theta process")
    qba_d = BASE.pos(
        process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"],
        "active bias discrete process",
    )
    sS2 = (sigma.lo * h * h * h) ** 2
    sa2 = sigma.lo * sigma.lo
    betaS = BASE.up(sS2 / BASE.rs_variance_lower(rs, sched))
    betaAcc = BASE.up(
        (fhi * fhi * qtheta + sa2 + (qba_d if mode == "A" else 0.0)) / ra
    )
    betaMag = BASE.up((mhi * mhi * qtheta) / rm)
    return BASE.up(BASE.up(betaS + betaAcc) + betaMag)


def _translation_direct_blocks(x: Interval, sigma: Interval, raw: dict,
                               beta: float, sched: dict):
    """Return conditioned (Omega,Sigma), process factor and norm diagnostics.

    The process family is represented as Q_scaled = x * Shape(x).  We keep the
    already-conditioned Shape(x)=Q_scaled/x interval family, transform it by C,
    and use x.lo only after the congruence.  For every concrete source value in
    the cell, x*Shape(x) >= x.lo*Shape(x) in Loewner order because Shape(x) is
    PSD.  This avoids reintroducing the small-x interval dependency eliminated
    by the exact RL^{-1} process certificate.
    """
    shape = FACTORED.process_shape_q(x)
    shape_norm = _matrix_abs_row_sum_upper(shape)
    qnorm = BASE.up(x.hi * shape_norm)
    factor = BASE.down(1.0 / BASE.up(1.0 + BASE.up(beta * qnorm)))

    C = FACTORED._C_INTERVAL
    Ct = BASE.matrix_transpose(C)
    shape_tilde = BASE.matrix_symmetric_hull(
        BASE.matrix_mul(BASE.matrix_mul(C, shape), Ct)
    )
    omega_scale = BASE.down(factor * x.lo)
    Omega_tilde = _scale_matrix(shape_tilde, omega_scale)

    h = sched["dt_s"]
    scales = [sigma.lo * h, sigma.lo * h * h, sigma.lo * h * h * h, sigma.lo]
    upper = raw["Sigma_diagonal_upper"]
    physical = [upper[6], upper[9], upper[12], upper[15]]
    Sigma = _diag_matrix([
        BASE.up(physical[i] / (scales[i] * scales[i])) for i in range(4)
    ])
    Sigma_tilde = BASE.matrix_symmetric_hull(
        BASE.matrix_mul(BASE.matrix_mul(C, Sigma), Ct)
    )
    return Omega_tilde, Sigma_tilde, factor, qnorm



# --- word-accumulated injected-noise floor -----------------------------------
#
# The blocks above compare the covariance at the *word* endpoint against the
# noise of a single IMU step.  That is a valid lower comparison but an extremely
# weak one: on the S channel alone the two differ by (T/h)^7.  The routines
# below build the accumulated floor the word actually injects, using only the
# shipping recursion, and the cell keeps whichever of the two certifies more.


@functools.lru_cache(maxsize=200000)
def _whitened_word_shape(lo: float, hi: float):
    """``C Q_scaled(X) C'`` on one X sub-cell, or ``None`` if it is not certifiable."""
    X = Interval.outward_bounds(lo, hi)
    if not (X.hi < WORD_SERIES_X_LIMIT):
        return None
    shape = FACTORED.process_shape_q(X)
    Q = [[Interval.outward_bounds(X.lo, X.lo) * shape[i][j] for j in range(4)]
         for i in range(4)]
    C = FACTORED._C_INTERVAL
    Oc = BASE.matrix_symmetric_hull(
        BASE.matrix_mul(BASE.matrix_mul(C, BASE.matrix_symmetric_hull(Q)),
                        BASE.matrix_transpose(C))
    )
    if not symmetric_positive_definite_ldlt(Oc)[0]:
        return None
    return tuple(tuple(Oc[i][j] for j in range(4)) for i in range(4))


@functools.lru_cache(maxsize=200000)
def _whitened_word_shape_inverse(block):
    """``(C Q_scaled(X) C')^-1``.  The congruence leaves a condition number near
    two, so this inverse is tight where the raw process matrix would not be."""
    if block is None:
        return None
    try:
        inv = matrix_inverse_gauss_jordan([list(row) for row in block])
    except Exception:
        return None
    return tuple(tuple(inv[i][j] for j in range(4)) for i in range(4))


def _x_subcells(lo: float, hi: float, count: int):
    if not (0.0 < lo <= hi):
        raise ValueError("positive x cell required")
    if hi == lo:
        return [(lo, hi)]
    ratio = (hi / lo) ** (1.0 / count)
    edges = [lo]
    for _ in range(count - 1):
        edges.append(edges[-1] * ratio)
    edges.append(hi)
    return [(BASE.down(edges[i]), BASE.up(edges[i + 1])) for i in range(count)]


def _translation_word_margin(x: Interval, sigma: Interval, rs: Interval, raw: dict,
                             live: dict, vector: dict, process: dict, sched: dict,
                             mode: str):
    """Certified translation margin against the word-accumulated noise floor.

    Over ``N`` pure prediction steps the shipping recursion injects exactly
    ``sum_k Phi^k Q Phi^k' = Q(N h)`` -- the *same* analytic family, evaluated at
    the word horizon.  Scaled by ``diag(sigma T, sigma T^2, sigma T^3, sigma)``
    the integrated-OU chain depends only on ``X = T/tau``, so the already
    certified shape family supplies the accumulated floor with no recursion and
    no interval amplification.

    ``C Q_scaled(X) C'`` stays well conditioned for every ``X`` the word reaches,
    but the congruence cancels three decades, so the enclosure is certifiable
    only on a narrow ``X`` cell.  The source cell is therefore refined and the
    margin is the minimum over its sub-cells.

    Every admissible horizon is scanned and the best is kept, because the trade
    is not monotone: a longer horizon injects more noise but also assimilates
    more measurement information, and on a source cell with a small ``R_S`` the
    second effect wins well before the full word.  ``N = 1`` reproduces the
    single-step comparison, so the scan can only improve on it.
    """
    h = sched["dt_s"]
    steps_cap = math.floor(BASE.down(float(raw["word_horizon_s_lower"]) / BASE.up(h)))
    if steps_cap < 1:
        raise RuntimeError("source word does not certainly contain one prediction step")
    doublings = int(math.floor(math.log2(steps_cap)))

    vc = vector["configured_measurement_bounds"]
    ra = BASE.down(BASE.pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rm = BASE.down(BASE.pos(vc["mag_measurement_std_uT"], "mag std") ** 2)
    fhi = BASE.pos(live["specific_force_norm_upper_mps2"], "force upper")
    mhi = BASE.pos(live["magnetic_vector_norm_upper_uT"], "mag upper")
    qtheta = BASE.pos(process["attitude_gyro_bias"]["theta_diagonal_lower"], "theta process")
    qba_d = BASE.pos(
        process["active_accelerometer_bias"]["Q_accel_bias_lambda_min_lower"],
        "active bias discrete process",
    )
    cadence_lo = float(raw["cadence_s"][0])
    upper = raw["Sigma_diagonal_upper"]
    physical = [upper[6], upper[9], upper[12], upper[15]]
    C = FACTORED._C_INTERVAL
    Ct = BASE.matrix_transpose(C)

    best = None
    for k in range(doublings, -1, -1):
        steps = 2 ** k
        horizon = BASE.down(steps * h)
        # The congruence certifies only on a narrow X cell and the width it
        # tolerates scales with X, so short horizons need far fewer sub-cells.
        span = BASE.up(steps * x.hi)
        count = max(8, min(WORD_X_SUBCELLS,
                           int(math.ceil(WORD_X_SUBCELLS * span / WORD_SUBCELL_REFERENCE_X))))
        cells = _x_subcells(BASE.down(steps * x.lo), span, count)
        # Probe the ends before paying for the whole refinement: when a horizon
        # is too long for the congruence it fails there first.
        if any(_whitened_word_shape(*cells[i]) is None for i in (0, -1)):
            continue
        blocks = [_whitened_word_shape(lo, hi) for lo, hi in cells]
        if any(b is None for b in blocks):
            continue

        word_scales = [
            BASE.up(sigma.lo * horizon),
            BASE.up(sigma.lo * horizon * horizon),
            BASE.up(sigma.lo * horizon * horizon * horizon),
            BASE.up(sigma.lo),
        ]
        # Sigma^{1/2} for the word-scaled translation coordinates.  Upper
        # rounding here only lowers the reported margin.
        sigma_root = [
            BASE.up(math.sqrt(BASE.up(physical[i])) / BASE.down(word_scales[i]))
            for i in range(4)
        ]

        # Measurement information the word can assimilate, per word-scaled
        # coordinate.  No shipping update has a v or p column: the pseudo-update
        # measures S, the accelerometer measures a_w with its attitude and bias
        # cross terms, the magnetometer measures attitude only.
        s_firings = BASE.up(BASE.up(horizon / BASE.down(cadence_lo)) + 1.0)
        info_S = BASE.up(
            s_firings * BASE.up(word_scales[2] * word_scales[2] / BASE.rs_variance_lower(rs, sched))
        )
        per_sample_aw = BASE.up(
            (fhi * fhi * qtheta + sigma.lo * sigma.lo + (qba_d if mode == "A" else 0.0)) / ra
        )
        per_sample_aw = BASE.up(per_sample_aw + BASE.up(mhi * mhi * qtheta / rm))
        info_aw = BASE.up(float(steps) * per_sample_aw)
        # G D G is diagonal because D is.
        scaled_information = [
            0.0,
            0.0,
            BASE.up(info_S * BASE.up(sigma_root[2] * sigma_root[2])),
            BASE.up(info_aw * BASE.up(sigma_root[3] * sigma_root[3])),
        ]

        margin = math.inf
        worst_index = 0
        for index, block in enumerate(blocks):
            information = _word_translation_information(block)
            if information is None:
                margin = 0.0
                break
            value = _word_translation_block_margin(
                information, sigma_root, scaled_information
            )
            if value <= 0.0:
                margin = 0.0
                break
            if value < margin:
                margin, worst_index = value, index
        if not margin > 0.0:
            continue
        candidate = (margin, blocks[worst_index], sigma_root, {
            "accumulated_prediction_steps": steps,
            "accumulated_horizon_s": horizon,
            "word_horizon_s_lower": float(raw["word_horizon_s_lower"]),
            "x_subcells": len(cells),
            "S_zero_firings_upper": s_firings,
            "S_channel_information_upper": info_S,
            "a_w_channel_information_upper": info_aw,
        })
        if best is None or candidate[0] > best[0]:
            best = candidate
        elif candidate[0] < best[0]:
            # The margin against horizon is unimodal: a longer horizon injects
            # more noise but also assimilates more measurement information.
            # Once a shorter horizon reports less, the peak is behind us.
            break
    if best is None:
        raise RuntimeError("no admissible word horizon certified the translation floor")
    return best


def _attitude_bias_word_blocks(raw: dict, live: dict, vector: dict, process: dict,
                               sched: dict):
    """Return (Omega_posterior, Sigma) for one (theta, gyro-bias) axis pair."""
    h = sched["dt_s"]
    doublings = WORDNOISE.word_step_doublings(float(raw["word_horizon_s_lower"]), h)
    steps = 2 ** doublings

    ab = process["attitude_gyro_bias"]
    qtheta = BASE.pos(ab["theta_diagonal_lower"], "theta process")
    qbias = BASE.pos(ab["gyro_bias_diagonal_lower"], "gyro bias process")
    cross = float(ab["cross_norm_upper"])
    rho = BASE.down(1.0 - BASE.up(cross / BASE.down(math.sqrt(qtheta * qbias))))
    if rho <= 0.0:
        raise RuntimeError("scaled attitude/bias process comparison lost positivity")

    coupling = BASE.up(h * BASE.up(math.sqrt(BASE.up(qbias / BASE.down(qtheta)))))
    vc = vector["configured_measurement_bounds"]
    ra = BASE.down(BASE.pos(vc["acc_measurement_std_mps2"], "acc std") ** 2)
    rm = BASE.down(BASE.pos(vc["mag_measurement_std_uT"], "mag std") ** 2)
    fhi = BASE.pos(live["specific_force_norm_upper_mps2"], "force upper")
    mhi = BASE.pos(live["magnetic_vector_norm_upper_uT"], "mag upper")
    # Attitude information per sample in the sqrt(q_theta) scaling.  The a_w and
    # S informations do not belong here: no shipping update measures the gyro
    # bias directly, and charging the accelerometer's a_w figure against these
    # coordinates is what collapsed the binding P3 channel.
    per_sample = BASE.up(BASE.up(fhi * fhi / ra + mhi * mhi / rm) * qtheta)
    Omega = WORDNOISE.attitude_bias_word_noise(
        rho, coupling, doublings, BASE.up(float(steps) * per_sample)
    )
    upper = raw["Sigma_diagonal_upper"]
    Sigma = _diag_matrix([
        BASE.up(upper[0] / qtheta),
        BASE.up(upper[3] / qbias),
    ])
    return Omega, Sigma, {
        "attitude_accumulated_prediction_steps": steps,
        "scaled_attitude_bias_process_floor_lower": rho,
        "scaled_bias_coupling_per_step": coupling,
        "attitude_information_per_sample_upper": per_sample,
    }

def _assemble_mode_matrices(mode: str, raw: dict, trans_Omega, trans_Sigma,
                            att_Omega=None, att_Sigma=None):
    n = 18 if mode == "H" else 21
    O = _zero_matrix(n)
    S = _zero_matrix(n)
    qpost = float(raw["post_measurement_scaled_Omega_lambda_min_lower"])
    scales2 = raw["comparison_scale_diagonal_squared"]
    upper = raw["Sigma_diagonal_upper"]

    placed = set()
    for axis in range(3):
        idx = [6 + axis, 9 + axis, 12 + axis, 15 + axis]
        placed.update(idx)
        for i in range(4):
            for j in range(4):
                O[idx[i]][idx[j]] = trans_Omega[i][j]
                S[idx[i]][idx[j]] = trans_Sigma[i][j]

    if att_Omega is not None:
        for axis in range(3):
            idx = [axis, 3 + axis]
            placed.update(idx)
            for i in range(2):
                for j in range(2):
                    O[idx[i]][idx[j]] = att_Omega[i][j]
                    S[idx[i]][idx[j]] = att_Sigma[i][j]

    for i in range(n):
        if i in placed:
            continue
        O[i][i] = Interval.outward_bounds(qpost, qpost)
        sui = BASE.up(upper[i] / scales2[i])
        S[i][i] = Interval.outward_bounds(sui, sui)
    return O, S


def mode_cell(mode: str, x: Interval, rho_trans: float, sigma: Interval,
              rs: Interval, live: dict, vector: dict, process: dict,
              sched: dict, alpha6: float) -> dict:
    raw = _ORIGINAL_MODE_CELL(
        mode, x, rho_trans, sigma, rs, live, vector, process, sched, alpha6
    )
    beta = _measurement_beta_upper(mode, sigma, rs, live, vector, process, sched)
    tO, tS, factor, qnorm = _translation_direct_blocks(x, sigma, raw, beta, sched)
    trans_delta = _certified_generalized_delta(tO, tS, MIN_USEFUL_DELTA)
    if trans_delta <= 0.0:
        raise RuntimeError(f"{mode} direct translation generalized inequality lost positivity")

    other = math.inf
    qpost = float(raw["post_measurement_scaled_Omega_lambda_min_lower"])
    qpost_interval = Interval.outward_bounds(qpost, qpost)
    trans_indices = {6+a for a in range(3)} | {9+a for a in range(3)} | {12+a for a in range(3)} | {15+a for a in range(3)}
    for i, (u, s2) in enumerate(zip(raw["Sigma_diagonal_upper"], raw["comparison_scale_diagonal_squared"])):
        if i in trans_indices:
            continue
        sui = BASE.up(u / s2)
        sint = Interval.outward_bounds(sui, sui)
        ratio = BASE.down(qpost_interval.lo / sint.hi)
        other = min(other, ratio)

    full_delta = BASE.down(min(trans_delta, other))
    if full_delta <= 0.0:
        raise RuntimeError(f"{mode} direct generalized matrix inequality lost positivity")
    fullO, fullS = _assemble_mode_matrices(mode, raw, tO, tS)
    shrink_count = 0
    full_ok = _spd_at_delta(fullO, fullS, full_delta)
    while not full_ok and shrink_count < 12:
        full_delta = BASE.down(0.5 * full_delta)
        shrink_count += 1
        full_ok = _spd_at_delta(fullO, fullS, full_delta)
    if not full_ok:
        raise RuntimeError(f"{mode} reported direct generalized delta did not re-certify")

    limiting = "translation_RL_inverse_block" if trans_delta <= other else "attitude_bias_or_active_ba_block"

    # Same inequality, with the word-accumulated injected-noise floor in place of
    # the single-step one.  Both are lower comparisons for Omega_word, so the
    # cell keeps whichever certifies more; a failure here can never lower the
    # single-step result.
    word = None
    try:
        wt_delta, wBlock, wSigmaRoot, wdiag = _translation_word_margin(
            x, sigma, rs, raw, live, vector, process, sched, mode
        )
        aO, aS, adiag = _attitude_bias_word_blocks(raw, live, vector, process, sched)
        wa_delta = _certified_generalized_delta(aO, aS, MIN_USEFUL_DELTA)
        if wt_delta <= 0.0 or wa_delta <= 0.0:
            raise RuntimeError("word-accumulated block lost positivity")
        wba = math.inf
        for i in (18, 19, 20):
            if i >= len(raw["Sigma_diagonal_upper"]):
                continue
            sui = BASE.up(raw["Sigma_diagonal_upper"][i] / raw["comparison_scale_diagonal_squared"][i])
            wba = min(wba, BASE.down(qpost / sui))
        # Each block is independently certified, so the mode margin is the
        # minimum over blocks of the *best* margin available for that block.
        # Bundling them would let one block's weaker route discard another
        # block's stronger one: on a short-tau source cell the translation
        # horizon is capped while the attitude floor still carries the full word.
        best_trans = max(wt_delta, trans_delta)
        best_att = max(wa_delta, other)
        word_delta = BASE.down(min(best_trans, best_att, wba))
        if not word_delta > 0.0:
            raise RuntimeError("word-accumulated margin lost positivity")

        # The word route places three 4x4 translation blocks on [v,p,S,a_w] of
        # each axis, three 2x2 blocks on (theta, gyro bias) of each axis, and the
        # remaining active-bias coordinates on the diagonal.  Those index sets
        # partition the mode, so the assembled comparison is block diagonal and
        # is positive definite exactly when every block is -- there is no full
        # assembly for the LDLT to see that the per-block tests do not.
        covered = set()
        for axis in range(3):
            covered.update({axis, 3 + axis, 6 + axis, 9 + axis, 12 + axis, 15 + axis})
        dimension = 18 if mode == "H" else 21
        covered.update(range(18, dimension))
        if covered != set(range(dimension)):
            raise RuntimeError("word-accumulated block partition does not cover the mode")
        def _blocks_hold(value: float) -> bool:
            translation = (
                value <= wt_delta if wt_delta >= trans_delta
                else _spd_at_delta(tO, tS, value)
            )
            attitude = (
                _spd_at_delta(aO, aS, value)
                if wa_delta >= other else value <= other
            )
            return translation and attitude

        wshrink = 0
        word_ok = _blocks_hold(word_delta)
        while not word_ok and wshrink < 12:
            word_delta = BASE.down(0.5 * word_delta)
            wshrink += 1
            word_ok = _blocks_hold(word_delta)
        if not word_ok:
            raise RuntimeError("word-accumulated generalized delta did not re-certify")
        word = {
            "translation_margin_lower": wt_delta,
            "attitude_bias_margin_lower": wa_delta,
            "best_translation_margin_lower": best_trans,
            "best_attitude_bias_margin_lower": best_att,
            "translation_route": "word" if wt_delta >= trans_delta else "single_step",
            "attitude_bias_route": "word" if wa_delta >= other else "single_step",
            "active_accelerometer_bias_margin_lower": None if wba is math.inf else wba,
            "margin_lower": word_delta,
            "recertified": word_ok,
            "block_diagonal_partition_covers_mode": True,
            "translation_inequality_form": "information: (1/delta) Sigma^-1 - (Omega_word^-1 + H' R^-1 H) is SPD",
            "rounding_boundary_downward_shrink_count": wshrink,
            "limiting_block": (
                "translation_RL_inverse_block" if best_trans <= min(best_att, wba)
                else "attitude_gyro_bias_block" if best_att <= wba
                else "active_accelerometer_bias_block"
            ),
            **wdiag,
            **adiag,
        }
    except Exception as exc:  # fail-safe: never worse than the single-step route
        word = {"margin_lower": 0.0, "unavailable_reason": str(exc)}

    word_delta = float(word["margin_lower"])
    if word_delta > full_delta:
        full_delta = word_delta
        limiting = word["limiting_block"]
        noise_route = "WORD_ACCUMULATED_INJECTED_NOISE_FLOOR"
    else:
        noise_route = "SINGLE_STEP_INJECTED_NOISE_FLOOR"

    out = dict(raw)
    out["word_accumulated_noise_floor"] = word
    out["injected_noise_floor_route"] = noise_route
    out["single_step_margin_lower"] = BASE.down(min(trans_delta, other))
    out["relative_Riccati_injection_margin_lower"] = full_delta
    out["direct_translation_generalized_margin_lower"] = trans_delta
    out["direct_nontranslation_margin_lower"] = other
    out["generalized_matrix_inequality"] = {
        "form": "Omega_word_RL_inverse_minus_delta_Sigma_upper_RL_inverse_is_SPD",
        "validated_interval_ldlt": True,
        "full_mode_dimension": len(fullO),
        "reported_delta_recertified": full_ok,
        "rounding_boundary_downward_shrink_count": shrink_count,
        "limiting_block": limiting,
        "measurement_information_beta_upper": beta,
        "translation_Qscaled_row_sum_norm_upper": qnorm,
        "translation_posterior_matrix_factor_lower": factor,
        "translation_process_representation": "x_lo*C*(Q_scaled/x)*C^T",
        "translation_congruence": "C=R*L_inverse applied to both Omega and Sigma",
        "full_margin_composition": "min(certified_translation_block, certified_nontranslation_interval_endpoint_blocks), followed by full H/A LDLT recertification; any recertification adjustment is downward only",
        "old_scalar_rho_over_max_scaled_upper_used": False,
    }
    return out


BASE.mode_cell = mode_cell
BASE._build_cached.cache_clear()


def build(domain_path: Path = DEFAULT_DOMAIN) -> dict:
    d = FACTORED.build(domain_path)
    out = dict(d)
    out["generalized_matrix_backend"] = "DIRECT_RL_INVERSE_CONGRUENCE_INTERVAL_LDLT"
    out["old_scalar_generalized_ratio_used"] = False
    return out


def validate(d: dict) -> list[str]:
    failures = FACTORED.validate(d)
    if d.get("generalized_matrix_backend") != "DIRECT_RL_INVERSE_CONGRUENCE_INTERVAL_LDLT":
        failures.append("direct generalized matrix backend missing")
    if d.get("old_scalar_generalized_ratio_used") is not False:
        failures.append("old scalar generalized ratio still active")
    for mode in ("H", "A"):
        m = d.get("modes", {}).get(mode, {}).get("matrix_comparison", {})
        g = m.get("generalized_matrix_inequality", {})
        if g.get("validated_interval_ldlt") is not True:
            failures.append(f"{mode} direct generalized inequality not LDLT certified")
        if g.get("reported_delta_recertified") is not True:
            failures.append(f"{mode} direct generalized delta not recertified")
        if g.get("old_scalar_rho_over_max_scaled_upper_used") is not False:
            failures.append(f"{mode} old scalar generalized ratio still active")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--domain", type=Path, default=DEFAULT_DOMAIN)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    d = build(args.domain)
    failures = validate(d)
    out = dict(d)
    out["validation_pass"] = not failures
    out["validation_failures"] = failures
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "backend": d["generalized_matrix_backend"],
        "cells": d["cell_partition"],
        "H": d["modes"]["H"],
        "A": d["modes"]["A"],
        "failures": failures,
    }, indent=2, sort_keys=True))
    return 0 if not failures else 2


if __name__ == "__main__":
    raise SystemExit(main())
