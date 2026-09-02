#!/usr/bin/env python3
"""Lean whole-word process-noise floor for the retained OU-III P3 route.

The source-uniform P3 covariance ceiling is a multi-second word quantity.  A
single 5 ms process covariance is a rigorous but useless lower comparison for
that ceiling, especially on the integrated displacement S coordinate where the
ratio scales as (T/h)^7.  This module supplies only the missing process algebra;
it contains no source scan, theorem promotion, usefulness threshold or P4 code.

Translation
===========
For the integrated-OU chain [v,p,S,a_w], N consecutive predictions inject
exactly Q(Nh), the same analytic family as one step.  In the physical word
scaling diag(sigma*T, sigma*T^2, sigma*T^3, sigma), that family depends only on
X=T/tau.  The dependency-preserving exact exponential series from
:mod:`ou3_p3_scaled_process` is therefore read at X.  To invert it tightly we
apply the exact rational congruence C=R L^-1 of the x->0 shape Gramian.  C is
only a conditioning transform; the resulting upper information matrix is
pulled back to the original word-scaled coordinates before use.

For every concrete X in an interval, Q_scaled(X)=X B(X) >= X_lo B(X) in
Loewner order because B(X) is positive definite.  We validate the congruenced
interval family before inversion, so the returned matrix is a rigorous upper
bound on Q_scaled(X)^-1.

Attitude / gyro bias
====================
In the existing P3 scaling the one-step process block dominates rho I and the
prediction is [[1,-c],[0,1]].  The exact doubling identity

  Omega_2N = Phi_N Omega_N Phi_N' + Omega_N

gives a whole-word lower floor without interval iteration over hundreds of
samples.  Measurement information is then applied through
(Omega^-1 + D)^-1 with D only on attitude; gyro bias is not measured directly.

All interval matrix operations remain outward rounded.  No floating
singular/eigen decomposition is used.
"""
from __future__ import annotations

from fractions import Fraction
import math

from ou3_interval import (
    Interval,
    hull,
    matrix_add,
    matrix_identity,
    matrix_mul,
    matrix_transpose,
    symmetric_gershgorin_upper,
    symmetric_positive_definite_ldlt,
)
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, matrix_symmetric_hull
import ou3_p3_scaled_process as SCALED

F = Fraction
WORD_EXACT_SERIES_MAX_X = 2.5


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def IF(q: Fraction) -> Interval:
    return I(float(q))


# Exact x->0 LDL conditioning transform from the integrated-OU process shape.
_L_INV = (
    (F(1), F(0), F(0), F(0)),
    (-F(3, 8), F(1), F(0), F(0)),
    (F(1, 15), -F(4, 9), F(1), F(0)),
    (-F(15, 2), F(30), -F(105, 2), F(1)),
)
_R = (F(1), F(10), F(100), F(2))
_CQ = tuple(tuple(_R[i] * _L_INV[i][j] for j in range(4)) for i in range(4))
C = [[IF(v) for v in row] for row in _CQ]
CT = matrix_transpose(C)


def _large_word_normalized_matrix(x: Interval):
    """Exact exponential-branch B(X)=Q_scaled(X)/X on 0.01<=X<=2.5.

    The per-step producer deliberately advertises only the deployed X<=0.25
    range, but its exact coefficient cancellation and range-reduced Lagrange
    remainder are not tied to that endpoint.  We reuse those audited scalar
    routines here and set a separate word-horizon cap.  Any interval that loses
    definiteness is rejected/subdivided by the consumer rather than silently
    widened into a proof.
    """
    if not (SCALED.BRANCH_X <= x.lo <= x.hi <= WORD_EXACT_SERIES_MAX_X):
        raise ValueError("word exact-series interval outside [0.01,2.5]")
    names = (
        ("vv", "vp", "vS", "va"),
        ("vp", "pp", "pS", "pa"),
        ("vS", "pS", "SS", "Sa"),
        ("va", "pa", "Sa", "aa"),
    )
    return [[SCALED._near_exact_normalized_entry(name, x) for name in row] for row in names]


def word_normalized_matrix(x: Interval):
    """Return an enclosure of B(X)=Q_scaled(X)/X for positive X<=2.5."""
    if not (0.0 < x.lo <= x.hi <= WORD_EXACT_SERIES_MAX_X):
        raise ValueError("word X interval outside audited positive range")
    if x.hi < SCALED.BRANCH_X:
        return SCALED.small_normalized_matrix(x)
    if x.lo >= SCALED.BRANCH_X:
        return _large_word_normalized_matrix(x)
    left = SCALED.small_normalized_matrix(
        Interval(x.lo, math.nextafter(SCALED.BRANCH_X, -math.inf))
    )
    right = _large_word_normalized_matrix(Interval(SCALED.BRANCH_X, x.hi))
    return [[hull(left[i][j], right[i][j]) for j in range(4)] for i in range(4)]


def translation_information_upper(x: Interval):
    """Upper bound on Q_scaled(X)^-1, or ``None`` if the interval is too wide."""
    B = word_normalized_matrix(x)
    xlo = I(x.lo)
    Qlo = matrix_symmetric_hull([[xlo * B[i][j] for j in range(4)] for i in range(4)])
    conditioned = matrix_symmetric_hull(matrix_mul(matrix_mul(C, Qlo), CT))
    if not symmetric_positive_definite_ldlt(conditioned)[0]:
        return None
    try:
        inv = matrix_inverse_gauss_jordan(conditioned)
    except Exception:
        return None
    info = matrix_symmetric_hull(matrix_mul(matrix_mul(CT, inv), C))
    return info


def translation_margin_from_information(
    information, sigma_root: list[float], measurement_information_diag: list[float]
) -> float:
    """Lower delta for (Omega^-1+D)^-1 >= delta*Sigma in Sigma-normalized form."""
    if len(sigma_root) != 4 or len(measurement_information_diag) != 4:
        raise ValueError("translation comparison requires four coordinates")
    if any(not (math.isfinite(v) and v > 0.0) for v in sigma_root):
        raise ValueError("Sigma roots must be finite positive")
    if any(not (math.isfinite(v) and v >= 0.0) for v in measurement_information_diag):
        raise ValueError("measurement information must be finite nonnegative")
    g = [I(v) for v in sigma_root]
    M = [[g[i] * information[i][j] * g[j] for j in range(4)] for i in range(4)]
    for i in range(4):
        M[i][i] = M[i][i] + I(measurement_information_diag[i])
    top = symmetric_gershgorin_upper(matrix_symmetric_hull(M))
    if not (math.isfinite(top) and top > 0.0):
        return 0.0
    return down(1.0 / top)


def word_step_doublings(word_horizon_lower_s: float, dt_s: float) -> int:
    if not (word_horizon_lower_s > 0.0 and dt_s > 0.0):
        raise ValueError("positive word horizon and dt required")
    steps = math.floor(down(word_horizon_lower_s / up(dt_s)))
    if steps < 1:
        raise RuntimeError("word does not certainly contain one prediction")
    return int(math.floor(math.log2(steps)))


def _diag(values: list[float]):
    n = len(values)
    return [[I(values[i]) if i == j else I(0.0) for j in range(n)] for i in range(n)]


def _measurement_posterior(Omega, information: list[float]):
    n = len(Omega)
    if len(information) != n:
        raise ValueError("information dimension mismatch")
    D = _diag(information)
    M = matrix_add(matrix_identity(n), matrix_mul(D, Omega))
    Minv = matrix_inverse_gauss_jordan(M)
    return matrix_symmetric_hull(matrix_mul(Omega, Minv))


def attitude_bias_word_noise(
    rho: float, coupling_per_step: float, doublings: int, attitude_information: float
):
    """Posterior lower floor for one scaled (theta,b_g) axis pair."""
    if not (0.0 < rho <= 1.0):
        raise ValueError("rho must lie in (0,1]")
    if doublings < 0 or not (math.isfinite(attitude_information) and attitude_information >= 0.0):
        raise ValueError("invalid doubling count or information")
    Phi = [[I(1.0), I(-abs(coupling_per_step))], [I(0.0), I(1.0)]]
    Omega = [[I(rho), I(0.0)], [I(0.0), I(rho)]]
    A = Phi
    for _ in range(doublings):
        Omega = matrix_symmetric_hull(
            matrix_add(matrix_mul(matrix_mul(A, Omega), matrix_transpose(A)), Omega)
        )
        A = matrix_mul(A, A)
    return _measurement_posterior(Omega, [attitude_information, 0.0])


def generalized_delta(Omega, Sigma, gate: float = 1.0e-18) -> float:
    """Locate a validated generalized lower margin Omega >= delta Sigma."""
    if gate <= 0.0:
        raise ValueError("positive search gate required")

    def holds(delta: float) -> bool:
        q = I(delta)
        A = matrix_symmetric_hull([
            [Omega[i][j] - q * Sigma[i][j] for j in range(len(Omega))]
            for i in range(len(Omega))
        ])
        return symmetric_positive_definite_ldlt(A)[0]

    if holds(gate):
        lo = gate
        hi = gate
        while hi < 1.0:
            trial = min(1.0, hi * 10.0)
            if not holds(trial):
                hi = trial
                break
            lo = trial
            hi = trial
            if trial == 1.0:
                return down(lo)
    else:
        hi = gate
        lo = gate
        for _ in range(80):
            lo /= 10.0
            if holds(lo):
                break
        else:
            return 0.0

    for _ in range(48):
        mid = math.sqrt(lo * hi)
        if holds(mid):
            lo = mid
        else:
            hi = mid
    return down(lo)
