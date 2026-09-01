#!/usr/bin/env python3
"""Correlation-preserving interval enclosure of the Kalman gain.

The coarse P4/P5 backends can lose essentially all useful information when they
form an interval box for ``S^-1`` independently of ``P H^T``.  This module avoids
that split.  For any finite point matrix K0, the exact shipping Kalman gain

    K = P H^T (H P H^T + R)^-1

satisfies

    (K-K0) S = B,
    B = P H^T - K0 S,
    S = H P H^T + R.

The same residual also has the algebraically equivalent representation

    B = (I-K0 H) P H^T - K0 R.

Both interval evaluations enclose the same exact B, so their entrywise
intersection is rigorous and often preserves important P/H correlations before
scalarization.  Since S >= R and R has a certified positive eigenvalue lower
r_min, every row obeys

    ||(K-K0)_i||_2 <= ||B_i||_2 / r_min.

Thus K can be enclosed without an interval inverse of S.  The point K0 only
improves tightness; it is not trusted as a proof value.

For Joseph-information consumers, the exact identity

    H K = I - R S^-1

also gives ``S^-1 = R^-1 (I-HK)`` from the same rigorous K enclosure.  This is
useful to #449 while the K enclosure itself is useful to #450's state-return
map.  No filter parameter or theorem-domain bound is changed here.
"""
from __future__ import annotations

import math
from typing import Sequence

from ou3_interval import (
    Interval,
    IntervalMatrix,
    down,
    up,
    matrix_add,
    matrix_identity,
    matrix_mul,
    matrix_point,
    matrix_sub,
    matrix_transpose,
    symmetric_gershgorin_lower,
)
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, shape


class CorrelatedGainFailure(RuntimeError):
    """Raised when the fail-closed assumptions for the gain enclosure are absent."""


def _mid(x: Interval) -> float:
    return 0.5 * x.lo + 0.5 * x.hi


def _point_mid_matrix(A: Sequence[Sequence[Interval]]) -> IntervalMatrix:
    return matrix_point([[_mid(x) for x in row] for row in A])


def _intersect(a: Interval, b: Interval) -> Interval:
    lo = max(a.lo, b.lo)
    hi = min(a.hi, b.hi)
    if lo > hi:
        raise CorrelatedGainFailure(
            f"algebraically identical gain residual enclosures became disjoint: {a} vs {b}"
        )
    return Interval(lo, hi)


def _matrix_intersection(A, B) -> IntervalMatrix:
    sa = shape(A)
    sb = shape(B)
    if sa != sb:
        raise ValueError(f"residual shape mismatch {sa} != {sb}")
    return [[_intersect(A[i][j], B[i][j]) for j in range(sa[1])] for i in range(sa[0])]


def _row_norm_upper(row: Sequence[Interval]) -> float:
    s = 0.0
    for x in row:
        a = x.abs_upper()
        s = up(s + up(a * a))
    return up(math.sqrt(s))


def _point_gain(P, H, R) -> list[list[float]]:
    """Return an untrusted finite point K0 used only as an enclosure center."""
    P0 = _point_mid_matrix(P)
    H0 = _point_mid_matrix(H)
    R0 = _point_mid_matrix(R)
    PHt0 = matrix_mul(P0, matrix_transpose(H0))
    S0 = matrix_add(matrix_mul(H0, PHt0), R0)
    S0inv = matrix_inverse_gauss_jordan(S0)
    K0box = matrix_mul(PHt0, S0inv)
    return [[_mid(x) for x in row] for row in K0box]


def gain_enclosure(P, H, R) -> dict:
    """Enclose K=P H^T S^-1 through its correlated gain equation.

    The only inverse-like inequality used for the correction is the certified
    spectral lower bound ``S >= R >= r_min I``.  No interval inverse of S is
    constructed.
    """
    n, np = shape(P)
    m, hn = shape(H)
    rm, rp = shape(R)
    if n != np or hn != n or rm != rp or rm != m:
        raise ValueError("P/H/R dimensions are inconsistent")

    r_min = symmetric_gershgorin_lower(R)
    if not (math.isfinite(r_min) and r_min > 0.0):
        raise CorrelatedGainFailure("R has no certified positive eigenvalue lower bound")

    PHt = matrix_mul(P, matrix_transpose(H))
    S = matrix_add(matrix_mul(H, PHt), R)
    K0f = _point_gain(P, H, R)
    K0 = matrix_point(K0f)

    b_direct = matrix_sub(PHt, matrix_mul(K0, S))
    IminusK0H = matrix_sub(matrix_identity(n), matrix_mul(K0, H))
    b_factored = matrix_sub(matrix_mul(IminusK0H, PHt), matrix_mul(K0, R))
    B = _matrix_intersection(b_direct, b_factored)

    row_residual_norm_upper = [_row_norm_upper(row) for row in B]
    row_gain_radius_upper = [up(x / r_min) for x in row_residual_norm_upper]
    K: IntervalMatrix = []
    for i in range(n):
        rad = row_gain_radius_upper[i]
        K.append([
            Interval(down(K0f[i][j] - rad), up(K0f[i][j] + rad))
            for j in range(m)
        ])

    return {
        "K": K,
        "K0": K0f,
        "PHt": PHt,
        "S": S,
        "residual_direct": b_direct,
        "residual_factored": b_factored,
        "residual_intersection": B,
        "R_eigenvalue_lower": r_min,
        "row_residual_norm_upper": row_residual_norm_upper,
        "row_gain_radius_upper": row_gain_radius_upper,
        "interval_S_inverse_formed": False,
        "correlated_gain_equation_used": True,
    }


def joseph_s_inverse_from_gain(H, R, K) -> IntervalMatrix:
    """Rigorous S^-1 enclosure from exact H K = I - R S^-1."""
    m, n = shape(H)
    nk, mk = shape(K)
    rm, rn = shape(R)
    if n != nk or m != mk or rm != rn or rm != m:
        raise ValueError("H/R/K dimensions are inconsistent")
    Rinv = matrix_inverse_gauss_jordan(R)
    return matrix_mul(Rinv, matrix_sub(matrix_identity(m), matrix_mul(H, K)))
