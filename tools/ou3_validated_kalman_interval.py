#!/usr/bin/env python3
"""Validated interval Kalman/Joseph primitives for OU-III source-word enclosure.

All matrix arithmetic is delegated to the outward-rounded interval backend.  The
innovation covariance is inverted only through the validated fixed-pivot solver in
``ou3_interval_linear_algebra``; no NumPy eigensolver/inverse participates in a proof
claim.

These routines are deliberately generic.  The later H/A word producer will supply the
actual implementation-derived prediction/correction/reset matrices and continuous
parameter boxes.  This module provides the rigorous covariance recursion needed to
propagate ``Sigma_KF(g)`` and its information metric through those words.
"""
from __future__ import annotations

from typing import Sequence

from ou3_interval import (
    Interval,
    IntervalMatrix,
    matrix_add,
    matrix_identity,
    matrix_mul,
    matrix_sub,
    matrix_transpose,
)
from ou3_interval_linear_algebra import (
    matrix_inverse_gauss_jordan,
    matrix_symmetric_hull,
    shape,
)


def _check_square(A: Sequence[Sequence[Interval]], label: str) -> int:
    n, m = shape(A)
    if n != m:
        raise ValueError(f"{label} must be square, got {(n, m)}")
    return n


def covariance_predict(
    F: Sequence[Sequence[Interval]],
    P: Sequence[Sequence[Interval]],
    Q: Sequence[Sequence[Interval]],
) -> IntervalMatrix:
    """Outward enclosure of P^- = F P F^T + Q."""
    n = _check_square(P, "P")
    if _check_square(F, "F") != n or _check_square(Q, "Q") != n:
        raise ValueError("prediction matrices must share one square state dimension")
    FP = matrix_mul(F, P)
    pred = matrix_add(matrix_mul(FP, matrix_transpose(F)), Q)
    return matrix_symmetric_hull(pred)


def innovation_covariance(
    H: Sequence[Sequence[Interval]],
    P: Sequence[Sequence[Interval]],
    R: Sequence[Sequence[Interval]],
) -> IntervalMatrix:
    """Outward enclosure of S = H P H^T + R."""
    n = _check_square(P, "P")
    r, c = shape(H)
    if c != n:
        raise ValueError("H column dimension must match P")
    if _check_square(R, "R") != r:
        raise ValueError("R dimension must match H row dimension")
    HP = matrix_mul(H, P)
    S = matrix_add(matrix_mul(HP, matrix_transpose(H)), R)
    return matrix_symmetric_hull(S)


def kalman_gain(
    P: Sequence[Sequence[Interval]],
    H: Sequence[Sequence[Interval]],
    R: Sequence[Sequence[Interval]],
) -> tuple[IntervalMatrix, IntervalMatrix]:
    """Return validated K=P H^T S^-1 and innovation covariance S."""
    S = innovation_covariance(H, P, R)
    Sinv = matrix_inverse_gauss_jordan(S)
    PHt = matrix_mul(P, matrix_transpose(H))
    K = matrix_mul(PHt, Sinv)
    return K, S


def joseph_measurement_update(
    P: Sequence[Sequence[Interval]],
    H: Sequence[Sequence[Interval]],
    R: Sequence[Sequence[Interval]],
) -> dict:
    """Validated Joseph covariance update and deterministic correction map.

    Returns interval enclosures of ``K``, ``S``, ``A_corr=I-KH`` and

        P^+ = A_corr P A_corr^T + K R K^T.

    The Joseph form is used because it preserves the exact implemented positive
    semidefinite decomposition and avoids an unsupported cancellation-sensitive
    ``(I-KH)P`` proof path.
    """
    n = _check_square(P, "P")
    K, S = kalman_gain(P, H, R)
    KH = matrix_mul(K, H)
    A = matrix_sub(matrix_identity(n), KH)

    AP = matrix_mul(A, P)
    left = matrix_mul(AP, matrix_transpose(A))
    KR = matrix_mul(K, R)
    right = matrix_mul(KR, matrix_transpose(K))
    Pplus = matrix_symmetric_hull(matrix_add(left, right))
    return {
        "K": K,
        "S": S,
        "A_correction": A,
        "P_plus": Pplus,
    }
