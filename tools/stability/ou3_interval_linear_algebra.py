#!/usr/bin/env python3
"""Auditable interval linear solves for the OU-III validated word backend.

The trusted result is always produced by outward-rounded interval arithmetic.
A round-to-nearest midpoint inverse may be used only as a *fixed numerical left
preconditioner*: if the resulting interval Gauss--Jordan pivots cannot exclude
zero, the solve is rejected.  The floating result is therefore never accepted
as an inverse or certificate by itself.

For exactly symmetric interval matrices, the generic inverse first attempts a
fixed-order interval LDL^T solve.  This is particularly important for theorem
matrices with known SPD provenance such as innovation covariance
S=H P H^T+R.  Every accepted LDL^T diagonal interval must remain strictly
positive.  If that condition cannot be certified, the implementation falls back
to the existing preconditioned/raw Gauss--Jordan paths; it never assumes SPD
merely from symmetry.
"""
from __future__ import annotations

import math
from typing import Sequence

from ou3_interval import (
    Interval,
    IntervalMatrix,
    hull,
    matrix_identity,
    matrix_mul,
)


class IntervalPivotError(RuntimeError):
    """Raised when a validated fixed-pivot elimination cannot exclude zero."""


def shape(A: Sequence[Sequence[Interval]]) -> tuple[int, int]:
    rows = len(A)
    cols = len(A[0]) if rows else 0
    if any(len(row) != cols for row in A):
        raise ValueError("ragged interval matrix")
    return rows, cols


def matrix_copy(A: Sequence[Sequence[Interval]]) -> IntervalMatrix:
    shape(A)
    return [[x for x in row] for row in A]


def matrix_symmetric_hull(A: Sequence[Sequence[Interval]]) -> IntervalMatrix:
    """Return one symmetric interval box containing A and A^T."""
    n, m = shape(A)
    if n != m:
        raise ValueError("symmetric hull requires a square matrix")
    out = matrix_copy(A)
    for i in range(n):
        for j in range(i + 1, n):
            x = hull(A[i][j], A[j][i])
            out[i][j] = x
            out[j][i] = x
    return out


def _exactly_symmetric(A: Sequence[Sequence[Interval]]) -> bool:
    n, m = shape(A)
    if n != m:
        return False
    return all(
        A[i][j].lo == A[j][i].lo and A[i][j].hi == A[j][i].hi
        for i in range(n) for j in range(i + 1, n)
    )


def matrix_solve_gauss_jordan(
    A: Sequence[Sequence[Interval]],
    B: Sequence[Sequence[Interval]],
) -> IntervalMatrix:
    """Enclose X=A^{-1}B using fixed-pivot interval Gauss--Jordan elimination."""
    n, m = shape(A)
    rb, cb = shape(B)
    if n != m:
        raise ValueError("interval solve requires a square coefficient matrix")
    if rb != n:
        raise ValueError(f"solve shape mismatch {(n, m)} and {(rb, cb)}")
    if n == 0:
        return []

    aug: IntervalMatrix = [
        [x for x in A[i]] + [x for x in B[i]]
        for i in range(n)
    ]
    width = n + cb
    zero = Interval.point(0.0)
    one = Interval.point(1.0)

    for k in range(n):
        pivot = aug[k][k]
        if pivot.lo <= 0.0 <= pivot.hi:
            raise IntervalPivotError(
                f"pivot {k} contains zero: [{pivot.lo!r}, {pivot.hi!r}]"
            )
        for j in range(width):
            if j != k:
                aug[k][j] = aug[k][j] / pivot
        aug[k][k] = one
        for i in range(n):
            if i == k:
                continue
            factor = aug[i][k]
            for j in range(width):
                if j != k:
                    aug[i][j] = aug[i][j] - factor * aug[k][j]
            aug[i][k] = zero

    return [[aug[i][n + j] for j in range(cb)] for i in range(n)]


def matrix_solve_spd_ldlt(
    A: Sequence[Sequence[Interval]],
    B: Sequence[Sequence[Interval]],
) -> IntervalMatrix:
    """Enclose X=A^{-1}B with a validated fixed-order interval LDL^T solve.

    Acceptance is self-certifying: the outward interval LDL^T factorization must
    produce strictly positive D intervals.  Therefore callers need not trust a
    floating factorization.  If interval dependency prevents strict positivity,
    ``IntervalPivotError`` is raised and no inverse claim is produced.
    """
    n, m = shape(A)
    rb, cb = shape(B)
    if n != m:
        raise ValueError("interval SPD solve requires a square coefficient matrix")
    if rb != n:
        raise ValueError(f"SPD solve shape mismatch {(n, m)} and {(rb, cb)}")
    if n == 0:
        return []

    S = matrix_symmetric_hull(A)
    zero = Interval.point(0.0)
    one = Interval.point(1.0)
    L: IntervalMatrix = [[zero for _ in range(n)] for _ in range(n)]
    D = [zero for _ in range(n)]
    for i in range(n):
        L[i][i] = one

    for j in range(n):
        dj = S[j][j]
        for k in range(j):
            dj = dj - L[j][k] * L[j][k] * D[k]
        if not (math.isfinite(dj.lo) and math.isfinite(dj.hi)) or dj.lo <= 0.0:
            raise IntervalPivotError(
                f"SPD LDLT diagonal {j} not strictly positive: [{dj.lo!r}, {dj.hi!r}]"
            )
        D[j] = dj
        for i in range(j + 1, n):
            num = S[i][j]
            for k in range(j):
                num = num - L[i][k] * L[j][k] * D[k]
            L[i][j] = num / D[j]

    Y: IntervalMatrix = [[zero for _ in range(cb)] for _ in range(n)]
    for i in range(n):
        for c in range(cb):
            y = B[i][c]
            for k in range(i):
                y = y - L[i][k] * Y[k][c]
            Y[i][c] = y

    Z: IntervalMatrix = [[Y[i][c] / D[i] for c in range(cb)] for i in range(n)]

    X: IntervalMatrix = [[zero for _ in range(cb)] for _ in range(n)]
    for i in reversed(range(n)):
        for c in range(cb):
            x = Z[i][c]
            for k in range(i + 1, n):
                x = x - L[k][i] * X[k][c]
            X[i][c] = x
    return X


def matrix_inverse_spd_ldlt(
    A: Sequence[Sequence[Interval]],
) -> IntervalMatrix:
    n, m = shape(A)
    if n != m:
        raise ValueError("interval SPD inverse requires a square matrix")
    return matrix_solve_spd_ldlt(A, matrix_identity(n))


def _midpoint_float_inverse(A: Sequence[Sequence[Interval]]) -> list[list[float]] | None:
    """Ordinary-float approximate inverse used only to choose a preconditioner."""
    n, m = shape(A)
    if n != m:
        raise ValueError("midpoint inverse requires square matrix")
    if n == 0:
        return []
    aug = []
    for i in range(n):
        row = [0.5 * (A[i][j].lo + A[i][j].hi) for j in range(n)]
        row += [1.0 if i == j else 0.0 for j in range(n)]
        if any(not math.isfinite(x) for x in row):
            return None
        aug.append(row)
    width = 2 * n
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(aug[i][k]))
        if not math.isfinite(aug[p][k]) or abs(aug[p][k]) <= 1e-300:
            return None
        if p != k:
            aug[k], aug[p] = aug[p], aug[k]
        pivot = aug[k][k]
        for j in range(width):
            aug[k][j] /= pivot
        for i in range(n):
            if i == k:
                continue
            f = aug[i][k]
            for j in range(width):
                aug[i][j] -= f * aug[k][j]
    C = [[aug[i][n + j] for j in range(n)] for i in range(n)]
    return C if all(math.isfinite(x) for row in C for x in row) else None


def _raw_inverse(A: Sequence[Sequence[Interval]]) -> IntervalMatrix:
    n, m = shape(A)
    if n != m:
        raise ValueError("interval inverse requires a square matrix")
    return matrix_solve_gauss_jordan(A, matrix_identity(n))


def matrix_inverse_gauss_jordan(
    A: Sequence[Sequence[Interval]],
) -> IntervalMatrix:
    """Enclose A^{-1}, preferring a validated symmetric LDL^T path when possible.

    Exact symmetry alone is never accepted as a proof of positive definiteness:
    the LDL^T path is used only when every outward D interval is strictly
    positive.  On failure it falls through to the prior midpoint-preconditioned
    and raw interval Gauss--Jordan implementation.
    """
    n, m = shape(A)
    if n != m:
        raise ValueError("interval inverse requires a square matrix")
    if n == 0:
        return []

    if _exactly_symmetric(A):
        try:
            return matrix_inverse_spd_ldlt(A)
        except IntervalPivotError:
            pass

    approx = _midpoint_float_inverse(A)
    if approx is not None:
        C: IntervalMatrix = [[Interval.point(x) for x in row] for row in approx]
        CA = matrix_mul(C, A)
        try:
            return matrix_solve_gauss_jordan(CA, C)
        except IntervalPivotError:
            pass
    return _raw_inverse(A)
