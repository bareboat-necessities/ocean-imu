#!/usr/bin/env python3
"""Verified inverse enclosure for small SPD innovation matrices.

The existing fixed-pivot interval Gauss-Jordan inverse can fail when a broad
interval family makes an intermediate pivot cross zero even though every real
innovation covariance is positive definite.  Falling back only to ``S >= R``
then destroys correlation and can make the Kalman-gain box enormous.

This module supplies a source-independent verified alternative.  For an
interval square matrix S, choose a binary64 point preconditioner C from the
inverse of the midpoint matrix and form

    E = I - C S.

If an outward-rounded infinity-norm bound q = ||E||_inf is strictly below one,
then I-E, C, and S are nonsingular and

    S^-1 = (I-E)^-1 C
         = C + E C + sum_{k>=2} E^k C.

The first two terms are evaluated with interval arithmetic.  The remaining
entrywise error is bounded by

    q^2 ||C||_inf / (1-q).

No floating-point inverse is trusted as an enclosure: C is only a fixed point
preconditioner.  Every acceptance claim comes from outward interval arithmetic
and the validated q<1 Neumann criterion.  The result is symmetrically
intersected because the intended innovation matrices are symmetric; callers
must still establish SPD/symmetry from their own source tuple.
"""
from __future__ import annotations

import math
from typing import Sequence

from ou3_interval import (
    Interval,
    matrix_abs_row_sum_upper,
    matrix_identity,
    matrix_mul,
    matrix_sub,
    up,
)
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, shape

SCHEMA = 1


class VerifiedInverseFailure(RuntimeError):
    pass


def _midpoint(x: Interval) -> float:
    return float(x.lo + 0.5 * (x.hi - x.lo))


def _point_preconditioner(S: Sequence[Sequence[Interval]]):
    n, m = shape(S)
    if n != m or n == 0:
        raise ValueError("verified inverse requires a nonempty square matrix")
    M = [[Interval.point(_midpoint(S[i][j])) for j in range(n)] for i in range(n)]
    try:
        Minv = matrix_inverse_gauss_jordan(M)
    except Exception as exc:
        raise VerifiedInverseFailure(f"midpoint matrix inverse failed: {type(exc).__name__}: {exc}") from exc
    C = [[Interval.point(_midpoint(Minv[i][j])) for j in range(n)] for i in range(n)]
    return C


def _intersect(a: Interval, b: Interval) -> Interval:
    lo = max(a.lo, b.lo)
    hi = min(a.hi, b.hi)
    if lo > hi:
        raise VerifiedInverseFailure("symmetric inverse enclosures became disjoint")
    return Interval(lo, hi)


def inverse_enclosure(S: Sequence[Sequence[Interval]]) -> tuple[list[list[Interval]], dict]:
    """Return a verified inverse enclosure and audit metadata.

    Raises ``VerifiedInverseFailure`` unless the midpoint-preconditioned Neumann
    criterion is strict.  It does not silently fall back to a weaker bound.
    """
    n, m = shape(S)
    if n != m or n == 0:
        raise ValueError("verified inverse requires a nonempty square matrix")
    C = _point_preconditioner(S)
    E = matrix_sub(matrix_identity(n), matrix_mul(C, S))
    q = matrix_abs_row_sum_upper(E)
    if not math.isfinite(q) or not q < 1.0:
        raise VerifiedInverseFailure(f"Neumann criterion failed: ||I-CS||_inf upper={q!r}")

    EC = matrix_mul(E, C)
    cnorm = matrix_abs_row_sum_upper(C)
    q2 = up(q * q)
    denom = 1.0 - q
    if not denom > 0.0:
        raise VerifiedInverseFailure("Neumann denominator is not positive")
    tail = up(up(q2 * cnorm) / denom)
    if not math.isfinite(tail):
        raise VerifiedInverseFailure("Neumann remainder bound is not finite")
    rem = Interval(-tail, tail)

    X = [[C[i][j] + EC[i][j] + rem for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            z = _intersect(X[i][j], X[j][i])
            X[i][j] = z
            X[j][i] = z

    return X, {
        "schema": SCHEMA,
        "qualification": "OU3_VERIFIED_MIDPOINT_NEUMANN_INVERSE",
        "dimension": n,
        "neumann_q_inf_upper": q,
        "preconditioner_inf_norm_upper": cnorm,
        "tail_entry_abs_upper": tail,
        "criterion_strict": True,
        "ordinary_float_inverse_used_as_enclosure": False,
        "fixed_point_preconditioner_only": True,
    }
