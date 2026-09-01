#!/usr/bin/env python3
"""Verified inverse enclosure for small certified-SPD innovation matrices.

For a source-certified symmetric positive-definite interval family ``S``, choose
a binary64 point preconditioner ``C`` from the inverse of a midpoint matrix and
form ``E = I-CS``.  If an outward-rounded bound ``||E||_inf < 1`` is certified,
the Neumann series gives a rigorous inverse enclosure.

Symmetry is never inferred from a generic square interval box.  Callers must
explicitly certify that their exact source family is symmetric and SPD.  Given
that certificate, paired interval entries are intersected before the solve;
this is rigorous because every exact source matrix has ``S_ij=S_ji``.  The same
source certificate justifies intersecting paired inverse entries afterward.
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

SCHEMA = 2


class VerifiedInverseFailure(RuntimeError):
    pass


def _midpoint(x: Interval) -> float:
    return float(x.lo + 0.5 * (x.hi - x.lo))


def _intersect(a: Interval, b: Interval) -> Interval:
    lo = max(a.lo, b.lo)
    hi = min(a.hi, b.hi)
    if lo > hi:
        raise VerifiedInverseFailure("certified-symmetric paired enclosures became disjoint")
    return Interval(lo, hi)


def _certified_symmetric_family(S: Sequence[Sequence[Interval]], *, symmetric_certified: bool):
    n, m = shape(S)
    if n != m or n == 0:
        raise ValueError("verified inverse requires a nonempty square matrix")
    if symmetric_certified is not True:
        raise VerifiedInverseFailure("source symmetry certificate is required")
    out = [[S[i][j] for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            z = _intersect(out[i][j], out[j][i])
            out[i][j] = z
            out[j][i] = z
    return out


def _point_preconditioner(S: Sequence[Sequence[Interval]]):
    n, m = shape(S)
    if n != m or n == 0:
        raise ValueError("verified inverse requires a nonempty square matrix")
    M = [[Interval.point(_midpoint(S[i][j])) for j in range(n)] for i in range(n)]
    try:
        Minv = matrix_inverse_gauss_jordan(M)
    except Exception as exc:
        raise VerifiedInverseFailure(
            f"midpoint matrix inverse failed: {type(exc).__name__}: {exc}"
        ) from exc
    return [[Interval.point(_midpoint(Minv[i][j])) for j in range(n)] for i in range(n)]


def inverse_enclosure(
    S: Sequence[Sequence[Interval]],
    *,
    symmetric_certified: bool = False,
    spd_certified: bool = False,
) -> tuple[list[list[Interval]], dict]:
    """Return a verified inverse enclosure for a source-certified SPD family.

    ``symmetric_certified`` and ``spd_certified`` are proof premises supplied by
    the source covariance construction.  This routine checks the interval
    algebra needed for inversion but does not manufacture either premise from
    a broad Cartesian matrix box.
    """
    if spd_certified is not True:
        raise VerifiedInverseFailure("source SPD certificate is required")
    Ssym = _certified_symmetric_family(S, symmetric_certified=symmetric_certified)
    n, _ = shape(Ssym)
    C = _point_preconditioner(Ssym)
    E = matrix_sub(matrix_identity(n), matrix_mul(C, Ssym))
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
        "source_symmetry_certified": True,
        "source_SPD_certified": True,
        "paired_input_intersection_used": True,
        "neumann_q_inf_upper": q,
        "preconditioner_inf_norm_upper": cnorm,
        "tail_entry_abs_upper": tail,
        "criterion_strict": True,
        "ordinary_float_inverse_used_as_enclosure": False,
        "fixed_point_preconditioner_only": True,
    }
