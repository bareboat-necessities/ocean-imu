#!/usr/bin/env python3
"""Auditable interval linear solves for the OU-III validated word backend.

The existing :mod:`ou3_interval` layer deliberately avoids ordinary floating-point
matrix factorizations.  This module extends that trusted layer with fixed-pivot
Gauss--Jordan elimination performed entirely with outward-rounded ``Interval``
operations.

For an interval matrix family A, every elimination pivot is itself enclosed by an
interval.  The solve is accepted only if each pivot interval excludes zero.  Under
that condition, the fixed elimination sequence is well-defined for every concrete
matrix in the family and the natural interval evaluation encloses every concrete
solve/inverse.  If a pivot touches zero the routine refuses the claim; it never
silently switches to a floating-point pivot or pseudoinverse.

This is intentionally conservative.  Branch subdivision/permutation may later be
used by the source-word backend when a wide box cannot certify a fixed pivot order.
"""
from __future__ import annotations

from typing import Sequence

from ou3_interval import Interval, IntervalMatrix, hull, matrix_identity


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


def matrix_solve_gauss_jordan(
    A: Sequence[Sequence[Interval]],
    B: Sequence[Sequence[Interval]],
) -> IntervalMatrix:
    """Enclose X=A^{-1}B using fixed-pivot interval Gauss--Jordan elimination.

    Every pivot interval must exclude zero.  Exact structural zeros/ones created by
    row normalization and elimination are written as point intervals; those entries
    are identities for every concrete elimination and therefore do not lose rigor.
    """
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

        # Normalize the pivot row.  The pivot entry itself is exactly one for every
        # concrete member after division by its own pivot.
        for j in range(width):
            if j != k:
                aug[k][j] = aug[k][j] / pivot
        aug[k][k] = one

        # Eliminate this column from every other row.  The eliminated entry is
        # exactly zero for every concrete member; other entries use outward interval
        # arithmetic and retain all correlation loss conservatively.
        for i in range(n):
            if i == k:
                continue
            factor = aug[i][k]
            for j in range(width):
                if j != k:
                    aug[i][j] = aug[i][j] - factor * aug[k][j]
            aug[i][k] = zero

    return [[aug[i][n + j] for j in range(cb)] for i in range(n)]


def matrix_inverse_gauss_jordan(
    A: Sequence[Sequence[Interval]],
) -> IntervalMatrix:
    """Enclose every inverse in an interval matrix family with certified pivots."""
    n, m = shape(A)
    if n != m:
        raise ValueError("interval inverse requires a square matrix")
    return matrix_solve_gauss_jordan(A, matrix_identity(n))
