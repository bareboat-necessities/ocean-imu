#!/usr/bin/env python3
"""Accumulated injected-noise algebra for the complete SEA3 P3 word.

This is algebra only: it does not generate a source, replay a trajectory, or
replace the complete SEA3 word.  It supplies exact prediction-noise accumulation
and conservative measurement-posterior operations used by the full H18/A21
restartable generalized-matrix comparison.

For N pure predictions,

    Omega_N = sum_{k=0}^{N-1} Phi^k Q (Phi^k)^T,

with the exact doubling identity

    Omega_2N = Phi^N Omega_N (Phi^N)^T + Omega_N.

For an upper measurement-information matrix D, assimilating that information
with the optimal gain gives the conservative lower noise covariance

    (Omega^-1 + D)^-1.

A consumer must prove that its D upper covers the shipping measurements and
must assemble the final full-state matrix itself.  This helper cannot promote
P3.
"""
from __future__ import annotations

import math

from ou3_interval import (
    Interval,
    matrix_add,
    matrix_identity,
    matrix_mul,
    matrix_transpose,
)
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan, matrix_symmetric_hull

SCHEMA = 1
QUALIFICATION = "OU3_COMPLETE_SEA3_WORD_NOISE_ACCUMULATION_ALGEBRA"


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def diagonal(values: list[float]):
    n = len(values)
    return [[I(values[i]) if i == j else I(0.0) for j in range(n)] for i in range(n)]


def accumulate_word_noise(Q, Phi, doublings: int, rebalance=None):
    if doublings < 0:
        raise ValueError("doubling count must be non-negative")
    Omega = matrix_symmetric_hull(Q)
    A = [list(row) for row in Phi]
    if rebalance is None:
        for _ in range(doublings):
            Omega = matrix_symmetric_hull(
                matrix_add(matrix_mul(matrix_mul(A, Omega), matrix_transpose(A)), Omega)
            )
            A = matrix_mul(A, A)
        return Omega
    if len(rebalance) != len(Omega) or any(not (v > 0.0) for v in rebalance):
        raise ValueError("rebalancing factors must be positive and match the matrix order")
    J = diagonal(list(rebalance))
    Jinv = diagonal([1.0 / v for v in rebalance])
    for _ in range(doublings):
        JA = matrix_mul(J, A)
        carried = matrix_mul(matrix_mul(J, Omega), matrix_transpose(J))
        Omega = matrix_symmetric_hull(
            matrix_add(matrix_mul(matrix_mul(JA, Omega), matrix_transpose(JA)), carried)
        )
        A = matrix_mul(matrix_mul(JA, A), Jinv)
    return Omega


def measurement_posterior(Omega, information):
    """Return conservative lower ``(Omega^-1 + D)^-1`` for PSD D upper."""
    n = len(Omega)
    if len(information) != n:
        raise ValueError("information diagonal does not match matrix order")
    if any(d < 0.0 or not math.isfinite(d) for d in information):
        raise ValueError("measurement information must be finite non-negative")
    D = diagonal(list(information))
    Oinv = matrix_inverse_gauss_jordan(matrix_symmetric_hull(Omega))
    M = matrix_symmetric_hull(matrix_add(Oinv, D))
    return matrix_symmetric_hull(matrix_inverse_gauss_jordan(M))


def attitude_bias_word_noise(
    rho_attitude: float,
    coupling_per_step: float,
    doublings: int,
    attitude_information: float,
):
    """Scaled one-axis (theta,b_g) accumulated/posterior lower.

    The scaled homogeneous transition is conservatively represented by
    ``[[1,-c],[0,1]]``.  Rotation is orthogonal and therefore does not reduce
    an isotropic process lower.  Only theta receives direct measurement
    information; b_g has no direct shipping measurement row.
    """
    if not (0.0 < rho_attitude <= 1.0):
        raise ValueError("scaled attitude/bias process floor must lie in (0,1]")
    Phi = [[I(1.0), I(-abs(coupling_per_step))], [I(0.0), I(1.0)]]
    Q = [[I(rho_attitude), I(0.0)], [I(0.0), I(rho_attitude)]]
    Omega = accumulate_word_noise(Q, Phi, doublings)
    return measurement_posterior(Omega, [attitude_information, 0.0])


def validate_algebra() -> list[str]:
    failures: list[str] = []
    try:
        q = [[I(1.0), I(0.0)], [I(0.0), I(1.0)]]
        p = [[I(1.0), I(0.1)], [I(0.0), I(1.0)]]
        o = accumulate_word_noise(q, p, 1)
        if len(o) != 2 or len(o[0]) != 2:
            failures.append("doubling algebra returned wrong shape")
        post = measurement_posterior(o, [1.0, 0.0])
        if len(post) != 2:
            failures.append("posterior algebra returned wrong shape")
    except Exception as exc:
        failures.append(f"algebra smoke failed: {exc}")
    return failures


if __name__ == "__main__":
    f = validate_algebra()
    print({"qualification": QUALIFICATION, "failures": f})
    raise SystemExit(0 if not f else 2)
