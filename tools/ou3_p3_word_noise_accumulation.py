#!/usr/bin/env python3
"""Word-accumulated injected-noise floor for the source-reachable P3 comparison.

P3 certifies the generalized endpoint inequality

    Omega_word - delta * Sigma_upper  >>  0,

where ``Omega_word`` is the process noise the deployed recursion injects over a
whole source word and ``Sigma_upper`` bounds the covariance at the word
endpoint.  The direct backend used a valid but extremely weak lower bound for
``Omega_word``: the *single* IMU step's ``Q``.  Because ``Sigma_upper`` is a
word-horizon quantity, that compares a 5 ms noise injection against a 3.1 s
covariance -- a ratio of ``(T/h)^7`` on the ``S`` channel alone, and it is what
held the P3 margin at ``3.8e-35``.

For pure prediction the source runs ``P <- Phi P Phi' + Q`` each step, so after
``N`` steps the injected floor is

    Omega_N = sum_{k=0}^{N-1} Phi^k Q (Phi^k)',

which obeys the exact doubling identity

    Omega_{2m} = Phi^m Omega_m (Phi^m)' + Omega_m,   Phi^{2m} = Phi^m Phi^m.

That identity is what this module supplies, for the (theta, gyro-bias) block.
The translation block does not need it: for the integrated-OU chain the same sum
is exactly ``Q(N h)`` -- the same analytic family the certificate already
validates, read at the word horizon rather than at one step -- so the direct
backend evaluates it there instead, with no recursion and no interval growth.

``N`` is taken from a *lower* bound on the word's step count, so the accumulated
matrix stays a lower bound on what the word actually injects while
``Sigma_upper`` stays an upper bound at the word's horizon upper.

The word's measurement updates shrink that floor.  For the Joseph update with
any implemented gain,

    (I-KH) Omega (I-KH)' + K R K'  >=  (Omega^-1 + H' R^-1 H)^-1,

so assimilating every admissible measurement of the word with the optimal gain
is a conservative lower comparison.  With an information upper bound ``D`` the
posterior floor is ``(Omega^-1 + D)^-1 = Omega (I + D Omega)^-1``, computed here
in that closed form for the small, well-conditioned attitude block.

``D`` is kept *block structured*.  A single scalar bound is what makes the
attitude/gyro-bias comparison collapse: the accelerometer's information lives in
the ``a_w`` coordinate and is worth four orders more than the attitude
information, and the gyro bias receives no direct measurement information at
all.  Charging the ``a_w`` figure against the bias coordinate erases the whole
accumulated bias floor.

Everything is outward rounded.  Every routine returns a lower comparison for the
same quantity the parent bounds, so a consumer takes the better of the two and
can never be widened by this module.
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

def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def _I(x: float) -> Interval:
    return Interval.outward_bounds(float(x), float(x))


def _diagonal(values: list[float]):
    n = len(values)
    return [[_I(values[i]) if i == j else _I(0.0) for j in range(n)] for i in range(n)]


def accumulate_word_noise(Q, Phi, doublings: int, rebalance=None):
    """``sum_{k<2^doublings} Phi^k Q (Phi^k)'`` by the exact doubling identity.

    Held in the raw comparison scaling the accumulated matrix spans fourteen
    decades after nine doublings -- the ``S`` row grows like ``N^7`` while the
    ``a_w`` row grows like ``N`` -- and interval arithmetic loses the
    correlation long before that: a relative input width of ``1e-6`` already
    leaves the enclosure indefinite.

    ``rebalance`` fixes that by carrying the recursion in the scaling of the
    *current* horizon.  With ``J = diag(f)`` mapping the ``T`` scaling to the
    ``2T`` scaling, the same exact identity reads

        Omega_{k+1} = (J Psi_k) Omega_k (J Psi_k)' + J Omega_k J',
        Psi_{k+1}   = J Psi_k Psi_k J^-1,

    and both stay order one at every level.  For the ``[v,p,S,a_w]`` chain
    ``f = (1/2, 1/4, 1/8, 1)`` -- exact binary fractions, so the rebalancing is
    rounding-free.  Pass ``None`` to keep the raw recursion.
    """
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
    J = _diagonal(list(rebalance))
    Jinv = _diagonal([1.0 / v for v in rebalance])
    for _ in range(doublings):
        JA = matrix_mul(J, A)
        JAt = matrix_transpose(JA)
        carried = matrix_mul(matrix_mul(J, Omega), matrix_transpose(J))
        Omega = matrix_symmetric_hull(
            matrix_add(matrix_mul(matrix_mul(JA, Omega), JAt), carried)
        )
        A = matrix_mul(matrix_mul(JA, A), Jinv)
    return Omega


def measurement_posterior(Omega, information: list[float]):
    """Lower comparison ``(Omega^-1 + D)^-1 = Omega (I + D Omega)^-1``, ``D=diag(d)``.

    ``information`` must be an upper bound on the measurement information the
    word can assimilate in these coordinates; a larger ``D`` only lowers the
    result, so an over-estimate stays conservative.
    """
    n = len(Omega)
    if len(information) != n:
        raise ValueError("information vector does not match matrix order")
    if any(d < 0.0 or not math.isfinite(d) for d in information):
        raise ValueError("measurement information must be finite non-negative")
    D = [[_I(information[i]) if i == j else _I(0.0) for j in range(n)] for i in range(n)]
    M = matrix_add(matrix_identity(n), matrix_mul(D, Omega))
    Minv = matrix_inverse_gauss_jordan(M)
    return matrix_symmetric_hull(matrix_mul(Omega, Minv))


def word_step_doublings(word_horizon_lower_s: float, dt_s: float) -> int:
    """Largest ``k`` with ``2^k`` prediction steps certainly inside the word."""
    if not (word_horizon_lower_s > 0.0 and dt_s > 0.0):
        raise ValueError("positive word horizon and step required")
    steps = math.floor(down(word_horizon_lower_s / up(dt_s)))
    if steps < 1:
        raise RuntimeError("source word does not certainly contain one prediction step")
    return int(math.floor(math.log2(steps)))


def attitude_bias_word_noise(rho_attitude: float, coupling_per_step: float,
                             doublings: int, attitude_information: float):
    """Accumulated (theta, gyro-bias) noise floor for one axis, scaled and posterior.

    In the comparison scaling ``diag(sqrt q_theta, sqrt q_bias)`` the source
    attitude/bias prediction is ``[[1, -c],[0, 1]]`` with
    ``c = h sqrt(q_bias/q_theta)``, and the scaled one-step process matrix
    dominates ``rho_attitude I`` -- the same scaled comparison the parent
    certificate already uses.  Only the attitude coordinate receives measurement
    information; the gyro bias is never measured directly.
    """
    if not (0.0 < rho_attitude <= 1.0):
        raise ValueError("scaled attitude/bias process floor must lie in (0,1]")
    Phi = [[_I(1.0), _I(-abs(coupling_per_step))], [_I(0.0), _I(1.0)]]
    Q = [[_I(rho_attitude), _I(0.0)], [_I(0.0), _I(rho_attitude)]]
    Omega = accumulate_word_noise(Q, Phi, doublings)
    return measurement_posterior(Omega, [attitude_information, 0.0])
