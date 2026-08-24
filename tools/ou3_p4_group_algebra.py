#!/usr/bin/env python3
"""Trusted exact-group algebra for the OU-III P4 source-word certificate.

This module mirrors the deployed attitude injection in KalmanOUCoreMath.h:

* for ||dtheta|| < 1e-2 it uses the same polynomial quaternion coefficients
  w = 1-theta^2/8+theta^4/384 and
  k = 1/2-theta^2/48+theta^4/3840, followed by quaternion normalization;
* otherwise it uses the exact axis-angle quaternion and normalization.

The result is converted to a rotation matrix with interval arithmetic.  The
module also provides the exact group energy V_R=tr(I-R)/2 and the Rodrigues
finite-correction energy identity used by the manuscript.  No small-angle
linearized attitude update is present here.

All algebraic operations use :mod:`ou3_interval`; trigonometric kernels use the
validated rational-Taylor layer in :mod:`ou3_validated_transcendentals`.  The
only square-root primitive is audited against exact binary64 rationals and
outward corrected with nextafter.
"""
from __future__ import annotations

from fractions import Fraction
import math
from typing import Sequence

from ou3_interval import (
    Interval,
    hull,
    matrix_add,
    matrix_identity,
    matrix_mul,
    matrix_transpose,
)
import ou3_validated_transcendentals as VT

SERIES_BRANCH_NORM = 1.0e-2


def down(x: float) -> float:
    return math.nextafter(float(x), -math.inf)


def up(x: float) -> float:
    return math.nextafter(float(x), math.inf)


def sqrt_point(x: float) -> Interval:
    """Outward enclosure of sqrt(x) for one nonnegative binary64 input."""
    x = float(x)
    if not math.isfinite(x) or x < 0.0:
        raise ValueError("sqrt_point requires finite x >= 0")
    if x == 0.0:
        return Interval.point(0.0)
    q = Fraction.from_float(x)
    f = math.sqrt(x)
    fq = Fraction.from_float(f)
    lo = f
    hi = f
    if fq * fq > q:
        lo = math.nextafter(lo, -math.inf)
    if fq * fq < q:
        hi = math.nextafter(hi, math.inf)
    # One extra representable value makes the containment independent of any
    # platform claim stronger than IEEE sqrt's usual near-correct rounding.
    return Interval(down(lo), up(hi))


def sqrt_interval(x: Interval) -> Interval:
    if x.lo < 0.0:
        raise ValueError("sqrt interval crosses negative values")
    return Interval(sqrt_point(x.lo).lo, sqrt_point(x.hi).hi)


def vector_norm_interval(v: Sequence[Interval]) -> Interval:
    total = Interval.point(0.0)
    for x in v:
        total = total + x.square()
    return sqrt_interval(total)


def skew(v: Sequence[Interval]):
    if len(v) != 3:
        raise ValueError("skew requires three components")
    z = Interval.point(0.0)
    x, y, zz = v
    return [[z, -zz, y], [zz, z, -x], [-y, x, z]]


def _matrix_scale(A, s: Interval):
    return [[s * A[i][j] for j in range(len(A[0]))] for i in range(len(A))]


def _matrix_hull(A, B):
    if len(A) != len(B) or (A and len(A[0]) != len(B[0])):
        raise ValueError("matrix hull shape mismatch")
    return [[hull(A[i][j], B[i][j]) for j in range(len(A[0]))] for i in range(len(A))]


def quaternion_rotation(q: Sequence[Interval]):
    """Rotation matrix of an already normalized interval quaternion [w,x,y,z]."""
    if len(q) != 4:
        raise ValueError("quaternion requires four components")
    w, x, y, z = q
    one = Interval.point(1.0)
    two = Interval.point(2.0)
    return [
        [one-two*(y.square()+z.square()), two*(x*y-z*w), two*(x*z+y*w)],
        [two*(x*y+z*w), one-two*(x.square()+z.square()), two*(y*z-x*w)],
        [two*(x*z-y*w), two*(y*z+x*w), one-two*(x.square()+y.square())],
    ]


def normalize_quaternion(q: Sequence[Interval]):
    n2 = Interval.point(0.0)
    for x in q:
        n2 = n2 + x.square()
    n = sqrt_interval(n2)
    if n.lo <= 0.0:
        raise RuntimeError("quaternion norm interval crosses zero")
    return [x / n for x in q]


def _series_quaternion(d: Sequence[Interval]):
    theta = vector_norm_interval(d)
    t2 = theta.square()
    t4 = t2.square()
    w = Interval.point(1.0) - Interval.point(1.0/8.0)*t2 + Interval.point(1.0/384.0)*t4
    k = Interval.point(0.5) - Interval.point(1.0/48.0)*t2 + Interval.point(1.0/3840.0)*t4
    return normalize_quaternion([w, k*d[0], k*d[1], k*d[2]])


def _axis_angle_quaternion(d: Sequence[Interval]):
    theta = vector_norm_interval(d)
    half = Interval(down(0.5*theta.lo), up(0.5*theta.hi))
    # The P4 chart/correction domain is kept below pi, so half-angle lies in a
    # monotone cos sector.  Use point enclosures at endpoints, with the lower
    # endpoint attained at the larger half-angle.
    c_lo = VT.cos_point(half.hi).lo
    c_hi = VT.cos_point(half.lo).hi
    w = Interval(c_lo, c_hi)
    sinc_half = VT.sinc_interval(half)
    k = Interval.point(0.5) * sinc_half
    return normalize_quaternion([w, k*d[0], k*d[1], k*d[2]])


def deployed_injection_rotation(d: Sequence[Interval]):
    """Outward enclosure of the exact deployed quat_from_delta_theta rotation."""
    theta = vector_norm_interval(d)
    if theta.hi < SERIES_BRANCH_NORM:
        return quaternion_rotation(_series_quaternion(d))
    if theta.lo >= SERIES_BRANCH_NORM:
        return quaternion_rotation(_axis_angle_quaternion(d))
    # A source box may straddle the C++ branch.  Both exact branch images are
    # enclosed and hulled; no branch is selected optimistically.
    return _matrix_hull(
        quaternion_rotation(_series_quaternion(d)),
        quaternion_rotation(_axis_angle_quaternion(d)),
    )


def rodrigues_rotation(rotation_vector: Sequence[Interval]):
    """Exact exp([r]x) enclosure using sinc/cosc, for comparison/audit."""
    theta = vector_norm_interval(rotation_vector)
    if theta.hi > VT.MAX_TRIG_ARGUMENT:
        raise ValueError("Rodrigues rotation exceeds audited trig range")
    K = skew(rotation_vector)
    K2 = matrix_mul(K, K)
    return matrix_add(
        matrix_add(matrix_identity(3), _matrix_scale(K, VT.sinc_interval(theta))),
        _matrix_scale(K2, VT.cosc_interval(theta)),
    )


def group_energy(R) -> Interval:
    """Exact V_R(R)=1/2 tr(I-R) interval."""
    if len(R) != 3 or any(len(row) != 3 for row in R):
        raise ValueError("group energy requires 3x3 matrix")
    tr = R[0][0] + R[1][1] + R[2][2]
    return Interval.point(0.5) * (Interval.point(3.0) - tr)


def corrected_group_energy(R_error, dtheta: Sequence[Interval]) -> Interval:
    return group_energy(matrix_mul(deployed_injection_rotation(dtheta), R_error))


def exact_energy_change_identity(R_error, dtheta: Sequence[Interval]) -> Interval:
    """Equation (widen-exact-energy-change), evaluated with validated intervals."""
    delta = vector_norm_interval(dtheta)
    S = VT.sinc_interval(delta)
    C = VT.cosc_interval(delta)
    # e_R = vex((R-R')/2)
    half = Interval.point(0.5)
    e = [
        half*(R_error[2][1]-R_error[1][2]),
        half*(R_error[0][2]-R_error[2][0]),
        half*(R_error[1][0]-R_error[0][1]),
    ]
    dot = Interval.point(0.0)
    for i in range(3):
        dot = dot + dtheta[i]*e[i]
    K2R = matrix_mul(matrix_mul(skew(dtheta), skew(dtheta)), R_error)
    tr = K2R[0][0] + K2R[1][1] + K2R[2][2]
    return S*dot - Interval.point(0.5)*C*tr


def point_vector(v: Sequence[float]) -> list[Interval]:
    return [Interval.point(float(x)) for x in v]


def point_matrix(A: Sequence[Sequence[float]]):
    return [[Interval.point(float(x)) for x in row] for row in A]


def transpose(R):
    return matrix_transpose(R)
