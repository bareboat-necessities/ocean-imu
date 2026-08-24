#!/usr/bin/env python3
"""Trusted exact-group algebra for the OU-III P4 source-word certificate.

The quantitative P4 chart uses the exact Cayley/Rodrigues coordinate

    c(R) = 2 tan(theta/2) u = 4 e_R / (1 + tr R),   theta < pi.

This chart is especially useful here because its local coordinate is exactly
``delta_theta`` and a left group product has the rational composition law

    c(R_d R) = (c_d + c + 0.5 c_d x c) /
               (1 - 0.25 c_d' c).

Thus the full source-varying Kalman information matrix from P3 can be lifted to
an exact SO(3) Lyapunov metric without discarding attitude--linear cross terms.
The antipodal set remains excluded exactly where the Cayley denominator
vanishes.

The deployed correction itself is mirrored exactly:

* for ||dtheta|| < 1e-2, the source polynomial quaternion coefficients are used
  and normalized;
* otherwise the axis-angle quaternion is used and normalized.

For Cayley composition the quaternion normalization cancels: if the unnormalized
source quaternion is (w,v), its exact Cayley vector is 2v/w.  This avoids a
spurious normalization dependency in the proof map.

All algebraic operations use :mod:`ou3_interval`; trigonometric kernels use the
validated rational-Taylor layer in :mod:`ou3_validated_transcendentals`.
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
CAYLEY_MONOTONE_NORM_MAX = 3.0


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


def dot(a: Sequence[Interval], b: Sequence[Interval]) -> Interval:
    if len(a) != len(b):
        raise ValueError("dot product dimension mismatch")
    total = Interval.point(0.0)
    for x, y in zip(a, b):
        total = total + x*y
    return total


def cross(a: Sequence[Interval], b: Sequence[Interval]) -> list[Interval]:
    if len(a) != 3 or len(b) != 3:
        raise ValueError("cross product requires three-vectors")
    return [
        a[1]*b[2]-a[2]*b[1],
        a[2]*b[0]-a[0]*b[2],
        a[0]*b[1]-a[1]*b[0],
    ]


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


def _vector_hull(a, b):
    if len(a) != len(b):
        raise ValueError("vector hull shape mismatch")
    return [hull(x, y) for x, y in zip(a, b)]


def quaternion_rotation(q: Sequence[Interval]):
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
    if theta.hi > CAYLEY_MONOTONE_NORM_MAX:
        raise ValueError("axis-angle quaternion proof box exceeds promoted P4 chart range")
    half = Interval(down(0.5*theta.lo), up(0.5*theta.hi))
    w = Interval(VT.cos_point(half.hi).lo, VT.cos_point(half.lo).hi)
    k = Interval.point(0.5) * VT.sinc_interval(half)
    return normalize_quaternion([w, k*d[0], k*d[1], k*d[2]])


def deployed_injection_rotation(d: Sequence[Interval]):
    theta = vector_norm_interval(d)
    if theta.hi < SERIES_BRANCH_NORM:
        return quaternion_rotation(_series_quaternion(d))
    if theta.lo >= SERIES_BRANCH_NORM:
        return quaternion_rotation(_axis_angle_quaternion(d))
    return _matrix_hull(
        quaternion_rotation(_series_quaternion(d)),
        quaternion_rotation(_axis_angle_quaternion(d)),
    )


def _series_injection_cayley(d: Sequence[Interval]) -> list[Interval]:
    """Exact Cayley vector of the normalized source series quaternion.

    Normalization cancels from 2v/w, so this is an exact rational expression
    in the same source polynomial coefficients.
    """
    theta = vector_norm_interval(d)
    t2 = theta.square()
    t4 = t2.square()
    w = Interval.point(1.0) - Interval.point(1.0/8.0)*t2 + Interval.point(1.0/384.0)*t4
    if w.lo <= 0.0:
        raise RuntimeError("series quaternion scalar part crosses zero")
    k = Interval.point(0.5) - Interval.point(1.0/48.0)*t2 + Interval.point(1.0/3840.0)*t4
    coeff = Interval.point(2.0) * k / w
    return [coeff*x for x in d]


def _axis_angle_injection_cayley(d: Sequence[Interval]) -> list[Interval]:
    theta = vector_norm_interval(d)
    if theta.hi > CAYLEY_MONOTONE_NORM_MAX:
        raise ValueError("correction norm exceeds promoted Cayley chart range")
    half = Interval(down(0.5*theta.lo), up(0.5*theta.hi))
    cos_half = Interval(VT.cos_point(half.hi).lo, VT.cos_point(half.lo).hi)
    if cos_half.lo <= 0.0:
        raise RuntimeError("correction reaches Cayley antipode")
    coeff = VT.sinc_interval(half) / cos_half
    return [coeff*x for x in d]


def deployed_injection_cayley(d: Sequence[Interval]) -> list[Interval]:
    """Exact Cayley vector of the deployed normalized quaternion correction."""
    theta = vector_norm_interval(d)
    if theta.hi < SERIES_BRANCH_NORM:
        return _series_injection_cayley(d)
    if theta.lo >= SERIES_BRANCH_NORM:
        return _axis_angle_injection_cayley(d)
    return _vector_hull(_series_injection_cayley(d), _axis_angle_injection_cayley(d))


def cayley_compose_left(c_left: Sequence[Interval], c_right: Sequence[Interval]) -> list[Interval]:
    """Exact Cayley coordinate of R(c_left) R(c_right)."""
    if len(c_left) != 3 or len(c_right) != 3:
        raise ValueError("Cayley composition requires three-vectors")
    denom = Interval.point(1.0) - Interval.point(0.25)*dot(c_left, c_right)
    if denom.lo <= 0.0 <= denom.hi:
        raise RuntimeError("Cayley composition denominator crosses antipodal singularity")
    cr = cross(c_left, c_right)
    num = [
        c_left[i] + c_right[i] + Interval.point(0.5)*cr[i]
        for i in range(3)
    ]
    return [x/denom for x in num]


def cayley_coordinate(R) -> list[Interval]:
    """Exact c(R)=4 e_R/(1+tr R), valid only off the antipodal set."""
    if len(R) != 3 or any(len(row) != 3 for row in R):
        raise ValueError("Cayley coordinate requires 3x3 matrix")
    tr = R[0][0] + R[1][1] + R[2][2]
    denom = Interval.point(1.0) + tr
    if denom.lo <= 0.0:
        raise RuntimeError("rotation interval reaches the Cayley antipodal set")
    e = [
        Interval.point(0.5)*(R[2][1]-R[1][2]),
        Interval.point(0.5)*(R[0][2]-R[2][0]),
        Interval.point(0.5)*(R[1][0]-R[0][1]),
    ]
    return [Interval.point(4.0)*x/denom for x in e]


def rotation_from_cayley(c: Sequence[Interval]):
    """Exact inverse Cayley map for c=2 tan(theta/2)u."""
    if len(c) != 3:
        raise ValueError("inverse Cayley map requires a three-vector")
    c2 = dot(c, c)
    denom = Interval.point(4.0) + c2
    K = skew(c)
    K2 = matrix_mul(K, K)
    return matrix_add(
        matrix_add(matrix_identity(3), _matrix_scale(K, Interval.point(4.0)/denom)),
        _matrix_scale(K2, Interval.point(2.0)/denom),
    )


def deployed_correct_cayley(c_error: Sequence[Interval], dtheta: Sequence[Interval]) -> list[Interval]:
    return cayley_compose_left(deployed_injection_cayley(dtheta), c_error)


def rodrigues_rotation(rotation_vector: Sequence[Interval]):
    theta = vector_norm_interval(rotation_vector)
    if theta.hi > CAYLEY_MONOTONE_NORM_MAX:
        raise ValueError("Rodrigues rotation exceeds promoted P4 chart range")
    K = skew(rotation_vector)
    K2 = matrix_mul(K, K)
    return matrix_add(
        matrix_add(matrix_identity(3), _matrix_scale(K, VT.sinc_interval(theta))),
        _matrix_scale(K2, VT.cosc_interval(theta)),
    )


def group_energy(R) -> Interval:
    if len(R) != 3 or any(len(row) != 3 for row in R):
        raise ValueError("group energy requires 3x3 matrix")
    tr = R[0][0] + R[1][1] + R[2][2]
    return Interval.point(0.5) * (Interval.point(3.0) - tr)


def corrected_group_energy(R_error, dtheta: Sequence[Interval]) -> Interval:
    return group_energy(matrix_mul(deployed_injection_rotation(dtheta), R_error))


def exact_energy_change_identity(R_error, dtheta: Sequence[Interval]) -> Interval:
    delta = vector_norm_interval(dtheta)
    S = VT.sinc_interval(delta)
    C = VT.cosc_interval(delta)
    e = [
        Interval.point(0.5)*(R_error[2][1]-R_error[1][2]),
        Interval.point(0.5)*(R_error[0][2]-R_error[2][0]),
        Interval.point(0.5)*(R_error[1][0]-R_error[0][1]),
    ]
    K = skew(dtheta)
    K2R = matrix_mul(matrix_mul(K, K), R_error)
    tr = K2R[0][0] + K2R[1][1] + K2R[2][2]
    return S*dot(dtheta, e) - Interval.point(0.5)*C*tr


def point_vector(v: Sequence[float]) -> list[Interval]:
    return [Interval.point(float(x)) for x in v]


def point_matrix(A: Sequence[Sequence[float]]):
    return [[Interval.point(float(x)) for x in row] for row in A]


def transpose(R):
    return matrix_transpose(R)
