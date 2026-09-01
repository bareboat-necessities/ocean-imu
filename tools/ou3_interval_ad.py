#!/usr/bin/env python3
"""First-order outward interval automatic differentiation for OU-III P4.

The path-dependent P4 route needs derivatives of the *complete nonlinear word*,
not scalar Lipschitz constants accumulated after the fact.  This module provides
an intentionally small interval-AD layer over :mod:`ou3_interval` and the exact
Cayley/deployed-quaternion algebra used by the shipping MEKF proof.

Each :class:`AD` value contains an outward interval value and one outward
interval derivative for every independent state coordinate.  Algebraic
operations use the natural interval extension.  The Cayley inverse rotation is
therefore differentiated exactly as a rational map.

The deployed correction quaternion needs one non-algebraic radial branch.  For
``r=||d|| <= 6`` the source axis-angle coefficients are

    w(u) = cos(sqrt(u)/2),
    k(u) = sin(sqrt(u)/2)/sqrt(u),       u=r^2,

with ``v=k d``.  On ``0 <= r <= 6 < 2*pi`` the integral representations give
rigorous derivative bounds

    -1/8 <= dw/du <= 0,
    -1/48 <= dk/du <= 0.

The second bound follows from
``k=1/2 int_0^1 cos(t r/2) dt`` and ``sin(x)<=x`` on the audited range.
Those bounds avoid differentiating a square root at the origin.  The strict
small-angle source-polynomial branch is differentiated algebraically.  When a
cell intersects the 1e-2 branch boundary, values and derivatives are hulled, so
the result is a generalized-Jacobian enclosure of both shipping branches.

This is a proof primitive.  It does not select trajectories, use finite
differences, or invoke NumPy floating-point eigensolvers.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Sequence

from ou3_interval import Interval, hull
import ou3_p5_deployed_quaternion_cayley_cell as QCOMP
import ou3_p4_group_algebra as GROUP


def I(x: float) -> Interval:
    """Return a point interval."""
    return Interval.point(float(x))


def _coerce_interval(x) -> Interval:
    if isinstance(x, Interval):
        return x
    if isinstance(x, (int, float)):
        return I(float(x))
    raise TypeError(f"cannot coerce {type(x)!r} to Interval")


@dataclass(frozen=True)
class AD:
    """One interval value with an interval first-derivative vector."""

    val: Interval
    der: tuple[Interval, ...]

    @property
    def n(self) -> int:
        return len(self.der)

    def _other(self, x) -> "AD":
        if isinstance(x, AD):
            if x.n != self.n:
                raise ValueError("AD derivative dimensions differ")
            return x
        return constant(x, self.n)

    def __add__(self, other):
        o = self._other(other)
        return AD(self.val + o.val, tuple(a + b for a, b in zip(self.der, o.der)))

    __radd__ = __add__

    def __neg__(self):
        return AD(-self.val, tuple(-x for x in self.der))

    def __sub__(self, other):
        return self + (-self._other(other))

    def __rsub__(self, other):
        return self._other(other) - self

    def __mul__(self, other):
        o = self._other(other)
        return AD(
            self.val * o.val,
            tuple(a * o.val + self.val * b for a, b in zip(self.der, o.der)),
        )

    __rmul__ = __mul__

    def reciprocal(self):
        if self.val.lo <= 0.0 <= self.val.hi:
            raise ZeroDivisionError("AD reciprocal interval crosses zero")
        inv = self.val.reciprocal()
        den = self.val.square()
        return AD(inv, tuple(-d / den for d in self.der))

    def __truediv__(self, other):
        return self * self._other(other).reciprocal()

    def __rtruediv__(self, other):
        return self._other(other) / self

    def square(self):
        """Range-aware square with the exact derivative ``2*x*x'``.

        Generic interval multiplication ``x*x`` loses the repeated-variable
        dependency and may produce a negative lower bound when ``x`` straddles
        zero.  That is especially damaging in Cayley/quaternion denominators,
        where sums of squares are known nonnegative.  Preserve the exact scalar
        square range while outward-enclosing its first derivative.
        """
        two = I(2.0)
        return AD(
            self.val.square(),
            tuple(two * self.val * d for d in self.der),
        )


def constant(x, n: int) -> AD:
    """Lift a scalar/interval to an AD constant."""
    z = I(0.0)
    return AD(_coerce_interval(x), tuple(z for _ in range(int(n))))


def independent(x, index: int, n: int) -> AD:
    """Lift one interval state coordinate with derivative e_index."""
    if not 0 <= int(index) < int(n):
        raise IndexError("independent AD coordinate outside derivative dimension")
    d = [I(0.0) for _ in range(int(n))]
    d[int(index)] = I(1.0)
    return AD(_coerce_interval(x), tuple(d))


def independent_vector(values: Sequence[Interval], n: int | None = None, offset: int = 0) -> list[AD]:
    """Create independent AD coordinates for one contiguous state slice."""
    nder = len(values) + int(offset) if n is None else int(n)
    return [independent(x, int(offset) + i, nder) for i, x in enumerate(values)]


def hull_ad(*xs: AD) -> AD:
    """Hull interval values and every derivative coordinate."""
    if not xs:
        raise ValueError("hull_ad requires at least one value")
    n = xs[0].n
    if any(x.n != n for x in xs):
        raise ValueError("AD hull derivative dimensions differ")
    return AD(
        hull(*(x.val for x in xs)),
        tuple(hull(*(x.der[j] for x in xs)) for j in range(n)),
    )


def dot(a: Sequence[AD], b: Sequence[AD]) -> AD:
    if len(a) != len(b) or not a:
        raise ValueError("AD dot requires equal nonempty vectors")
    y = constant(0.0, a[0].n)
    for x, z in zip(a, b):
        y = y + x * z
    return y


def squared_norm(a: Sequence[AD]) -> AD:
    """Return ``sum_i a_i^2`` without repeated-variable interval loss."""
    if not a:
        raise ValueError("AD squared_norm requires a nonempty vector")
    y = constant(0.0, a[0].n)
    for x in a:
        y = y + x.square()
    return y


def cross(a: Sequence[AD], b: Sequence[AD]) -> list[AD]:
    if len(a) != 3 or len(b) != 3:
        raise ValueError("AD cross requires three-vectors")
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]


def matmul(A: Sequence[Sequence[AD]], B: Sequence[Sequence[AD]]) -> list[list[AD]]:
    if not A or not B or len(A[0]) != len(B):
        raise ValueError("AD matrix multiply shape mismatch")
    n = A[0][0].n
    out = []
    for i in range(len(A)):
        row = []
        for j in range(len(B[0])):
            y = constant(0.0, n)
            for k in range(len(B)):
                y = y + A[i][k] * B[k][j]
            row.append(y)
        out.append(row)
    return out


def matvec(A: Sequence[Sequence[AD]], x: Sequence[AD]) -> list[AD]:
    if not A or len(A[0]) != len(x):
        raise ValueError("AD matrix/vector shape mismatch")
    n = x[0].n
    out = []
    for row in A:
        y = constant(0.0, n)
        for a, b in zip(row, x):
            y = y + a * b
        out.append(y)
    return out


def rotation_from_cayley(c: Sequence[AD]) -> list[list[AD]]:
    """Exact AD inverse Cayley map for c=2 tan(theta/2) u."""
    if len(c) != 3:
        raise ValueError("Cayley AD vector must have length three")
    n = c[0].n
    c2 = squared_norm(c)
    den = constant(4.0, n) + c2
    four_over = constant(4.0, n) / den
    two_over = constant(2.0, n) / den
    z = constant(0.0, n)
    x, y, zz = c
    K = [[z, -zz, y], [zz, z, -x], [-y, x, z]]
    K2 = matmul(K, K)
    R = []
    for i in range(3):
        row = []
        for j in range(3):
            base = constant(1.0 if i == j else 0.0, n)
            row.append(base + four_over * K[i][j] + two_over * K2[i][j])
        R.append(row)
    return R


def _norm_upper(values: Sequence[Interval]) -> float:
    s = 0.0
    for x in values:
        a = x.abs_upper()
        s = math.nextafter(s + math.nextafter(a * a, math.inf), math.inf)
    return math.nextafter(math.sqrt(max(0.0, s)), math.inf)


def _small_quaternion_ad(d: Sequence[AD]) -> tuple[AD, list[AD]]:
    n = d[0].n
    u = squared_norm(d)
    u2 = u.square()
    w = constant(1.0, n) - constant(1.0 / 8.0, n) * u + constant(1.0 / 384.0, n) * u2
    k = constant(0.5, n) - constant(1.0 / 48.0, n) * u + constant(1.0 / 3840.0, n) * u2
    return w, [k * x for x in d]


def _axis_quaternion_ad(d: Sequence[AD], norm_upper: float) -> tuple[AD, list[AD]]:
    """Axis-angle quaternion AD using rigorous radial derivative bounds."""
    vals = [x.val for x in d]
    q = QCOMP._axis_homogeneous(vals, norm_upper)
    if q is None:
        raise RuntimeError("axis quaternion requested outside axis branch")
    wv, vv = q
    n = d[0].n
    # w_u and k_u are non-positive on r/2 in [0,3].  The interval integral
    # bounds are global on the promoted correction range and include r->0.
    w_u = Interval(-1.0 / 8.0, 0.0)
    k_u = Interval(-1.0 / 48.0, 0.0)
    # k itself is needed for dv/dd.  QCOMP returns v=k*d but division by a
    # component cell would be unsafe; use the source-wide coefficient range.
    half_lo = 0.5 * GROUP.SERIES_BRANCH_NORM
    half_hi = 0.5 * float(norm_upper)
    # sinc is positive and decreasing on [0,3].  QCOMP already validates this
    # range; its homogeneous vector coefficient is k=0.5*sinc(r/2).
    import ou3_validated_transcendentals as VT
    k_val = I(0.5) * VT.sinc_interval(Interval(half_lo, half_hi))

    du = []
    for j in range(n):
        y = I(0.0)
        for x in d:
            y = y + I(2.0) * x.val * x.der[j]
        du.append(y)
    w = AD(wv, tuple(w_u * qj for qj in du))
    out = []
    for i in range(3):
        der = []
        for j in range(n):
            der.append(k_val * d[i].der[j] + d[i].val * k_u * du[j])
        out.append(AD(vv[i], tuple(der)))
    return w, out


def deployed_quaternion_ad(d: Sequence[AD]) -> tuple[AD, list[AD], float]:
    """Generalized-Jacobian enclosure of the shipping correction quaternion."""
    if len(d) != 3:
        raise ValueError("deployed correction AD requires a three-vector")
    vals = [x.val for x in d]
    dn = _norm_upper(vals)
    if dn > QCOMP.MAX_CORRECTION_NORM:
        raise ValueError("correction cell exceeds validated deployed quaternion range")
    parts: list[tuple[AD, list[AD]]] = []
    # The component box may include small-branch points whenever its norm upper
    # reaches the origin; the algebraic polynomial enclosure is valid on the
    # strict-small intersection and using the full component box only widens it.
    parts.append(_small_quaternion_ad(d))
    if dn >= GROUP.SERIES_BRANCH_NORM:
        parts.append(_axis_quaternion_ad(d, max(dn, GROUP.SERIES_BRANCH_NORM)))
    if len(parts) == 1:
        return parts[0][0], parts[0][1], dn
    w = hull_ad(*(p[0] for p in parts))
    v = [hull_ad(*(p[1][i] for p in parts)) for i in range(3)]
    return w, v, dn


def deployed_correct_cayley(c: Sequence[AD], d: Sequence[AD]) -> list[AD]:
    """Differentiate exact homogeneous shipping-quaternion/Cayley composition."""
    if len(c) != 3 or len(d) != 3:
        raise ValueError("deployed Cayley correction requires three-vectors")
    n = c[0].n
    w, v, _dn = deployed_quaternion_ad(d)
    W = constant(2.0, n) * w - dot(v, c)
    if W.val.lo <= 0.0 <= W.val.hi:
        raise RuntimeError("AD deployed correction can reach resulting Cayley antipode")
    vxc = cross(v, c)
    V = [w * c[i] + constant(2.0, n) * v[i] + vxc[i] for i in range(3)]
    return [constant(2.0, n) * x / W for x in V]


def jacobian(y: Sequence[AD]) -> list[list[Interval]]:
    """Extract the outward interval Jacobian of an AD output vector."""
    if not y:
        return []
    n = y[0].n
    if any(x.n != n for x in y):
        raise ValueError("AD output derivative dimensions differ")
    return [[x.der[j] for j in range(n)] for x in y]


def values(y: Sequence[AD]) -> list[Interval]:
    return [x.val for x in y]


def interval_matrix_op2_upper(A: Sequence[Sequence[Interval]]) -> float:
    """Rigorous ||A||_2 upper via sqrt(||A||_1 ||A||_inf)."""
    if not A:
        return 0.0
    rows = []
    for row in A:
        s = 0.0
        for x in row:
            s = math.nextafter(s + x.abs_upper(), math.inf)
        rows.append(s)
    cols = []
    for j in range(len(A[0])):
        s = 0.0
        for i in range(len(A)):
            s = math.nextafter(s + A[i][j].abs_upper(), math.inf)
        cols.append(s)
    return math.nextafter(math.sqrt(math.nextafter(max(rows) * max(cols), math.inf)), math.inf)
