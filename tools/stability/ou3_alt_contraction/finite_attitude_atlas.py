"""Exact four-chart attitude graph for the finite joint24 runtime.

The chart with largest homogeneous quaternion component is always defined.
Chart changes are proof coordinates only: they change neither the physical
reference, nominal quaternion, nor the shipping error covariance. See
docs/ou3-alt-attitude-atlas.md for the universal covering and transport proof.
No chart origin, other than chart zero's, denotes zero attitude error.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F

from . import finite_measurement_graph as M
from . import finite_prediction_graph as P


def slots(chart):
    if type(chart) is not int or chart not in range(4):
        raise ValueError('attitude chart must be an integer in [0,3]')
    return tuple(i for i in range(4) if i != chart)


def homogeneous(c, chart=0):
    indices = slots(chart)
    q = [F(0)]*4
    q[chart] = F(2)
    for i, x in zip(indices, M.vec(c, 3)):
        q[i] = x
    return tuple(q)


@dataclass(frozen=True)
class Point:
    chart: int
    coordinates: tuple

    def __post_init__(self):
        slots(self.chart)
        object.__setattr__(self, 'coordinates', tuple(M.vec(self.coordinates, 3)))

    @property
    def quaternion(self):
        return homogeneous(self.coordinates, self.chart)


def encode(q, *, chart=None):
    q = M.vec(q, 4)
    if not any(q):
        raise ValueError('zero homogeneous quaternion')
    if chart is None:
        # Stable tie rule; invariant under every nonzero projective scaling.
        chart = max(range(4), key=lambda i: abs(q[i]))
    indices = slots(chart)
    if q[chart] == 0:
        raise ValueError('requested attitude chart has zero denominator')
    return Point(chart, tuple(2*q[i]/q[chart] for i in indices))


def rotation_numerator(q):
    w, x, y, z = q
    return [[w*w+x*x-y*y-z*z, 2*(x*y-w*z), 2*(x*z+w*y)],
            [2*(x*y+w*z), w*w-x*x+y*y-z*z, 2*(y*z-w*x)],
            [2*(x*z-w*y), 2*(y*z+w*x), w*w-x*x-y*y+z*z]]


def rotation(point):
    q = point.quaternion
    return M.scaled(rotation_numerator(q), 1/M.dot(q, q))


@dataclass(frozen=True)
class Transport:
    before: Point
    after: Point
    denominator: F
    # Exact c_after = affine * (c_before,1), with the SAME hard quaternion graph.
    affine: tuple


def transport(point, *, left=(1,0,0,0), right=(1,0,0,0), chart_after=None):
    """E_next = L E R; ``right`` is already conjugated if required.

    No differentiability, small angle, independent reset state or denominator
    bound is assumed. The maximum-component rule supplies the denominator.
    """
    if not isinstance(point, Point):
        raise TypeError('attitude atlas Point required')
    left, right = M.vec(left, 4), M.vec(right, 4)
    if not any(left) or not any(right):
        raise ValueError('nonzero left and right quaternion factors required')
    columns = [P.quat_mul(P.quat_mul(left, [F(i == j) for i in range(4)]), right)
               for j in range(4)]
    A = M.transpose(columns)
    qp = M.mv(A, point.quaternion)
    after = encode(qp, chart=chart_after)
    den = qp[after.chart]
    affine = tuple(tuple([2*A[i][j]/den for j in slots(point.chart)] +
                         [4*A[i][point.chart]/den]) for i in slots(after.chart))
    if tuple(M.mv(affine, (*point.coordinates, F(1)))) != after.coordinates:
        raise AssertionError('finite atlas transport identity failed')
    return Transport(point, after, den, affine)


def residual_factors(kind, point, *, f_hat=None, R_hat=None, m_hat=None):
    """Actual shipping H, finite Hbar and SAME-state attitude offset.

    r = Hbar*z + offset(c,chart) + nu. For charts 1..3 the offset is
    (R(q(c))-I)*f_hat, a hard rational state function, NEVER a free supply.
    In chart zero the existing exact Cayley secant has zero offset.
    """
    if point.chart == 0:
        H, Hb = M.residual_factors(kind, point.coordinates,
                                  f_hat=f_hat, R_hat=R_hat, m_hat=m_hat)
        return H, Hb, (F(0),)*3
    H, Hb = M.zeros(3, 21), M.zeros(3, 24)
    offset = (F(0),)*3
    if kind in ('accelerometer', 'magnetometer'):
        f = M.vec(f_hat if kind == 'accelerometer' else m_hat, 3)
        E = rotation(point)
        offset = tuple(M.mv(M.plus(E, M.eye(3), -1), f))
        Htheta = M.scaled(M.skew(f), -1)
        for i in range(3):
            H[i][:3] = Htheta[i]
        if kind == 'accelerometer':
            R = M.mat(R_hat, 3, 3)
            ER = M.mm(E, R)
            for i in range(3):
                H[i][15:18], Hb[i][15:18] = R[i], ER[i]
                H[i][18+i] = Hb[i][18+i] = F(1)
    elif kind == 'S_zero':
        for i in range(3):
            H[i][12+i] = Hb[i][12+i] = F(1)
    else:
        raise ValueError('unsupported measurement kind')
    return H, Hb, offset


@dataclass(frozen=True)
class Descriptor:
    equality: list
    before: list
    after: list
    transport: Transport

    @property
    def denominator(self):
        return self.transport.denominator


def descriptor(point, d, w, k, N, S, Hbar, offset, *, alpha):
    """Inverse-free finite measurement graph, including exact chart transport.

    chi=(z24,solve3,nu3,1), with c=z[:3], d=N_theta*solve and the SAME
    computed offset. Projection retains physical beta in the final three slots.
    """
    d, offset = M.vec(d, 3), M.vec(offset, 3)
    w, k = M.rational(w), M.rational(k)
    N, S, Hb = M.mat(N, 21, 3), M.mat(S, 3, 3), M.mat(Hbar, 3, 24)
    tr = transport(point, right=P.quat_conj([w, *[k*x for x in d]]))
    before = [M.eye(24)[i]+[F(0)]*7 for i in range(24)]
    after = [row[:] for row in before]
    for i, row in enumerate(tr.affine):
        after[i] = list(row[:3])+[F(0)]*27+[row[3]]
    for i in range(3, 21):
        after[i][24:27] = [-x for x in N[i]]
    after = M.mm(M.projection_matrix(alpha), after)
    eq = [[-x for x in Hb[i]]+S[i]+[-F(i == j) for j in range(3)]+[-offset[i]]
          for i in range(3)]
    for i in range(3):
        eq.append([F(0)]*24+N[i]+[F(0)]*3+[-d[i]])
        eq.append([F(i == j) for j in range(24)]+[F(0)]*6+[-point.coordinates[i]])
    return Descriptor(eq, before, after, tr)


def readiness():
    return {
        'all_nonzero_relative_quaternions_covered': True,
        'selected_unit_quaternion_component_abs_lower': F(1,2),
        'selected_coordinate_abs_upper': F(2),
        'selected_coordinate_squared_norm_upper': F(12),
        'finite_left_right_and_chart_transport_identity': True,
        'inverse_free_measurement_descriptor_with_chart_transport': True,
        'chart_origin_is_zero_attitude_error_for_every_chart': False,
        'compatible_coercive_storage_closed': False,
        'source_uniform_deployment_word_closed': False,
        'storage_search_allowed': False,
    }
