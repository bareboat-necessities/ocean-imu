"""Conditional finite joint24/21-covariance core, not a source-uniform word.

Every successor consumes the preceding finite state AND full shipping covariance.
No Jacobian products, independent gain snapshots or per-word S origins occur.
The algebraic composition theorem and outstanding graph obligations are in
``docs/ou3-alt-finite-core-composition.md``. Rational executions test identities;
they do not admit BRMM histories, frontend states or floating-point branches.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYSICAL
from tools.stability.ou3_alt_contraction.finite_prediction_deployed_step import step_quaternion_polynomial


def quaternion(q):
    q = P.vec(q, 4)
    if not any(q):
        raise ValueError('zero homogeneous quaternion')
    return tuple(q)


def rotation(q):
    """Rotation of a nonzero homogeneous quaternion, without a square root."""
    w, *v = quaternion(q)
    n2 = w*w + M.dot(v, v)
    V = M.skew(v)
    return M.plus(M.eye(3), M.scaled(M.plus(M.scaled(V, w), M.mm(V, V)), 2/n2))


def cayley(q):
    q = quaternion(q)
    if q[0] == 0:
        raise ValueError('finite physical Cayley chart pole')
    return tuple(2*x/q[0] for x in q[1:])


@dataclass(frozen=True)
class Reference(PHYSICAL.PhysicalKinematics):
    """Existing physical endpoint with persistent, non-qualifying source labels.

    Inherit the physical time/Live-origin/kinematic contract instead of creating
    a competing source model. Labels never establish BRMM/BIAS admission.
    """
    history_id: str
    bias_root: str
    bias_family: str

    def __post_init__(self):
        super().__post_init__()
        if any(not isinstance(x, str) or not x for x in (self.history_id, self.bias_root)):
            raise ValueError('persistent physical/source ancestry required')
        if self.bias_family not in ('BIAS0', 'BIAS1', 'BIAS2'):
            raise ValueError('one declared BIAS family required')


@dataclass(frozen=True)
class State:
    mode: str
    z: tuple
    covariance: tuple
    q_hat: tuple
    reference: Reference

    def __post_init__(self):
        if self.mode not in ('H', 'A'):
            raise ValueError('mode must be H or A')
        if not isinstance(self.reference, Reference):
            raise TypeError('physical reference required')
        z = tuple(P.vec(self.z, 24))
        cov = M.mat(self.covariance, 21, 21)
        if cov != M.transpose(cov):
            raise ValueError('full symmetric 21-state covariance required, also in H')
        qh = quaternion(self.q_hat)
        if z[:3] != cayley(P.quat_mul(self.reference.q_world_to_body, P.quat_conj(qh))):
            raise ValueError('finite attitude error is not the physical/nominal relative rotation')
        if z[21:24] != self.reference.beta:
            raise ValueError('joint24 beta must be the SAME physical bias')
        object.__setattr__(self, 'z', z)
        object.__setattr__(self, 'covariance', tuple(tuple(r) for r in cov))
        object.__setattr__(self, 'q_hat', qh)


def prediction(state, segment, *, gyro_body, axis_coefficients, F21, Q21, phi_hat=None):
    """Existing finite PhysicalSegment predictor plus conditional covariance.

    F21/Q21 must still be attached to the literal shipping coefficient/repair
    graph by the source provider. Supplying matrices does NOT close that graph.
    No floor, due-S event, frontend update or mode release is implicit here.
    """
    if not isinstance(segment, PHYSICAL.PhysicalSegment):
        raise TypeError('existing same-history PhysicalSegment required')
    before, after = segment.before, segment.after
    if before != state.reference or not isinstance(after, Reference):
        raise ValueError('physical segment must extend this exact predecessor')
    for name in ('history_id', 'bias_root', 'bias_family'):
        if getattr(after, name) != getattr(before, name):
            raise ValueError(name + ' cannot restart at a word boundary')
    # PhysicalSegment already checks duration, one-time origin, q15 primitive
    # recurrence and the single shared bias driver; do not duplicate its model.
    g = P.vec(gyro_body, 3)
    omega_hat = [g[i] - (before.gyro_bias[i]-state.z[3+i]) for i in range(3)]
    qn = step_quaternion_polynomial([-segment.h*x for x in omega_hat])
    z = PHYSICAL.prediction(state.z, segment, nominal_step=qn,
                            axis_coefficients=axis_coefficients,
                            active_bias=state.mode == 'A', phi_hat=phi_hat)
    A, Q = M.mat(F21, 21, 21), M.mat(Q21, 21, 21)
    if Q != M.transpose(Q):
        raise ValueError('symmetric real process covariance required')
    cov = M.plus(M.mm(M.mm(A, state.covariance), M.transpose(A)), Q)
    return State(state.mode, tuple(z), cov, P.quat_mul(qn, state.q_hat), after)


def solve3(S, rhs):
    """Exact evaluation of S*x=rhs; NOT an emulation of Eigen's LDLT decisions."""
    S, rhs = M.mat(S, 3, 3), P.vec(rhs, 3)
    a = [row+[x] for row, x in zip(S, rhs)]
    for j in range(3):
        pivot = next((i for i in range(j, 3) if a[i][j]), None)
        if pivot is None:
            raise ValueError('singular innovation: no accepted nonsingular graph witness')
        a[j], a[pivot] = a[pivot], a[j]
        divisor = a[j][j]
        a[j] = [x/divisor for x in a[j]]
        for i in range(3):
            if i != j:
                multiplier = a[i][j]
                a[i] = [x-multiplier*y for x, y in zip(a[i], a[j])]
    return [row[3] for row in a]


def covariance_reset(Pj, d):
    """Literal G=I+[d]x/2 covariance law, NOT the finite mean reset map.

    Only attitude rows/columns change; retain all 21 covariance coordinates.
    """
    Pj, d = M.mat(Pj, 21, 21), P.vec(d, 3)
    if Pj != M.transpose(Pj):
        raise ValueError('symmetric real covariance required for reset')
    G = M.plus(M.eye(3), M.scaled(M.skew(d), F(1, 2)))
    out = [r[:] for r in Pj]
    aa = M.mm(M.mm(G, [r[:3] for r in Pj[:3]]), M.transpose(G))
    cross = M.mm(G, [r[3:] for r in Pj[:3]])
    for i in range(3):
        out[i][:3] = aa[i]
        out[i][3:] = cross[i]
        for j in range(3, 21):
            out[j][i] = cross[i][j-3]
    return out


@dataclass(frozen=True)
class Accepted:
    state: State
    innovation: tuple
    numerator: tuple
    gain: tuple
    innovation_solution: tuple
    correction: tuple
    cayley_denominator: F


def measurement(state, kind, *, R, observed=None, magnetic_reference=None,
                gravity=F(196133, 20000), innovation_shift=0, alpha=1,
                radius=F(2, 5)):
    """Execute one CONDITIONAL accepted, finite, polynomial-quaternion graph.

    Sensors are in the common deheeled frame with modeled temperature removed.
    Shift is the SAME safe-LDLT retry shift in the gain AND Joseph graph. Its
    literal retry predicate/value, finite sensor guards and applied R still
    require the provider; no acceptance is inferred from this solver's success.
    Outside the polynomial injection branch this rational evaluator raises:
    the all-w,k finite identity remains available in finite_measurement_graph.
    """
    ref, z = state.reference, state.z
    Rh = rotation(state.q_hat)
    if kind == 'accelerometer':
        g = P.rational(gravity)
        aw_hat = [ref.acceleration[i]-z[15+i] for i in range(3)]
        fhat = M.mv(Rh, [aw_hat[0], aw_hat[1], aw_hat[2]-g])
        y = P.vec(observed, 3)
        residual = [y[i]-fhat[i]-(ref.beta[i]-z[18+i]) for i in range(3)]
        H, Hb = M.residual_factors(kind, z[:3], f_hat=fhat, R_hat=Rh)
        physical = M.mv(rotation(ref.q_world_to_body), [ref.acceleration[0], ref.acceleration[1], ref.acceleration[2]-g])
        nu = [y[i]-physical[i]-ref.beta[i] for i in range(3)]
    elif kind == 'magnetometer':
        mr, y = P.vec(magnetic_reference, 3), P.vec(observed, 3)
        mhat = M.mv(Rh, mr)
        residual = [y[i]-mhat[i] for i in range(3)]
        H, Hb = M.residual_factors(kind, z[:3], m_hat=mhat)
        physical = M.mv(rotation(ref.q_world_to_body), mr)
        nu = [y[i]-physical[i] for i in range(3)]
    elif kind == 'S_zero':
        if observed is not None or magnetic_reference is not None:
            raise ValueError('S source comes only from the persistent physical primitive')
        residual = [z[12+i]-ref.centered_S[i] for i in range(3)]
        H, Hb = M.residual_factors(kind, z[:3])
        nu = [-x for x in ref.centered_S]
    else:
        raise ValueError('unsupported accepted measurement')
    if residual != [x+y for x, y in zip(M.mv(Hb, z), nu)]:
        raise AssertionError('physical finite residual identity failed')
    N, S = M.measurement_operands(state.covariance, R, H, active_bias=state.mode == 'A')
    shift = P.rational(innovation_shift)
    if shift < 0:
        raise ValueError('safe-LDLT diagonal shift cannot be negative')
    S = M.plus(S, M.scaled(M.eye(3), shift))
    q = solve3(S, residual)
    K = [solve3(S, row) for row in N]
    d = M.mv(N, q)
    if d != M.mv(K, residual) or M.mm(K, S) != N:
        raise AssertionError('inverse-free gain/correction identities failed')
    w, k = M.small_quaternion_parameters(d[:3])
    epre = [z[18+i]-d[18+i] for i in range(3)]
    if not M.projection_graph_holds(epre, ref.beta, alpha, radius):
        raise ValueError('radial projection factor detached from the SAME physical beta')
    graph = M.descriptor(z[:3], d[:3], w, k, N, S, Hb, alpha=alpha)
    chi = list(z)+q+nu+[F(1)]
    if any(M.mv(graph.equality, chi)):
        raise AssertionError('finite descriptor equality failed')
    znext = M.mv(graph.after, chi)
    # Reuse the checked exact rank-three lemma from PR #523, including the
    # SAME repaired innovation. Nonzero solve defects are not discarded.
    Pj = M.solved_joseph_covariance(state.covariance, N, S, K)
    cov = covariance_reset(Pj, d[:3])
    qhat = P.quat_mul([w, *[k*x for x in d[:3]]], state.q_hat)
    nxt = State(state.mode, tuple(znext), cov, qhat, ref)
    return Accepted(nxt, tuple(map(tuple, S)), tuple(map(tuple, N)), tuple(map(tuple, K)),
                    tuple(q), tuple(d), graph.denominator)


def compose_accepted(state, events):
    """All finite accepted-event tuples, including empty ones; retain every prefix.

    This is conditional composition, not admission of the event schedule. No
    replay/root fitting, copied P/K, new physical reference or metric is used.
    """
    prefixes = [state]
    for event in events:
        state = measurement(state, **event).state
        prefixes.append(state)
    return tuple(prefixes)
