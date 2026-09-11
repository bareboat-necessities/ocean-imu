"""Finite physical prediction identities; no estimator-copy reference motion.

All functions below are real-arithmetic identities.  Rational inputs give exact
regressions.  The polynomial routines also accept a polynomial-ring scalar so
coefficient proofs need not sample states.  Neither a supplied physical step nor
its bias family name is a certificate of COMPLETE-BRMM/BIAS admission.

The physical attitude step is Q_true(t+h) conjugate(Q_true(t)), up to a
nonzero projective scale.  In particular it is NOT the shipping quaternion
routine applied to ``(-omega_hat+e_bg)*h``.  The physical/sampled-step defect is
retained explicitly by ``physical_step_from_defect``.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from typing import Sequence

from . import finite_measurement_graph as M
from . import finite_prediction_graph as P


def quaternion_product(a, b):
    """Hamilton product, polynomial over all eight scalar indeterminates."""
    if len(a) != 4 or len(b) != 4:
        raise ValueError("two quaternions are required")
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return [aw*bw-ax*bx-ay*by-az*bz,
            aw*bx+ax*bw+ay*bz-az*by,
            aw*by-ax*bz+ay*bw+az*bx,
            aw*bz+ax*by-ay*bx+az*bw]


def conjugate(q):
    if len(q) != 4:
        raise ValueError("a quaternion is required")
    return [q[0], -q[1], -q[2], -q[3]]


def physical_step_from_defect(defect, sampled_physical_step):
    """Q_phys = D_omega Q_sample.  No bound on D_omega is asserted here.

    Q_sample must use the physical sampled angular rate, not omega_hat.  The
    source relation must retain omega_hat = omega_sample + e_bg + n_g and
    D_omega = Q_phys conjugate(Q_sample) (projectively).  Treating a missing
    D_omega as identity silently removes angular integration/model forcing.
    """
    return quaternion_product(defect, sampled_physical_step)


def attitude_polynomials(c, physical_step, nominal_step):
    """Return (denominator, numerator) for denom*c_next = numerator.

    Normalization of either quaternion cancels.  This is an exact finite
    fractional map, with a generally NONZERO value at c=0.  No derivative is
    formed or multiplied.
    """
    if len(c) != 3:
        raise ValueError("three Cayley coordinates required")
    q = quaternion_product(
        quaternion_product(physical_step, [2, *c]), conjugate(nominal_step))
    return q[0], [2*x for x in q[1:]]


def attitude(c, physical_step, nominal_step):
    c = M.vec(c, 3)
    physical_step, nominal_step = M.vec(physical_step, 4), M.vec(nominal_step, 4)
    if not M.dot(physical_step, physical_step) or not M.dot(nominal_step, nominal_step):
        raise ValueError("zero step quaternion")
    denominator, numerator = attitude_polynomials(c, physical_step, nominal_step)
    if not denominator:
        raise ValueError("physical prediction crosses the Cayley pole")
    return [x/denominator for x in numerator]


def linear_prediction_polynomials(e, a0, a1, J0, J1, J2, *, h, va, pa, Sa, alpha):
    """One axis, [e_v,e_p,e_S,e_aw], with its SAME physical moments.

    J0 = integral_0^h a(s) ds, J1 = integral_0^h (h-s) a(s) ds,
    J2 = integral_0^h (h-s)^2 a(s)/2 ds.  These are not three independently
    bounded inputs.  The caller must retain their joint physical history.
    """
    ev, ep, eS, ea = e
    return [ev + va*ea + J0-va*a0,
            ep + h*ev + pa*ea + J1-pa*a0,
            eS + h*ep + F(1, 2)*h*h*ev + Sa*ea + J2-Sa*a0,
            alpha*ea + a1-alpha*a0]


@dataclass(frozen=True)
class PhysicalKinematics:
    """A physical endpoint, not a copy of the filter's OU state."""
    time: F
    q_world_to_body: tuple[F, ...]
    velocity: tuple[F, ...]
    position: tuple[F, ...]
    centered_S: tuple[F, ...]
    acceleration: tuple[F, ...]
    gyro_bias: tuple[F, ...]
    beta: tuple[F, ...]
    live_origin: F

    def __post_init__(self):
        object.__setattr__(self, "time", M.rational(self.time))
        object.__setattr__(self, "live_origin", M.rational(self.live_origin))
        for key in ("q_world_to_body", "velocity", "position", "centered_S",
                    "acceleration", "gyro_bias", "beta"):
            n = 4 if key == "q_world_to_body" else 3
            object.__setattr__(self, key, tuple(M.vec(getattr(self, key), n)))
        if not M.dot(self.q_world_to_body, self.q_world_to_body):
            raise ValueError("zero physical quaternion")
        if self.time < self.live_origin:
            raise ValueError("physical endpoint precedes its Live origin")
        if self.time == self.live_origin and any(self.centered_S):
            raise ValueError("fresh centered physical S must be zero")


@dataclass(frozen=True)
class PhysicalSegment:
    before: PhysicalKinematics
    after: PhysicalKinematics
    J0: tuple[F, ...]
    J1: tuple[F, ...]
    J2: tuple[F, ...]
    phi_true: F
    bias_driver: tuple[F, ...]

    def __post_init__(self):
        for key in ("J0", "J1", "J2", "bias_driver"):
            object.__setattr__(self, key, tuple(M.vec(getattr(self, key), 3)))
        object.__setattr__(self, "phi_true", M.rational(self.phi_true))
        if self.h <= 0 or self.before.live_origin != self.after.live_origin:
            raise ValueError("positive duration and one persistent Live origin required")
        if not 0 < self.phi_true <= 1:
            raise ValueError("physical bias factor must lie in (0,1]")
        h, b, a = self.h, self.before, self.after
        for i in range(3):
            expected = (b.velocity[i]+self.J0[i],
                        b.position[i]+h*b.velocity[i]+self.J1[i],
                        b.centered_S[i]+h*b.position[i]+h*h*b.velocity[i]/2+self.J2[i],
                        self.phi_true*b.beta[i]+self.bias_driver[i])
            actual = (a.velocity[i], a.position[i], a.centered_S[i], a.beta[i])
            if expected != actual:
                raise ValueError("physical moment/bias recurrence is not one history")

    @property
    def h(self):
        return self.after.time-self.before.time

    @property
    def physical_rotation_step(self):
        # Dividing by ||q_before||^2 is unnecessary in projective coordinates.
        return quaternion_product(self.after.q_world_to_body,
                                  conjugate(self.before.q_world_to_body))

    @property
    def q15(self):
        return (self.before.acceleration+self.after.acceleration+
                self.J0+self.J1+self.J2)


def prediction(z, segment: PhysicalSegment, *, nominal_step, axis_coefficients,
               active_bias: bool, phi_hat=None):
    """Exact conditional joint24 prediction, with no detached truth slots.

    Coefficients must be the actual shipping mean coefficients for this same
    event (including safe_phi_A_coeffs' polynomial branch).  This algebra
    routine does not certify their runtime/tuner ancestry.
    """
    z = M.vec(z, 24)
    if z[21:24] != list(segment.before.beta):
        raise ValueError("joint24 truth is not this physical predecessor")
    if not isinstance(active_bias, bool):
        raise TypeError("actual H18/A21 mode required")
    if len(axis_coefficients) != 3:
        raise ValueError("three actual axis coefficient tuples required")
    out = list(z)
    out[:3] = attitude(z[:3], segment.physical_rotation_step, nominal_step)
    for i in range(3):
        out[3+i] = z[3+i] + segment.after.gyro_bias[i]-segment.before.gyro_bias[i]
        va, pa, Sa, alpha, h = M.vec(axis_coefficients[i], 5)
        if h != segment.h:
            raise ValueError("mean coefficients and physical segment use different h")
        v = linear_prediction_polynomials(
            [z[6+i], z[9+i], z[12+i], z[15+i]],
            segment.before.acceleration[i], segment.after.acceleration[i],
            segment.J0[i], segment.J1[i], segment.J2[i],
            h=h, va=va, pa=pa, Sa=Sa, alpha=alpha)
        for j in range(4):
            out[6+3*j+i] = v[j]
    ph = M.rational(phi_hat) if active_bias else F(1)
    if not 0 < ph <= 1:
        raise ValueError("invalid active estimator bias factor")
    eb, beta = P.bias_prediction(z[18:21], z[21:24], segment.bias_driver,
                                 ph, segment.phi_true)
    out[18:21], out[21:24] = eb, beta
    if out[21:24] != list(segment.after.beta):
        raise AssertionError("shared physical driver lost")
    return out
