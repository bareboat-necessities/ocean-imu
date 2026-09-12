"""Exact necessary BRMM source constraints carried by finite physical prefixes.

These predicates are an OUTER inclusion, never a source-admission oracle.
In particular, they do not certify a bounded wave-generator realization,
Q/O recurrence, a shared angular-rate function, or a BIAS generating function.

For acceleration moments J=(J0,J1,J2), the unnormalized Gramian is G(h).
Concatenation shifts the first segment's kernels by the second duration:

    J_ab = T(b) J_a + J_b,
    G(a+b) = T(b) G(a) T(b)^T + G(b).

The projection identity in docs/ou3-alt-source-continuation.md proves
E(J_ab,a+b) <= E(J_a,a)+E(J_b,b), where E=h*x^T G(1)^-1*x,
x=(J0/h,J1/h^2,J2/h^3), summed over all THREE spatial axes. The prefix
budget below is derived from the actual segment moments, not a free supply.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
import json
from pathlib import Path

import ou3_brmm_acceleration_moment_iqc as IQC
from . import finite_measurement_graph as M
from . import finite_physical_prediction as PHYS

R = M.rational
ZERO = (F(0), F(0), F(0))
DOMAIN = Path(__file__).resolve().parents[3] / 'tools/stability/ou3_proof_operating_domain.json'


@dataclass(frozen=True)
class Bounds:
    acceleration: F
    velocity: F
    position: F
    centered_S: F
    angular_rate_upper: F

    def __post_init__(self):
        for name in self.__dataclass_fields__:
            value = R(getattr(self, name))
            if value <= 0:
                raise ValueError('positive physical source bound required: '+name)
            object.__setattr__(self, name, value)


def declared_bounds():
    d = json.loads(DOMAIN.read_text(), parse_float=F)
    e = d['complete_brmm_physical_envelope']
    # pi < 22/7 follows from integral_0^1 x^4(1-x)^4/(1+x^2) dx = 22/7-pi > 0.
    # This is an outward rate enclosure, not a tighter/redefined physical cap.
    rate = R(e['body_rate_norm_upper_deg_s'])*F(22, 7)/180
    return Bounds(e['wave_acceleration_norm_upper_mps2'], e['wave_velocity_norm_upper_mps'],
                  e['wave_position_norm_upper_m'], e['centered_primitive_D_S_upper_m_s'], rate)


BOUNDS = declared_bounds()


def gram(h):
    """Polynomial in h; also accepts an exact polynomial-ring scalar."""
    h2 = h*h; h3 = h2*h; h4 = h3*h; h5 = h4*h
    return ((h, F(1, 2)*h2, F(1, 6)*h3),
            (F(1, 2)*h2, F(1, 3)*h3, F(1, 8)*h4),
            (F(1, 6)*h3, F(1, 8)*h4, F(1, 20)*h5))


def shift(h):
    return ((1, 0, 0), (h, 1, 0), (F(1, 2)*h*h, h, 1))


def concatenate(Ja, Jb, second_duration):
    """Exact kernel translation, not an entrywise interval hull."""
    b = second_duration
    return (tuple(Ja[0][i]+Jb[0][i] for i in range(3)),
            tuple(Ja[1][i]+b*Ja[0][i]+Jb[1][i] for i in range(3)),
            tuple(Ja[2][i]+b*Ja[1][i]+F(1, 2)*b*b*Ja[0][i]+Jb[2][i]
                  for i in range(3)))


def energy(duration, moments):
    """Minimum integral ||a||^2 compatible with all nine acceleration moments."""
    h = R(duration)
    if h <= 0 or len(moments) != 3:
        raise ValueError('positive duration and J0/J1/J2 required')
    J = tuple(tuple(M.vec(j, 3)) for j in moments)
    axes = tuple((J[0][i]/h, J[1][i]/(h*h), J[2][i]/(h*h*h)) for i in range(3))
    return h*IQC.joint_normalized_moment_energy(axes)


def norm2(v):
    return sum((x*x for x in v), F(0))


def check_endpoint(endpoint: PHYS.PhysicalKinematics):
    if not isinstance(endpoint, PHYS.PhysicalKinematics):
        raise TypeError('finite physical endpoint required')
    for name, cap in (('acceleration', BOUNDS.acceleration), ('velocity', BOUNDS.velocity),
                      ('position', BOUNDS.position), ('centered_S', BOUNDS.centered_S)):
        if norm2(getattr(endpoint, name)) > cap*cap:
            raise ValueError('BRMM physical '+name+' vector cap exceeded')


def check_segment(segment: PHYS.PhysicalSegment):
    """Necessary vector/moment/rotation constraints on the SAME physical step."""
    if not isinstance(segment, PHYS.PhysicalSegment):
        raise TypeError('finite physical segment required')
    check_endpoint(segment.before); check_endpoint(segment.after)
    J = (segment.J0, segment.J1, segment.J2)
    E = energy(segment.h, J)
    if E > segment.h*BOUNDS.acceleration*BOUNDS.acceleration:
        raise ValueError('BRMM coupled three-axis acceleration moment IQC violated')
    q0, q1 = segment.before.q_world_to_body, segment.after.q_world_to_body
    dot = sum((a*b for a, b in zip(q0, q1)), F(0))
    # 4 sin^2(theta/2) <= theta^2 <= (Omega*h)^2. Both quaternion signs
    # and all nonzero projective scalings represent the same physical rotation.
    chord2 = 4*(1-dot*dot/(norm2(q0)*norm2(q1)))
    if chord2 > (BOUNDS.angular_rate_upper*segment.h)**2:
        raise ValueError('BRMM consecutive physical rotation chord cap exceeded')
    return E


@dataclass(frozen=True)
class Prefix:
    duration: F = F(0)
    J0: tuple = ZERO
    J1: tuple = ZERO
    J2: tuple = ZERO
    segment_energy_sum: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, 'duration', R(self.duration))
        object.__setattr__(self, 'segment_energy_sum', R(self.segment_energy_sum))
        for name in ('J0', 'J1', 'J2'):
            object.__setattr__(self, name, tuple(M.vec(getattr(self, name), 3)))
        if self.duration < 0 or self.segment_energy_sum < 0:
            raise ValueError('nonnegative prefix duration and derived energy required')
        if not self.duration:
            if self.segment_energy_sum or any(any(j) for j in self.moments):
                raise ValueError('empty prefix has zero moments and zero energy')
        elif not self.minimum_energy <= self.segment_energy_sum <= self.duration*BOUNDS.acceleration**2:
            raise ValueError('correlated prefix acceleration energy budget violated')

    @property
    def moments(self):
        return self.J0, self.J1, self.J2

    @property
    def minimum_energy(self):
        return energy(self.duration, self.moments) if self.duration else F(0)


def append(prefix: Prefix, segment: PHYS.PhysicalSegment):
    if not isinstance(prefix, Prefix):
        raise TypeError('carried physical moment prefix required')
    E = check_segment(segment)
    J = concatenate(prefix.moments, (segment.J0, segment.J1, segment.J2), segment.h)
    return Prefix(prefix.duration+segment.h, *J, prefix.segment_energy_sum+E)
