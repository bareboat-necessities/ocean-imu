#!/usr/bin/env python3
"""Exact physical-generator potential binding for the retained primitive graph.

This is a necessary correlated OUTER transition, not a numerical generator
solver and not source admission from sampled points. Every actual bounded
shaping history satisfies it. It carries the SAME generator state, primitive
moments and one-time Live potential across every word. In particular the new
physical bound cannot be attached as an independent additive S disturbance.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from typing import Sequence

import ou3_brmm_physical_wave_source as WAVE
import ou3_brmm_centered_primitive_transition as PRIMITIVE


def _state_and_potential(cert: dict, state: Sequence) -> tuple[tuple[F, ...], tuple[F, F, F]]:
    if WAVE.verify_certificate(cert):
        raise ValueError("invalid physical generator certificate")
    if cert["kind"] != "bounded_physical_shaping_potential":
        raise ValueError("this finite-state binding requires a shaping realization")
    x = tuple(WAVE.rational(v) for v in state)
    P = WAVE.matrix(cert["P"]); L = WAVE.matrix(cert["L"])
    if len(x) != len(P):
        raise ValueError("generator state dimension changed")
    energy = sum((x[i]*P[i][j]*x[j] for i in range(len(x)) for j in range(len(x))), F(0))
    if energy > F(cert["invariant_state_radius"])**2:
        raise ValueError("physical shaping state left its derived invariant")
    phi = tuple(sum((a*b for a, b in zip(row, x)), F(0)) for row in L)
    return x, phi


@dataclass(frozen=True)
class Boundary:
    source_id: str
    generator_certificate_id: str
    generator_state: tuple[F, ...]
    live_potential: tuple[F, F, F]
    primitive: PRIMITIVE.PrimitiveState


def enter_live(*, certificate: dict, source_id: str, origin_id: str,
               state_id: str, generator_state: Sequence, p: Sequence, v: Sequence) -> Boundary:
    """Center S once, without resetting physical p, v or the generator."""
    if not source_id or not origin_id or not state_id:
        raise ValueError("source, state and actual Live-origin identities required")
    x, phi = _state_and_potential(certificate, generator_state)
    return Boundary(source_id, WAVE.certificate_id(certificate), x, phi,
                    PRIMITIVE.PrimitiveState(state_id, origin_id,
                        PRIMITIVE._v3(v), PRIMITIVE._v3(p), (F(0), F(0), F(0))))


def validate_boundary(boundary: Boundary, certificate: dict) -> None:
    if boundary.generator_certificate_id != WAVE.certificate_id(certificate):
        raise ValueError("physical generator changed without a qualified switch")
    _, phi = _state_and_potential(certificate, boundary.generator_state)
    if len(boundary.live_potential) != 3:
        raise ValueError("Live potential dimension changed")
    # A constructed Boundary must also retain the consequence of the generator
    # invariant at the original Live time, not merely a self-consistent S offset.
    cap = F(certificate["potential_norm_upper_m_s"])
    if sum((WAVE.rational(x)**2 for x in boundary.live_potential), F(0)) > cap**2:
        raise ValueError("Live potential is outside the derived generator image bound")
    if tuple(phi[i]-boundary.live_potential[i] for i in range(3)) != boundary.primitive.S_L:
        raise ValueError("centered S detached from the same bounded generator potential")


def advance(boundary: Boundary, certificate: dict, moments: PRIMITIVE.MomentWitness,
            *, next_generator_state: Sequence, next_state_id: str, h: F) -> Boundary:
    """Intersect the existing J0/J1/J2 recurrence with the generator image.

    Moment feasibility, continuous generator dynamics, output p/v/a and IMU
    attachment remain in the full source provider; this function does not
    assert them from endpoints. It enforces the exact necessary equality that
    was missing from the old indefinite source graph.
    """
    validate_boundary(boundary, certificate)
    x, phi = _state_and_potential(certificate, next_generator_state)
    primitive = PRIMITIVE.advance(boundary.primitive, moments,
                                  next_primitive_id=next_state_id, h=h)
    if primitive.S_L != tuple(phi[i]-boundary.live_potential[i] for i in range(3)):
        raise ValueError("acceleration moments/position and generator potential are not the same history")
    result = Boundary(boundary.source_id, boundary.generator_certificate_id, x,
                      boundary.live_potential, primitive)
    validate_boundary(result, certificate)
    return result
