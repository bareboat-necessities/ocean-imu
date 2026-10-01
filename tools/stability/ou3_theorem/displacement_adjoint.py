"""Exact displacement-chord contraction against a four-S spline.

The result is intentionally fail-closed: the standard four-S spline has zero
exterior jets, so it cannot by itself represent a displacement chord, whose
integration-by-parts functional has a boundary-velocity term.  This module
states the missing boundary row exactly instead of hiding it in a norm box.
"""
from __future__ import annotations
from fractions import Fraction as F
from .signed_temporal import atoms, jet

def four_s_boundary_jets(knots):
    k=tuple(F(x) for x in knots)
    if len(k)!=4: raise ValueError("four S knots required")
    return {
      "left":tuple(jet(k,k[0],q,"left") for q in range(3)),
      "right":tuple(jet(k,k[-1],q,"right") for q in range(3)),
    }

def displacement_chord_identity(duration, displacement_chord):
    """q.(p1-p0)=T q.v0 + int_0^T (T-t) q.a dt."""
    T,P=F(duration),F(displacement_chord)
    if T<=0 or P<0: raise ValueError("T>0 and P>=0 required")
    return {"duration":T,"chord":P,
            "left_velocity_coefficient":T,
            "acceleration_kernel":"T-t",
            "right_velocity_form_coefficient":T}

def spline_can_annihilate_displacement_boundary(knots):
    jets=four_s_boundary_jets(knots)
    # psi,psi',psi'' vanish outside by construction.  The displacement chord
    # needs a nonzero coefficient on v at an endpoint.  No scalar multiple of
    # this spline can supply it.
    return jets["left"]==(F(0),F(0),F(0)) and jets["right"]==(F(0),F(0),F(0))

def certificate():
    knots=(F(0),F(1),F(2),F(3))
    jets=four_s_boundary_jets(knots)
    return {
      "qualification":"OU3_DISPLACEMENT_FOUR_S_ADJOINT_V1",
      "four_S_exterior_jets":{k:list(map(str,v)) for k,v in jets.items()},
      "four_S_boundary_jets_vanish":True,
      "displacement_chord_requires_boundary_velocity":True,
      "displacement_chord_in_standard_four_S_row_space":False,
      "exact_missing_functional":"T*q^T*v(t0) (or T*q^T*v(t1) in the reversed identity)",
      "consequence":"four-S+accelerometer interior adjoint alone cannot consume Delta_p; retain a boundary velocity row or augment the proof multiplier",
      "taking_norm_before_boundary_completion_invalid":True,
      "source_uniform_outer_compatibility_closed":False,
      "theorem_closed":False,
    }
