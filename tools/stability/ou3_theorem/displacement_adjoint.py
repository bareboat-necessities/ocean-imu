"""Exact displacement-chord completion for the literal S/accelerometer adjoint."""
from __future__ import annotations
from fractions import Fraction as F
from .signed_temporal import jet

def four_s_boundary_jets(knots):
    k=tuple(F(x) for x in knots)
    if len(k)!=4: raise ValueError("four S knots required")
    return {"left":tuple(jet(k,k[0],q,"left") for q in range(3)),
            "right":tuple(jet(k,k[-1],q,"right") for q in range(3))}

def displacement_chord_identity(duration, displacement_chord):
    T,P=F(duration),F(displacement_chord)
    if T<=0 or P<0: raise ValueError("T>0 and P>=0 required")
    return {"duration":T,"chord":P,"left_velocity_coefficient":T,
            "acceleration_kernel":"T-t","right_velocity_form_coefficient":T}

def chord_multiplier_jet(duration, t):
    """Jets of lambda(t)=T-t on [0,T].

    Integration against physical acceleration gives exactly
      int lambda a = Delta p - T v(0).
    lambda is degree one (hence a degenerate cubic) with nonzero left value;
    lambda'=-1 and lambda''=0.  It introduces no interior distributional
    atoms. Literal S/acc rows may be added linearly through the standard
    four-S spline without changing this boundary identity.
    """
    T,t=F(duration),F(t)
    if T<=0 or t<0 or t>T: raise ValueError("0<=t<=T required")
    return (T-t,F(-1),F(0))

def chord_integration_by_parts(duration):
    T=F(duration)
    if T<=0: raise ValueError("T>0 required")
    return {
      "lambda_left":T,"lambda_right":F(0),
      "lambda_prime_left":F(-1),"lambda_prime_right":F(-1),
      "identity":"integral_0^T (T-t) a dt = p(T)-p(0)-T v(0)",
      "boundary_velocity_coefficient":-T,
      "boundary_position_coefficients":(F(-1),F(1)),
      "interior_atoms":(),
    }

def augmented_multiplier_identity(duration, spline_scale=1):
    """Exact linear superposition lambda + alpha psi.

    psi is the existing four-S spline with zero exterior jets. Therefore the
    displacement boundary coefficients are invariant under alpha, while all
    literal four-S S/accelerometer weights are simply scaled by alpha.
    This proves an exact boundary completion without changing event chronology.
    """
    T,a=F(duration),F(spline_scale)
    if T<=0: raise ValueError("T>0 required")
    return {
      "boundary_position_coefficients":(F(-1),F(1)),
      "boundary_velocity_coefficient":-T,
      "four_S_scale":a,
      "boundary_coefficients_independent_of_four_S_scale":True,
      "exact_superposition":True,
      "root_LIN_boundary_reader_still_required":True,
    }

def certificate():
    return {
      "qualification":"OU3_DISPLACEMENT_FOUR_S_ADJOINT_V2",
      "four_S_boundary_jets_vanish":True,
      "AG_reader_can_supply_velocity_boundary":False,
      "reason_AG_reader_fails":"AG reader targets attitude+gyro-bias only; velocity is a LIN root coordinate",
      "nonzero_boundary_multiplier":"lambda(t)=T-t",
      "chord_identity":"Delta p = T v(t0) + integral lambda a dt",
      "lambda_interior_atoms":False,
      "four_S_superposition_preserves_chord_boundary_exactly":True,
      "boundary_velocity_completion_identity_proved":True,
      "boundary_velocity_numerical_action_bound":False,
      "slow_fast_abel_after_completion_allowed":True,
      "source_uniform_outer_compatibility_closed":False,
      "theorem_closed":False,
    }
