"""Conditional magnetic capture algebra, with an explicit full-frame premise.

Let A_i = R_tilt_hat_i R_true_WB_i and G_L = Rz(-psi_true_at_handoff).
If every accepted sample satisfies ||A_i-G_L||_2 <= 2*s, then its rotated
field is G_L B + e_i with ||e_i|| <= 2*Bmax*s + bHImax + nmax.
The arithmetic mean satisfies the same bound, without a statistical gain.

A gravity-direction/tilt bound does NOT supply this full-frame bound: stripping
observer yaw leaves the physical vessel heading in A_i. Heading excursion over
acquisition and between the accepted samples and handoff must be retained.
In particular the declared 0.02-rad gravity bound cannot be substituted for s.

The conditional numerical example s=0.01 and handoff tilt <=0.02 gives E<=8.5,
sin(yaw_error)<=17/30, yaw error<0.61 and total angle<0.63<pi/4. Its arithmetic
is valid; the source-uniform frame and handoff premises are unproved.

For the finite atlas word, no small-frame premise is necessary: proper
rotations satisfy ||A-G||_2<=2, and derive(1) supplies the unrestricted mean
image bound. This does not establish nonvanishing north or a storage basin.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_mag_source_qualification as Q


def R(x): return M.rational(x)
DECLARED_STARTUP_TILT_RAD_MAX=F(1,50) # 0.02 rad
YAW_CERT_RAD=F(61,100)
FULL_ATTITUDE_CERT_RAD=F(63,100)

@dataclass(frozen=True)
class Bound:
    half_tilt_sin_max:F
    earth_rotation_error_max:F
    mean_perturbation_max:F
    horizontal_true_min:F
    sin_yaw_error_max:F
    capture_nonzero:bool
    def __post_init__(self):
        for n in ('half_tilt_sin_max','earth_rotation_error_max','mean_perturbation_max',
                  'horizontal_true_min','sin_yaw_error_max'):
            object.__setattr__(self,n,R(getattr(self,n)))

@dataclass(frozen=True)
class FreshAttitudeCertificate:
    magnetic:Bound
    yaw_rad_upper:F
    tilt_rad_upper:F
    full_attitude_rad_upper:F
    sin_yaw_test_lower:F
    below_pi_over_4:bool


def derive(half_tilt_sin_max,envelope:Q.Envelope|None=None):
    e=Q.default_envelope() if envelope is None else envelope
    if not isinstance(e,Q.Envelope) or e.assumption_id!=Q.ASSUMPTION_ID:
        raise ValueError('MAG-BMM150-DET-v1 envelope required')
    s=R(half_tilt_sin_max)
    if s<0 or s>1: raise ValueError('sin(delta_tilt/2) bound must lie in [0,1]')
    rot=2*e.world_field_norm_max*s
    total=e.hard_iron_norm_max+e.residual_norm_max+rot
    H=e.world_field_horizontal_min
    capture=total<H
    ratio=total/H if H>0 else F(10**9)
    return Bound(s,rot,total,H,ratio,capture)


def require_capture(bound:Bound):
    if not isinstance(bound,Bound): raise TypeError('startup magnetic capture bound required')
    if not bound.capture_nonzero or not (bound.sin_yaw_error_max<1):
        raise ValueError('deterministic startup magnetic perturbation can erase horizontal north')
    return True


def default_tilt_margin_half_sin():
    e=Q.default_envelope()
    numerator=e.world_field_horizontal_min-e.hard_iron_norm_max-e.residual_norm_max
    denominator=2*e.world_field_norm_max
    if numerator<=0: return F(0)
    return numerator/denominator


def unrestricted_frame_bound(envelope:Q.Envelope|None=None):
    """For every proper A,G and v, ||(A-G)v||<=||Av||+||Gv||=2||v||.

    Apply this to each SAME physical field sample and then the positive-weight
    mean. No heading excursion, handoff lag or observer-accuracy cap is added.
    Exact finite transport retains the rotations and discrepancy themselves;
    this image bound never replaces attitude error by independent forcing.
    """
    return derive(F(1),envelope)


def conditional_fresh_attitude_certificate(*, half_frame_sin_max, handoff_tilt_rad_max):
    """Prove only the implication from explicitly supplied full-frame bounds."""
    b=derive(half_frame_sin_max); require_capture(b)
    tilt=R(handoff_tilt_rad_max)
    if tilt<0: raise ValueError('nonnegative handoff tilt bound required')
    x=YAW_CERT_RAD
    # Alternating Taylor truncation: sin(x) >= x-x^3/6 for 0<=x<=1.
    sin_lower=x-x*x*x/F(6)
    if not sin_lower>b.sin_yaw_error_max:
        raise AssertionError('0.61-rad yaw certificate does not dominate magnetic sine bound')
    full=x+tilt
    # pi>3 is sufficient: pi/4>3/4, and 63/100<3/4.
    below=full<F(3,4)
    if not below: raise AssertionError('startup attitude certificate does not fit pi/4')
    return FreshAttitudeCertificate(b,x,tilt,full,sin_lower,True)


def readiness():
    c=conditional_fresh_attitude_certificate(half_frame_sin_max=F(1,100),
                                             handoff_tilt_rad_max=F(1,50))
    coarse=unrestricted_frame_bound()
    return {
      'deterministic_average_does_not_claim_sqrtN_improvement':True,
      'earth_field_rotation_chord_bound_materialized':True,
      'hard_iron_and_residual_combined_without_independence_assumption':True,
      'horizontal_nonvanishing_capture_condition_materialized':True,
      'yaw_sine_error_supply_bound_materialized':True,
      'current_envelope_half_tilt_sin_margin_is_4_over_75':default_tilt_margin_half_sin()==F(4,75),
      'conditional_full_frame_0p02_and_handoff_tilt_0p02_imply_pi_over_4_entry':c.below_pi_over_4,
      'unrestricted_full_frame_chord_bound_closed':True,
      'unrestricted_frame_operator_norm_upper':F(2),
      'unrestricted_mean_perturbation_upper_uT':coarse.mean_perturbation_max,
      'small_frame_accuracy_required_by_finite_atlas_word':False,
      'full_accumulation_to_handoff_frame_bound_source_qualified':False,
      'declared_startup_tilt_0p02_rad_attached_via_sin_x_le_x':False,
      'real_arithmetic_yaw_lt_0p61_rad_certified':False,
      'real_arithmetic_full_SO3_lt_pi_over_4_certified':False,
      'atan2_AngleAxis_binary32_correspondence_attached':False,
      'startup_capture_closed':False,
      'ALT_STARTUP_PASS':False,
    }
