"""Deterministic startup north/yaw capture bound under MAG-BMM150-DET-v1.

No stochastic averaging gain is claimed. If each accepted de-yawed magnetic
sample is generated from one fixed world field B plus bounded body hard iron,
bounded residual and a bounded tilt-frame rotation error, then every sample and
therefore its arithmetic mean lies in the same deterministic perturbation ball.

For true field magnitude <= Bmax and relative tilt rotation angle delta,

    ||(R_hat R_true' - I) B|| <= 2 Bmax sin(delta/2).

Hence the startup mean is B + e with

    ||e|| <= bHImax + nmax + 2 Bmax sin(delta/2) = E.

The horizontal component obeys the same norm bound. If Htrue >= Hmin and
E < Hmin, the learned horizontal vector cannot vanish and the yaw-direction
error theta satisfies |sin(theta)| <= E/Hmin.

The canonical operating domain already declares startup world-averaged gravity
direction error <= 0.02 rad. Since sin(x)<=x for x>=0, this supplies
sin(delta/2)<=0.01 without any floating transcendental evaluation. Under the
commissioned 5 uT hard-iron / 2 uT residual envelope this gives E<=8.5 uT and
|sin(theta)|<=17/30. The final atan2/binary32 and full-SO(3) 45-degree entrance
comparison remain separate obligations.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_mag_source_qualification as Q


def R(x): return M.rational(x)
DECLARED_STARTUP_TILT_RAD_MAX=F(1,50) # 0.02 rad canonical operating-domain assumption

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
    """Open half-angle-sine threshold permitted by the commissioned envelope."""
    e=Q.default_envelope()
    numerator=e.world_field_horizontal_min-e.hard_iron_norm_max-e.residual_norm_max
    denominator=2*e.world_field_norm_max
    if numerator<=0: return F(0)
    return numerator/denominator


def declared_startup_tilt_capture():
    # delta <= 1/50 => delta/2 <= 1/100 and sin(delta/2)<=delta/2<=1/100.
    b=derive(F(1,100))
    require_capture(b)
    return b


def readiness():
    b=declared_startup_tilt_capture()
    return {
      'deterministic_average_does_not_claim_sqrtN_improvement':True,
      'earth_field_rotation_chord_bound_materialized':True,
      'hard_iron_and_residual_combined_without_independence_assumption':True,
      'horizontal_nonvanishing_capture_condition_materialized':True,
      'yaw_sine_error_supply_bound_materialized':True,
      'current_envelope_half_tilt_sin_margin_is_4_over_75':default_tilt_margin_half_sin()==F(4,75),
      'declared_startup_tilt_0p02_rad_attached_via_sin_x_le_x':True,
      'declared_tilt_plus_magnetic_envelope_yields_E_8p5_uT':b.mean_perturbation_max==F(17,2),
      'declared_tilt_plus_magnetic_envelope_yields_sin_yaw_le_17_over_30':b.sin_yaw_error_max==F(17,30),
      'atan2_binary32_yaw_bound_attached':False,
      'full_SO3_45deg_startup_entrance_closed':False,
      'startup_capture_closed':False,
      'ALT_STARTUP_PASS':False,
    }
