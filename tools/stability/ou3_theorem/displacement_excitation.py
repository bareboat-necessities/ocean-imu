"""Kinematic consequences of MARINE displacement-span excitation.

No instantaneous acceleration floor is inferred.  The lemmas expose signed
velocity/acceleration moments suitable for the existing same-history Abel/S
machinery.
"""
from __future__ import annotations
import math

def chord_moment_lower(duration_s, displacement_span_m, velocity_bound_mps):
    """Necessary signed acceleration moment after choosing a span chord.

    For t0<t1, T=t1-t0 and unit u along p(t1)-p(t0),
      u.(p1-p0)=T u.v0 + integral_0^T (T-s) u.a(s) ds.
    Hence |integral (T-s) u.a ds| >= max(0,P_E-T Vmax).
    Reversing time gives the analogous endpoint-v1 identity.  This bound can
    be zero for slow displacement; it must not be promoted to pointwise a.
    """
    vals=(duration_s,displacement_span_m,velocity_bound_mps)
    if not all(math.isfinite(x) and x>=0 for x in vals) or duration_s<=0:
        raise ValueError("finite T>0 and nonnegative span/velocity required")
    return max(0.0, displacement_span_m-duration_s*velocity_bound_mps)

def velocity_chord_lower(duration_s, displacement_span_m):
    """Mean-value lower bound: some |u.v| >= P_E/T on the chord interval."""
    if not (math.isfinite(duration_s) and duration_s>0 and
            math.isfinite(displacement_span_m) and displacement_span_m>=0):
        raise ValueError("finite T>0 and nonnegative span required")
    return displacement_span_m/duration_s

def zero_translation_excluded(displacement_e_m):
    if not math.isfinite(displacement_e_m) or displacement_e_m<=0:
        raise ValueError("qualified P_E>0 required")
    return True

def certificate():
    return {
      "qualification":"OU3_MARINE_DISPLACEMENT_EXCITATION_V1",
      "whole_window_displacement_diameter_required":True,
      "zero_translation_moving_alias_excluded_conditionally_on_P_E_positive":True,
      "some_velocity_projection_lower":"P_E/T_P",
      "signed_acceleration_moment_lower":"max(0,P_E-T_P*V_max)",
      "instantaneous_acceleration_floor_inferred":False,
      "numerical_T_P_and_P_E_qualification":"OPEN",
      "four_S_and_bias_abel_link_required":True,
      "full_moving_compatibility_excluded":False,
      "theorem_closed":False,
    }
