"""Derive the actual fresh H18 joint24 error from startup Live entry.

No fresh-entry error box is postulated.  Given the SAME physical ``Reference``
and the real-arithmetic shipping estimator state after ``goLive/enterLive_``, the
joint24 coordinates are defined exactly as physical truth minus estimator:

  c      = Cayley(q_true_WB * conjugate(q_hat_WB))
  e_bg   = b_g,true - b_g,hat
  e_v    = v_true - v_hat
  e_p    = p_true - p_hat
  e_S    = S_true(centered at the one Live origin) - S_hat
  e_aw   = a_true - a_w,hat
  e_ba   = beta_true - b_a,hat
  beta   = beta_true

The resulting object is the existing full ``finite_core.State`` in H mode with
the full 21-state covariance.  Thus startup does not rely on covariance
consistency or an independently selected entry-error radius.

This is exact real arithmetic.  Source qualification of the physical endpoint,
deployment normalization/clock roundoff, and proof that the resulting values lie
inside the later storage basin remain separate obligations.
"""
from __future__ import annotations

from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_startup_live_entry as LIVE
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE


def _sub(a,b): return tuple(x-y for x,y in zip(a,b))


def build(entry:LIVE.Result,reference:CORE.Reference,*,scope:SCOPE.Scope):
    if not isinstance(entry,LIVE.Result): raise TypeError('fresh startup Live-entry result required')
    if not isinstance(reference,CORE.Reference): raise TypeError('same-history physical Reference required')
    SCOPE.assert_certified_scope(scope)
    if not entry.zero_wind_heel or entry.startup_stage!='Live' or entry.acc_bias_updates_enabled:
        raise ValueError('fresh certified entry must be zero-heel Live H18')
    if reference.time != reference.live_origin:
        raise ValueError('fresh entry Reference must be at the one-time Live origin')
    if any(reference.centered_S):
        raise ValueError('fresh physical centered S must be zero at Live origin')

    x=entry.x
    # Shipping NX layout with gyro bias enabled and accelerometer bias present:
    # att[0:3], bg[3:6], v[6:9], p[9:12], S[12:15], aw[15:18], ba[18:21].
    relative=P.quat_mul(reference.q_world_to_body,P.quat_conj(entry.q_hat))
    c=tuple(CORE.cayley(relative))
    z=(c+
       _sub(reference.gyro_bias,x[3:6])+
       _sub(reference.velocity,x[6:9])+
       _sub(reference.position,x[9:12])+
       _sub(reference.centered_S,x[12:15])+
       _sub(reference.acceleration,x[15:18])+
       _sub(reference.beta,x[18:21])+
       tuple(reference.beta))
    if len(z)!=24: raise AssertionError('joint24 layout error')
    return CORE.State('H',z,entry.P,entry.q_hat,reference)


def readiness():
    return {
      'fresh_entry_joint24_derived_from_truth_minus_estimator':True,
      'fresh_attitude_error_from_same_true_and_nominal_quaternions':True,
      'one_time_Live_origin_enforced_at_entry':True,
      'fresh_centered_physical_S_zero_enforced':True,
      'full_21_state_covariance_retained_in_H18':True,
      'covariance_consistency_not_used_as_entry_assumption':True,
      'independent_fresh_entry_error_box_removed':True,
      'deployment_binary32_correspondence_closed':False,
      'fresh_entry_inside_storage_basin_proved':False,
      'complete_same_history_startup_to_Live_word':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
