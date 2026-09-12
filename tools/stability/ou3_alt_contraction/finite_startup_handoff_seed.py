"""Startup proxy-to-MEKF handoff seed/control relation for ALT.

Shipping handOffToMekf_ uses the persistent startup proxy quaternion.  If a
magnetic yaw gauge was learned before Live, boatQuatWithAbsoluteYaw_ removes the
proxy's current heading and writes the pending absolute yaw; otherwise the proxy
quaternion is used directly.  goLive then receives:

  tilt sigma = 0.035 rad,
  yaw sigma  = 0.087 rad when north-gauged, else 1.5708 rad,
  allow_acc_bias = false.

This module proves the seed/control choice only.  The exact
initialize_from_attitude covariance/state mutation is a separate obligation.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as VERT
from tools.stability.ou3_alt_contraction import finite_mag_startup_ready as READY
from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as TILT

TILT_SIGMA=F(35,1000)
YAW_SIGMA_GAUGED=F(87,1000)
YAW_SIGMA_FREE=F(15708,10000)

@dataclass(frozen=True)
class Result:
    q_seed:tuple
    tilt_sigma:F
    yaw_sigma:F
    allow_acc_bias:bool
    gauged:bool


def seed(proxy:VERT.State,pending:READY.PendingYaw|None,*,
         proxy_q_norm:TILT.SqrtWitness|None=None,
         proxy_yaw_half:TILT.YawHalfWitness|None=None):
    if not isinstance(proxy,VERT.State): raise TypeError('persistent startup proxy state required')
    if pending is None:
        if proxy_q_norm is not None or proxy_yaw_half is not None:
            raise ValueError('ungauged handoff uses proxy quaternion directly and consumes no yaw-strip witnesses')
        q=tuple(M.vec(proxy.q,4))
        if M.dot(q,q)==0: raise ValueError('zero startup proxy quaternion')
        # Shipping checks allFinite on q_seed before goLive; normalization is
        # performed by initialize_from_attitude.  Keep this exact predecessor.
        return Result(q,TILT_SIGMA,YAW_SIGMA_FREE,False,False)
    if not isinstance(pending,READY.PendingYaw): raise TypeError('pending startup yaw lock required')
    if not isinstance(proxy_q_norm,TILT.SqrtWitness):
        raise TypeError('gauged handoff requires proxy normalization witness')
    tilt=TILT.yaw_removed(proxy.q,q_norm=proxy_q_norm,yaw_half=proxy_yaw_half)
    q=tuple(P.quat_mul(pending.q_abs,tilt.q_tilt))
    if M.dot(q,q)!=1: raise AssertionError('absolute-yaw times proxy-tilt seed lost unit norm')
    return Result(q,TILT_SIGMA,YAW_SIGMA_GAUGED,False,True)


def readiness():
    return {
      'ungauged_handoff_uses_proxy_quaternion_directly':True,
      'gauged_handoff_strips_proxy_yaw_before_pending_absolute_yaw':True,
      'gauged_yaw_sigma_0087_materialized':True,
      'free_yaw_sigma_15708_materialized':True,
      'tilt_sigma_0035_materialized':True,
      'handoff_allow_acc_bias_false':True,
      'initialize_from_attitude_state_covariance_attached':False,
      'handoff_transcendental_binary32_attached':False,
      'startup_capture_closed':False,
      'ALT_STARTUP_PASS':False,
    }
