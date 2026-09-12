"""Finite startup ``goLive``/``enterLive_`` composition for ALT.

Scope: certified zero-wind-heel deployment and the magnetically gauged handoff.
The wrapper calls inner ``goLive(..., allow_acc_bias=false)``; inner goLive first
executes ``initialize_from_attitude`` and then ``enterLive_``.  ``enterLive_``
commits the already-qualified tuner operating point, seats the a_w covariance on
the committed stationary covariance, clears all a_w cross-covariances, keeps
accelerometer-bias learning held, installs Live R_S and enters StartupStage::Live.

This module captures the exact real-arithmetic state/covariance mutation after
the handoff reset.  Binary32 setter/normalization correspondence and wrapper
clock/event arithmetic remain separate obligations.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_startup_handoff_init as INIT
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE

OFF_AW=15
OFF_BA=18
N=21

@dataclass(frozen=True)
class Result:
    x:tuple
    P:tuple
    active:ACTIVE.ActiveParameters
    startup_stage:str
    startup_stage_t:F
    acc_bias_updates_enabled:bool
    gauged:bool
    zero_wind_heel:bool


def enter_live(handoff:INIT.Result,active:ACTIVE.ActiveParameters,*,scope:SCOPE.Scope,
               accel_bias_locked=True,acc_bias_hold=False):
    if not isinstance(handoff,INIT.Result): raise TypeError('initialized startup handoff required')
    if not handoff.gauged: raise ValueError('this startup Live-entry certificate requires gauged handoff')
    if not isinstance(active,ACTIVE.ActiveParameters): raise TypeError('committed active tuner parameters required')
    SCOPE.assert_certified_scope(scope)
    if not accel_bias_locked:
        raise ValueError('fresh certified H18 entry requires startup accelerometer-bias lock retained')

    x=list(M.vec(handoff.x,N)); P=[list(r) for r in M.mat(handoff.P,N,N)]

    # shipping reset_aw_covariance_to_stationary(): clear every AW cross block
    # and replace P_aw,aw with the same committed Sigma_aw.
    for i in range(N):
        if i<OFF_AW or i>=OFF_AW+3:
            for a in range(3):
                P[OFF_AW+a][i]=F(0); P[i][OFF_AW+a]=F(0)
    for a in range(3):
        for b in range(3): P[OFF_AW+a][OFF_AW+b]=active.Sigma_aw[a][b]

    # goLive is called with allow_acc_bias=false and the startup lock remains
    # engaged, hence enterLive_ keeps the shipping MEKF on H18.
    allow_bias=(not accel_bias_locked) and (not acc_bias_hold)
    if allow_bias: raise AssertionError('fresh startup unexpectedly released accelerometer bias')

    if active.R_S is None:
        raise ValueError('fresh Live entry requires committed Live R_S')

    return Result(tuple(x),tuple(tuple(r) for r in P),active,'Live',F(0),False,
                  True,True)


def readiness():
    return {
      'zero_wind_heel_scope_consumed':True,
      'goLive_initialize_then_enterLive_order_attached':True,
      'fresh_aw_covariance_seated_on_same_committed_Sigma_aw':True,
      'fresh_aw_cross_covariances_zeroed':True,
      'fresh_accelerometer_bias_gate_remains_H18':True,
      'committed_Live_RS_required':True,
      'fresh_real_arithmetic_H18_entry_composed':True,
      'wrapper_live_clock_binary32_attached':False,
      'setter_and_normalization_binary32_attached':False,
      'complete_same_history_startup_to_Live_word':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
