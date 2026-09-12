"""Finite asynchronous magnetometer gate and accelerometer-bias mode release.

``SeaStateFusionFilter_OU_III::updateMag`` is an external event, not part of each
IMU sample.  If magnetometry is enabled and wrapper time has reached the delay,
it calls the MEKF magnetic measurement, then increments ``mag_updates_applied``
regardless of whether the MEKF accepted that correction.  The first attempted
update timestamps ``first_mag_update_time``.  Only after all of

  * Live stage,
  * count >= configured unlock count (250 by default),
  * time-first_mag_time > 1 s,

may ``accel_bias_locked`` clear.  Clearing the lock enables BA estimation only
when no external hold is active.

This module proves that control recurrence and the literal covariance effects of
``set_acc_bias_updates_enabled``.  It deliberately does not invent a magnetic
sensor/source packet or an accepted magnetic update; those are separate finite
word obligations.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P


def R(x): return P.rational(x)
OFF_BA=18

@dataclass(frozen=True)
class Config:
    with_mag:bool=True
    mag_delay:F=F(7)
    unlock_count:int=250
    initial_ba_std:F=F(4,1000)
    def __post_init__(self):
        if not isinstance(self.with_mag,bool): raise TypeError('literal with_mag branch required')
        d,s=R(self.mag_delay),R(self.initial_ba_std)
        if d<0 or s<0 or not isinstance(self.unlock_count,int) or self.unlock_count<0:
            raise ValueError('valid mag delay/unlock count/BA std required')
        object.__setattr__(self,'mag_delay',d); object.__setattr__(self,'initial_ba_std',s)

@dataclass(frozen=True)
class State:
    updates:int=0
    first_time:F|None=None
    locked:bool=True
    hold:bool=False
    def __post_init__(self):
        if not isinstance(self.updates,int) or self.updates<0: raise ValueError('nonnegative mag update count required')
        if self.first_time is not None:
            t=R(self.first_time)
            if t<0: raise ValueError('nonnegative first-mag time required')
            object.__setattr__(self,'first_time',t)
        if not isinstance(self.locked,bool) or not isinstance(self.hold,bool):
            raise TypeError('literal bias lock/hold flags required')

@dataclass(frozen=True)
class GateResult:
    state:State
    filter_state:CORE.State
    attempted:bool
    unlocked_now:bool
    enabled_bias_now:bool


def _active(state:CORE.State,initial_std):
    """Literal set_acc_bias_updates_enabled(true) H->A covariance effect."""
    if state.mode=='A': return state
    floor=R(initial_std)**2
    cov=[list(r) for r in state.covariance]
    for i in range(3): cov[OFF_BA+i][OFF_BA+i]=max(cov[OFF_BA+i][OFF_BA+i],floor)
    return CORE.State('A',state.z,tuple(tuple(r) for r in cov),state.q_hat,state.reference)


def _held(state:CORE.State):
    """Literal set_acc_bias_updates_enabled(false) A->H covariance effect."""
    if state.mode=='H': return state
    cov=[list(r) for r in state.covariance]
    for i in range(3):
        k=OFF_BA+i
        for j in range(21):
            if OFF_BA <= j < OFF_BA+3: continue
            cov[k][j]=F(0); cov[j][k]=F(0)
    return CORE.State('H',state.z,tuple(tuple(r) for r in cov),state.q_hat,state.reference)


def update_mag_call(control:State,filter_state:CORE.State,cfg:Config,*,time,live,
                    measurement_state:CORE.State|None=None):
    """Apply wrapper control after one external updateMag call.

    ``measurement_state`` is the MEKF state returned by the separately proved
    magnetic measurement branch.  It is required iff the wrapper attempts the
    measurement.  Acceptance is intentionally not an input: wrapper count/time
    bookkeeping occurs after the call even when SafeLDLT rejects internally.
    """
    if not isinstance(control,State) or not isinstance(filter_state,CORE.State) or not isinstance(cfg,Config):
        raise TypeError('mag control/filter/config required')
    if not isinstance(live,bool): raise TypeError('literal Live branch required')
    t=R(time)
    if t<0: raise ValueError('nonnegative wrapper time required')
    attempted=cfg.with_mag and t>=cfg.mag_delay
    if not attempted:
        if measurement_state is not None: raise ValueError('gated-off mag call consumes no measurement successor')
        return GateResult(control,filter_state,False,False,False)
    if not isinstance(measurement_state,CORE.State):
        raise TypeError('attempted mag call requires separately composed MEKF measurement successor')
    if measurement_state.reference != filter_state.reference:
        raise ValueError('mag measurement successor detached from same physical endpoint')

    n=control.updates+1
    first=t if control.first_time is None else control.first_time
    locked=control.locked
    unlocked=False
    enabled=False
    fs=measurement_state
    if locked and live and n>=cfg.unlock_count and (t-first)>1:
        locked=False; unlocked=True
        if not control.hold:
            fs=_active(fs,cfg.initial_ba_std); enabled=(filter_state.mode=='H' and fs.mode=='A') or fs.mode=='A'
    nxt=State(n,first,locked,control.hold)
    return GateResult(nxt,fs,True,unlocked,enabled)


def set_hold(control:State,filter_state:CORE.State,cfg:Config,*,hold,live):
    """Literal setAccBiasHold control effect outside updateMag."""
    if not isinstance(hold,bool) or not isinstance(live,bool): raise TypeError('literal hold/Live branches required')
    if hold==control.hold: return GateResult(control,filter_state,False,False,False)
    fs=filter_state
    enabled=False
    if hold:
        fs=_held(fs)
    elif not control.locked and live:
        fs=_active(fs,cfg.initial_ba_std); enabled=True
    return GateResult(State(control.updates,control.first_time,control.locked,hold),fs,False,False,enabled)


def readiness():
    return {
      'mag_delay_and_with_mag_gate_materialized':True,
      'attempt_count_independent_of_measurement_acceptance':True,
      'first_attempt_timestamp_materialized':True,
      'unlock_count_and_strict_one_second_guard_materialized':True,
      'external_hold_blocks_unlock_enable_but_not_lock_clear':True,
      'H_to_A_BA_variance_floor_materialized':True,
      'A_to_H_zero_BA_cross_covariance_materialized':True,
      'magnetic_measurement_same_source_packet_attached':False,
      'Rmag_runtime_ancestry_attached':False,
      'mag_delay_clock_binary64_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
