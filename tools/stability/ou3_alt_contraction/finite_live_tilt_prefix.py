"""Compose the persistent Live tilt watchdog with the default Live IMU prefix.

Shipping evaluates the watchdog immediately after the accelerometer correction
and before the measurement-only tuner/direction suffix.  ``finite_live_imu_prefix``
already proves that the represented tuner suffix and periodic a_w-sync state are
disjoint from the MEKF hard-reset state.  We therefore execute the ordinary
prefix with its old no-reset branch, evaluate the exact watchdog recurrence from
the same post-accelerometer sample, and apply the preserve-yaw hard reset to the
MEKF component of the product successor when it fires.  This commutation changes
no coupled state and avoids duplicating frontend/WPE algebra.

The tilt angle itself is a same-expression runtime operand whose acos/binary32
ancestry remains open.  A firing reset requires ``PreserveYawWitness``; a
nonfiring branch is forbidden from consuming one.
"""
from __future__ import annotations
from dataclasses import dataclass, replace

from tools.stability.ou3_alt_contraction import finite_live_imu_prefix as LIVE
from tools.stability.ou3_alt_contraction import finite_tilt_watchdog as WATCH
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR


@dataclass(frozen=True)
class State:
    live: LIVE.State
    watchdog: WATCH.State
    def __post_init__(self):
        if not isinstance(self.live,LIVE.State) or not isinstance(self.watchdog,WATCH.State):
            raise TypeError('Live IMU and watchdog states required')


@dataclass(frozen=True)
class Result:
    state: State
    live: LIVE.Result
    watchdog: WATCH.StepResult
    reset_applied: bool


def step(state:State,raw:SENSOR.RawImuSample,segment,*,tilt_deg,
         reset_witness:WATCH.PreserveYawWitness|None=None,
         reset_tilt_sigma=None,reset_yaw_sigma=None,**live_kwargs):
    if not isinstance(state,State): raise TypeError('Live+watchdog state required')
    if 'tilt_reset_due' in live_kwargs:
        raise TypeError('tilt-reset branch is owned by persistent watchdog composer')
    out=LIVE.step(state.live,raw,segment,tilt_reset_due=False,**live_kwargs)
    wd=WATCH.step(state.watchdog,dt=segment.h,tilt_deg=tilt_deg)
    mekf=out.state.mekf
    if wd.fired:
        if reset_witness is None:
            raise ValueError('firing Live tilt watchdog requires preserve-yaw reset witness')
        rkw={}
        if reset_tilt_sigma is not None: rkw['tilt_sigma']=reset_tilt_sigma
        if reset_yaw_sigma is not None: rkw['yaw_sigma']=reset_yaw_sigma
        mekf=WATCH.preserve_yaw_reset(mekf,out.guarded,reset_witness,**rkw)
    else:
        if reset_witness is not None or reset_tilt_sigma is not None or reset_yaw_sigma is not None:
            raise ValueError('nonfiring Live tilt watchdog consumes no reset operands')
    live_state=replace(out.state,mekf=mekf)
    return Result(State(live_state,wd.state),out,wd,wd.fired)


def readiness():
    return {
      'persistent_tilt_overlimit_and_cooldown_carried_in_Live_word':True,
      'watchdog_evaluated_from_post_accel_sample':True,
      'watchdog_reset_commutes_with_measurement_only_tuner_suffix':True,
      'firing_edge_applies_preserve_yaw_hard_reset_structure':True,
      'nonfiring_edge_consumes_no_reset_witness':True,
      'watchdog_tilt_acos_binary32_attached':False,
      'preserve_yaw_atan_asin_angleaxis_binary32_attached':False,
      'async_magnetometer_branch_attached':False,
      'source_uniform_COMPLETE_BRMM_bounds_attached':False,
      'deployment_finite_precision_closed':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
