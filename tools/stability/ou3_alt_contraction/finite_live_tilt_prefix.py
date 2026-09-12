"""Compose the persistent Live tilt watchdog with the default Live IMU prefix.

Shipping evaluates the watchdog immediately after the accelerometer correction
and before the measurement-only tuner/direction suffix. ``finite_live_imu_prefix``
already proves that the represented tuner suffix and periodic a_w-sync state are
disjoint from the MEKF hard-reset state. We therefore execute the ordinary
prefix with its old no-reset branch, evaluate the watchdog from the same
post-accelerometer state, and apply the preserve-yaw hard reset to the MEKF
component of the product successor when it fires.

``step`` is retained as a low-level conditional relation. The theorem-facing
``step_from_shipping_operands`` removes the free tilt scalar and final reset
quaternion by deriving both from the actual state and guarded sample.
"""
from __future__ import annotations
from dataclasses import dataclass, replace

from tools.stability.ou3_alt_contraction import finite_live_imu_prefix as LIVE
from tools.stability.ou3_alt_contraction import finite_tilt_watchdog as WATCH
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_tilt_reset_runtime as RESET


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


def step_from_shipping_operands(state:State,raw:SENSOR.RawImuSample,segment,*,
         old_q_norm=None,old_yaw_half=None,acc_tilt=None,pitch_cos=None,
         pitch_half=None,roll_half=None,reset_tilt_sigma=None,reset_yaw_sigma=None,
         **live_kwargs):
    """Theorem-facing Live step with no free tilt scalar or reset quaternion.

    The ordinary prefix is evaluated through the post-accelerometer state, the
    >70deg predicate is derived from that exact nominal attitude, and a firing
    reset derives its preserve-yaw quaternion/down axis from the same predecessor
    and guarded accelerometer. Deployment threshold/libm boundary cases remain
    fail-closed in ``finite_tilt_reset_runtime``.
    """
    if not isinstance(state,State): raise TypeError('Live+watchdog state required')
    if 'tilt_reset_due' in live_kwargs:
        raise TypeError('tilt-reset branch is owned by persistent watchdog composer')
    out=LIVE.step(state.live,raw,segment,tilt_reset_due=False,**live_kwargs)
    over=RESET.watchdog_over_limit(out.accelerometer.state)
    wd=WATCH.step_over_limit(state.watchdog,dt=segment.h,over_limit=over)
    mekf=out.state.mekf
    operands=(old_q_norm,old_yaw_half,acc_tilt,pitch_cos,pitch_half,roll_half,
              reset_tilt_sigma,reset_yaw_sigma)
    if wd.fired:
        if old_q_norm is None or acc_tilt is None or pitch_cos is None or pitch_half is None:
            raise ValueError('firing preserve-yaw reset requires same-operand arithmetic witnesses')
        witness=RESET.preserve_yaw_witness(
            out.accelerometer.state,out.guarded,old_q_norm=old_q_norm,
            old_yaw_half=old_yaw_half,acc_tilt=acc_tilt,pitch_cos=pitch_cos,
            pitch_half=pitch_half,roll_half=roll_half)
        rkw={}
        if reset_tilt_sigma is not None: rkw['tilt_sigma']=reset_tilt_sigma
        if reset_yaw_sigma is not None: rkw['yaw_sigma']=reset_yaw_sigma
        mekf=WATCH.preserve_yaw_reset(mekf,out.guarded,witness,**rkw)
    elif any(x is not None for x in operands):
        raise ValueError('nonfiring Live watchdog consumes no reset arithmetic operands')
    live_state=replace(out.state,mekf=mekf)
    return Result(State(live_state,wd.state),out,wd,wd.fired)


def readiness():
    return {
      'persistent_tilt_overlimit_and_cooldown_carried_in_Live_word':True,
      'watchdog_evaluated_from_post_accel_sample':True,
      'watchdog_reset_commutes_with_measurement_only_tuner_suffix':True,
      'firing_edge_applies_preserve_yaw_hard_reset_structure':True,
      'theorem_entry_derives_watchdog_predicate_from_post_accel_attitude':True,
      'theorem_entry_derives_reset_quaternion_from_same_guarded_sample':True,
      'nonfiring_edge_consumes_no_reset_witness':True,
      # Retain the older fail-closed key for downstream ledgers/tests while the
      # more precise key below records that only threshold/libm correspondence
      # remains open, not the state dependency itself.
      'watchdog_tilt_acos_binary32_attached':False,
      'watchdog_threshold_boundary_and_binary32_attached':False,
      'preserve_yaw_atan_asin_angleaxis_binary32_attached':False,
      'async_magnetometer_branch_attached':False,
      'source_uniform_COMPLETE_BRMM_bounds_attached':False,
      'deployment_finite_precision_closed':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
