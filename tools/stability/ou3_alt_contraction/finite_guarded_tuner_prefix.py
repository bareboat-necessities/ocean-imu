"""Guard-persistent shipping frontend/tuner prefix for ALT.

This is the shipping-level wrapper around ``finite_tuner_frontend_prefix``.
One raw IMU packet first advances the persistent AccelVibrationGuard.  The exact
guard successor is wrapped as ``GuardedImuSample`` and then consumed by private
Mahony / band / stillness / tuner / later-WPE logic.  Guard state is carried to
the next sample, so no sample may restart the LP cascade, detector RMS, or
engagement slew.

The same ``GuardedImuSample`` object is intentionally returned for the MEKF
accelerometer event.  The full Live word can therefore interleave prediction and
accelerometer correction between Mahony and tuner processing without inventing a
second ``acc_in``.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as GUARD
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_tuner_frontend_prefix as TUNER


@dataclass(frozen=True)
class State:
    guard: GUARD.State
    tuner: TUNER.State
    def __post_init__(self):
        if not isinstance(self.guard,GUARD.State) or not isinstance(self.tuner,TUNER.State):
            raise TypeError('guard and tuner-prefix states required')


@dataclass(frozen=True)
class Result:
    state: State
    guarded_sample: SENSOR.GuardedImuSample
    tuner: TUNER.Result


def step(state:State,raw:SENSOR.RawImuSample,*,dt,guard_cfg:GUARD.Config,
         guard_decay:GUARD.DecayWitness|None=None,guard_rms:GUARD.RmsWitness|None=None,
         **tuner_kwargs):
    if not isinstance(state,State) or not isinstance(raw,SENSOR.RawImuSample):
        raise TypeError('guarded-prefix state and raw packet required')
    guarded=SENSOR.guarded_sample(raw,state.guard,guard_cfg,dt=dt,
                                  decay=guard_decay,rms=guard_rms)
    if 'sample' in tuner_kwargs or 'dt' in tuner_kwargs:
        raise TypeError('sample and dt are owned by guard-persistent prefix')
    tout=TUNER.step(state.tuner,guarded,dt=dt,**tuner_kwargs)
    return Result(State(guarded.guard.state,tout.state),guarded,tout)


def readiness():
    return {
      'raw_accel_guard_state_persists_across_samples':True,
      'one_guard_successor_feeds_private_vertical_and_future_MEKF_event':True,
      'guarded_sample_preserves_raw_COMPLETE_BRMM_sensor_ancestry':True,
      'tuner_before_WPE_order_retained_under_guard_wrapper':True,
      'guard_exp_sqrt_binary32_ancestry_attached':False,
      'Racc_inflation_from_same_guard_excess_attached':False,
      'finite_MEKF_event_interleaved_in_same_sample_word':False,
      'sensor_residual_source_bounds_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
