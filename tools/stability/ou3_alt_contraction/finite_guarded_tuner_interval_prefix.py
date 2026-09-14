"""Guard-persistent shipping frontend with interval SpectralMSE tuner state.

This is the startup/source-facing wrapper for
``finite_tuner_frontend_interval_prefix``.  One raw IMU packet first advances
the literal persistent AccelVibrationGuard.  The exact guarded packet is then
consumed by the private Mahony / band / stillness / interval-tuner / later-WPE
chain.  No second acceleration sample can be injected between the guard and the
tuner frontend.

The wrapper closes same-history guard persistence for the exact-real interval
tuner relation.  It does not yet close guard exp/sqrt binary32 ancestry, Racc
inflation, MEKF event interleaving, or the binary32 SpectralMSE sqrt/pow/exp
correspondence needed to choose the actual committed R_S float.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as GUARD
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_tuner_frontend_interval_prefix as TUNER

QUALIFICATION='OU3_ALT_GUARDED_INTERVAL_TUNER_PREFIX_V1'


@dataclass(frozen=True)
class State:
    guard:GUARD.State
    tuner:TUNER.State
    def __post_init__(self):
        if not isinstance(self.guard,GUARD.State) or not isinstance(self.tuner,TUNER.State):
            raise TypeError('guard and interval tuner-prefix states required')


@dataclass(frozen=True)
class Result:
    state:State
    guarded_sample:SENSOR.GuardedImuSample
    tuner:TUNER.Result


def step(state:State,raw:SENSOR.RawImuSample,*,dt,guard_cfg:GUARD.Config,
         guard_decay:GUARD.DecayWitness|None=None,guard_rms:GUARD.RmsWitness|None=None,
         **tuner_kwargs):
    if not isinstance(state,State) or not isinstance(raw,SENSOR.RawImuSample):
        raise TypeError('guarded interval-prefix state and raw packet required')
    guarded=SENSOR.guarded_sample(raw,state.guard,guard_cfg,dt=dt,
                                  decay=guard_decay,rms=guard_rms)
    if 'sample' in tuner_kwargs or 'dt' in tuner_kwargs:
        raise TypeError('sample and dt are owned by guard-persistent interval prefix')
    tout=TUNER.step(state.tuner,guarded,dt=dt,**tuner_kwargs)
    return Result(State(guarded.guard.state,tout.state),guarded,tout)


def readiness():
    t=TUNER.readiness()
    return {
      'qualification':QUALIFICATION,
      'raw_accel_guard_state_persists_across_interval_tuner_samples':True,
      'one_guard_successor_feeds_private_vertical_band_tuner_and_WPE':True,
      'guarded_sample_preserves_raw_COMPLETE_BRMM_sensor_ancestry':True,
      'postCold_general_SpectralMSE_interval_candidate_composed':t['postCold_general_SpectralMSE_interval_candidate_composed'],
      'actual_0p2Hz_prior_postCold_candidate_representable_on_guarded_path':t['actual_0p2Hz_prior_postCold_candidate_representable'],
      'guard_exp_sqrt_binary32_ancestry_attached':False,
      'Racc_inflation_from_same_guard_excess_attached':False,
      'finite_MEKF_event_interleaved_in_same_sample_word':False,
      'binary32_sqrt_pow_exp_correspondence_closed':False,
      'goLive_interval_RS_to_actual_MEKF_commit_closed':False,
      'source_uniform_complete_startup_reachability_closed':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
