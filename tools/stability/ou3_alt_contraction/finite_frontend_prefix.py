"""Finite same-history raw-IMU -> frontend prefix for ALT.

One ``RawImuSample`` owns the body-frame gyro/accelerometer packet.  The same
packet advances the private Mahony vertical observer.  Shipping then uses that
vertical successor in two causally different ways: tracker/stillness and the
sigma tuner run first using a READ-ONLY view of the WPE state carried into the
sample; the current vertical sample advances WPE only later.  This module keeps
that exact ordering while retaining the full StillnessAdapter tracker-frequency
state for correspondence tests.

No source bound, tracker algorithm, transcendental binary32 result, or stability
property is invented here.  The stability-side tuner prefix projects out the
tracker-frequency state entirely.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as RAW
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_frontend_runtime as FRONT
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as WPE
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as STILL
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P


def R(x): return P.rational(x)


@dataclass(frozen=True)
class State:
    vertical: V.State
    wpe: WPE.WPEState
    band: BAND.BandState
    stats: BAND.StatsState
    tracker_lpf: FRONT.LPFState
    stillness: STILL.State
    sample_index: int = 0
    time: F = F(0)
    def __post_init__(self):
        if not isinstance(self.vertical,V.State): raise TypeError('vertical state required')
        if not isinstance(self.wpe,WPE.WPEState): raise TypeError('WPE state required')
        if not isinstance(self.band,BAND.BandState): raise TypeError('band state required')
        if not isinstance(self.stats,BAND.StatsState): raise TypeError('statistics state required')
        if not isinstance(self.tracker_lpf,FRONT.LPFState): raise TypeError('tracker LPF state required')
        if not isinstance(self.stillness,STILL.State): raise TypeError('stillness state required')
        if not isinstance(self.sample_index,int) or self.sample_index<0: raise ValueError('nonnegative sample index required')
        object.__setattr__(self,'time',R(self.time))
        if self.time<0: raise ValueError('nonnegative runtime time required')


@dataclass(frozen=True)
class Result:
    state: State
    raw_sample: RAW.RawImuSample
    vertical: V.Result
    preupdate_wpe: BAND.WPEFrequencyView
    band: BAND.FrontendResult
    tracker_lpf: FRONT.LPFResult
    stillness: STILL.Result
    tracker: FRONT.TrackerOutputWitness
    wpe: WPE.UpdateResult


def step(state:State, sample:RAW.RawImuSample, *, dt,
         vertical_cfg:V.Config, accel_invnorm:V.InvSqrtWitness|None,
         quat_invnorm:V.InvSqrtWitness|None, seed:V.SeedWitness|None,
         wpe_cfg:WPE.WPEConfig, wpe_decay:WPE.ExpWitness,
         wpe_moment_decay:WPE.ExpWitness|None=None,
         wpe_period_witness:WPE.PeriodWitness|None=None,
         wpe_log_witness:WPE.LogUpdateWitness|None=None,
         wpe_current_period=None, wpe_current_frequency=None,
         wpe_post_output:WPE.CanonicalOutputWitness|None=None,
         band_cfg:BAND.BandConfig=None, stats_cfg:BAND.StatsConfig=None,
         band_decay:BAND.BandDecayWitness=None,
         variance_decay:BAND.VarianceDecayWitness=None,
         bench_noise_sigma=F(0), noise_sqrt:BAND.NoiseSqrtWitness|None=None,
         tracker_lpf_decay:FRONT.LPFDecayWitness=None,
         tracker:FRONT.TrackerOutputWitness=None,
         still_cfg:STILL.Config=None,
         still_relax:STILL.RelaxWitness|None=None,
         still_attenuation:STILL.AttenuationWitness|None=None):
    if not isinstance(state,State) or not isinstance(sample,RAW.RawImuSample):
        raise TypeError('finite frontend state and RawImuSample required')
    dt=R(dt)
    if dt<=0: raise ValueError('positive dt required')
    if any(x is None for x in (band_cfg,stats_cfg,band_decay,variance_decay,tracker_lpf_decay,tracker,still_cfg)):
        raise TypeError('all runtime configs/branch witnesses required')

    vertical=RAW.vertical_step_from_raw(state.vertical,vertical_cfg,sample,dt=dt,
        accel_invnorm=accel_invnorm,quat_invnorm=quat_invnorm,seed=seed)
    view=FRONT.wpe_view(state.wpe,
        current_period=wpe_current_period if state.wpe.usable_period else None,
        current_frequency=wpe_current_frequency if state.wpe.usable_period else None)
    lpf=FRONT.tracker_lpf_step(state.tracker_lpf,vertical,decay=tracker_lpf_decay)
    still=FRONT.stillness_step_from_vertical(state.stillness,still_cfg,lpf,tracker,dt=dt,
        relax=still_relax,attenuation=still_attenuation)
    band=FRONT.band_step_from_vertical_view(state.band,state.stats,view,vertical,dt=dt,
        band_cfg=band_cfg,stats_cfg=stats_cfg,band_decay=band_decay,
        variance_decay=variance_decay,bench_noise_sigma=bench_noise_sigma,
        noise_sqrt=noise_sqrt)
    wpe=FRONT.wpe_step_from_vertical(state.wpe,wpe_cfg,vertical,dt=dt,decay=wpe_decay,
        moment_decay=wpe_moment_decay,period_witness=wpe_period_witness,
        log_witness=wpe_log_witness,current_period=wpe_current_period,
        current_frequency=wpe_current_frequency,post_output=wpe_post_output)

    nxt=State(vertical.state,wpe.state,band.band_state,band.stats_state,
              lpf.state,still.state,state.sample_index+1,state.time+dt)
    return Result(nxt,sample,vertical,view,band,lpf,still,tracker,wpe)


def assert_prediction_packet(result:Result, mekf_state, gyro_body_raw):
    if not isinstance(result,Result): raise TypeError('frontend prefix result required')
    return RAW.assert_prediction_gyro(result.raw_sample,mekf_state,gyro_body_raw)


def assert_measurement_packet(result:Result, accel_body_raw):
    if not isinstance(result,Result): raise TypeError('frontend prefix result required')
    return RAW.assert_acc_measurement_input(result.raw_sample,accel_body_raw)


def readiness():
    return {
      'same_raw_packet_private_vertical_and_MEKF_bindable':True,
      'same_vertical_successor_tuner_LPF_and_later_WPE':True,
      'tuner_uses_preupdate_WPE_state':True,
      'WPE_prior_takeover_band_stats_temporal_state_carried':True,
      'stillness_temporal_state_carried':True,
      'successive_frontend_samples_compose_without_state_restart':True,
      'tracker_algorithm_attached':False,
      'sensor_residual_source_bounds_attached':False,
      'frontend_transcendental_binary32_attached':False,
      'finite_MEKF_event_composed_in_same_function':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
