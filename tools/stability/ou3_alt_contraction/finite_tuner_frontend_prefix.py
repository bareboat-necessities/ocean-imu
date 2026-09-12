"""Finite raw-IMU -> tracker-free tuner-candidate sample prefix.

This is the adaptation-side sample composer for ALT.  It keeps persistent
private-Mahony, WPE, adaptive-band/statistics, tracker-input LPF, projected
stillness, and TuneState memory.  A single RawImuSample supplies the private
vertical observer.  That same vertical successor supplies WPE, sigma-band and
LPF.  The tuner-relevant stillness projection is computed from the LPF output
without any dominant-frequency tracker state, then the same band/statistics and
stillness successors generate the tau/sigma/R_S candidate.

The staged commit at the *next* IMU boundary is intentionally not performed in
this module.  This module proves the sample-k candidate ancestry; the existing
finite_tuner_commit module proves the subsequent boundary transaction.  The
next composition obligation is to connect those two states without inserting
an independent TuneState or band-noise-floor value.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as RAW
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V
from tools.stability.ou3_alt_contraction import finite_frontend_runtime as FRONT
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as WPE
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as BAND
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as STILL_FULL
from tools.stability.ou3_alt_contraction import finite_tuner_stillness_projection as STILL
from tools.stability.ou3_alt_contraction import finite_tuner_projection_bridge as BRIDGE
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction.finite_tuner_commit import TuneState
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
    tune: TuneState
    last_adapt_time: F = F(0)
    sample_index: int = 0
    time: F = F(0)
    def __post_init__(self):
        types=((self.vertical,V.State),(self.wpe,WPE.WPEState),(self.band,BAND.BandState),
               (self.stats,BAND.StatsState),(self.tracker_lpf,FRONT.LPFState),
               (self.stillness,STILL.State),(self.tune,TuneState))
        if any(not isinstance(x,t) for x,t in types): raise TypeError('invalid tuner-prefix state component')
        object.__setattr__(self,'last_adapt_time',R(self.last_adapt_time)); object.__setattr__(self,'time',R(self.time))
        if self.last_adapt_time<0 or self.time<self.last_adapt_time: raise ValueError('invalid tuner-prefix clocks')
        if not isinstance(self.sample_index,int) or self.sample_index<0: raise ValueError('nonnegative sample index required')


@dataclass(frozen=True)
class Result:
    state: State
    raw_sample: RAW.RawImuSample
    vertical: V.Result
    wpe: WPE.UpdateResult
    band: BAND.FrontendResult
    tracker_lpf: FRONT.LPFResult
    stillness: STILL.Result
    candidate: CAND.CandidateResult


def step(state:State,sample:RAW.RawImuSample,*,dt,
         vertical_cfg:V.Config,accel_invnorm:V.InvSqrtWitness|None,
         quat_invnorm:V.InvSqrtWitness|None,seed:V.SeedWitness|None,
         wpe_cfg:WPE.WPEConfig,wpe_decay:WPE.ExpWitness,
         wpe_moment_decay:WPE.ExpWitness|None=None,
         wpe_period_witness:WPE.PeriodWitness|None=None,
         wpe_log_witness:WPE.LogUpdateWitness|None=None,
         wpe_current_period=None,wpe_current_frequency=None,
         wpe_post_output:WPE.CanonicalOutputWitness|None=None,
         band_cfg:BAND.BandConfig=None,stats_cfg:BAND.StatsConfig=None,
         band_decay:BAND.BandDecayWitness=None,
         variance_decay:BAND.VarianceDecayWitness=None,
         bench_noise_sigma=F(0),noise_sqrt:BAND.NoiseSqrtWitness|None=None,
         tracker_lpf_decay:FRONT.LPFDecayWitness=None,
         still_cfg:STILL_FULL.Config=None,
         still_attenuation:STILL_FULL.AttenuationWitness|None=None,
         candidate_cfg:CAND.CandidateConfig=None,sigma_wave_sqrt=None,
         spectral:CAND.SpectralWitness=None,ema:CAND.EmaWitness=None):
    if not isinstance(state,State) or not isinstance(sample,RAW.RawImuSample):
        raise TypeError('finite tuner-prefix state and RawImuSample required')
    dt=R(dt)
    if dt<=0: raise ValueError('positive dt required')
    required=(band_cfg,stats_cfg,band_decay,variance_decay,tracker_lpf_decay,still_cfg,candidate_cfg,sigma_wave_sqrt,spectral,ema)
    if any(x is None for x in required): raise TypeError('all represented runtime configs/witnesses required')

    vertical=RAW.vertical_step_from_raw(state.vertical,vertical_cfg,sample,dt=dt,
        accel_invnorm=accel_invnorm,quat_invnorm=quat_invnorm,seed=seed)
    wpe=FRONT.wpe_step_from_vertical(state.wpe,wpe_cfg,vertical,dt=dt,decay=wpe_decay,
        moment_decay=wpe_moment_decay,period_witness=wpe_period_witness,
        log_witness=wpe_log_witness,current_period=wpe_current_period,
        current_frequency=wpe_current_frequency,post_output=wpe_post_output)
    band=FRONT.band_step_from_vertical(state.band,state.stats,wpe,vertical,dt=dt,
        band_cfg=band_cfg,stats_cfg=stats_cfg,band_decay=band_decay,
        variance_decay=variance_decay,bench_noise_sigma=bench_noise_sigma,noise_sqrt=noise_sqrt)
    lpf=FRONT.tracker_lpf_step(state.tracker_lpf,vertical,decay=tracker_lpf_decay)
    still=STILL.step(state.stillness,still_cfg,a_vert_up_lp=lpf.output,dt=dt,
                     attenuation=still_attenuation)
    now=state.time+dt
    cand=BRIDGE.step(state.tune,band,still,candidate_cfg,sigma_wave_sqrt=sigma_wave_sqrt,
                     dt=dt,time=now,last_adapt_time=state.last_adapt_time,
                     spectral=spectral,ema=ema)
    nxt=State(vertical.state,wpe.state,band.band_state,band.stats_state,lpf.state,
              still.state,cand.tune_next,cand.last_adapt_time_after,
              state.sample_index+1,now)
    return Result(nxt,sample,vertical,wpe,band,lpf,still,cand)


def readiness():
    return {
      'raw_IMU_to_private_vertical_same_packet':True,
      'same_vertical_WPE_band_tracker_LPF':True,
      'tracker_free_tuner_stillness_projection':True,
      'band_statistics_and_stillness_generate_candidate_same_sample':True,
      'TuneState_and_adapt_clock_persist_across_samples':True,
      'dominant_frequency_tracker_absent_from_OU_tuner_prefix':True,
      'next_boundary_staged_commit_composed':False,
      'sensor_residual_source_bounds_attached':False,
      'transcendental_binary32_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
