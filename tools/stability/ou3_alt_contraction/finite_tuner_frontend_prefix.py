"""Finite raw-IMU -> tracker-free tuner sample prefix.

Persistent state carries the private Mahony observer, WPE, adaptive band and
statistics, tracker-input LPF, tuner-relevant StillnessAdapter projection,
TuneState, adaptation clock/pending bit, and the shipping startup stage clock.

Literal shipping order is enforced:
  1. current raw packet advances private vertical;
  2. tracker LPF/stillness and adaptive band/statistics advance, with the tuner
     reading a READ-ONLY view of the WPE state carried into the sample;
  3. Cold/TunerWarm/TunerReady/Live stage logic decides whether the candidate
     EMA executes; Cold always returns after the frontend update, while
     TunerWarm executes adaptation as soon as tuner frequency is ready and
     promotes to TunerReady only when variance is ready and the PRE-UPDATE WPE
     usable latch is already true;
  4. only after that does the current vertical sample advance WPE for the next
     IMU sample.

A pending online candidate must be committed at the next IMU boundary before
this function may consume another sample.  goLive()/attitude handoff remains an
external hybrid transition; this module represents the startup tuner stages on
either side of that handoff but does not invent its attitude qualification.
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
STAGES=('Cold','TunerWarm','TunerReady','Live')


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
    pending: bool = False
    sample_index: int = 0
    time: F = F(0)
    stage: str = 'Live'
    stage_time: F = F(0)
    warmup_sec: F = F(5)
    def __post_init__(self):
        types=((self.vertical,V.State),(self.wpe,WPE.WPEState),(self.band,BAND.BandState),
               (self.stats,BAND.StatsState),(self.tracker_lpf,FRONT.LPFState),
               (self.stillness,STILL.State),(self.tune,TuneState))
        if any(not isinstance(x,t) for x,t in types): raise TypeError('invalid tuner-prefix state component')
        object.__setattr__(self,'last_adapt_time',R(self.last_adapt_time)); object.__setattr__(self,'time',R(self.time))
        object.__setattr__(self,'stage_time',R(self.stage_time)); object.__setattr__(self,'warmup_sec',R(self.warmup_sec))
        if self.last_adapt_time<0 or self.time<self.last_adapt_time: raise ValueError('invalid tuner-prefix clocks')
        if self.stage not in STAGES or self.stage_time<0 or self.warmup_sec<0: raise ValueError('invalid startup stage/clock')
        if not isinstance(self.pending,bool): raise TypeError('literal pending bit required')
        if not isinstance(self.sample_index,int) or self.sample_index<0: raise ValueError('nonnegative sample index required')


@dataclass(frozen=True)
class Result:
    state: State
    raw_sample: RAW.RawImuSample
    vertical: V.Result
    preupdate_wpe: BAND.WPEFrequencyView
    band: BAND.FrontendResult
    tracker_lpf: FRONT.LPFResult
    stillness: STILL.Result
    candidate: CAND.CandidateResult | None
    wpe: WPE.UpdateResult
    stage_before: str
    stage_after: str


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
         candidate_cfg:CAND.CandidateConfig|None=None,sigma_wave_sqrt=None,
         spectral:CAND.SpectralWitness|None=None,ema:CAND.EmaWitness|None=None):
    if not isinstance(state,State) or not isinstance(sample,RAW.RawImuSample):
        raise TypeError('finite tuner-prefix state and RawImuSample required')
    if state.pending:
        raise ValueError('pending tuner state must be committed at next IMU boundary before another sample')
    dt=R(dt)
    if dt<=0: raise ValueError('positive dt required')
    required_front=(band_cfg,stats_cfg,band_decay,variance_decay,tracker_lpf_decay,still_cfg)
    if any(x is None for x in required_front): raise TypeError('all represented frontend configs/witnesses required')

    vertical=RAW.vertical_step_from_raw(state.vertical,vertical_cfg,sample,dt=dt,
        accel_invnorm=accel_invnorm,quat_invnorm=quat_invnorm,seed=seed)
    view=FRONT.wpe_view(state.wpe,
        current_period=wpe_current_period if state.wpe.usable_period else None,
        current_frequency=wpe_current_frequency if state.wpe.usable_period else None)
    lpf=FRONT.tracker_lpf_step(state.tracker_lpf,vertical,decay=tracker_lpf_decay)
    still=STILL.step(state.stillness,still_cfg,a_vert_up_lp=lpf.output,dt=dt,
                     attenuation=still_attenuation)
    band=FRONT.band_step_from_vertical_view(state.band,state.stats,view,vertical,dt=dt,
        band_cfg=band_cfg,stats_cfg=stats_cfg,band_decay=band_decay,
        variance_decay=variance_decay,bench_noise_sigma=bench_noise_sigma,noise_sqrt=noise_sqrt)

    now=state.time+dt
    stage_clock=state.stage_time+dt
    stage_after=state.stage
    next_stage_clock=stage_clock
    cand=None

    # update_tuner() has already updated band/statistics above before entering
    # this switch, exactly as shipping does.
    if state.stage == 'Cold':
        if stage_clock >= state.warmup_sec:
            stage_after='TunerWarm'; next_stage_clock=F(0)
        if any(x is not None for x in (candidate_cfg,sigma_wave_sqrt,spectral,ema)):
            raise ValueError('Cold tuner branch returns before candidate operands are consumed')
    else:
        # For the represented valid-frequency sample, SeaStateAutoTuner::update
        # has just stored a positive bounded frequency, so isFreqReady is true.
        if state.stage == 'TunerWarm' and band.stats_state.var_ready and view.state.usable_period:
            stage_after='TunerReady'; next_stage_clock=F(0)
        if any(x is None for x in (candidate_cfg,sigma_wave_sqrt,spectral,ema)):
            raise TypeError('post-Cold tuner branch requires candidate witnesses')
        cand=BRIDGE.step(state.tune,band,still,candidate_cfg,sigma_wave_sqrt=sigma_wave_sqrt,
                         dt=dt,time=now,last_adapt_time=state.last_adapt_time,
                         spectral=spectral,ema=ema)

    # Current sample updates WPE after tuner/stage logic; its new usable latch is
    # therefore invisible to the TunerWarm promotion above until next sample.
    wpe=FRONT.wpe_step_from_vertical(state.wpe,wpe_cfg,vertical,dt=dt,decay=wpe_decay,
        moment_decay=wpe_moment_decay,period_witness=wpe_period_witness,
        log_witness=wpe_log_witness,current_period=wpe_current_period,
        current_frequency=wpe_current_frequency,post_output=wpe_post_output)

    tune_next=state.tune if cand is None else cand.tune_next
    last_adapt=state.last_adapt_time if cand is None else cand.last_adapt_time_after
    pending=False if cand is None else cand.pending_after
    nxt=State(vertical.state,wpe.state,band.band_state,band.stats_state,lpf.state,
              still.state,tune_next,last_adapt,pending,state.sample_index+1,now,
              stage_after,next_stage_clock,state.warmup_sec)
    return Result(nxt,sample,vertical,view,band,lpf,still,cand,wpe,state.stage,stage_after)


def readiness():
    return {
      'raw_IMU_to_private_vertical_same_packet':True,
      'same_vertical_band_tracker_LPF_and_later_WPE':True,
      'tuner_uses_WPE_state_from_sample_entry':True,
      'current_sample_WPE_update_cannot_feed_own_candidate':True,
      'tracker_free_tuner_stillness_projection':True,
      'Cold_frontend_update_then_early_return_materialized':True,
      'Cold_to_TunerWarm_warmup_transition_materialized':True,
      'TunerWarm_candidate_before_ready_materialized':True,
      'TunerWarm_to_TunerReady_requires_variance_and_preupdate_WPE_usable':True,
      'TunerReady_and_Live_candidate_recurrence_materialized':True,
      'TuneState_adapt_clock_and_pending_persist_across_samples':True,
      'pending_boundary_cannot_be_skipped':True,
      'dominant_frequency_tracker_absent_from_OU_tuner_prefix':True,
      'goLive_attitude_handoff_qualification_attached':False,
      'next_boundary_staged_commit_composed':False,
      'sensor_residual_source_bounds_attached':False,
      'transcendental_binary32_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
