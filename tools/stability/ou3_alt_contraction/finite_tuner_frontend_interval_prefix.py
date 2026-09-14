"""Tracker-free tuner frontend with an exact-real interval SpectralMSE state.

This is the theorem-facing counterpart of ``finite_tuner_frontend_prefix`` for
ordinary shipping SpectralMSE cells whose sqrt/pow roots are irrational.  It
reuses the same lower runtime primitives and literal within-sample ordering:

  private Mahony successor
    -> sample-entry WPE read-only view
    -> tracker LPF / stillness + adaptive band/statistics
    -> startup stage + tuner candidate
    -> current-sample WPE update for the NEXT sample.

Only the tuner state type changes: tau and sigma remain exact rationals while
R_S is carried as the rigorous exact-real interval from
``finite_tuner_candidate_interval``.  No representative R_S is selected.
Consequently this module closes real-arithmetic frontend composition but does
NOT yet authorize goLive/application to the MEKF; the actual binary32
sqrt/pow/exp result must first be related to this interval.
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
from tools.stability.ou3_alt_contraction import finite_tuner_candidate_interval as ICAND
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P


def R(x): return P.rational(x)
STAGES=('Cold','TunerWarm','TunerReady','Live')
QUALIFICATION='OU3_ALT_INTERVAL_TUNER_FRONTEND_PREFIX_V1'


@dataclass(frozen=True)
class State:
    vertical:V.State
    wpe:WPE.WPEState
    band:BAND.BandState
    stats:BAND.StatsState
    tracker_lpf:FRONT.LPFState
    stillness:STILL.State
    tune:ICAND.IntervalTuneState
    last_adapt_time:F=F(0)
    pending:bool=False
    sample_index:int=0
    time:F=F(0)
    stage:str='Live'
    stage_time:F=F(0)
    warmup_sec:F=F(5)
    def __post_init__(self):
        types=((self.vertical,V.State),(self.wpe,WPE.WPEState),(self.band,BAND.BandState),
               (self.stats,BAND.StatsState),(self.tracker_lpf,FRONT.LPFState),
               (self.stillness,STILL.State),(self.tune,ICAND.IntervalTuneState))
        if any(not isinstance(x,t) for x,t in types): raise TypeError('invalid interval tuner-prefix state component')
        for n in ('last_adapt_time','time','stage_time','warmup_sec'):
            object.__setattr__(self,n,R(getattr(self,n)))
        if self.last_adapt_time<0 or self.time<self.last_adapt_time or self.stage_time<0 or self.warmup_sec<0:
            raise ValueError('invalid interval tuner clocks')
        if self.stage not in STAGES: raise ValueError('invalid startup stage')
        if not isinstance(self.pending,bool): raise TypeError('literal pending bit required')
        if not isinstance(self.sample_index,int) or self.sample_index<0: raise ValueError('nonnegative sample index required')


@dataclass(frozen=True)
class Result:
    state:State
    raw_sample:RAW.RawImuSample|RAW.GuardedImuSample
    vertical:V.Result
    preupdate_wpe:BAND.WPEFrequencyView
    band:BAND.FrontendResult
    tracker_lpf:FRONT.LPFResult
    stillness:STILL.Result
    candidate:ICAND.CandidateIntervalResult|None
    wpe:WPE.UpdateResult
    stage_before:str
    stage_after:str


def continue_after_vertical(state:State,sample,vertical:V.Result,*,dt,
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
                            ema:CAND.EmaWitness|None=None,spectral_bits:int=96):
    if not isinstance(state,State) or not isinstance(sample,(RAW.RawImuSample,RAW.GuardedImuSample)) or not isinstance(vertical,V.Result):
        raise TypeError('interval state, raw/guarded sample and private vertical successor required')
    if state.pending: raise ValueError('pending interval tuner state must be committed before another sample')
    dt=R(dt)
    if dt<=0: raise ValueError('positive dt required')
    if any(x is None for x in (band_cfg,stats_cfg,band_decay,variance_decay,tracker_lpf_decay,still_cfg)):
        raise TypeError('all represented frontend configs/witnesses required')

    view=FRONT.wpe_view(state.wpe,
        current_period=wpe_current_period if state.wpe.usable_period else None,
        current_frequency=wpe_current_frequency if state.wpe.usable_period else None)
    lpf=FRONT.tracker_lpf_step(state.tracker_lpf,vertical,decay=tracker_lpf_decay)
    still=STILL.step(state.stillness,still_cfg,a_vert_up_lp=lpf.output,dt=dt,
                     attenuation=still_attenuation)
    band=FRONT.band_step_from_vertical_view(state.band,state.stats,view,vertical,dt=dt,
        band_cfg=band_cfg,stats_cfg=stats_cfg,band_decay=band_decay,
        variance_decay=variance_decay,bench_noise_sigma=bench_noise_sigma,noise_sqrt=noise_sqrt)

    now=state.time+dt; stage_clock=state.stage_time+dt
    stage_after=state.stage; next_stage_clock=stage_clock; cand=None
    if state.stage=='Cold':
        if stage_clock>=state.warmup_sec:
            stage_after='TunerWarm'; next_stage_clock=F(0)
        if any(x is not None for x in (candidate_cfg,sigma_wave_sqrt,ema)):
            raise ValueError('Cold tuner branch returns before interval candidate operands are consumed')
    else:
        if state.stage=='TunerWarm' and band.stats_state.var_ready and view.state.usable_period:
            stage_after='TunerReady'; next_stage_clock=F(0)
        if any(x is None for x in (candidate_cfg,sigma_wave_sqrt,ema)):
            raise TypeError('post-Cold interval tuner branch requires candidate operands')
        cand=BRIDGE.step_interval(state.tune,band,still,candidate_cfg,
            sigma_wave_sqrt=sigma_wave_sqrt,dt=dt,time=now,last_adapt_time=state.last_adapt_time,
            ema=ema,bits=spectral_bits)

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


def step(state:State,sample,*,dt,vertical_cfg:V.Config,accel_invnorm:V.InvSqrtWitness|None,
         quat_invnorm:V.InvSqrtWitness|None,seed:V.SeedWitness|None,**suffix_kwargs):
    if not isinstance(state,State) or not isinstance(sample,(RAW.RawImuSample,RAW.GuardedImuSample)):
        raise TypeError('interval tuner-prefix state and raw/guarded sample required')
    if state.pending: raise ValueError('pending interval tuner state must be committed before another sample')
    dt=R(dt)
    if isinstance(sample,RAW.GuardedImuSample):
        vertical=RAW.vertical_step_from_guarded(state.vertical,vertical_cfg,sample,dt=dt,
            accel_invnorm=accel_invnorm,quat_invnorm=quat_invnorm,seed=seed)
    else:
        vertical=RAW.vertical_step_from_raw(state.vertical,vertical_cfg,sample,dt=dt,
            accel_invnorm=accel_invnorm,quat_invnorm=quat_invnorm,seed=seed)
    return continue_after_vertical(state,sample,vertical,dt=dt,**suffix_kwargs)


def readiness():
    b=BRIDGE.readiness()
    return {
      'qualification':QUALIFICATION,
      'same_private_vertical_WPE_band_stillness_order_as_shipping_frontend':True,
      'sample_entry_WPE_view_precedes_current_sample_WPE_update':True,
      'Cold_frontend_and_stage_transition_materialized':True,
      'postCold_general_SpectralMSE_interval_candidate_composed':b['general_SpectralMSE_interval_candidate_projection_available'],
      'tau_sigma_frontend_state_remains_exact':True,
      'RS_frontend_state_is_rigorous_exact_real_interval':True,
      'actual_0p2Hz_prior_postCold_candidate_representable':True,
      'point_RS_representative_selected':False,
      'goLive_interval_RS_to_actual_MEKF_commit_closed':False,
      'binary32_sqrt_pow_exp_correspondence_closed':False,
      'guard_state_persisted_here':False,
      'source_uniform_complete_startup_reachability_closed':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
