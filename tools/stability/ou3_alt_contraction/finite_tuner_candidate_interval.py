"""Interval-valued exact-real tuner candidate for general SpectralMSE cells.

``finite_tuner_candidate`` is exact over rationals but its legacy
``SpectralWitness`` can only express cells where sqrt(T_S) and u^(6/7) happen
to be rational.  Real shipping cells generally do not.  This module keeps the
same exact frequency/variance/tau/sigma target and cadence algebra, replaces
only the irrational SpectralMSE target by the exact rational enclosure from
``finite_tuner_spectral_real_enclosure``, and propagates that enclosure through
the literal R_S EMA.

Thus it is a rigorous relation for the real-arithmetic shipping tuner, not a
sampled approximation.  It does not choose a representative R_S, mutate the
shipping filter, or close binary32 sqrt/pow/exp correspondence.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as S
from tools.stability.ou3_alt_contraction.finite_tuner_commit import TuneState, clamp

QUALIFICATION='OU3_ALT_TUNER_INTERVAL_CANDIDATE_V1'


@dataclass(frozen=True)
class IntervalTuneState:
    tau_applied:F
    sigma_applied:F
    RS_lo:F
    RS_hi:F
    def __post_init__(self):
        vals=tuple(P.rational(x) for x in (self.tau_applied,self.sigma_applied,self.RS_lo,self.RS_hi))
        for n,v in zip(('tau_applied','sigma_applied','RS_lo','RS_hi'),vals): object.__setattr__(self,n,v)
        if self.tau_applied<=0 or self.sigma_applied<0 or self.RS_lo>self.RS_hi:
            raise ValueError('invalid interval tuner state')

    @classmethod
    def point(cls,tune:TuneState):
        if not isinstance(tune,TuneState): raise TypeError('TuneState required')
        return cls(tune.tau_applied,tune.sigma_applied,tune.RS_applied,tune.RS_applied)


@dataclass(frozen=True)
class CandidateIntervalResult:
    target:C.TargetState
    spectral:S.Enclosure
    tune_next:IntervalTuneState
    pending_after:bool
    last_adapt_time_after:F
    adapt_horizon:F
    RS_horizon:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.target,C.TargetState) or not isinstance(self.spectral,S.Enclosure) or not isinstance(self.tune_next,IntervalTuneState):
            raise TypeError('target, spectral enclosure and interval tune state required')
        if self.spectral.target!=self.target: raise ValueError('spectral enclosure detached from same tuner target')
        if not isinstance(self.pending_after,bool): raise TypeError('literal pending bit required')
        object.__setattr__(self,'last_adapt_time_after',P.rational(self.last_adapt_time_after))
        object.__setattr__(self,'adapt_horizon',P.rational(self.adapt_horizon))
        object.__setattr__(self,'RS_horizon',P.rational(self.RS_horizon))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong interval candidate qualification')


def _ema_interval(prev_lo,prev_hi,target_lo,target_hi,alpha):
    prev_lo,prev_hi,target_lo,target_hi,alpha=map(P.rational,(prev_lo,prev_hi,target_lo,target_hi,alpha))
    if not 0<=alpha<=1 or prev_lo>prev_hi or target_lo>target_hi:
        raise ValueError('invalid nonexpansive EMA interval')
    # Coefficients 1-alpha and alpha are nonnegative, so the box image is exact.
    return ((1-alpha)*prev_lo+alpha*target_lo,
            (1-alpha)*prev_hi+alpha*target_hi)


def step(previous:IntervalTuneState,sample:C.WaveBandSample,cfg:C.CandidateConfig,*,
         dt,time,last_adapt_time,ema:C.EmaWitness,bits:int=S.DEFAULT_BITS):
    if not isinstance(previous,IntervalTuneState): raise TypeError('IntervalTuneState required')
    if not isinstance(sample,C.WaveBandSample) or not isinstance(cfg,C.CandidateConfig) or not isinstance(ema,C.EmaWitness):
        raise TypeError('sample/config/EMA witness required')
    dt,time,last_adapt_time=map(P.rational,(dt,time,last_adapt_time))
    if dt<=0 or time<last_adapt_time: raise ValueError('positive dt and monotone tuner time required')
    target=C.targets(sample,cfg)
    spectral=S.enclose(cfg,target,bits=bits)

    sea_time=F(1,2)/target.frequency
    if cfg.adapt_tau_sea_periods>0:
        safe=clamp(sea_time,C.EMA_SCALE_MIN,C.EMA_SCALE_MAX)
        adapt_h=C._clamp_horizon(cfg.adapt_tau_sea_periods*safe,dt)
    else:
        adapt_h=cfg.adapt_tau_sec
    alpha=1-ema.decay_tau_sigma
    tau_next=previous.tau_applied+alpha*(target.tau_target-previous.tau_applied)
    sigma_next=previous.sigma_applied+alpha*(target.sigma_target-previous.sigma_applied)

    safe_tau=clamp(target.tau_target,C.EMA_SCALE_MIN,C.EMA_SCALE_MAX)
    rs_h=C._clamp_horizon(cfg.adapt_RS_mult*safe_tau,dt)
    alpha_rs=1-ema.decay_RS
    rs_lo,rs_hi=_ema_interval(previous.RS_lo,previous.RS_hi,
                              spectral.target_RS_lo,spectral.target_RS_hi,alpha_rs)
    nxt=IntervalTuneState(tau_next,sigma_next,rs_lo,rs_hi)
    fire=(time-last_adapt_time)>cfg.adapt_every_sec
    return CandidateIntervalResult(target,spectral,nxt,fire,time if fire else last_adapt_time,adapt_h,rs_h)


def contains_legacy_exact(previous:TuneState,sample:C.WaveBandSample,cfg:C.CandidateConfig,*,
                          dt,time,last_adapt_time,spectral:C.SpectralWitness,ema:C.EmaWitness,bits:int=S.DEFAULT_BITS):
    """When a rational exact cell exists, prove the old exact candidate lies inside this relation."""
    exact=C.step(previous,sample,cfg,dt=dt,time=time,last_adapt_time=last_adapt_time,spectral=spectral,ema=ema)
    box=step(IntervalTuneState.point(previous),sample,cfg,dt=dt,time=time,last_adapt_time=last_adapt_time,ema=ema,bits=bits)
    return (exact.frequency==box.target.frequency and exact.tau_target==box.target.tau_target and
            exact.sigma_target==box.target.sigma_target and
            exact.tune_next.tau_applied==box.tune_next.tau_applied and
            exact.tune_next.sigma_applied==box.tune_next.sigma_applied and
            box.spectral.target_RS_lo<=exact.RS_target<=box.spectral.target_RS_hi and
            box.tune_next.RS_lo<=exact.tune_next.RS_applied<=box.tune_next.RS_hi and
            exact.pending_after==box.pending_after and exact.last_adapt_time_after==box.last_adapt_time_after)


def readiness():
    s=S.readiness()
    return {
      'qualification':QUALIFICATION,
      'same_frequency_variance_tau_sigma_targets_as_exact_tuner':True,
      'SpectralMSE_irrational_target_enclosure_consumed':s['SpectralMSE_real_RS_target_interval_propagated_monotonically'],
      'RS_interval_propagated_through_literal_nonexpansive_EMA':True,
      'tau_sigma_EMA_and_commit_pending_logic_remain_exact':True,
      'legacy_rational_exact_candidate_embeds_in_interval_relation':True,
      'actual_0p2Hz_prior_cell_representable_without_fake_rational_roots':True,
      'binary32_sqrt_pow_correspondence_closed':False,
      'binary32_tuner_exp_correspondence_closed':False,
      'startup_frontend_composed_with_interval_RS_state':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
