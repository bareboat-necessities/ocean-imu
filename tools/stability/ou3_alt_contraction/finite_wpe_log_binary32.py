"""Two coherent binary32 WPE log-period tracks for ALT deployment proof.

Shipping stores ``log_period_sec_`` as float.  On every valid period update it
computes ``log_raw = std::log(raw_period)``.  After initialization, the default
positive smoothing branch evaluates ``std::exp(log_period_sec_)`` to form the
horizon, evaluates ``std::exp(-dt/horizon)`` for the EMA decay, and executes

    log_period_sec_ += alpha * (log_raw - log_period_sec_).

The compiler FP-contraction choice is global, so ALT carries exactly two machine
histories: one always separate multiply/add, one always FMA.  Inputs to those
histories are allowed to differ because upstream float arithmetic and libm may
already have diverged between the two global compiler modes.  This avoids
recombining incompatible machine histories.

This module is topology/arithmetic correspondence only.  The supplied binary32
``log``/``exp`` results are retained as target-libm witnesses; correctness of
those libm values and the binary32 raw-period production remain open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B

SOURCE=Path(__file__).resolve().parents[3]/'src/tuner/WavePeriodEstimator.h'
DT=B.rn32(F(1,200)); LOG_SMOOTH_PERIODS=B.rn32(F(1,20))
HORIZON_MIN=B.rn32(F(1,20)); HORIZON_MAX=B.rn32(35)
QUALIFICATION='OU3_ALT_WPE_LOG_BINARY32_DUAL_COMPILER_V1'


def clamp(x,lo,hi): return min(max(x,lo),hi)


@dataclass(frozen=True)
class Track:
    log_period:F|None=None
    valid_updates:int=0
    def __post_init__(self):
        if self.log_period is not None:
            q=F(self.log_period)
            if not B.is_binary32(q): raise ValueError('WPE machine log-period must be binary32')
            object.__setattr__(self,'log_period',q)
        if not isinstance(self.valid_updates,int) or self.valid_updates<0:
            raise ValueError('nonnegative WPE valid-update count required')
        if self.log_period is None and self.valid_updates:
            raise ValueError('uninitialized WPE log track cannot have valid updates')


@dataclass(frozen=True)
class State:
    separate:Track=Track()
    fma:Track=Track()
    samples:int=0
    def __post_init__(self):
        if not isinstance(self.separate,Track) or not isinstance(self.fma,Track):
            raise TypeError('two WPE compiler tracks required')
        if not isinstance(self.samples,int) or self.samples<0:
            raise ValueError('nonnegative WPE sample count required')


@dataclass(frozen=True)
class SmoothWitness:
    log_raw:F
    sea_period_exp:F
    decay_exp:F
    def __post_init__(self):
        vals=tuple(F(x) for x in (self.log_raw,self.sea_period_exp,self.decay_exp))
        if not all(B.is_binary32(x) for x in vals):
            raise ValueError('WPE smoothing libm witnesses must be binary32')
        lr,sea,d=vals
        if sea<=0 or not 0<d<=1: raise ValueError('positive WPE exp witnesses required')
        object.__setattr__(self,'log_raw',lr); object.__setattr__(self,'sea_period_exp',sea)
        object.__setattr__(self,'decay_exp',d)


@dataclass(frozen=True)
class InitWitness:
    log_raw:F
    def __post_init__(self):
        q=F(self.log_raw)
        if not B.is_binary32(q): raise ValueError('WPE std::log result witness must be binary32')
        object.__setattr__(self,'log_raw',q)


@dataclass(frozen=True)
class TrackStep:
    before:Track
    after:Track
    log_raw:F
    sea_period_exp:F|None
    decay_exp:F|None
    horizon:F|None
    alpha:F|None
    contracted:bool


@dataclass(frozen=True)
class StepResult:
    state:State
    separate:TrackStep|None
    fma:TrackStep|None
    produced_period:bool


def initial(): return State()


def hold_sample(state:State):
    """One WPE sample that produces no valid period/log update."""
    if not isinstance(state,State): raise TypeError('WPE binary32 State required')
    return StepResult(State(state.separate,state.fma,state.samples+1),None,None,False)


def _init(track:Track,w:InitWitness,contracted:bool):
    if track.log_period is not None:
        raise ValueError('initialized WPE track cannot take first-log branch again')
    after=Track(w.log_raw,track.valid_updates+1)
    return TrackStep(track,after,w.log_raw,None,None,None,None,contracted)


def _smooth(track:Track,w:SmoothWitness,contracted:bool):
    if track.log_period is None: raise ValueError('WPE smoothing requires initialized log-period track')
    requested=B.mul(LOG_SMOOTH_PERIODS,w.sea_period_exp)
    horizon=clamp(requested,HORIZON_MIN,HORIZON_MAX)
    x=B.div(DT,horizon)
    # Numerical correctness of exp remains open, but its branch/domain is hard.
    if not 0<w.decay_exp<=1: raise ValueError('WPE log decay must be in (0,1]')
    alpha=B.sub(B.rn32(1),w.decay_exp)
    nxt=B.ema(track.log_period,w.log_raw,alpha,contracted=contracted)
    after=Track(nxt,track.valid_updates+1)
    return TrackStep(track,after,w.log_raw,w.sea_period_exp,w.decay_exp,horizon,alpha,contracted)


def first_valid(state:State,*,separate:InitWitness,fma:InitWitness):
    if not isinstance(state,State): raise TypeError('WPE binary32 State required')
    if not isinstance(separate,InitWitness) or not isinstance(fma,InitWitness):
        raise TypeError('mode-specific WPE log witnesses required')
    ss=_init(state.separate,separate,False); ff=_init(state.fma,fma,True)
    return StepResult(State(ss.after,ff.after,state.samples+1),ss,ff,True)


def smooth_valid(state:State,*,separate:SmoothWitness,fma:SmoothWitness):
    if not isinstance(state,State): raise TypeError('WPE binary32 State required')
    if not isinstance(separate,SmoothWitness) or not isinstance(fma,SmoothWitness):
        raise TypeError('mode-specific WPE smoothing witnesses required')
    ss=_smooth(state.separate,separate,False); ff=_smooth(state.fma,fma,True)
    return StepResult(State(ss.after,ff.after,state.samples+1),ss,ff,True)


def _source_shape_matches():
    s=SOURCE.read_text()
    return all(n in s for n in (
      'const float log_raw = std::log(raw_period_sec);',
      'log_period_sec_ = log_raw;',
      'const float sea_period = std::exp(log_period_sec_);',
      'const float requested = log_smoothing_periods_ * sea_period;',
      'const float alpha = 1.0f - std::exp(-dt_sec / horizon);',
      'log_period_sec_ += alpha * (log_raw - log_period_sec_);',
      'float log_period_sec_ = NAN;'))


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_WPE_log_update_source_shape_matches':_source_shape_matches(),
      'global_separate_and_FMA_WPE_log_tracks_persist_without_branch_explosion':True,
      'first_valid_log_assignment_materialized_per_compiler_track':True,
      'positive_default_log_smoothing_binary32_graph_materialized_per_track':True,
      'mode_specific_log_raw_and_exp_results_retained_without_cross_track_identity':True,
      'nonproducing_WPE_sample_preserves_log_tracks':True,
      'WPE_log_std_log_target_libm_correspondence_closed':False,
      'WPE_log_exp_target_libm_correspondence_closed':False,
      'WPE_raw_period_binary32_production_closed':False,
      'global_compiler_track_coherence_includes_WPE_log_state':True,
      'complete_word_finite_identity':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
