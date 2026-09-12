"""Default unweighted MagAutoTuner accepted-window recurrence.

Shipping's default ``enable_quality_weighting`` is false, so every sample that
has already passed MagAutoTuner's acceptance gates contributes weight one. This
module closes the resulting sufficient-statistic recurrence and the minimum
count/window finalize gate. It deliberately does not decide whether a raw
sample is accepted; those finite sensor/quality gates remain a separate source
obligation.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M


def R(x): return M.rational(x)

@dataclass(frozen=True)
class Config:
    min_samples:int=250
    min_window:F=F(10)
    max_window:F=F(0)
    sample_dt:F=F(1,200)
    quality_weighting:bool=False
    def __post_init__(self):
        if not isinstance(self.min_samples,int): raise TypeError('integer min_samples required')
        if self.min_samples<0: raise ValueError('nonnegative min_samples required')
        for n in ('min_window','max_window','sample_dt'):
            v=R(getattr(self,n))
            if v<0: raise ValueError('nonnegative magnetic window times required')
            object.__setattr__(self,n,v)
        if not isinstance(self.quality_weighting,bool): raise TypeError('literal quality-weighting branch required')

@dataclass(frozen=True)
class State:
    world_sum:tuple=(F(0),F(0),F(0))
    norm_sum:F=F(0)
    accepted_count:int=0
    accepted_window:F=F(0)
    weight_sum:F=F(0)
    def __post_init__(self):
        object.__setattr__(self,'world_sum',tuple(M.vec(self.world_sum,3)))
        for n in ('norm_sum','accepted_window','weight_sum'):
            v=R(getattr(self,n))
            if v<0: raise ValueError('nonnegative accumulator scalar required')
            object.__setattr__(self,n,v)
        if not isinstance(self.accepted_count,int) or self.accepted_count<0:
            raise ValueError('nonnegative accepted count required')

@dataclass(frozen=True)
class AcceptedResult:
    state:State
    finalize_due:bool


def accepted_sample(state:State,cfg:Config,*,world_sample,world_norm,dt):
    if not isinstance(state,State) or not isinstance(cfg,Config):
        raise TypeError('mag accumulator state/config required')
    if cfg.quality_weighting:
        raise NotImplementedError('weighted MagAutoTuner branch is not the deployed default proof branch')
    m=tuple(M.vec(world_sample,3)); n=R(world_norm); d=R(dt)
    if n<0: raise ValueError('nonnegative magnetic norm required')
    # Shipping uses supplied positive finite dt, otherwise configured sample_dt.
    d_use=d if d>0 else cfg.sample_dt
    nxt=State(tuple(state.world_sum[i]+m[i] for i in range(3)),
              state.norm_sum+n,state.accepted_count+1,
              state.accepted_window+d_use,state.weight_sum+1)
    count_ok=nxt.accepted_count>=max(1,cfg.min_samples)
    timed_out=cfg.max_window>0 and nxt.accepted_window>=cfg.max_window
    window_ok=(cfg.min_window<=0 or nxt.accepted_window>=cfg.min_window or timed_out)
    return AcceptedResult(nxt,count_ok and window_ok)


def mean(state:State):
    if not isinstance(state,State): raise TypeError('mag accumulator state required')
    if state.accepted_count<=0 or state.weight_sum<=F(1,10**6):
        raise ValueError('MagAutoTuner mean unavailable before positive accepted weight')
    return tuple(x/state.weight_sum for x in state.world_sum)


def readiness():
    return {
      'default_quality_weighting_disabled_branch_materialized':True,
      'accepted_world_sum_count_time_weight_recurrence_materialized':True,
      'dt_fallback_to_configured_sample_period_materialized':True,
      'min_count_min_window_and_timeout_finalize_gate_materialized':True,
      'weighted_mean_from_same_accumulator_materialized':True,
      'raw_mag_sample_acceptance_gates_attached':False,
      'world_sample_tilt_quaternion_binary32_attached':False,
      'norm_and_dt_finite_binary32_attached':False,
      'optional_quality_weighted_branch_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
