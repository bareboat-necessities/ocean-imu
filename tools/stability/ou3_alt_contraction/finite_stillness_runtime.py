"""Finite StillnessAdapter recurrence used by the OU-III tuner.

The dominant-frequency tracker does not choose the OU wave frequency, but its
existing StillnessAdapter contributes one boolean/time state to the tuner: when
still, the wave-band variance is multiplied by exp(-still_time/1s).  This module
retains that narrow path explicitly.

Tracker output and leveled tracker-input acceleration remain upstream inputs;
exp binary32 ancestry remains an explicit witness obligation.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P


def R(x): return P.rational(x)
def clamp(x,lo,hi): return min(max(x,lo),hi)


@dataclass(frozen=True)
class Config:
    gravity:F=F(980665,100000)
    energy_alpha:F=F(1,20)
    energy_thresh:F=F(8,10000)
    still_thresh:F=F(2)
    relax_tau:F=F(1)
    target_freq:F=F(1,5)
    def __post_init__(self):
        for n in ('gravity','energy_alpha','energy_thresh','still_thresh','relax_tau','target_freq'):
            object.__setattr__(self,n,R(getattr(self,n)))
        if self.gravity<=0 or not 0<self.energy_alpha<=1 or self.energy_thresh<0 or self.still_thresh<0 or self.relax_tau<=0 or self.target_freq<=0:
            raise ValueError('invalid stillness configuration')


@dataclass(frozen=True)
class State:
    energy_ema:F=F(0)
    still_time:F=F(0)
    freq_init:bool=False
    freq_state:F=F(3,10)
    last_is_still:bool=False
    def __post_init__(self):
        object.__setattr__(self,'energy_ema',R(self.energy_ema)); object.__setattr__(self,'still_time',R(self.still_time)); object.__setattr__(self,'freq_state',R(self.freq_state))
        if not isinstance(self.freq_init,bool) or not isinstance(self.last_is_still,bool): raise TypeError('literal stillness flags required')
        if self.energy_ema<0 or self.still_time<0: raise ValueError('invalid stillness state')


@dataclass(frozen=True)
class RelaxWitness:
    decay:F
    def __post_init__(self):
        d=R(self.decay)
        if not 0<d<=1: raise ValueError('relax exp decay must be in (0,1]')
        object.__setattr__(self,'decay',d)


@dataclass(frozen=True)
class AttenuationWitness:
    value:F
    def __post_init__(self):
        v=R(self.value)
        if not 0<v<=1: raise ValueError('stillness attenuation must lie in (0,1]')
        object.__setattr__(self,'value',v)


@dataclass(frozen=True)
class Result:
    state:State
    tracker_frequency_out:F
    variance_attenuation:F


def step(s:State,cfg:Config,*,a_vert_up_lp,dt,tracker_frequency,
         relax:RelaxWitness|None=None,attenuation:AttenuationWitness|None=None):
    a,dt,f=map(R,(a_vert_up_lp,dt,tracker_frequency))
    if dt<=0: raise ValueError('positive stillness dt required')
    fs=f if (not s.freq_init) else s.freq_state
    a_norm=a/cfg.gravity; inst=a_norm*a_norm
    energy=(1-cfg.energy_alpha)*s.energy_ema+cfg.energy_alpha*inst
    is_still=energy<cfg.energy_thresh
    if is_still:
        st=min(F(60),s.still_time+dt)
        if st>cfg.still_thresh:
            if relax is None: raise ValueError('relaxing stillness branch requires exp witness')
            alpha=1-relax.decay
            fs=fs+alpha*(cfg.target_freq-fs)
        else:
            if relax is not None: raise ValueError('pre-hold stillness branch consumes no relax witness')
            fs=f
        if attenuation is None: raise ValueError('still variance branch requires attenuation witness')
        att=attenuation.value
    else:
        if relax is not None or attenuation is not None: raise ValueError('moving branch consumes no stillness witnesses')
        st=F(0); fs=f; att=F(1)
    return Result(State(energy,st,True,fs,is_still),fs,att)


def readiness():
    return {
      'energy_ema_and_stillness_threshold_materialized':True,
      'still_time_hold_and_60s_cap_materialized':True,
      'tracker_frequency_relaxation_materialized':True,
      'variance_attenuation_branch_materialized':True,
      'tracker_frequency_and_vertical_input_same_history_attached':False,
      'relax_and_attenuation_exp_binary32_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
