"""Tracker-free projection of StillnessAdapter relevant to OU-III tuning.

Shipping's OU tuner reads only ``isStill()`` and ``getStillTime()`` from
StillnessAdapter.  Those coordinates are functions of the LPF vertical sample,
dt, and prior energy/still-time state.  The dominant-frequency tracker output
changes only StillnessAdapter's frequency state/output, which is used by the
separate direction/reporting path and is not read by ``update_tuner``.

Therefore the ALT stability word may retain the exact adaptation-relevant
projection without carrying any frequency-tracker algorithm state.  This is an
exact state projection, not an assumption that the tracker behaves nicely.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_stillness_runtime as FULL


def R(x): return P.rational(x)


@dataclass(frozen=True)
class State:
    energy_ema:F=F(0)
    still_time:F=F(0)
    last_is_still:bool=False
    def __post_init__(self):
        object.__setattr__(self,'energy_ema',R(self.energy_ema))
        object.__setattr__(self,'still_time',R(self.still_time))
        if self.energy_ema<0 or self.still_time<0: raise ValueError('invalid projected stillness state')
        if not isinstance(self.last_is_still,bool): raise TypeError('literal stillness branch required')


@dataclass(frozen=True)
class Result:
    state:State
    variance_attenuation:F


def project(full:FULL.State):
    if not isinstance(full,FULL.State): raise TypeError('full StillnessAdapter state required')
    return State(full.energy_ema,full.still_time,full.last_is_still)


def step(s:State,cfg:FULL.Config,*,a_vert_up_lp,dt,
         attenuation:FULL.AttenuationWitness|None=None):
    if not isinstance(s,State) or not isinstance(cfg,FULL.Config):
        raise TypeError('projected state and shipping stillness config required')
    a,dt=R(a_vert_up_lp),R(dt)
    if dt<=0: raise ValueError('positive stillness dt required')
    inst=(a/cfg.gravity)**2
    energy=(1-cfg.energy_alpha)*s.energy_ema+cfg.energy_alpha*inst
    is_still=energy<cfg.energy_thresh
    if is_still:
        st=min(F(60),s.still_time+dt)
        if attenuation is None: raise ValueError('still tuner branch requires attenuation witness')
        att=attenuation.value
    else:
        if attenuation is not None: raise ValueError('moving tuner branch consumes no attenuation witness')
        st=F(0); att=F(1)
    return Result(State(energy,st,is_still),att)


def assert_matches_full(projected:Result,full:FULL.Result):
    """Check exact adaptation-relevant equality for any tracker-frequency path."""
    if not isinstance(projected,Result) or not isinstance(full,FULL.Result):
        raise TypeError('projected and full stillness results required')
    if projected.state != project(full.state):
        raise ValueError('tracker-free stillness projection detached from shipping state')
    if projected.variance_attenuation != full.variance_attenuation:
        raise ValueError('variance attenuation detached from shipping stillness branch')
    return True


def readiness():
    return {
      'tuner_reads_only_isStill_and_stillTime':True,
      'energy_still_time_projection_independent_of_tracker_frequency':True,
      'tracker_algorithm_not_required_for_OU_tuner_stability_word':True,
      'stillness_attenuation_exp_binary32_attached':False,
      'direction_tracker_separate_from_filter_stability_word':True,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
