"""Finite AccelVibrationGuard runtime for the ALT shipping word.

Shipping conditions the raw body-frame accelerometer BEFORE private Mahony,
MEKF acceleration, and the tilt watchdog.  This module is the literal finite
real-arithmetic state machine for the valid-input/default-shape guard:

* optional exact identity when disabled;
* first-sample seeding with exact identity output;
* configurable 1..4-pole low-pass cascade;
* two-pole detector high-pass cascade;
* detector mean-square EMA and RMS witness;
* engagement target, slew and 1e-4 rail snapping;
* exact blend from raw acceleration to low-passed acceleration;
* excessRms() used by the later R_acc inflation branch.

The exp/sqrt evaluations are explicit witnesses.  Their deployed binary32
ancestry and source bounds remain open; no numerical approximation is silently
substituted.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M


def R(x): return M.rational(x)
def V(x): return tuple(M.vec(x,3))
def add(a,b): return tuple(a[i]+b[i] for i in range(3))
def sub(a,b): return tuple(a[i]-b[i] for i in range(3))
def scale(s,a): s=R(s); return tuple(s*a[i] for i in range(3))
def sq(a): return tuple(x*x for x in a)
def clamp(x,lo,hi): return min(max(x,lo),hi)
ZERO=(F(0),F(0),F(0))
WEIGHT_EPS=F(1,10000)


@dataclass(frozen=True)
class Config:
    cutoff_hz:F=F(14)
    detect_hz:F=F(25)
    engage_lo:F=F(3,100)
    engage_hi:F=F(8,100)
    slew_tau:F=F(5)
    poles:int=2
    removed_rms_hz:F=F(1,20)
    def __post_init__(self):
        for n in ('cutoff_hz','detect_hz','engage_lo','engage_hi','slew_tau','removed_rms_hz'):
            object.__setattr__(self,n,R(getattr(self,n)))
        if self.cutoff_hz<0 or self.detect_hz<=0 or self.engage_lo<0 or self.engage_hi<0 or self.slew_tau<=0 or self.removed_rms_hz<=0:
            raise ValueError('invalid vibration-guard configuration')
        if not isinstance(self.poles,int) or not 1<=self.poles<=4:
            raise ValueError('shipping guard poles must be 1..4')
    @property
    def enabled(self): return self.cutoff_hz>0


@dataclass(frozen=True)
class State:
    stages:tuple=(ZERO,ZERO,ZERO,ZERO)
    detect_stages:tuple=(ZERO,ZERO)
    removed_ms:tuple=ZERO
    weight:F=F(0)
    initialized:bool=False
    def __post_init__(self):
        if len(self.stages)!=4 or len(self.detect_stages)!=2: raise ValueError('guard state shape mismatch')
        object.__setattr__(self,'stages',tuple(V(x) for x in self.stages))
        object.__setattr__(self,'detect_stages',tuple(V(x) for x in self.detect_stages))
        rms=V(self.removed_ms)
        if any(x<0 for x in rms): raise ValueError('removed mean-square state must be nonnegative')
        object.__setattr__(self,'removed_ms',rms)
        w=R(self.weight)
        if not 0<=w<=1: raise ValueError('guard engagement must lie in [0,1]')
        object.__setattr__(self,'weight',w)
        if not isinstance(self.initialized,bool): raise TypeError('literal initialization branch required')


@dataclass(frozen=True)
class DecayWitness:
    alpha:F
    gamma:F
    beta:F
    slew:F
    def __post_init__(self):
        for n in ('alpha','gamma','beta','slew'):
            v=R(getattr(self,n)); object.__setattr__(self,n,v)
            if not 0<=v<=1: raise ValueError('guard decay/slew witness outside [0,1]')


@dataclass(frozen=True)
class RmsWitness:
    value:F
    def __post_init__(self):
        v=R(self.value)
        if v<0: raise ValueError('RMS witness must be nonnegative')
        object.__setattr__(self,'value',v)


@dataclass(frozen=True)
class Result:
    state:State
    output:tuple
    low_passed:tuple
    detector_high_passed:tuple
    removed_rms:F
    excess_rms:F
    target:F
    raw_weight:F


def _rms(removed_ms,rms:RmsWitness):
    s=sum(removed_ms)
    if rms.value*rms.value != s:
        raise ValueError('guard RMS witness detached from same detector mean-square state')
    return rms.value


def step(state:State,cfg:Config,acc,dt,*,decay:DecayWitness|None=None,rms:RmsWitness|None=None):
    if not isinstance(state,State) or not isinstance(cfg,Config): raise TypeError('guard state/config required')
    a=V(acc); dt=R(dt)
    if dt<=0: raise ValueError('positive finite dt branch required')

    if not cfg.enabled:
        if decay is not None or rms is not None: raise ValueError('disabled guard consumes no exp/sqrt witnesses')
        return Result(state,a,a,ZERO,F(0),F(0),F(0),state.weight)

    if not state.initialized:
        if decay is not None or rms is not None: raise ValueError('first-sample guard seed consumes no exp/sqrt witnesses')
        stages=(a,a,a,a)
        det=(a,ZERO)
        nxt=State(stages,det,ZERO,F(0),True)
        return Result(nxt,a,a,ZERO,F(0),F(0),F(0),F(0))

    if decay is None or rms is None: raise ValueError('initialized enabled guard requires decay and RMS witnesses')
    # Header names alpha/gamma as the retained old-state coefficients exp(-...).
    alpha,gamma,beta,slew=decay.alpha,decay.gamma,decay.beta,decay.slew

    stages=list(state.stages)
    low=a
    for i in range(cfg.poles):
        stages[i]=add(scale(1-alpha,low),scale(alpha,stages[i]))
        low=stages[i]

    detectors=list(state.detect_stages)
    hp=a
    for i in range(2):
        detectors[i]=add(scale(1-gamma,hp),scale(gamma,detectors[i]))
        hp=sub(hp,detectors[i])

    hpsq=sq(hp)
    removed=tuple(state.removed_ms[i]+beta*(hpsq[i]-state.removed_ms[i]) for i in range(3))
    rr=_rms(removed,rms)

    if not (cfg.engage_lo>0):
        target=F(1)
    elif cfg.engage_hi>cfg.engage_lo:
        target=(rr-cfg.engage_lo)/(cfg.engage_hi-cfg.engage_lo)
    else:
        target=F(1) if rr>=cfg.engage_lo else F(0)
    target=clamp(target,F(0),F(1))

    raw_weight=state.weight+slew*(target-state.weight)
    weight=raw_weight
    if weight<WEIGHT_EPS: weight=F(0)
    elif weight>1-WEIGHT_EPS: weight=F(1)

    if weight<=0: out=a
    elif weight>=1: out=low
    else: out=add(a,scale(weight,sub(low,a)))

    nxt=State(tuple(stages),tuple(detectors),removed,weight,True)
    excess=max(F(0),rr-cfg.engage_lo)
    return Result(nxt,out,low,hp,rr,excess,target,raw_weight)


def readiness():
    return {
      'guard_disabled_identity_branch':True,
      'guard_first_sample_seed_identity_branch':True,
      'guard_lowpass_cascade_materialized':True,
      'guard_two_pole_detector_materialized':True,
      'guard_removed_mean_square_and_RMS_materialized':True,
      'guard_engagement_slew_and_rail_snapping_materialized':True,
      'guard_conditioned_accel_and_excess_RMS_materialized':True,
      'guard_exp_sqrt_binary32_ancestry_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
