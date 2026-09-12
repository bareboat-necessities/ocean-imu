"""Finite WavePeriodEstimator runtime recurrence for the ALT same-history word.

This is a literal finite-state representation of the shipping WPE update. It
keeps the two shared high-pass stages, leaky velocity/elevation proxies,
period-scaled moment horizon, weighted moments, moment-ratio period extraction,
canonical log-period smoothing and the one-way usable-period latch on one
history.

The deployed exp/log/sqrt evaluations are represented by explicit runtime
witnesses. Algebraic relations that can be checked exactly are checked here;
transcendental ancestry and binary32 roundoff remain open obligations. This is
therefore finite runtime attachment, not source admission or a stability gate.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_prediction_graph as P


def R(x): return P.rational(x)
def clamp(x, lo, hi): return min(max(x, lo), hi)


@dataclass(frozen=True)
class WPEConfig:
    lambda_: F
    moment_horizon_periods: F
    log_smoothing_periods: F
    min_horizon_sec: F
    max_horizon_sec: F
    def __post_init__(self):
        vals=[R(x) for x in (self.lambda_,self.moment_horizon_periods,self.log_smoothing_periods,
                              self.min_horizon_sec,self.max_horizon_sec)]
        lam,mhp,lsp,hmin,hmax=vals
        if lam<=0 or mhp<1 or lsp<0 or hmin<=0 or hmax<hmin:
            raise ValueError('invalid WPE configuration')
        for n,v in zip(('lambda_','moment_horizon_periods','log_smoothing_periods','min_horizon_sec','max_horizon_sec'),vals):
            object.__setattr__(self,n,v)


@dataclass(frozen=True)
class WPEState:
    accel_prev: F=F(0)
    hp1: F=F(0)
    hp1_prev: F=F(0)
    hp2: F=F(0)
    velocity: F=F(0)
    elevation: F=F(0)
    velocity_mean: F=F(0)
    velocity_sq: F=F(0)
    elevation_mean: F=F(0)
    elevation_sq: F=F(0)
    weight: F=F(0)
    elapsed: F=F(0)
    raw_period: F|None=None
    log_period: F|None=None
    usable_period: bool=False
    last_moment_horizon: F=F(20)
    last_log_horizon: F=F(0)
    def __post_init__(self):
        for n in ('accel_prev','hp1','hp1_prev','hp2','velocity','elevation','velocity_mean','velocity_sq',
                  'elevation_mean','elevation_sq','weight','elapsed','last_moment_horizon','last_log_horizon'):
            object.__setattr__(self,n,R(getattr(self,n)))
        if self.raw_period is not None: object.__setattr__(self,'raw_period',R(self.raw_period))
        if self.log_period is not None: object.__setattr__(self,'log_period',R(self.log_period))
        if not isinstance(self.usable_period,bool): raise TypeError('usable latch must be literal bool')
        if self.weight<0 or self.elapsed<0: raise ValueError('invalid WPE state')


@dataclass(frozen=True)
class ExpWitness:
    value: F
    def __post_init__(self):
        v=R(self.value)
        if v<=0: raise ValueError('exp witness must be positive')
        object.__setattr__(self,'value',v)


@dataclass(frozen=True)
class PeriodWitness:
    sqrt_omega_sq: F
    raw_period: F
    log_raw_period: F
    def __post_init__(self):
        for n in ('sqrt_omega_sq','raw_period','log_raw_period'):
            object.__setattr__(self,n,R(getattr(self,n)))
        if self.sqrt_omega_sq<=0 or self.raw_period<=0: raise ValueError('positive period witnesses required')


@dataclass(frozen=True)
class LogUpdateWitness:
    current_period: F
    decay: F
    next_log_period: F
    def __post_init__(self):
        for n in ('current_period','decay','next_log_period'):
            object.__setattr__(self,n,R(getattr(self,n)))
        if self.current_period<=0 or not 0<self.decay<=1: raise ValueError('invalid log smoothing witnesses')


@dataclass(frozen=True)
class CanonicalOutputWitness:
    """Post-update exp(log_period) and exp(-log_period) outputs."""
    period: F
    frequency: F
    def __post_init__(self):
        p,f=R(self.period),R(self.frequency)
        if p<=0 or f<=0 or p*f!=1:
            raise ValueError('canonical WPE period/frequency must be positive reciprocals')
        object.__setattr__(self,'period',p); object.__setattr__(self,'frequency',f)


@dataclass(frozen=True)
class UpdateResult:
    state: WPEState
    produced_period: bool
    frequency: F|None
    period: F|None
    moment_horizon: F|None
    moment_alpha: F|None


def _period_for_horizon(current_period:F|None):
    return current_period if current_period is not None and current_period>0 else F(6)


def update(s:WPEState,cfg:WPEConfig,*,dt,vertical_accel,decay:ExpWitness,
           moment_decay:ExpWitness|None=None,period_witness:PeriodWitness|None=None,
           log_witness:LogUpdateWitness|None=None,current_period:F|None=None,
           current_frequency:F|None=None,post_output:CanonicalOutputWitness|None=None):
    """One literal valid-input shipping WPE update.

    ``current_*`` are outputs of the canonical log state on entry. ``post_output``
    is the output of the post-update log state and therefore the value consumed
    by ``hasUsablePeriod`` and downstream tuning on this sample. Their exp/log
    ancestry remains an explicit finite-precision obligation.
    """
    dt=R(dt); x=R(vertical_accel)
    if dt<=0: raise ValueError('shipping-valid positive dt required')
    d=decay.value
    if d>1: raise ValueError('shipping decay exp(-lambda*dt) must be <=1')
    if s.log_period is None:
        if current_period is not None or current_frequency is not None:
            raise ValueError('no visible period before canonical log state exists')
        cp=cf=None
    else:
        if current_period is None or current_frequency is None:
            raise ValueError('finite log state requires period/frequency witnesses')
        cp,cf=R(current_period),R(current_frequency)
        if cp<=0 or cf<=0 or cp*cf!=1:
            raise ValueError('entry WPE period/frequency must be reciprocal')

    gain=(1-d)/cfg.lambda_
    stage1=d*(s.hp1+x-s.accel_prev)
    stage2=d*(s.hp2+stage1-s.hp1_prev)
    vel=d*s.velocity+gain*stage2
    elev=d*s.elevation+gain*vel
    elapsed=s.elapsed+dt
    base=replace(s,accel_prev=x,hp1=stage1,hp1_prev=stage1,hp2=stage2,
                 velocity=vel,elevation=elev,elapsed=elapsed)

    moment_start=F(3)/cfg.lambda_
    if elapsed < moment_start:
        if any(w is not None for w in (moment_decay,period_witness,log_witness,post_output)):
            raise ValueError('pre-moment-start branch consumes no moment/period witnesses')
        return UpdateResult(base,False,cf,cp,None,None)

    requested=cfg.moment_horizon_periods*_period_for_horizon(cp)
    horizon=clamp(requested,cfg.min_horizon_sec,cfg.max_horizon_sec)
    if moment_decay is None: raise ValueError('moment branch requires exp(-dt/horizon) witness')
    md=moment_decay.value
    if md>1: raise ValueError('moment decay must be <=1')
    alpha=1-md
    weight=md*s.weight+alpha
    vm=md*s.velocity_mean+alpha*vel
    vs=md*s.velocity_sq+alpha*vel*vel
    em=md*s.elevation_mean+alpha*elev
    es=md*s.elevation_sq+alpha*elev*elev
    base=replace(base,weight=weight,velocity_mean=vm,velocity_sq=vs,
                 elevation_mean=em,elevation_sq=es,last_moment_horizon=horizon)
    if weight <= F(1,1000):
        if any(w is not None for w in (period_witness,log_witness,post_output)):
            raise ValueError('insufficient-weight branch consumes no period witness')
        return UpdateResult(base,False,cf,cp,horizon,alpha)

    vmean=vm/weight; emean=em/weight
    vvar=max(F(0),vs/weight-vmean*vmean)
    evar=max(F(0),es/weight-emean*emean)
    if evar<=F(1,10**12) or vvar<=F(1,10**12):
        if any(w is not None for w in (period_witness,log_witness,post_output)):
            raise ValueError('degenerate-moment branch consumes no period witness')
        return UpdateResult(base,False,cf,cp,horizon,alpha)
    omega_sq=vvar/evar-cfg.lambda_*cfg.lambda_
    if omega_sq<=F(1,10**8):
        if any(w is not None for w in (period_witness,log_witness,post_output)):
            raise ValueError('nonpositive omega branch consumes no period witness')
        return UpdateResult(base,False,cf,cp,horizon,alpha)
    if period_witness is None or post_output is None:
        raise ValueError('valid moment ratio requires period and post-canonical witnesses')
    pw=period_witness
    if pw.sqrt_omega_sq*pw.sqrt_omega_sq != omega_sq:
        raise ValueError('sqrt(omega_sq) witness detached from same WPE moments')
    raw=pw.raw_period

    if s.log_period is None:
        if log_witness is not None: raise ValueError('first valid period initializes log state directly')
        logp=pw.log_raw_period; log_h=F(0)
    elif cfg.log_smoothing_periods<=0:
        if log_witness is not None: raise ValueError('zero smoothing consumes no EMA witness')
        logp=pw.log_raw_period; log_h=dt
    else:
        if log_witness is None: raise ValueError('log smoothing branch requires current-period/decay witness')
        lw=log_witness
        if lw.current_period != cp: raise ValueError('log horizon period detached from same canonical state')
        requested_log=cfg.log_smoothing_periods*cp
        log_h=clamp(requested_log,max(F(1,20),min(dt,F(35))),F(35))
        la=1-lw.decay
        expected=s.log_period+la*(pw.log_raw_period-s.log_period)
        if lw.next_log_period != expected: raise ValueError('log-period EMA successor detached from same state/raw period')
        logp=expected

    pp,pf=post_output.period,post_output.frequency
    usable=s.usable_period
    if not usable:
        usable_floor=F(4)/cfg.lambda_
        history=elapsed-moment_start
        if elapsed>=usable_floor and history>=pp: usable=True
    nxt=replace(base,raw_period=raw,log_period=logp,usable_period=usable,last_log_horizon=log_h)
    return UpdateResult(nxt,True,pf,pp,horizon,alpha)


def readiness():
    return {
      'two_high_pass_stages_materialized':True,
      'leaky_velocity_elevation_materialized':True,
      'period_scaled_weighted_moments_materialized':True,
      'moment_ratio_period_branch_materialized':True,
      'canonical_log_period_state_materialized':True,
      'post_log_output_and_one_way_usable_latch_materialized':True,
      'exp_log_sqrt_pi_binary32_ancestry_attached':False,
      'vertical_accel_frontend_same_history_attached':False,
      'adaptive_band_variance_state_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
