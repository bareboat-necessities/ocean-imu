"""Binary32 WPE moment/raw-period machine recurrence for ALT finite master.

This is the missing deployed arithmetic between the machine vertical sample and
``std::log(raw_period)``.  It carries the two high-pass stages, leaky
velocity/elevation states, weighted moments and every shipping early-return
branch.  A valid branch derives ``raw_period`` from the SAME machine moments;
there is no independent raw-period input port.

Transcendentals are retained as deployment witnesses.  ``exp`` witnesses must
land in the binary32 RNE cell of a rigorous enclosure of the same rounded
argument.  ``sqrt`` is checked against the correctly-rounded binary32 square-root
cell.  These relations do not yet certify a particular target libm/compiler.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as EXP
from tools.stability.ou3_alt_contraction import finite_binary32_mahony_startup as SQRT
from tools.stability.ou3_alt_contraction import finite_wpe_runtime as SHADOW

SOURCE=Path(__file__).resolve().parents[3]/'src/tuner/WavePeriodEstimator.h'
QUALIFICATION='OU3_ALT_WPE_MOMENT_BINARY32_V1'
ZERO=B.rn32(0); ONE=B.rn32(1); TWO=B.rn32(2); THREE=B.rn32(3); SIX=B.rn32(6)
PI_F=SQRT.M.value(0x40490FDB)
TWO_PI=B.mul(TWO,PI_F)
WEIGHT_GATE=B.rn32(F(1,1000)); VAR_GATE=B.rn32(F(1,10**12)); OMEGA_GATE=B.rn32(F(1,10**8))


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(name+' must be actual binary32')
    return q


def _sum_products(a,b,c,d,*,successor=None):
    """Allowed scalar compiler outcomes for ``a*b + c*d``."""
    ab=B.mul(a,b); cd=B.mul(c,d)
    vals=tuple(sorted(set((B.add(ab,cd),B.fma(a,b,cd),B.fma(c,d,ab)))))
    if successor is None: return vals
    y=_q(successor,'WPE stored sum')
    if y not in vals: raise ValueError('WPE stored sum outside legal contraction set')
    return y


def _ema(previous,target,alpha,*,successor=None):
    om=B.sub(ONE,alpha)
    vals=_sum_products(om,previous,alpha,target)
    if successor is None: return vals
    y=_q(successor,'WPE EMA successor')
    if y not in vals: raise ValueError('WPE EMA successor outside legal contraction set')
    return y


def _exp_minus(arg,witness,name):
    a=_q(arg,name+' argument'); w=_q(witness,name+' result')
    if a<0 or a>60: raise ValueError(name+' argument outside certified exp enclosure')
    lo,hi=EXP.exp_minus_enclosure(a)
    if not EXP._interval_hits_rne_cell(lo,hi,w):
        raise ValueError(name+' result detached from same rounded exp argument RNE cell')
    return w


@dataclass(frozen=True)
class Config:
    lambda_:F
    moment_horizon_periods:F
    min_horizon_sec:F
    max_horizon_sec:F
    def __post_init__(self):
        for n in ('lambda_','moment_horizon_periods','min_horizon_sec','max_horizon_sec'):
            object.__setattr__(self,n,_q(getattr(self,n),'WPE '+n))
        if self.lambda_<=0 or self.moment_horizon_periods<1 or self.min_horizon_sec<=0 or self.max_horizon_sec<self.min_horizon_sec:
            raise ValueError('invalid machine WPE configuration')
    @classmethod
    def from_shadow(cls,cfg:SHADOW.WPEConfig):
        if not isinstance(cfg,SHADOW.WPEConfig): raise TypeError('shadow WPEConfig required')
        return cls(B.rn32(cfg.lambda_),B.rn32(cfg.moment_horizon_periods),
                   B.rn32(cfg.min_horizon_sec),B.rn32(cfg.max_horizon_sec))


@dataclass(frozen=True)
class State:
    accel_prev:F=ZERO
    hp1:F=ZERO
    hp1_prev:F=ZERO
    hp2:F=ZERO
    velocity:F=ZERO
    elevation:F=ZERO
    velocity_mean:F=ZERO
    velocity_sq:F=ZERO
    elevation_mean:F=ZERO
    elevation_sq:F=ZERO
    weight:F=ZERO
    elapsed:F=ZERO
    raw_period:F|None=None
    last_moment_horizon:F|None=None
    samples:int=0
    def __post_init__(self):
        for n in ('accel_prev','hp1','hp1_prev','hp2','velocity','elevation','velocity_mean','velocity_sq',
                  'elevation_mean','elevation_sq','weight','elapsed'):
            object.__setattr__(self,n,_q(getattr(self,n),'WPE '+n))
        if self.raw_period is not None: object.__setattr__(self,'raw_period',_q(self.raw_period,'WPE raw period'))
        if self.last_moment_horizon is not None: object.__setattr__(self,'last_moment_horizon',_q(self.last_moment_horizon,'WPE horizon'))
        if not isinstance(self.samples,int) or self.samples<0: raise ValueError('nonnegative WPE sample count required')


@dataclass(frozen=True)
class MomentSuccessors:
    weight:F
    velocity_mean:F
    velocity_sq:F
    elevation_mean:F
    elevation_sq:F
    def __post_init__(self):
        for n in ('weight','velocity_mean','velocity_sq','elevation_mean','elevation_sq'):
            object.__setattr__(self,n,_q(getattr(self,n),'WPE '+n+' successor'))


@dataclass(frozen=True)
class StepResult:
    before:State
    state:State
    vertical_accel:F
    dt:F
    decay:F
    gain:F
    horizon:F|None
    moment_decay:F|None
    alpha:F|None
    branch:str
    velocity_var:F|None=None
    elevation_var:F|None=None
    omega_sq:F|None=None
    sqrt_omega:F|None=None
    raw_period:F|None=None
    def __post_init__(self):
        if self.state.samples!=self.before.samples+1: raise ValueError('WPE sample count did not advance')
        if self.raw_period!=self.state.raw_period and self.branch=='valid-period':
            raise ValueError('WPE valid raw period detached from successor state')


def _clamp(x,lo,hi): return min(max(x,lo),hi)


def step(state:State,cfg:Config,*,dt,vertical_accel,decay_exp,
         canonical_period=None,moment_decay_exp=None,moment_successors:MomentSuccessors|None=None,
         velocity_successor=None,elevation_successor=None,velocity_var_successor=None,elevation_var_successor=None,sqrt_omega=None):
    if not isinstance(state,State) or not isinstance(cfg,Config): raise TypeError('machine WPE State/Config required')
    h=_q(dt,'WPE dt'); x=_q(vertical_accel,'WPE vertical acceleration')
    if h<=0: raise ValueError('positive WPE dt required')
    lamdt=B.mul(cfg.lambda_,h); decay=_exp_minus(lamdt,decay_exp,'WPE leak decay')
    gain=B.div(B.sub(ONE,decay),cfg.lambda_) if cfg.lambda_>B.rn32(F(1,10**9)) else h

    t1=B.add(state.hp1,x); t1=B.sub(t1,state.accel_prev); hp1=B.mul(decay,t1)
    t2=B.add(state.hp2,hp1); t2=B.sub(t2,state.hp1_prev); hp2=B.mul(decay,t2)
    vvals=_sum_products(decay,state.velocity,gain,hp2)
    if velocity_successor is None: raise TypeError('WPE velocity stored successor required')
    vel=_q(velocity_successor,'WPE velocity successor')
    if vel not in vvals: raise ValueError('WPE velocity successor outside legal contraction set')
    evals=_sum_products(decay,state.elevation,gain,vel)
    if elevation_successor is None: raise TypeError('WPE elevation stored successor required')
    elev=_q(elevation_successor,'WPE elevation successor')
    if elev not in evals: raise ValueError('WPE elevation successor outside legal contraction set')
    elapsed=B.add(state.elapsed,h)
    base=replace(state,accel_prev=x,hp1=hp1,hp1_prev=hp1,hp2=hp2,velocity=vel,elevation=elev,
                 elapsed=elapsed,samples=state.samples+1)

    moment_start=B.div(THREE,cfg.lambda_)
    if elapsed<moment_start:
        if any(w is not None for w in (canonical_period,moment_decay_exp,moment_successors,sqrt_omega)):
            raise ValueError('pre-moment-start WPE branch consumes no moment/period witnesses')
        return StepResult(state,base,x,h,decay,gain,None,None,None,'pre-moment-start')

    if canonical_period is None:
        period=SIX
    else:
        period=_q(canonical_period,'WPE canonical period for moment horizon')
        if period<=0: raise ValueError('positive canonical period required')
    requested=B.mul(cfg.moment_horizon_periods,period)
    horizon=_clamp(requested,cfg.min_horizon_sec,cfg.max_horizon_sec)
    if moment_decay_exp is None: raise TypeError('moment branch requires exp result')
    mdarg=B.div(h,horizon); md=_exp_minus(mdarg,moment_decay_exp,'WPE moment decay')
    alpha=B.sub(ONE,md)
    if not isinstance(moment_successors,MomentSuccessors): raise TypeError('moment branch requires stored moment successors')
    ms=moment_successors
    if ms.weight not in _ema(state.weight,ONE,alpha): raise ValueError('WPE weight successor detached')
    if ms.velocity_mean not in _ema(state.velocity_mean,vel,alpha): raise ValueError('WPE velocity mean successor detached')
    vel2=B.mul(vel,vel)
    if ms.velocity_sq not in _ema(state.velocity_sq,vel2,alpha): raise ValueError('WPE velocity square successor detached')
    if ms.elevation_mean not in _ema(state.elevation_mean,elev,alpha): raise ValueError('WPE elevation mean successor detached')
    elev2=B.mul(elev,elev)
    if ms.elevation_sq not in _ema(state.elevation_sq,elev2,alpha): raise ValueError('WPE elevation square successor detached')
    base=replace(base,weight=ms.weight,velocity_mean=ms.velocity_mean,velocity_sq=ms.velocity_sq,
                 elevation_mean=ms.elevation_mean,elevation_sq=ms.elevation_sq,last_moment_horizon=horizon)
    if ms.weight<=WEIGHT_GATE:
        if sqrt_omega is not None: raise ValueError('insufficient-weight WPE branch consumes no sqrt')
        return StepResult(state,base,x,h,decay,gain,horizon,md,alpha,'insufficient-weight')

    vm=B.div(ms.velocity_mean,ms.weight); em=B.div(ms.elevation_mean,ms.weight)
    vsecond=B.div(ms.velocity_sq,ms.weight); esecond=B.div(ms.elevation_sq,ms.weight)
    vchoices=tuple(sorted(set((max(ZERO,B.sub(vsecond,B.mul(vm,vm))),
                               max(ZERO,B.fma(-vm,vm,vsecond))))))
    echoices=tuple(sorted(set((max(ZERO,B.sub(esecond,B.mul(em,em))),
                               max(ZERO,B.fma(-em,em,esecond))))))
    if velocity_var_successor is None or elevation_var_successor is None:
        raise TypeError('post-weight WPE branch requires stored variance results')
    vvar=_q(velocity_var_successor,'WPE velocity variance'); evar=_q(elevation_var_successor,'WPE elevation variance')
    if vvar not in vchoices or evar not in echoices:
        raise ValueError('WPE stored variance outside legal compiler arithmetic set')
    if evar<=VAR_GATE or vvar<=VAR_GATE:
        if sqrt_omega is not None: raise ValueError('degenerate WPE branch consumes no sqrt')
        return StepResult(state,base,x,h,decay,gain,horizon,md,alpha,'degenerate-moments',vvar,evar)
    ratio=B.div(vvar,evar); lam2=B.mul(cfg.lambda_,cfg.lambda_); omega=B.sub(ratio,lam2)
    if omega<=OMEGA_GATE:
        if sqrt_omega is not None: raise ValueError('nonpositive-omega WPE branch consumes no sqrt')
        return StepResult(state,base,x,h,decay,gain,horizon,md,alpha,'nonpositive-omega',vvar,evar,omega)
    if sqrt_omega is None: raise TypeError('valid WPE moment ratio requires sqrt result')
    root=_q(sqrt_omega,'WPE sqrt omega')
    if root!=SQRT.sqrt32(omega):
        raise ValueError('WPE sqrt result detached from same binary32 omega_sq RNE cell')
    raw=B.div(TWO_PI,root)
    nxt=replace(base,raw_period=raw)
    return StepResult(state,nxt,x,h,decay,gain,horizon,md,alpha,'valid-period',vvar,evar,omega,root,raw)


def _source_shape_matches():
    s=SOURCE.read_text()
    return all(x in s for x in (
      'const float decay = std::exp(-lambda_ * dt_sec);',
      'const float stage1 = decay * (high_pass_1_ + vertical_accel_ms2 - accel_prev_);',
      'velocity_ = decay * velocity_ + gain * stage2;',
      'weight_ = (1.0f - alpha) * weight_ + alpha;',
      'const float ratio_sq = velocity_var / elevation_var;',
      'const float omega_sq = ratio_sq - lambda_ * lambda_;',
      '2.0f * 3.14159265358979323846f / std::sqrt(omega_sq);'))


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_WPE_moment_source_shape_matches':_source_shape_matches(),
      'two_high_pass_and_leaky_integrator_binary32_graph_materialized':True,
      'weighted_moment_binary32_graph_materialized':True,
      'all_raw_period_early_return_branches_materialized':True,
      'raw_period_has_no_independent_input_port':True,
      'raw_period_derived_from_same_machine_moments_and_sqrt_result':True,
      'exp_witnesses_bound_to_same_rounded_arguments_RNE_cells':True,
      'sqrt_witness_bound_to_same_binary32_omega_RNE_cell':True,
      'target_exp_sqrt_libm_correspondence_closed':False,
      'compiler_profile_selection_closed':False,
      'canonical_period_horizon_ancestry_closed':False,
      'same_machine_vertical_input_ancestry_closed':False,
      'dual_compiler_persistent_moment_history_attached':False,
      'source_uniform_raw_period_supply_bound_closed':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,
    }
