"""Binary32 relation for the shipping AccelVibrationGuard input path.

The shipping filter's public ``update`` boundary receives ``float dt`` and
``Eigen::Vector3f`` gyro/acceleration.  This module first certifies those stored
binary32 values as the round-to-nearest-even images of the admitted finite-real
source values, then executes the actual AccelVibrationGuard branch structure.
It never rounds or reuses the exact-model guard successor.

The state includes the conditioning cascade, two-pole detector, detector RMS
EMA and engagement weight.  Local multiply/add contraction alternatives are
retained.  exp/sqrt witnesses are tied to the same rounded arguments by rigorous
real enclosures.  Native target libm, Eigen expression evaluation and mutable
runtime guard-configuration ancestry remain explicit open obligations.
"""
from __future__ import annotations

from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as HORIZON
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_binary32_mahony as MAH
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT

SRC=Path(__file__).resolve().parents[3]/'src/tuner/AccelVibrationGuard.h'
QUALIFICATION='OU3_ALT_MACHINE_ACCEL_GUARD_BINARY32_V1'
ZERO=B.rn32(0); ONE=B.rn32(1); TWO=B.rn32(2); PI_F=MAH.value(0x40490FDB)
MAX_SAMPLES=HORIZON.MAX_STEPS  # shared conditional timeout-crossing + 600 budget; not reachability
Vec=tuple[F,F,F]


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must be actual binary32')
    return q

def _v(v,name):
    if len(v)!=3: raise ValueError(f'{name} must be vec3')
    return tuple(_q(x,name) for x in v)

def _vzero(): return (ZERO,ZERO,ZERO)
def _sub(a,b): return tuple(B.sub(x,y) for x,y in zip(a,b))
def _abs2(a): return tuple(B.mul(x,x) for x in a)
def _uniq(xs): return tuple(sorted(set(F(x) for x in xs)))

def _sum3_values(v):
    vals=[]
    for p in permutations(v): vals.append(B.add(B.add(p[0],p[1]),p[2]))
    return _uniq(vals)

def _lin_values(alpha,x,old):
    om=B.sub(ONE,alpha); ax=B.mul(om,x); ao=B.mul(alpha,old)
    return _uniq((B.add(ax,ao),B.fma(om,x,ao),B.fma(alpha,old,ax)))
def _vec_lin(alpha,x,old,successor):
    y=_v(successor,'LP successor')
    for i in range(3):
        if y[i] not in _lin_values(alpha,x[i],old[i]): raise ValueError('guard LP successor outside local contraction set')
    return y

def _exp_minus(mag,witness,name):
    m=_q(mag,name+' argument'); w=_q(witness,name+' result')
    if m<0 or m>60: raise ValueError(name+' argument outside retained exp enclosure')
    lo,hi=EXP.exp_minus_enclosure(m)
    if not EXP._interval_hits_rne_cell(lo,hi,w): raise ValueError(name+' result detached from SAME rounded argument RNE cell')
    return w

def _sqrt(v,witness,name):
    x=_q(v,name+' radicand'); w=_q(witness,name+' result')
    if x<0: raise ValueError(name+' negative radicand')
    lo,hi=ROOT.sqrt_enclosure(x)
    if not EXP._interval_hits_rne_cell(lo,hi,w): raise ValueError(name+' result detached from SAME radicand RNE cell')
    return w

def api_vec(exact,stored,name):
    """Certify the actual float API store, not an exact-shadow successor."""
    s=_v(stored,name)
    if tuple(B.rn32(F(x)) for x in exact)!=s: raise ValueError(name+' detached from binary32 API rounding')
    return s

def api_scalar(exact,stored,name):
    s=_q(stored,name)
    if B.rn32(F(exact))!=s: raise ValueError(name+' detached from binary32 API rounding')
    return s

@dataclass(frozen=True)
class Config:
    cutoff_hz:F=B.rn32(14)
    poles:int=2
    detect_hz:F=B.rn32(25)
    engage_lo:F=B.rn32(F(3,100))
    engage_hi:F=B.rn32(F(8,100))
    slew_tau:F=B.rn32(5)
    removed_rms_hz:F=B.rn32(F(1,20))
    weight_epsilon:F=B.rn32(F(1,10000))
    qualification:str=QUALIFICATION
    def __post_init__(self):
        for n in ('cutoff_hz','detect_hz','engage_lo','engage_hi','slew_tau','removed_rms_hz','weight_epsilon'):
            object.__setattr__(self,n,_q(getattr(self,n),n))
        if not isinstance(self.poles,int) or not 1<=self.poles<=4: raise ValueError('guard poles outside shipping range')
        if self.detect_hz<=0 or self.slew_tau<=0 or self.removed_rms_hz<=0: raise ValueError('positive guard constants required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong guard config qualification')

@dataclass(frozen=True)
class State:
    stages:tuple[Vec,Vec,Vec,Vec]=(_vzero(),_vzero(),_vzero(),_vzero())
    detect_stages:tuple[Vec,Vec]=(_vzero(),_vzero())
    removed_ms:Vec=_vzero()
    weight:F=ZERO
    initialized:bool=False
    samples:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        object.__setattr__(self,'stages',tuple(_v(x,'guard stage') for x in self.stages))
        object.__setattr__(self,'detect_stages',tuple(_v(x,'detector stage') for x in self.detect_stages))
        object.__setattr__(self,'removed_ms',_v(self.removed_ms,'removed_ms'))
        object.__setattr__(self,'weight',_q(self.weight,'guard weight'))
        if not ZERO<=self.weight<=ONE: raise ValueError('guard weight outside [0,1]')
        if not isinstance(self.initialized,bool) or not isinstance(self.samples,int) or not 0<=self.samples<=MAX_SAMPLES: raise ValueError('invalid guard persistent state')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong guard state qualification')

@dataclass(frozen=True)
class Result:
    before:State
    state:State
    raw_gyro:Vec
    raw_acc:Vec
    conditioned_acc:Vec
    removed_rms:F
    excess_rms:F
    enabled:bool
    seeded:bool
    qualification:str=QUALIFICATION
    def __post_init__(self):
        object.__setattr__(self,'raw_gyro',_v(self.raw_gyro,'raw gyro'))
        object.__setattr__(self,'raw_acc',_v(self.raw_acc,'raw accel'))
        object.__setattr__(self,'conditioned_acc',_v(self.conditioned_acc,'conditioned accel'))
        object.__setattr__(self,'removed_rms',_q(self.removed_rms,'removed RMS'))
        object.__setattr__(self,'excess_rms',_q(self.excess_rms,'excess RMS'))
        if self.state.samples!=self.before.samples+1: raise ValueError('guard sample count did not advance')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong guard result qualification')


def step(state:State,cfg:Config,*,raw_gyro,raw_acc,dt,
         lp_alpha_exp=None,lp_successors=None,
         detect_gamma_exp=None,detect_successors=None,
         removed_beta_exp=None,removed_ms_successor=None,
         removed_rms_sqrt=None,slew_exp=None,weight_successor=None,
         output_successor=None):
    if not isinstance(state,State) or not isinstance(cfg,Config): raise TypeError('machine guard State/Config required')
    if state.samples>=MAX_SAMPLES: raise ValueError('machine guard exceeded bounded proof horizon')
    gyro=_v(raw_gyro,'raw gyro'); acc=_v(raw_acc,'raw accel'); h=_q(dt,'guard dt')
    if h<=0: raise ValueError('positive finite guard dt required')
    enabled=cfg.cutoff_hz>0
    if not enabled:
        nxt=State(state.stages,state.detect_stages,state.removed_ms,state.weight,state.initialized,state.samples+1)
        return Result(state,nxt,gyro,acc,acc,ZERO,ZERO,False,False)
    if not state.initialized:
        stages=(acc,acc,acc,acc); det=(acc,_vzero())
        nxt=State(stages,det,_vzero(),ZERO,True,state.samples+1)
        return Result(state,nxt,gyro,acc,acc,ZERO,ZERO,True,True)
    if None in (lp_alpha_exp,detect_gamma_exp,removed_beta_exp,removed_rms_sqrt,slew_exp,weight_successor):
        raise TypeError('initialized guard requires all exp/RMS/weight witnesses')
    # conditioning cascade
    lpmag=B.mul(B.mul(B.mul(TWO,PI_F),cfg.cutoff_hz),h); alpha=_exp_minus(lpmag,lp_alpha_exp,'guard LP exp')
    if lp_successors is None or len(lp_successors)!=cfg.poles: raise TypeError('one conditioning successor per active pole required')
    stages=list(state.stages); low=acc
    for i in range(cfg.poles):
        stages[i]=_vec_lin(alpha,low,state.stages[i],lp_successors[i]); low=stages[i]
    # two-pole detector; high_passed -= updated stage after each pole
    dmag=B.mul(B.mul(B.mul(TWO,PI_F),cfg.detect_hz),h); gamma=_exp_minus(dmag,detect_gamma_exp,'guard detector exp')
    if detect_successors is None or len(detect_successors)!=2: raise TypeError('two detector-stage successors required')
    det=list(state.detect_stages); hp=acc
    for i in range(2):
        det[i]=_vec_lin(gamma,hp,state.detect_stages[i],detect_successors[i]); hp=_sub(hp,det[i])
    # beta = 1-exp(-2*pi*kRemovedRmsHz*dt)
    bmag=B.mul(B.mul(B.mul(TWO,PI_F),cfg.removed_rms_hz),h); e=_exp_minus(bmag,removed_beta_exp,'guard RMS exp'); beta=B.sub(ONE,e)
    hp2=_abs2(hp)
    if removed_ms_successor is None: raise TypeError('removed-ms successor required')
    rmsms=_v(removed_ms_successor,'removed-ms successor')
    for i in range(3):
        delta=B.sub(hp2[i],state.removed_ms[i]); prod=B.mul(beta,delta)
        vals=_uniq((B.add(state.removed_ms[i],prod),B.fma(beta,delta,state.removed_ms[i])))
        if rmsms[i] not in vals: raise ValueError('removed-ms successor outside shipping contraction set')
    sums=_sum3_values(rmsms); rad=max(ZERO,sums[0])
    rms=_sqrt(rad,removed_rms_sqrt,'guard removed RMS')
    # target and engagement slew
    if not (cfg.engage_lo>0): target=ONE
    elif cfg.engage_hi>cfg.engage_lo: target=B.div(B.sub(rms,cfg.engage_lo),B.sub(cfg.engage_hi,cfg.engage_lo))
    else: target=ONE if rms>=cfg.engage_lo else ZERO
    target=max(ZERO,min(ONE,target))
    smag=B.div(h,cfg.slew_tau); se=_exp_minus(smag,slew_exp,'guard slew exp'); slew=B.sub(ONE,se)
    delta=B.sub(target,state.weight); prod=B.mul(slew,delta)
    vals=_uniq((B.add(state.weight,prod),B.fma(slew,delta,state.weight)))
    w=_q(weight_successor,'guard weight successor')
    if w not in vals: raise ValueError('guard weight successor outside shipping contraction set')
    if w<cfg.weight_epsilon: w=ZERO
    elif w>B.sub(ONE,cfg.weight_epsilon): w=ONE
    if w<=0: out=acc
    elif w>=1: out=low
    else:
        if output_successor is None: raise TypeError('partially engaged guard requires output successor')
        out=_v(output_successor,'guard mixed output')
        for i in range(3):
            d=B.sub(low[i],acc[i]); p=B.mul(w,d); ov=_uniq((B.add(acc[i],p),B.fma(w,d,acc[i])))
            if out[i] not in ov: raise ValueError('guard mixed output outside shipping contraction set')
    excess=max(ZERO,B.sub(rms,cfg.engage_lo))
    nxt=State(tuple(stages),tuple(det),rmsms,w,True,state.samples+1)
    return Result(state,nxt,gyro,acc,out,rms,excess,True,False)


def _source_shape_matches():
    s=SRC.read_text()
    return all(x in s for x in (
      'std::exp(-2.0f * static_cast<float>(M_PI) * cutoff_hz_ * dt)',
      'stages_[i] = (1.0f - alpha) * low_passed + alpha * stages_[i];',
      'std::exp(-2.0f * static_cast<float>(M_PI) * detect_hz_ * dt)',
      'removed_ms_ += beta * (high_passed.cwiseAbs2() - removed_ms_);',
      'weight_ += slew * (target - weight_);',
      'return acc + weight_ * (low_passed - acc);'))

def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_guard_source_shape_matches':_source_shape_matches(),
      'binary32_update_API_boundary_materialized':True,
      'first_sample_guard_seed_and_exact_transparency_materialized':True,
      'conditioning_and_detector_cascades_materialized':True,
      'removed_RMS_EMA_and_sqrt_materialized':True,
      'engagement_slew_rails_and_mixed_output_materialized':True,
      'exp_and_sqrt_same_argument_RNE_cells_required':True,
      'target_libm_and_Eigen_expression_correspondence_closed':False,
      'mutable_guard_runtime_config_ancestry_closed':False,
      'source_uniform_guard_supply_bound_closed':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
