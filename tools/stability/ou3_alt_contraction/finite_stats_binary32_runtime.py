"""Binary32 SeaStateAutoTuner moment/variance runtime for ALT.

This covers the shipping statistics path after the adaptive-band output is a
machine float.  It materializes frequency/horizon arithmetic, binds the single
``alpha_var`` exp call to its rounded argument, and carries one actual machine
DebiasedEMA state whose four successors must belong to the finite local
multiply-add contraction sets.  ``getAccelVariance`` is then evaluated with
literal binary32 divisions, multiply, subtraction and zero floor.

No exact-real variance is substituted.  Upstream band-machine membership and
target expf correspondence remain explicit obligations.
"""
from __future__ import annotations

from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as HORIZON
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as FREQ
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as R
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as E

SOURCE=Path(__file__).resolve().parents[3]/'src/tuner/SeaStateAutoTuner.h'
ZERO=B.rn32(0); ONE=B.rn32(1); TWO=B.rn32(2); HALF=B.rn32(F(1,2))
TIME_MIN=B.rn32(F(1,2)); TIME_MAX=B.rn32(6); HORIZON_MIN=B.rn32(F(1,20)); HORIZON_MAX=B.rn32(35)
READY_WEIGHT=B.rn32(F(1,10**6))
QUALIFICATION='OU3_ALT_STATS_BINARY32_RUNTIME_V1'; MAX_SAMPLES=HORIZON.MAX_STEPS


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must be actual binary32')
    return q

def _uniq(v): return tuple(sorted(set(F(x) for x in v)))

def _sum_products(a,b,c,d):
    ab=B.mul(a,b); cd=B.mul(c,d)
    return _uniq((B.add(ab,cd),B.fma(a,b,cd),B.fma(c,d,ab)))

def _scale_plus(a,b,c):
    ab=B.mul(a,b)
    return _uniq((B.add(ab,c),B.fma(a,b,c)))


@dataclass(frozen=True)
class State:
    frequency:F|None=None
    tau_var:F=ZERO
    mean_value:F=ZERO
    mean_weight:F=ZERO
    sq_value:F=ZERO
    sq_weight:F=ZERO
    samples:int=0
    def __post_init__(self):
        if self.frequency is not None: object.__setattr__(self,'frequency',_q(self.frequency,'stats frequency'))
        for n in ('tau_var','mean_value','mean_weight','sq_value','sq_weight'):
            object.__setattr__(self,n,_q(getattr(self,n),n))
        if not isinstance(self.samples,int) or not 0<=self.samples<=MAX_SAMPLES: raise ValueError('stats sample count outside bounded startup+word horizon')
        if self.mean_weight<0 or self.sq_weight<0: raise ValueError('stats weights must be nonnegative')
    @property
    def var_ready(self): return self.mean_weight>READY_WEIGHT and self.sq_weight>READY_WEIGHT


@dataclass(frozen=True)
class Coefficients:
    input_frequency:F
    frequency:F
    dt:F
    sea_time:F
    T_eff:F
    tau_requested:F
    tau_var:F
    exp_decay:F
    alpha:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        f=F(self.input_frequency)
        if not FREQ.is_finite_input(f): raise ValueError('binary32 input frequency required')
        object.__setattr__(self,'input_frequency',f)
        for n in ('frequency','dt','sea_time','T_eff','tau_requested','tau_var','exp_decay','alpha'):
            object.__setattr__(self,n,_q(getattr(self,n),n))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong stats coefficient qualification')
        if min(self.input_frequency,self.frequency,self.dt,self.tau_var,self.exp_decay)<=0 or not 0<=self.alpha<=1:
            raise ValueError('invalid stats coefficient domain')


@dataclass(frozen=True)
class Envelope:
    before:State
    coefficients:Coefficients
    accel:F
    accel_sq:F
    decay:F
    mean_values:tuple
    mean_weights:tuple
    sq_values:tuple
    sq_weights:tuple
    def __post_init__(self):
        if not isinstance(self.before,State) or not isinstance(self.coefficients,Coefficients): raise TypeError('stats state/coefficient relation required')
        object.__setattr__(self,'accel',_q(self.accel,'stats accel')); object.__setattr__(self,'accel_sq',_q(self.accel_sq,'stats accel square')); object.__setattr__(self,'decay',_q(self.decay,'stats decay'))
        for n in ('mean_values','mean_weights','sq_values','sq_weights'):
            vals=_uniq(getattr(self,n)); object.__setattr__(self,n,vals)
            if not vals or not all(B.is_binary32(x) for x in vals): raise ValueError('stats contraction set must contain binary32 values')

    def accepts(self,state:State):
        return (isinstance(state,State) and state.samples==self.before.samples+1
                and state.frequency==self.coefficients.frequency and state.tau_var==self.coefficients.tau_var
                and state.mean_value in self.mean_values and state.mean_weight in self.mean_weights
                and state.sq_value in self.sq_values and state.sq_weight in self.sq_weights)


def coefficients(cfg:R.StatsConfig,*,frequency,dt,exp_decay):
    if not isinstance(cfg,R.StatsConfig): raise TypeError('StatsConfig required')
    f=F(frequency); h=_q(dt,'stats dt')
    if not FREQ.is_finite_input(f): raise ValueError('binary32 stats input frequency required')
    if f<=0 or h<=0: raise ValueError('positive stats frequency/dt required')
    fmin=B.rn32(cfg.f_min); fmax=B.rn32(cfg.f_max); fe=min(max(f,fmin),fmax)
    sea=min(max(B.div(HALF,fe),TIME_MIN),TIME_MAX); teff=B.mul(TWO,sea)
    k=B.rn32(cfg.K_periods); tmin=B.rn32(cfg.tau_var_min); tmax=B.rn32(cfg.tau_var_max)
    req=min(max(B.mul(k,teff),tmin),tmax)
    lo=HORIZON_MIN
    if h>lo: lo=min(h,HORIZON_MAX)
    tau=min(max(req,lo),HORIZON_MAX)
    x=B.div(h,tau); e=_q(exp_decay,'stats exp decay')
    elo,ehi=E.exp_minus_enclosure(x)
    if not E._interval_hits_rne_cell(elo,ehi,e): raise ValueError('stats exp witness detached from SAME rounded dt/tau argument RNE cell')
    a=B.sub(ONE,e)
    return Coefficients(f,fe,h,sea,teff,req,tau,e,a)


def envelope(state:State,coef:Coefficients,*,accel):
    if not isinstance(state,State) or not isinstance(coef,Coefficients): raise TypeError('stats State and Coefficients required')
    if state.samples>=MAX_SAMPLES: raise ValueError('stats runtime exceeded bounded startup+word horizon')
    x=_q(accel,'stats accel'); d=B.sub(ONE,coef.alpha); x2=B.mul(x,x)
    mv=_sum_products(d,state.mean_value,coef.alpha,x)
    mw=_scale_plus(d,state.mean_weight,coef.alpha)
    sv=_sum_products(d,state.sq_value,coef.alpha,x2)
    sw=_scale_plus(d,state.sq_weight,coef.alpha)
    return Envelope(state,coef,x,x2,d,mv,mw,sv,sw)


def step(state:State,coef:Coefficients,*,accel,successor:State):
    env=envelope(state,coef,accel=accel)
    if not isinstance(successor,State) or not env.accepts(successor):
        raise ValueError('actual machine stats successor outside same-step contraction set')
    return successor,env


def variance(state:State):
    if not isinstance(state,State): raise TypeError('stats machine State required')
    if not state.var_ready: return ZERO
    mu=B.div(state.mean_value,state.mean_weight)
    second=B.div(state.sq_value,state.sq_weight)
    return max(ZERO,B.sub(second,B.mul(mu,mu)))


def variance_outcomes(state:State):
    """No-reassociation results of A_sq.get() - mu*mu, including one FMA."""
    if not isinstance(state,State): raise TypeError('stats machine State required')
    if not state.var_ready: return (ZERO,)
    mu=B.div(state.mean_value,state.mean_weight)
    second=B.div(state.sq_value,state.sq_weight)
    return _uniq((max(ZERO,B.sub(second,B.mul(mu,mu))),
                  max(ZERO,B.fma(-mu,mu,second))))


def _source_shape_matches():
    s=SOURCE.read_text()
    needles=('frequency_hz = f_eff;','const float sea_time_sec =',
      'const float T_eff = 2.0f * sea_time_sec;','const float alpha_var = 1.0f - std::exp(-dt_s / tau_var_sec);',
      'value  = (1.0f - alpha) * value + alpha * x;','weight = (1.0f - alpha) * weight + alpha;',
      'A_mean.update(accel, alpha_var);','A_sq.update(accel * accel, alpha_var);',
      'const float mu = A_mean.get();','return std::max(0.0f, A_sq.get() - mu * mu);')
    return all(x in s for x in needles)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_tuner_stats_source_shape_matches':_source_shape_matches(),
      'frequency_horizon_and_alpha_binary32_graph_materialized':True,
      'stats_exp_result_bound_to_same_dt_tau_argument_by_RNE_cell':True,
      'actual_machine_mean_and_square_EMA_successors_bound_to_local_contraction_sets':True,
      'mean_and_square_weights_use_same_alpha':True,
      'binary32_debiased_variance_readout_materialized':True,
      'no_exact_real_accel_variance_substituted':True,
      'target_exp_libm_correspondence_closed':False,
      'upstream_machine_band_output_attached':False,
      'target_compiler_contraction_membership_closed':False,
      'startup_frontend_machine_history_attached':False,
      'Live_600_step_machine_history_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
