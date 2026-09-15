"""Uniform machine WPE supply bounds from a bounded vertical input.

Exact rational inductive inequalities, not a trace or a contraction proof.
The bounds preserve the actual recurrence and all local FMA alternatives.
Correctly rounded transcendental witnesses define the named arithmetic model;
target libm correspondence and upstream observer totality remain separate.
"""
from __future__ import annotations
from fractions import Fraction as F
from functools import lru_cache
from copy import deepcopy

from tools.stability.ou3_alt_contraction import finite_wpe_moment_binary32 as M
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as LOG
from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as EXP
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT
from tools.stability.ou3_alt_contraction import finite_wpe_usable_binary32 as USABLE

QUALIFICATION='OU3_ALT_BOUNDED_INPUT_WPE_UNIFORM_SUPPLIES_V1'
U=F(1,2**24); ETA=F(1,2**150)
DT=B.rn32(F(1,200))
LAMBDA=B.mul(M.TWO_PI,B.rn32(F(1,50)))
INPUT=F(32)
HP1=F(2**17); HP2=F(2**29); V=F(2**33); E=F(2**37)
RAW_LO=F(1,2**57); RAW_HI=F(2**16)
LOG_RAW_LO=F(-40); LOG_RAW_HI=F(12); LOG_ABS=F(48)
STATE_BOUNDS={'accel_prev':INPUT,'hp1':HP1,'hp1_prev':HP1,'hp2':HP2,
              'velocity':V,'elevation':E,'velocity_mean':2*V,'velocity_sq':2*V*V,
              'elevation_mean':2*E,'elevation_sq':2*E*E,'weight':F(2),'elapsed':F(2**17)}


def rounded_abs(x):
    """Absolute RNE envelope valid for normal/subnormal results before overflow."""
    x=F(x)
    if x<0: raise ValueError('nonnegative absolute bound required')
    return (1+U)*x+ETA


@lru_cache(maxsize=1024)
def exp_interval(x):
    x=F(x)
    lo,hi=EXP.exp_minus_enclosure(abs(x))
    return (lo,hi) if x<=0 else (1/hi,1/lo)


def _log_unit(x,terms=32):
    """log(x)=2 atanh((x-1)/(x+1)), x in [1,2], positive tail."""
    z=(x-1)/(x+1)
    if not 0<=z<=F(1,3): raise ValueError('log range reduction outside [1,2]')
    total=sum((2*z**(2*k+1)/F(2*k+1) for k in range(terms)),F(0))
    tail=2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
    return total,total+tail


@lru_cache(maxsize=1024)
def log_interval(x):
    x=F(x)
    if not RAW_LO<=x<=RAW_HI: raise ValueError('WPE log input outside proved period bounds')
    y=x; exponent=0
    while y<1: y*=2; exponent-=1
    while y>=2: y/=2; exponent+=1
    lo,hi=_log_unit(y); l2,h2=_log_unit(F(2))
    return (lo+exponent*l2,hi+exponent*h2) if exponent>=0 else (lo+exponent*h2,hi+exponent*l2)


def _hits_cell(lo,hi,value):
    if value is None or not B.is_binary32(F(value)): return False
    value=F(value)
    # Monotonic RNE of an enclosing interval supplies a conservative witness
    # relation, including a possible boundary cell. No libm accuracy is inferred.
    return B.rn32(lo)<=value<=B.rn32(hi)


def check_exp(argument,result):
    argument=F(argument)
    if not -60<=argument<=60: raise ValueError('WPE exp argument outside proved domain')
    if not _hits_cell(*exp_interval(argument),result):
        raise ValueError('WPE exp witness detached from correctly rounded same argument')


def check_log(argument,result):
    if not _hits_cell(*log_interval(F(argument)),result):
        raise ValueError('WPE log witness detached from correctly rounded same raw period')


def unique_exp_minus(x):
    lo,hi=exp_interval(-F(x)); a,b=B.rn32(lo),B.rn32(hi)
    if a!=b: raise ArithmeticError('fixed WPE coefficient crosses an RNE cell')
    return a


def check_config(cfg):
    if (cfg.lambda_,cfg.moment_horizon_periods,cfg.min_horizon_sec,cfg.max_horizon_sec)!=(LAMBDA,F(4),F(20),F(180)):
        raise ValueError('uniform WPE bounds require unchanged default construction constants')


@lru_cache(maxsize=1)
def _build():
    USABLE.readiness()  # exact shipping header hash, including operation order
    R=rounded_abs
    decay=unique_exp_minus(B.mul(LAMBDA,DT))
    gain=B.div(1-decay,LAMBDA)
    alphas=tuple(B.sub(1,unique_exp_minus(B.div(DT,F(t)))) for t in (20,180))
    # Each envelope is affine in alpha. Endpoints therefore cover EVERY
    # rounded horizon/exp result in this interval, including branch changes.
    after={'hp1':R(decay*R(R(HP1+INPUT)+INPUT)),
           'hp2':R(decay*R(R(HP2+HP1)+HP1)),
           'velocity':R(R(decay*V)+R(gain*HP2)),
           'elevation':R(R(decay*E)+R(gain*V))}
    for name,target,squared in (('weight',F(1),False),('velocity_mean',V,False),
            ('elevation_mean',E,False),('velocity_sq',V,True),('elevation_sq',E,True)):
        vals=[]
        for alpha in alphas:
            forcing=R(R(alpha*target)*target) if squared else R(alpha*target)
            vals.append(R(R(R(1-alpha)*STATE_BOUNDS[name])+forcing))
        after[name]=max(vals)
    margins={name:STATE_BOUNDS[name]-value for name,value in after.items()}
    # Above the literal weight gate, every debiasing division and square is
    # finite. Clamping a nonnegative second moment minus a square cannot
    # increase it, for either separate operations or fused subtraction.
    vvar=R(2*V*V/M.WEIGHT_GATE)
    evar=R(2*E*E/M.WEIGHT_GATE)
    ratio=R(vvar/M.VAR_GATE)
    mean_square=R(R(2*E/M.WEIGHT_GATE)**2)
    sqrt_upper=(1+U)*ROOT.sqrt_enclosure(ratio)[1]+ETA
    sqrt_lower=(1-U)*ROOT.sqrt_enclosure(M.OMEGA_GATE)[0]-ETA
    raw_lower=(1-U)*M.TWO_PI/sqrt_upper-ETA
    raw_upper=R(M.TWO_PI/sqrt_lower)
    log_alphas=tuple(B.sub(1,unique_exp_minus(B.div(DT,t))) for t in (LOG.HORIZON_MIN,LOG.HORIZON_MAX))
    log_after=[]
    for a in log_alphas:
        # Keep convex-combination cancellation; charge the separately rounded
        # difference/product before the final add (also dominates the FMA path).
        difference=F(40)+LOG_ABS
        error=a*(U*difference+ETA)+U*a*((1+U)*difference+ETA)+ETA
        log_after.append(R((1-a)*LOG_ABS+a*40+error))
    log_margin=LOG_ABS-max(log_after)
    exp48_hi=exp_interval(LOG_ABS)[1]
    max_intermediate=2*max(mean_square,ratio,evar,R(4*R(exp48_hi)),*after.values())
    assert all(x>0 for x in margins.values()) and log_margin>0
    assert RAW_LO<raw_lower<raw_upper<RAW_HI
    assert exp_interval(LOG_RAW_LO)[1]<RAW_LO and exp_interval(LOG_RAW_HI)[0]>RAW_HI
    assert max_intermediate<F(2**120)  # far below binary32 overflow, not a useful ISS gain
    return {'qualification':QUALIFICATION,'vertical_input_abs_upper':INPUT,
            'state_bounds':STATE_BOUNDS,'induction_margins':margins,
            'moment_alpha_interval':(min(alphas),max(alphas)),
            'raw_period_interval':(RAW_LO,RAW_HI),'raw_period_derived_lower':raw_lower,
            'raw_period_derived_upper':raw_upper,'log_raw_interval':(LOG_RAW_LO,LOG_RAW_HI),
            'log_state_abs_upper':LOG_ABS,'log_induction_margin':log_margin,
            'all_intermediates_abs_below':F(2**120),
            'reset_to_every_finite_prefix_bounded_input_induction_closed':True,
            'source_uniform_raw_period_bound_under_input_contract_closed':True,
            'no_positive_variance_or_period_production_assumed':True,
            'both_local_FMA_and_separate_moment_results_covered':True,
            'elapsed_bound_follows_monotone_RNE_and_2pow17_fixed_point':B.add(F(2**17),DT)==F(2**17),
            'upstream_observer_totality_or_startup_reachability_closed':False,
            'target_libm_and_compiler_correspondence_closed':False,
            'storage_search_allowed':False}


def build():
    # Callers may mutate reports while auditing false promotions. Keep the
    # validator's reference independent of those mutations.
    return deepcopy(_build())


def check_state(state):
    check_config(state.cfg)
    for mom in (state.separate,state.fma):
        for name,bound in STATE_BOUNDS.items():
            if abs(getattr(mom,name))>bound: raise ValueError('WPE state exceeds uniform '+name+' bound')
        if min(mom.weight,mom.velocity_sq,mom.elevation_sq)<0:
            raise ValueError('WPE nonnegative accumulator invariant violated')
        if mom.raw_period is not None and not RAW_LO<=mom.raw_period<=RAW_HI:
            raise ValueError('WPE retained raw period exceeds uniform bound')
    for track in (state.logs.separate,state.logs.fma):
        if track.log_period is not None and abs(track.log_period)>LOG_ABS:
            raise ValueError('WPE log state exceeds uniform bound')


def check_mode_before(track,witness):
    md=witness.moment.get('moment_decay_exp')
    if md is not None:
        alpha=B.sub(1,md); lo,hi=build()['moment_alpha_interval']
        if not lo<=alpha<=hi: raise ValueError('WPE moment alpha outside uniform default-horizon bound')
    if witness.horizon is not None: check_exp(witness.horizon.log_period,witness.horizon.period_exp)
    if witness.raw_log is not None: check_log(witness.raw_log.raw_period,witness.raw_log.log_raw)
    if isinstance(witness.log,LOG.SmoothWitness):
        check_exp(track.log_period,witness.log.sea_period_exp)
        horizon=LOG.clamp(B.mul(LOG.LOG_SMOOTH_PERIODS,witness.log.sea_period_exp),LOG.HORIZON_MIN,LOG.HORIZON_MAX)
        check_exp(-B.div(DT,horizon),witness.log.decay_exp)
    if witness.usable is not None:
        if witness.usable.log_period is None: raise ValueError('bounded WPE cannot produce a nonfinite log')
        check_exp(witness.usable.log_period,witness.usable.period_exp)


def validate(report):
    return [k+' differs from bounded-input WPE induction' for k,v in build().items() if report.get(k)!=v]
