"""Uniform machine WPE supply bounds from a bounded vertical input.

Exact rational inductive inequalities, not a trace or a contraction proof.
The bounds preserve the actual recurrence and all local FMA alternatives.
Both strict RNE and a named error-inclusive libm relation are available. The
declared MEMS domain supplies the vertical input; target instruction/link
correspondence and startup capture remain separate.
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
RNE_PROFILE='rne'
ERROR_PROFILE='exp-sqrt-relative-2^-20-log-absolute-2^-14'
EXP_SQRT_REL_ERROR=F(1,2**20)
LOG_ABS_ERROR=F(1,2**14)
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
def log_interval(x,scale_exponent=0):
    x=F(x)
    scale=_scale(scale_exponent)
    if not RAW_LO/scale<=x<=RAW_HI: raise ValueError('WPE log input outside proved period bounds')
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


def _profile_errors(profile):
    if profile==RNE_PROFILE: return F(0),F(0)
    if profile==ERROR_PROFILE: return EXP_SQRT_REL_ERROR,LOG_ABS_ERROR
    raise ValueError('unknown WPE transcendental relation')


def exp_result_interval(argument,libm_profile=RNE_PROFILE):
    relative,_=_profile_errors(libm_profile)
    lo,hi=exp_interval(F(argument))
    if relative: lo,hi=(1-relative)*lo-ETA,(1+relative)*hi+ETA
    return B.rn32(lo),B.rn32(hi)


def check_exp(argument,result,libm_profile=RNE_PROFILE):
    argument=F(argument)
    if not -60<=argument<=60: raise ValueError('WPE exp argument outside proved domain')
    lo,hi=exp_result_interval(argument,libm_profile)
    if result is None or not B.is_binary32(F(result)) or not lo<=F(result)<=hi:
        raise ValueError('WPE exp witness detached from named same-argument relation')


def check_log(argument,result,scale_exponent=0,libm_profile=RNE_PROFILE):
    _,absolute=_profile_errors(libm_profile)
    lo,hi=log_interval(F(argument),scale_exponent)
    if not _hits_cell(lo-absolute,hi+absolute,result):
        raise ValueError('WPE log witness detached from named same-raw-period relation')


def check_sqrt(argument,result,libm_profile=RNE_PROFILE):
    relative,_=_profile_errors(libm_profile)
    if not relative:
        if result!=M.SQRT.sqrt32(argument):
            raise ValueError('WPE sqrt result detached from same binary32 omega_sq RNE cell')
        return
    argument=F(argument)
    if argument<=0: raise ValueError('positive WPE sqrt argument required')
    lo,hi=ROOT.sqrt_enclosure(argument)
    if not _hits_cell((1-relative)*lo-ETA,(1+relative)*hi+ETA,result):
        raise ValueError('WPE sqrt result detached from named same-argument relation')


def unique_exp_minus(x):
    lo,hi=exp_interval(-F(x)); a,b=B.rn32(lo),B.rn32(hi)
    if a!=b: raise ArithmeticError('fixed WPE coefficient crosses an RNE cell')
    return a


def check_config(cfg):
    if (cfg.lambda_,cfg.moment_horizon_periods,cfg.min_horizon_sec,cfg.max_horizon_sec)!=(LAMBDA,F(4),F(20),F(180)):
        raise ValueError('uniform WPE bounds require unchanged default construction constants')


def _scale(exponent):
    if type(exponent) is not int or exponent<0:
        raise ValueError('nonnegative integer WPE supply scale exponent required')
    return F(2**exponent)


def scale_for_input(input_abs_upper,libm_profile=RNE_PROFILE):
    """Choose a containing amplitude envelope; never infer source admission."""
    x=F(input_abs_upper)
    if x<0: raise ValueError('nonnegative WPE input norm bound required')
    exponent=0
    while INPUT*2**exponent<x: exponent+=1
    build(exponent,libm_profile)  # fail closed if THIS enclosure cannot prove totality
    return exponent


@lru_cache(maxsize=32)
def _build(scale_exponent=0,libm_profile=RNE_PROFILE):
    USABLE.readiness()  # exact shipping header hash, including operation order
    R=rounded_abs
    relative_error,log_error=_profile_errors(libm_profile)
    scale=_scale(scale_exponent)
    # Bound the exponent before constructing large rational powers. At this
    # scale the variance division enclosure already exceeds binary32 range.
    if scale_exponent>5:
        raise ValueError('WPE amplitude enclosure cannot prove binary32 totality')
    bounds={name:bound*(scale*scale if name in ('velocity_sq','elevation_sq')
                       else 1 if name in ('weight','elapsed') else scale)
            for name,bound in STATE_BOUNDS.items()}
    input_bound=INPUT*scale; hp1=HP1*scale; hp2=HP2*scale
    velocity=V*scale; elevation=E*scale
    raw_lo=RAW_LO/scale; log_raw_lo=LOG_RAW_LO-scale_exponent
    decay_lo,decay=exp_result_interval(-B.mul(LAMBDA,DT),libm_profile)
    gain=B.div(1-decay_lo,LAMBDA)
    alphas=tuple(B.sub(1,d) for t in (20,180)
                 for d in exp_result_interval(-B.div(DT,F(t)),libm_profile))
    assert 0<decay_lo<=decay<1 and 0<min(alphas)<=max(alphas)<1
    # Each envelope is affine in alpha. Endpoints therefore cover EVERY
    # rounded horizon/exp result in this interval, including branch changes.
    after={'hp1':R(decay*R(R(hp1+input_bound)+input_bound)),
           'hp2':R(decay*R(R(hp2+hp1)+hp1)),
           'velocity':R(R(decay*velocity)+R(gain*hp2)),
           'elevation':R(R(decay*elevation)+R(gain*velocity))}
    for name,target,squared in (('weight',F(1),False),('velocity_mean',velocity,False),
            ('elevation_mean',elevation,False),('velocity_sq',velocity,True),('elevation_sq',elevation,True)):
        vals=[]
        for alpha in alphas:
            forcing=R(R(alpha*target)*target) if squared else R(alpha*target)
            vals.append(R(R(R(1-alpha)*bounds[name])+forcing))
        after[name]=max(vals)
    margins={name:bounds[name]-value for name,value in after.items()}
    # Above the literal weight gate, every debiasing division and square is
    # finite. Clamping a nonnegative second moment minus a square cannot
    # increase it, for either separate operations or fused subtraction.
    vvar=R(2*velocity*velocity/M.WEIGHT_GATE)
    evar=R(2*elevation*elevation/M.WEIGHT_GATE)
    ratio=R(vvar/M.VAR_GATE)
    mean_square=R(R(2*elevation/M.WEIGHT_GATE)**2)
    sqrt_upper=R((1+relative_error)*ROOT.sqrt_enclosure(ratio)[1]+ETA)
    sqrt_lower=(1-U)*((1-relative_error)*ROOT.sqrt_enclosure(M.OMEGA_GATE)[0]-ETA)-ETA
    raw_lower=(1-U)*M.TWO_PI/sqrt_upper-ETA
    raw_upper=R(M.TWO_PI/sqrt_lower)
    log_alphas=tuple(B.sub(1,d) for t in (LOG.HORIZON_MIN,LOG.HORIZON_MAX)
                     for d in exp_result_interval(-B.div(DT,t),libm_profile))
    assert 0<min(log_alphas)<=max(log_alphas)<1
    log_target=R(abs(log_raw_lo)+log_error)
    log_after=[]
    for a in log_alphas:
        # Keep convex-combination cancellation; charge the separately rounded
        # difference/product before the final add (also dominates the FMA path).
        difference=log_target+LOG_ABS
        error=a*(U*difference+ETA)+U*a*((1+U)*difference+ETA)+ETA
        log_after.append(R((1-a)*LOG_ABS+a*log_target+error))
    log_margin=LOG_ABS-max(log_after)
    exp48_hi=exp_result_interval(LOG_ABS,libm_profile)[1]
    max_intermediate=2*max(mean_square,ratio,evar,R(4*R(exp48_hi)),*after.values())
    assert all(x>0 for x in margins.values()) and log_margin>0
    assert raw_lo<raw_lower<raw_upper<RAW_HI
    assert exp_interval(log_raw_lo)[1]<raw_lo and exp_interval(LOG_RAW_HI)[0]>RAW_HI
    finite_ceiling=F(2**128-2**104)
    if max_intermediate>=finite_ceiling:
        raise ValueError('WPE amplitude enclosure cannot prove binary32 totality')
    coarse_ceiling=F(2**(120+2*scale_exponent))
    assert max_intermediate<coarse_ceiling
    return {'qualification':QUALIFICATION,'vertical_input_abs_upper':input_bound,
            'scale_exponent':scale_exponent,
            'libm_profile':libm_profile,
            'exp_and_sqrt_relative_error_allowance_before_RNE_enclosure':relative_error,
            'log_absolute_error_allowance_before_RNE_enclosure':log_error,
            'leak_decay_interval':(decay_lo,decay),
            'state_bounds':bounds,'induction_margins':margins,
            'moment_alpha_interval':(min(alphas),max(alphas)),
            'raw_period_interval':(raw_lo,RAW_HI),'raw_period_derived_lower':raw_lower,
            'raw_period_derived_upper':raw_upper,'log_raw_interval':(log_raw_lo,LOG_RAW_HI),
            'log_state_abs_upper':LOG_ABS,'log_induction_margin':log_margin,
            'all_intermediates_abs_below':coarse_ceiling,
            'derived_intermediate_bound':max_intermediate,
            'finite_binary32_intermediate_margin':finite_ceiling-max_intermediate,
            'reset_to_every_finite_prefix_bounded_input_induction_closed':True,
            'source_uniform_raw_period_bound_under_input_contract_closed':True,
            'no_positive_variance_or_period_production_assumed':True,
            'both_local_FMA_and_separate_moment_results_covered':True,
            'elapsed_bound_follows_monotone_RNE_and_2pow17_fixed_point':B.add(F(2**17),DT)==F(2**17),
            'upstream_observer_totality_or_startup_reachability_closed':False,
            'target_libm_and_compiler_correspondence_closed':False,
            'storage_search_allowed':False}


def build(scale_exponent=0,libm_profile=RNE_PROFILE):
    # Callers may mutate reports while auditing false promotions. Keep the
    # validator's reference independent of those mutations.
    return deepcopy(_build(scale_exponent,libm_profile))


def check_state(state):
    check_config(state.cfg)
    report=build(state.supply_scale_exponent,state.libm_profile)
    for mom in (state.separate,state.fma):
        for name,bound in report['state_bounds'].items():
            if abs(getattr(mom,name))>bound: raise ValueError('WPE state exceeds uniform '+name+' bound')
        if min(mom.weight,mom.velocity_sq,mom.elevation_sq)<0:
            raise ValueError('WPE nonnegative accumulator invariant violated')
        if mom.raw_period is not None and not report['raw_period_interval'][0]<=mom.raw_period<=RAW_HI:
            raise ValueError('WPE retained raw period exceeds uniform bound')
    for track in (state.logs.separate,state.logs.fma):
        if track.log_period is not None and abs(track.log_period)>LOG_ABS:
            raise ValueError('WPE log state exceeds uniform bound')


def check_mode_before(track,witness,scale_exponent=0,libm_profile=RNE_PROFILE):
    md=witness.moment.get('moment_decay_exp')
    if md is not None:
        alpha=B.sub(1,md); lo,hi=build(scale_exponent,libm_profile)['moment_alpha_interval']
        if not lo<=alpha<=hi: raise ValueError('WPE moment alpha outside uniform default-horizon bound')
    if witness.horizon is not None: check_exp(witness.horizon.log_period,witness.horizon.period_exp,libm_profile)
    if witness.raw_log is not None: check_log(witness.raw_log.raw_period,witness.raw_log.log_raw,scale_exponent,libm_profile)
    if isinstance(witness.log,LOG.SmoothWitness):
        check_exp(track.log_period,witness.log.sea_period_exp,libm_profile)
        horizon=LOG.clamp(B.mul(LOG.LOG_SMOOTH_PERIODS,witness.log.sea_period_exp),LOG.HORIZON_MIN,LOG.HORIZON_MAX)
        check_exp(-B.div(DT,horizon),witness.log.decay_exp,libm_profile)
    if witness.usable is not None:
        if witness.usable.log_period is None: raise ValueError('bounded WPE cannot produce a nonfinite log')
        check_exp(witness.usable.log_period,witness.usable.period_exp,libm_profile)


def validate(report):
    try: expected=build(report.get('scale_exponent',0),report.get('libm_profile',RNE_PROFILE))
    except (TypeError,ValueError): return ['invalid WPE amplitude supply envelope']
    return [k+' differs from bounded-input WPE induction' for k,v in expected.items() if report.get(k)!=v]


@lru_cache(maxsize=64)
def _live_iss_certificate(supply_norm_upper):
    """Conditional Live supply from the EXISTING nine-coordinate ISS norm.

    In the dormant guard scope, raw accelerometer residual = nu_acc + thermal.
    Their squared norms are two blocks of ||w||², hence the residual norm is
    at most sqrt(2) W < 3W/2. No startup sensor amplitude is substituted for W.
    Observer totality, target correspondence and source membership are premises.
    """
    from tools.stability.ou3_alt_contraction import finite_startup_sensor_contract as SENSOR
    w=F(supply_norm_upper)
    if w<0: raise ValueError('nonnegative symbolic Live ISS norm bound required')
    SENSOR.build()  # source-locked normalization/projection lemma, not startup capture
    residual=F(3,2)*w
    raw=SENSOR.G+SENSOR.A+SENSOR.BIAS_NORM+residual
    stored=(1+U)*raw+2*ETA
    q2=F(139,125)
    derr=4*q2*((1+U)**8-1)+16*ETA
    value=(q2+2*derr)*stored
    for _ in range(8): value=rounded_abs(value)
    value=rounded_abs(value+B.rn32(SENSOR.G))
    exponent=scale_for_input(value)
    return {'ISS_norm_upper':w,'raw_accel_residual_norm_upper':residual,
            'raw_accel_norm_upper':raw,'stored_accel_norm_upper':stored,
            'vertical_abs_upper_after_defined_Mahony_update':value,
            'scale_exponent':exponent,
            'vertical_input_abs_upper':INPUT*2**exponent,
            'raw_residual_derived_from_same_nu_acc_and_thermal_blocks':True,
            'startup_sensor_residual_bound_reused_for_Live':False,
            'dormant_guard_required':True,
            'observer_totality_assumed_not_proved':True,
            'unrestricted_Live_ISS_amplitude_domain_closed':False}


def live_iss_certificate(supply_norm_upper):
    return deepcopy(_live_iss_certificate(F(supply_norm_upper)))


@lru_cache(maxsize=2)
def _physical_input_certificate(libm_profile):
    """Compose the declared MEMS packet domain with the actual WPE recurrence.

    Source qualification happens before execution in the shared Live input
    contract. The observer's all-time lattice induction carries construction
    integral history and establishes the vertical bound on every finite prefix.
    Neither induction assumes a startup deadline or a particular input phase.
    """
    from tools.stability.ou3_alt_contraction import finite_live_input_contract as INPUT_DOMAIN
    sensor=INPUT_DOMAIN.build(); prefix=sensor['prefix']
    if not (prefix['all_initialized_finite_prefixes_totality_closed']
            and prefix['vertical_input_abs_512_closed_under_prefix_premises']):
        raise RuntimeError('MEMS input contract lost initialized observer totality')
    exponent=scale_for_input(prefix['vertical_abs_upper'],libm_profile)
    supply=build(exponent,libm_profile)
    assert exponent==4 and supply['vertical_input_abs_upper']==512
    return {'input_qualification':sensor['qualification'],
            'input_profile':sensor['profile'],'libm_profile':libm_profile,
            'raw_accel_component_upper':sensor['raw_accel_component_upper'],
            'raw_gyro_component_upper':sensor['raw_gyro_component_upper'],
            'all_initialized_finite_prefixes_totality_closed':True,
            'vertical_abs_upper':prefix['vertical_abs_upper'],
            'vertical_input_abs_upper':supply['vertical_input_abs_upper'],
            'scale_exponent':exponent,
            'source_uniform_WPE_supplies_under_declared_MEMS_prefix_closed':True,
            'startup_deadline_required_for_WPE_supply_bound':False,
            'supply_does_not_require_a_numeric_ISS_W_limit':True,
            'startup_integral_not_reset_at_Live':prefix['startup_integral_not_reset_at_Live'],
            'input_admission_independent_of_execution_success':True,
            'initialized_observer_totality_under_seed_and_prefix_premises_closed':True,
            'startup_seed_capture_and_target_qualification_proved_here':False,
            'target_libm_correspondence_qualified':False,
            'storage_search_allowed':False}


def physical_input_certificate(libm_profile=RNE_PROFILE):
    return deepcopy(_physical_input_certificate(libm_profile))


def target_error_profile_certificate():
    """Join target exp/log approximation proofs to the actual WPE argument sets.

    Common scalar/division semantics, sqrt and final firmware object selection
    are qualified separately. This bridge never assumes ideal libm rounding.
    """
    from tools.stability.ou3_alt_contraction import target_wpe_libm as TARGET
    proof=TARGET.certificate(); correspondence=TARGET.profile_correspondence()
    exp=proof['exp']; log=proof['log']; supply=build(4,ERROR_PROFILE)
    if not (exp['all_input_exp_error_bound_proved_under_scalar_premises']
            and log['all_input_log_error_bound_proved_under_scalar_premises']):
        raise RuntimeError('target WPE library approximation proof missing')
    exp_lo,exp_hi=exp['argument_interval']
    log_lo,log_hi=log['argument_interval']
    raw_lo,raw_hi=supply['raw_period_interval']
    if not (exp_lo<=-LOG_ABS<=LOG_ABS<=exp_hi and log_lo<=raw_lo<=raw_hi<=log_hi):
        raise ValueError('target WPE library domain does not cover same source word')
    if not (exp['total_relative_exp_error']<EXP_SQRT_REL_ERROR
            and log['total_absolute_log_error']<LOG_ABS_ERROR):
        raise ValueError('target WPE library error exceeds stable supply allowance')
    return {'target_qualification':proof['qualification'],'libm_profile':ERROR_PROFILE,
            'exp_relative_error_upper':exp['total_relative_exp_error'],
            'log_absolute_error_upper':log['total_absolute_log_error'],
            'actual_WPE_arguments_covered_by_target_library_proof':True,
            'target_exp_log_satisfy_WPE_error_profile_under_scalar_and_link_premises':True,
            'pinned_WPE_exp_log_approximation_correspondence_closed':correspondence[
                'pinned_WPE_exp_log_approximation_correspondence_closed'],
            'pinned_WPE_sqrt_approximation_correspondence_closed':correspondence[
                'pinned_WPE_sqrt_approximation_correspondence_closed'],
            'transcendental_correct_rounding_assumed':False,
            'scalar_sqrt_ROM_and_firmware_link_qualification_supplied_here':False}
