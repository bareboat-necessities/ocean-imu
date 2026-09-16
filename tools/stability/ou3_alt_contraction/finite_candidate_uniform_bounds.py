"""All-time finite candidate and pending-commit arithmetic for configured ALT.

The same source-owned sigma target and tau frequency drive SpectralMSE, both
EMA rates, all three persistent candidates and the next pending transaction.
This is a finite-supply induction, not a startup deadline or contraction proof.
The named pow profile retains the compiled exponents and a positive finite
output range; exact-root shadow residuals remain distinct and are not erased.
"""
from __future__ import annotations

from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_frontend_uniform_bounds as FRONT
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_qeff_cache_binary32 as Q
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as LIBM
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_binary32 as SPEC
from tools.stability.ou3_alt_contraction import finite_wpe_uniform_bounds as W
from tools.stability.ou3_alt_contraction import target_powf_range as POW

QUALIFICATION='OU3_ALT_CONFIGURED_CANDIDATE_UNIFORM_FINITE_SUPPLIES_V1'
POW_PROFILE='positive-powf-output-2^-45-to-2^18'
POW_BASE=(F(1,2**44),F(2**17))
POW_OUTPUT=(F(1,2**45),F(2**18))
COMMON_ALPHA=(F(1,2**10),F(1,16))
RS_ALPHA=(F(1,2**12),F(1,64))
STATE_BOUNDS={'tau':(F(1,4),F(25,2)), 'sigma':(F(1,2**12),F(8)),
              'rs':(F(1,16),F(200))}
TARGET_BOUNDS={'tau':(B.div(F(1,2),D.MAX_FREQ),D.MAX_TAU),
               'sigma':(F(1,2**11),F(4)), 'rs':(D.MIN_RS,D.MAX_RS)}
GAMMA=(1+W.U)**8-1
TAIL=F(2**16)*W.ETA


def _ema_interval(previous,target,alphas):
    """RNE(p + alpha*RNE(target-p)), separate or fused outer update.

    Error is bounded by gamma8*(p+a*(p+t))+tail for nonnegative p,t.
    Both the upper and lower affine coefficients in p,t stay positive for the
    declared alpha ranges, so actual box extrema occur at their endpoints.
    This retains the cancellation of p in target-p before bounding errors.
    """
    lower=[]; upper=[]
    for a in alphas:
        assert 1-a-GAMMA*(1+a)>0 and a*(1-GAMMA)>0
        p,t=previous[0],target[0]
        lower.append((1-a)*p+a*t-GAMMA*(p+a*(p+t))-TAIL)
        p,t=previous[1],target[1]
        upper.append((1-a)*p+a*t+GAMMA*(p+a*(p+t))+TAIL)
    return min(lower),max(upper)


def require_config(cfg):
    if not isinstance(cfg,D.DeploymentConfig):
        raise TypeError('source-owned candidate deployment configuration required')
    constants={'min_freq':D.MIN_FREQ,'max_freq':D.MAX_FREQ,'tau_coeff':D.TAU_COEFF,
        'sigma_coeff':D.SIGMA_COEFF,'min_tau':D.MIN_TAU,'max_tau':D.MAX_TAU,
        'max_sigma':D.MAX_SIGMA,'pseudo_tau_ratio':D.PSEUDO_RATIO,
        'pseudo_min':D.PSEUDO_MIN,'pseudo_max':D.PSEUDO_MAX,'min_RS':D.MIN_RS,
        'max_RS':D.MAX_RS,'rs_mse_coeff':D.RS_MSE_COEFF,
        'adapt_tau_sec':D.ADAPT_TAU_SEC,'adapt_tau_sea_periods':D.ADAPT_TAU_SEA_PERIODS,
        'adapt_RS_mult':D.ADAPT_RS_MULT,'adapt_RS_slew_log':D.ADAPT_RS_SLEW_LOG,
        'adapt_every_sec':D.ADAPT_EVERY_SEC}
    if not cfg.clamp_enabled or any(getattr(cfg,name)!=value for name,value in constants.items()):
        raise ValueError('candidate finite supplies require actual shipping defaults')
    if cfg.accel_noise_density!=Q.default_r_a_binary32():
        raise ValueError('qeff cache input detached from configured source density')
    require_pow(cfg.qeff_cache.qeff_argument,cfg.qeff_cache.exponent,cfg.qeff_cache.result)


def require_commit_config(cfg):
    from tools.stability.ou3_alt_contraction import finite_tuner_commit as COMMIT
    if not isinstance(cfg,COMMIT.CommitConfig): raise TypeError('source commit configuration required')
    expected={'pseudo_tau_ratio':D.PSEUDO_RATIO,'pseudo_period_min':D.PSEUDO_MIN,
        'pseudo_period_max':D.PSEUDO_MAX,'pseudo_fixed_period':D.PSEUDO_NOMINAL,
        'min_R_S':D.MIN_RS,'max_R_S':D.MAX_RS,'S_factor':F(1),
        'R_S_x_factor':B.rn32(F(72,100)),'R_S_y_factor':B.rn32(F(72,100))}
    if not cfg.tau_scaled_cadence or cfg.cubic_rs_law:
        raise ValueError('candidate supplies require configured SpectralMSE cadence')
    if any(B.rn32(getattr(cfg,name))!=value for name,value in expected.items()):
        raise ValueError('candidate commit detached from configured scalar factors')


def require_pow(argument,exponent,result):
    argument,exponent,result=map(F,(argument,exponent,result))
    if not all(B.is_binary32(x) for x in (argument,exponent,result)):
        raise ValueError('actual binary32 pow operands and result required')
    if exponent==Q.EXPONENT:
        valid=argument==B.mul(F(2),Q.default_r_a_binary32())
    elif exponent==SPEC.POW_EXPONENT:
        valid=POW_BASE[0]<=argument<=POW_BASE[1]
    else:
        valid=False
    if not valid or not POW_OUTPUT[0]<=result<=POW_OUTPUT[1]:
        raise ValueError('pow witness outside named compiled-exponent target range')


@lru_cache(maxsize=1)
def _build():
    frontend=FRONT.build()
    pow_profile=POW.profile_correspondence()
    assert D._source_shape_matches() and Q._source_shape_matches() and SPEC._source_shape_matches()
    from tools.stability.ou3_alt_contraction import finite_tuner_tau_binary32 as TAU
    from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_binary32 as RSA
    from tools.stability.ou3_alt_contraction import finite_tuner_machine_boundary_commit as COMMIT
    assert TAU._source_shape_matches() and RSA._source_shape_matches() and COMMIT._source_shape_matches()
    # The source sigma floor is evaluated before the coefficient and max clamp.
    sigma_min=B.mul(D.SIGMA_COEFF,W.M.SQRT.sqrt32(LIBM.VAR_FLOOR))
    assert sigma_min>TARGET_BOUNDS['sigma'][0]
    assert frontend['sigma_target_upper']==TARGET_BOUNDS['sigma'][1]
    common_horizon=tuple(B.mul(D.ADAPT_TAU_SEA_PERIODS,x) for x in (F(1,2),F(6)))
    rs_horizon=tuple(B.mul(D.ADAPT_RS_MULT,x) for x in (F(1,2),F(6)))
    alpha_bounds={}
    for name,horizons,box in (('common',common_horizon,COMMON_ALPHA),('rs',rs_horizon,RS_ALPHA)):
        lo,hi=FRONT._alpha_interval(B.div(W.DT,horizons[1]),B.div(W.DT,horizons[0]))
        assert box[0]<=lo<=hi<=box[1]
        alpha_bounds[name]=(lo,hi)
    after={}; margins={}
    for name,previous in STATE_BOUNDS.items():
        image=_ema_interval(previous,TARGET_BOUNDS[name],RS_ALPHA if name=='rs' else COMMON_ALPHA)
        after[name]=image
        margins[name]=(image[0]-previous[0],previous[1]-image[1])
    assert all(lo>0 and hi>0 for lo,hi in margins.values())
    # Existing tau roundoff certification therefore remains defined forever.
    assert STATE_BOUNDS['tau'][1]<13
    for name,value in (('tau',B.rn32(F(11,10))),('sigma',B.rn32(F(1,100))),('rs',F(1,2))):
        assert STATE_BOUNDS[name][0]<value<STATE_BOUNDS[name][1]
    t2lo=B.mul(D.MIN_TAU,D.MIN_TAU); t2hi=B.mul(D.MAX_TAU,D.MAX_TAU)
    sablo=SPEC.SIGMA_AB_MIN; sabhi=B.div(F(4),D.SIGMA_COEFF)
    ulo=B.mul(B.mul(sablo,t2lo),t2lo); uhi=B.mul(B.mul(sabhi,t2hi),t2hi)
    assert POW_BASE[0]<ulo<=uhi<POW_BASE[1]
    # The actual newlib powf object is separately pinned at both compiled
    # exponents.  This establishes a finite positive output range for the
    # candidate graph; it does not claim the whole firmware call graph.
    assert pow_profile['pinned_tuner_powf_positive_finite_range_closed']
    assert tuple(pow_profile['actual_exponents']) == (SPEC.POW_EXPONENT, Q.EXPONENT)
    sqrtlo=W.M.SQRT.sqrt32(D.PSEUDO_MIN); sqrthi=W.M.SQRT.sqrt32(D.PSEUDO_MAX)
    assert F(1,16)<sqrtlo<=sqrthi<F(1,2)
    coeff=W.rounded_abs(D.RS_MSE_COEFF*POW_OUTPUT[1])
    powered=W.rounded_abs(coeff*POW_OUTPUT[1])
    rawrs=W.rounded_abs(powered/sqrtlo)
    # tau^3 is eagerly calculated by rs_target_from_law_ even in SpectralMSE.
    eager_tau3=W.rounded_abs(W.rounded_abs(D.MAX_TAU**2)*D.MAX_TAU)
    largest=max(rawrs,eager_tau3,frontend['largest_ordinary_intermediate_abs_upper'])
    assert rawrs<2**36 and largest<F((2**24-1)*2**104)
    return {'qualification':QUALIFICATION,'pow_profile':POW_PROFILE,
        'spectral_pow_base_interval':POW_BASE,'derived_spectral_pow_base_interval':(ulo,uhi),
        'spectral_pow_compiled_exponent':SPEC.POW_EXPONENT,
        'qeff_pow_argument':B.mul(F(2),Q.default_r_a_binary32()),
        'qeff_pow_compiled_exponent':Q.EXPONENT,'pow_output_interval':POW_OUTPUT,
        'sqrt_TS_interval':(sqrtlo,sqrthi),'actual_alpha_intervals':alpha_bounds,
        'candidate_state_intervals':dict(STATE_BOUNDS),'candidate_target_intervals':dict(TARGET_BOUNDS),
        'candidate_inductive_images':after,'strict_inductive_margins':margins,
        'raw_RS_before_clamp_upper':rawrs,'eager_tau3_upper':eager_tau3,
        'committed_tau_interval':STATE_BOUNDS['tau'],
        'committed_pseudo_period_interval':(D.PSEUDO_MIN,D.PSEUDO_MAX),
        'committed_stationary_aw_diagonal_upper':F(64),
        'committed_RS_diagonal_upper':F(10000),
        'candidate_arithmetic_totality_under_named_pow_range_closed':True,
        'pinned_powf_positive_finite_range_under_scalar_profile_closed':
            pow_profile['pinned_tuner_powf_positive_finite_range_closed'],
        'compiled_pow_exponents_identified_with_exact_rational_exponents':True,
        'all_finite_prefixes_from_literal_candidate_seeds_covered':True,
        'existing_tau_roundoff_predecessor_domain_inductively_closed':True,
        'startup_deadline_required':False,'exact_shadow_pow_residuals_discarded':False,
        'full_target_candidate_compiler_and_pow_profile_correspondence_closed':False,
        'startup_capture_and_MEKF_covariance_totality_proved_here':False,'storage_search_allowed':False}


def build(): return deepcopy(_build())


def require_state(state):
    from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as PRODUCT
    if not isinstance(state,PRODUCT.State): raise TypeError('persistent whole machine TuneState required')
    for name,(lo,hi) in STATE_BOUNDS.items():
        track=getattr(state,name)
        if not all(lo<=getattr(track,mode)<=hi for mode in ('separate','fma')):
            raise ValueError('carried '+name+' outside all-time candidate interval')


def require_inputs(cfg,*,tau,sigma,pow_result,sqrt_result):
    """Check source-produced operands before spectral arithmetic executes."""
    require_config(cfg)
    if not TARGET_BOUNDS['tau'][0]<=tau<=D.MAX_TAU or not TARGET_BOUNDS['sigma'][0]<=sigma<=4:
        raise ValueError('spectral operands detached from configured tau/sigma targets')
    tau2=B.mul(tau,tau); u=B.mul(B.mul(max(B.div(sigma,D.SIGMA_COEFF),SPEC.SIGMA_AB_MIN),tau2),tau2)
    ts=min(max(B.mul(D.PSEUDO_RATIO,tau),D.PSEUDO_MIN),D.PSEUDO_MAX)
    require_pow(u,SPEC.POW_EXPONENT,pow_result)
    if sqrt_result!=W.M.SQRT.sqrt32(ts):
        raise ValueError('spectral sqrt detached from same rounded cadence argument')


def require_result(previous,result):
    from tools.stability.ou3_alt_contraction import finite_tuner_machine_candidate_step as CAND
    if not isinstance(result,CAND.Result): raise TypeError('source-executed whole candidate result required')
    require_state(previous); require_state(result.product.state)
    for mode,alpha_range in (('separate',COMMON_ALPHA),('fma',COMMON_ALPHA)):
        m=getattr(result,mode); tau=m.common_alpha.step; spec=m.spectral_input.spectral.machine
        if tau.exp_profile!=LIBM.EXP_PROFILE or m.rs_alpha.machine.exp_profile!=LIBM.EXP_PROFILE:
            raise ValueError('qualified candidate lost target exponential profile')
        if not alpha_range[0]<=tau.alpha<=alpha_range[1] or not RS_ALPHA[0]<=m.rs_alpha.machine.alpha<=RS_ALPHA[1]:
            raise ValueError('same-source candidate alpha outside all-input enclosure')
        require_inputs(m.sigma_join.deployment_cfg,tau=tau.tau_target,
            sigma=m.sigma_join.machine.sigma_target,pow_result=spec.pow_witness.result,
            sqrt_result=spec.sqrt_witness.result)
    return build()
