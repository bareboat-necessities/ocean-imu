"""Binary32 deployment graph for the shipping SpectralMSE R_S target.

Shipping evaluates the cadence, sigma division/floor, tau products, one per-
candidate pow, one sqrt, and the final coefficient products/division in
binary32.  Ordinary operations are modeled exactly with the shared RNE kernel;
libm results remain explicit witnesses with exact same-argument error intervals.

The deployment entry consumes ``DeploymentConfig`` and therefore uses the
actual cached binary32 q_eff pow value.  The exact-real theorem path uses the
separate mathematical q_eff root interval; those two values are intentionally
not identified here.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C
from tools.stability.ou3_alt_contraction import finite_tuner_deployment_config as D
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
SIGMA_AB_MIN=B.rn32(F(1,10**6))
SIX=B.rn32(6); SEVEN=B.rn32(7); POW_EXPONENT=B.div(SIX,SEVEN)
QUALIFICATION='OU3_ALT_SPECTRAL_BINARY32_GRAPH_V2'


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must be an actual binary32 value')
    return q


def _positive(x,name):
    q=_q(x,name)
    if q<=0: raise ValueError(f'{name} must be positive')
    return q


@dataclass(frozen=True)
class LibmWitness:
    operation:str
    argument:F
    result:F
    exponent:F|None=None
    def __post_init__(self):
        if self.operation not in ('sqrt','pow'): raise ValueError('unsupported SpectralMSE libm operation')
        a=_positive(self.argument,'libm argument'); r=_positive(self.result,'libm result')
        object.__setattr__(self,'argument',a); object.__setattr__(self,'result',r)
        if self.operation=='sqrt':
            if self.exponent is not None: raise ValueError('sqrt witness consumes no exponent')
        else:
            if self.exponent is None: raise TypeError('pow witness requires compiled binary32 exponent')
            e=_positive(self.exponent,'pow exponent')
            if e!=POW_EXPONENT: raise ValueError('pow exponent detached from shipping 6.0f/7.0f evaluation')
            object.__setattr__(self,'exponent',e)


@dataclass(frozen=True)
class RootSupply:
    witness:LibmWitness
    true_lo:F
    true_hi:F
    error_lo:F
    error_hi:F
    def __post_init__(self):
        lo,hi,elo,ehi=map(F,(self.true_lo,self.true_hi,self.error_lo,self.error_hi))
        if lo<=0 or lo>hi or elo>ehi: raise ValueError('invalid root supply interval')
        if elo!=self.witness.result-hi or ehi!=self.witness.result-lo:
            raise ValueError('root error interval detached from same machine witness/root enclosure')
        object.__setattr__(self,'true_lo',lo); object.__setattr__(self,'true_hi',hi)
        object.__setattr__(self,'error_lo',elo); object.__setattr__(self,'error_hi',ehi)


@dataclass(frozen=True)
class Step:
    tau:F; sigma:F
    pseudo_ratio:F; pseudo_min:F; pseudo_max:F
    c_sigma:F; rs_mse_coeff:F; rs_qeff_pow:F
    requested_TS:F; TS:F; sigma_div:F; sigma_aB:F; tau2:F; u:F
    pow_witness:LibmWitness; sqrt_witness:LibmWitness
    pow_supply:RootSupply; sqrt_supply:RootSupply
    coeff_product:F; powered_product:F; raw_RS:F
    def __post_init__(self):
        for n in ('tau','sigma','pseudo_ratio','pseudo_min','pseudo_max','c_sigma','rs_mse_coeff','rs_qeff_pow',
                  'requested_TS','TS','sigma_div','sigma_aB','tau2','u','coeff_product','powered_product','raw_RS'):
            object.__setattr__(self,n,F(getattr(self,n)))
        if self.pow_witness.argument!=self.u or self.sqrt_witness.argument!=self.TS:
            raise ValueError('SpectralMSE libm witnesses detached from machine arithmetic arguments')
        if self.pow_supply.witness!=self.pow_witness or self.sqrt_supply.witness!=self.sqrt_witness:
            raise ValueError('SpectralMSE root supplies detached from witnesses')


def _root_supply(w:LibmWitness,*,bits:int=ROOT.DEFAULT_BITS):
    if w.operation=='sqrt': lo,hi=ROOT.sqrt_enclosure(w.argument,bits)
    else: lo,hi=ROOT.pow_6_7_enclosure(w.argument,bits)
    return RootSupply(w,lo,hi,w.result-hi,w.result-lo)


def step(*,tau,sigma,pseudo_ratio,pseudo_min,pseudo_max,c_sigma,
         rs_mse_coeff,rs_qeff_pow,pow_result,sqrt_result,bits:int=ROOT.DEFAULT_BITS):
    vals={n:_positive(v,n) for n,v in {
        'tau':tau,'pseudo_ratio':pseudo_ratio,'pseudo_min':pseudo_min,'pseudo_max':pseudo_max,
        'c_sigma':c_sigma,'rs_mse_coeff':rs_mse_coeff,'rs_qeff_pow':rs_qeff_pow}.items()}
    sig=_q(sigma,'sigma')
    if sig<0: raise ValueError('sigma must be nonnegative')
    if vals['pseudo_max']<vals['pseudo_min']: raise ValueError('invalid pseudo cadence bounds')
    requested=B.mul(vals['pseudo_ratio'],vals['tau'])
    TS=min(max(requested,vals['pseudo_min']),vals['pseudo_max'])
    sdiv=B.div(sig,vals['c_sigma']); sab=max(sdiv,SIGMA_AB_MIN)
    tau2=B.mul(vals['tau'],vals['tau'])
    u=B.mul(B.mul(sab,tau2),tau2)
    pw=LibmWitness('pow',u,F(pow_result),POW_EXPONENT)
    sw=LibmWitness('sqrt',TS,F(sqrt_result))
    ps=_root_supply(pw,bits=bits); ss=_root_supply(sw,bits=bits)
    cp=B.mul(vals['rs_mse_coeff'],vals['rs_qeff_pow'])
    pp=B.mul(cp,pw.result); rs=B.div(pp,sw.result)
    return Step(vals['tau'],sig,vals['pseudo_ratio'],vals['pseudo_min'],vals['pseudo_max'],
                vals['c_sigma'],vals['rs_mse_coeff'],vals['rs_qeff_pow'],requested,TS,sdiv,sab,tau2,u,
                pw,sw,ps,ss,cp,pp,rs)


def step_from_config(cfg:C.CandidateConfig,*,tau,sigma,pow_result,sqrt_result,bits:int=ROOT.DEFAULT_BITS):
    """Legacy component entry; not the theorem representation of shipping defaults."""
    if not isinstance(cfg,C.CandidateConfig): raise TypeError('CandidateConfig required')
    return step(tau=tau,sigma=sigma,pseudo_ratio=cfg.pseudo_tau_ratio,
                pseudo_min=cfg.pseudo_min,pseudo_max=cfg.pseudo_max,c_sigma=cfg.sigma_coeff,
                rs_mse_coeff=cfg.rs_mse_coeff,rs_qeff_pow=cfg.qeff_pow,
                pow_result=pow_result,sqrt_result=sqrt_result,bits=bits)


def step_from_deployment_config(cfg:D.DeploymentConfig,*,tau,sigma,pow_result,sqrt_result,bits:int=ROOT.DEFAULT_BITS):
    """Shipping machine entry using the explicit produced binary32 qeff cache."""
    if not isinstance(cfg,D.DeploymentConfig): raise TypeError('DeploymentConfig required')
    return step(tau=tau,sigma=sigma,pseudo_ratio=cfg.pseudo_tau_ratio,
                pseudo_min=cfg.pseudo_min,pseudo_max=cfg.pseudo_max,c_sigma=cfg.sigma_coeff,
                rs_mse_coeff=cfg.rs_mse_coeff,rs_qeff_pow=cfg.qeff_pow,
                pow_result=pow_result,sqrt_result=sqrt_result,bits=bits)


def _source_shape_matches():
    s=SOURCE.read_text()
    needles=(
      'return std::min(std::max(pseudo_update_tau_ratio_ * tau,',
      'const float sigma_aB = std::max(sigma / c_sigma, 1e-6f);',
      'const float tau2 = tau * tau;',
      'const float u = sigma_aB * tau2 * tau2;',
      'std::pow(u, 6.0f / 7.0f)', '/ std::sqrt(TS);',
      'rs_qeff_pow_ = std::pow(2.0f * r_a, 1.0f / 14.0f);')
    return all(x in s for x in needles)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_SpectralMSE_source_shape_matches':_source_shape_matches(),
      'pseudo_cadence_binary32_mul_and_clamp_materialized':True,
      'sigma_div_floor_tau2_u_binary32_graph_materialized':True,
      'pow_and_sqrt_kept_as_distinct_same_argument_binary32_witnesses':True,
      'final_coefficient_multiply_multiply_divide_graph_materialized':True,
      'deployment_entry_consumes_explicit_produced_qeff_binary32_cache':True,
      'deployment_machine_qeff_cache_not_identified_with_exact_real_qeff_root':True,
      'pow_witness_error_interval_against_exact_root_of_same_machine_u_exposed':True,
      'sqrt_witness_error_interval_against_exact_root_of_same_machine_TS_exposed':True,
      'target_libm_pow_correspondence_closed':False,
      'target_libm_sqrt_correspondence_closed':False,
      'cached_qeff_pow_target_libm_correspondence_closed':False,
      'source_uniform_root_supply_bounds_closed':False,
      'binary32_SpectralMSE_target_correspondence_closed':False,
      'source_uniform_complete_startup_reachability_closed':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
