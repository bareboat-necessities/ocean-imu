"""Exact rational enclosures for the real SpectralMSE nonlinear roots.

The legacy exact tuner uses rational algebraic witnesses for sqrt(T_S),
u^(6/7), and indirectly q_eff^(1/14). Ordinary shipping cells make these
roots irrational, while deployed q_eff is also cached as a binary32 libm result.
This module keeps the mathematical roots and the machine cache distinct.

``enclose`` remains the exact-rational component path. ``enclose_deployment``
consumes the deployment-qualified interface and uses its rigorous interval for
the exact (2*r_a)^(1/14) root.  The higher candidate layer enforces the concrete
DeploymentConfig class; this lower root primitive checks its qualification token
to avoid a root->deployment->cache->root import cycle.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C

DEFAULT_BITS=96
QUALIFICATION='OU3_ALT_SPECTRAL_REAL_ROOT_ENCLOSURE_V2'
DEPLOYMENT_QUALIFICATION='OU3_ALT_TUNER_DEPLOYMENT_CONFIG_V1'


def _root_enclosure(x,n:int,bits:int=DEFAULT_BITS):
    x=F(x)
    if x<0 or not isinstance(n,int) or n<=0: raise ValueError('nonnegative rational and positive root degree required')
    if not isinstance(bits,int) or bits<16: raise ValueError('at least 16 bisection bits required')
    if x==0: return F(0),F(0)
    lo=F(0); hi=F(1)
    while hi**n < x: hi*=2
    for _ in range(bits):
        mid=(lo+hi)/2
        if mid**n <= x: lo=mid
        else: hi=mid
    if not lo**n<=x<=hi**n: raise AssertionError('root enclosure lost exact containment')
    return lo,hi


def sqrt_enclosure(x,bits:int=DEFAULT_BITS): return _root_enclosure(x,2,bits)
def pow_6_7_enclosure(u,bits:int=DEFAULT_BITS):
    u=F(u)
    if u<0: raise ValueError('SpectralMSE u must be nonnegative')
    return _root_enclosure(u**6,7,bits)


@dataclass(frozen=True)
class Enclosure:
    target:C.TargetState
    TS:F; u:F
    sqrt_TS_lo:F; sqrt_TS_hi:F
    u_pow_6_7_lo:F; u_pow_6_7_hi:F
    raw_RS_lo:F; raw_RS_hi:F
    target_RS_lo:F; target_RS_hi:F
    qeff_root_lo:F|None=None
    qeff_root_hi:F|None=None
    bits:int=DEFAULT_BITS
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.target,C.TargetState): raise TypeError('TargetState required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong spectral enclosure qualification')
        for n in ('TS','u','sqrt_TS_lo','sqrt_TS_hi','u_pow_6_7_lo','u_pow_6_7_hi','raw_RS_lo','raw_RS_hi','target_RS_lo','target_RS_hi'):
            object.__setattr__(self,n,F(getattr(self,n)))
        if self.qeff_root_lo is not None: object.__setattr__(self,'qeff_root_lo',F(self.qeff_root_lo))
        if self.qeff_root_hi is not None: object.__setattr__(self,'qeff_root_hi',F(self.qeff_root_hi))
        if not (0<self.sqrt_TS_lo<=self.sqrt_TS_hi and 0<=self.u_pow_6_7_lo<=self.u_pow_6_7_hi):
            raise ValueError('invalid root enclosure')
        if not self.sqrt_TS_lo**2<=self.TS<=self.sqrt_TS_hi**2: raise ValueError('sqrt(T_S) enclosure detached')
        if not self.u_pow_6_7_lo**7<=self.u**6<=self.u_pow_6_7_hi**7: raise ValueError('u^(6/7) enclosure detached')
        if (self.qeff_root_lo is None)!=(self.qeff_root_hi is None): raise ValueError('qeff root interval must be carried as a pair')
        if self.qeff_root_lo is not None and not 0<self.qeff_root_lo<=self.qeff_root_hi: raise ValueError('invalid qeff root interval')
        if not 0<=self.raw_RS_lo<=self.raw_RS_hi: raise ValueError('invalid raw R_S enclosure')
        if not self.target_RS_lo<=self.target_RS_hi: raise ValueError('invalid clamped R_S enclosure')


def _geometry(cfg,target,bits):
    TS=C.clamp(cfg.pseudo_tau_ratio*target.tau_target,cfg.pseudo_min,cfg.pseudo_max)
    sigma_aB=max(target.sigma_target/cfg.sigma_coeff,C.SIGMA_AB_MIN)
    u=sigma_aB*target.tau_target**4
    slo,shi=sqrt_enclosure(TS,bits); plo,phi=pow_6_7_enclosure(u,bits)
    return TS,u,slo,shi,plo,phi


def enclose(cfg:C.CandidateConfig,target:C.TargetState,*,bits:int=DEFAULT_BITS):
    if not isinstance(cfg,C.CandidateConfig) or not isinstance(target,C.TargetState): raise TypeError('CandidateConfig and TargetState required')
    TS,u,slo,shi,plo,phi=_geometry(cfg,target,bits)
    k=cfg.rs_mse_coeff*cfg.qeff_pow
    raw_lo=k*plo/shi; raw_hi=k*phi/slo
    if cfg.clamp_enabled:
        out_lo=C.clamp(raw_lo,cfg.min_RS,cfg.max_RS); out_hi=C.clamp(raw_hi,cfg.min_RS,cfg.max_RS)
    else: out_lo,out_hi=raw_lo,raw_hi
    return Enclosure(target,TS,u,slo,shi,plo,phi,raw_lo,raw_hi,out_lo,out_hi,None,None,bits)


def enclose_deployment(cfg,target:C.TargetState,*,bits:int=DEFAULT_BITS):
    """Deployment exact-real relation; cached binary32 qeff is not used here."""
    if getattr(cfg,'qualification',None)!=DEPLOYMENT_QUALIFICATION or not isinstance(target,C.TargetState):
        raise TypeError('deployment-qualified tuner config and TargetState required')
    try: qlo,qhi=cfg.qeff_exact_root_interval
    except (AttributeError,TypeError,ValueError) as exc:
        raise TypeError('deployment config lost exact qeff root interval') from exc
    qlo,qhi=F(qlo),F(qhi)
    if not 0<qlo<=qhi: raise ValueError('invalid deployment qeff root interval')
    TS,u,slo,shi,plo,phi=_geometry(cfg,target,bits)
    raw_lo=cfg.rs_mse_coeff*qlo*plo/shi
    raw_hi=cfg.rs_mse_coeff*qhi*phi/slo
    if cfg.clamp_enabled:
        out_lo=C.clamp(raw_lo,cfg.min_RS,cfg.max_RS); out_hi=C.clamp(raw_hi,cfg.min_RS,cfg.max_RS)
    else: out_lo,out_hi=raw_lo,raw_hi
    return Enclosure(target,TS,u,slo,shi,plo,phi,raw_lo,raw_hi,out_lo,out_hi,qlo,qhi,bits)


def contains_legacy_exact_witness(cfg:C.CandidateConfig,target:C.TargetState,w:C.SpectralWitness,box:Enclosure):
    if not isinstance(w,C.SpectralWitness) or not isinstance(box,Enclosure): raise TypeError('spectral witness/enclosure required')
    if box.target!=target: raise ValueError('enclosure detached from target')
    exact=C.spectral_RS(cfg,target,w)
    return (box.sqrt_TS_lo<=w.sqrt_TS<=box.sqrt_TS_hi and
            box.u_pow_6_7_lo<=w.u_pow_6_7<=box.u_pow_6_7_hi and
            box.target_RS_lo<=exact<=box.target_RS_hi)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'sqrt_TS_real_root_enclosed_by_exact_rationals':True,
      'u_pow_6_7_real_root_enclosed_by_exact_rationals':True,
      'SpectralMSE_real_RS_target_interval_propagated_monotonically':True,
      'deployment_exact_qeff_fourteenth_root_interval_propagated_monotonically':True,
      'deployment_exact_qeff_root_not_replaced_by_cached_binary32_value':True,
      'irrational_roots_do_not_require_fake_rational_equalities':True,
      'legacy_exact_rational_spectral_cells_embed_in_interval_theorem':True,
      'binary32_qeff_pow_target_libm_correspondence_closed':False,
      'binary32_sqrt_target_libm_correspondence_closed':False,
      'binary32_pow_target_libm_correspondence_closed':False,
      'finite_tuner_candidate_interval_recurrence_composed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
