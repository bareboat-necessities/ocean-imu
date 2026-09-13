"""Exact rational enclosures for the real SpectralMSE sqrt/pow roots.

The older finite tuner shadow encoded the two nonlinear SpectralMSE roots by
requiring rational witnesses ``r`` with

    r^2 = T_S
    p^7 = u^6.

That is useful for hand-picked algebra fixtures but is not a universal proof:
ordinary shipping cells (including the 0.2 Hz startup prior) make these roots
irrational.  This module removes that artificial restriction at the REAL
arithmetic layer.  It brackets both roots by exact dyadic rationals using only
integer/Fraction comparisons, then propagates those intervals monotonically to
an exact-real R_S target enclosure.

No libm result is accepted here and no float value is substituted for an
irrational root.  Binary32 ``sqrt``/``pow`` correspondence remains a separate
deployment obligation.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_tuner_candidate as C

DEFAULT_BITS=96
QUALIFICATION='OU3_ALT_SPECTRAL_REAL_ROOT_ENCLOSURE_V1'


def _root_enclosure(x,n:int,bits:int=DEFAULT_BITS):
    x=F(x)
    if x<0 or not isinstance(n,int) or n<=0: raise ValueError('nonnegative rational and positive root degree required')
    if not isinstance(bits,int) or bits<16: raise ValueError('at least 16 bisection bits required')
    if x==0: return F(0),F(0)
    # First find a compact dyadic upper bracket.  Starting at x itself is sound
    # but destroys useful bisection bits for seventh roots of large x=u^6.
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
    TS:F
    u:F
    sqrt_TS_lo:F
    sqrt_TS_hi:F
    u_pow_6_7_lo:F
    u_pow_6_7_hi:F
    raw_RS_lo:F
    raw_RS_hi:F
    target_RS_lo:F
    target_RS_hi:F
    bits:int=DEFAULT_BITS
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.target,C.TargetState): raise TypeError('TargetState required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong spectral enclosure qualification')
        for n in ('TS','u','sqrt_TS_lo','sqrt_TS_hi','u_pow_6_7_lo','u_pow_6_7_hi','raw_RS_lo','raw_RS_hi','target_RS_lo','target_RS_hi'):
            object.__setattr__(self,n,F(getattr(self,n)))
        if not (0<self.sqrt_TS_lo<=self.sqrt_TS_hi and 0<=self.u_pow_6_7_lo<=self.u_pow_6_7_hi):
            raise ValueError('invalid root enclosure')
        if not self.sqrt_TS_lo**2<=self.TS<=self.sqrt_TS_hi**2:
            raise ValueError('sqrt(T_S) enclosure detached')
        if not self.u_pow_6_7_lo**7<=self.u**6<=self.u_pow_6_7_hi**7:
            raise ValueError('u^(6/7) enclosure detached')
        if not 0<=self.raw_RS_lo<=self.raw_RS_hi:
            raise ValueError('invalid raw R_S enclosure')
        if not self.target_RS_lo<=self.target_RS_hi:
            raise ValueError('invalid clamped R_S enclosure')


def enclose(cfg:C.CandidateConfig,target:C.TargetState,*,bits:int=DEFAULT_BITS):
    if not isinstance(cfg,C.CandidateConfig) or not isinstance(target,C.TargetState):
        raise TypeError('CandidateConfig and TargetState required')
    TS=C.clamp(cfg.pseudo_tau_ratio*target.tau_target,cfg.pseudo_min,cfg.pseudo_max)
    sigma_aB=max(target.sigma_target/cfg.sigma_coeff,C.SIGMA_AB_MIN)
    u=sigma_aB*target.tau_target**4
    slo,shi=sqrt_enclosure(TS,bits); plo,phi=pow_6_7_enclosure(u,bits)
    k=cfg.rs_mse_coeff*cfg.qeff_pow
    raw_lo=k*plo/shi; raw_hi=k*phi/slo
    if cfg.clamp_enabled:
        out_lo=C.clamp(raw_lo,cfg.min_RS,cfg.max_RS)
        out_hi=C.clamp(raw_hi,cfg.min_RS,cfg.max_RS)
    else:
        out_lo,out_hi=raw_lo,raw_hi
    return Enclosure(target,TS,u,slo,shi,plo,phi,raw_lo,raw_hi,out_lo,out_hi,bits)


def contains_legacy_exact_witness(cfg:C.CandidateConfig,target:C.TargetState,w:C.SpectralWitness,box:Enclosure):
    """Regression bridge: exact-rational legacy witnesses are contained when they exist."""
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
      'irrational_roots_do_not_require_fake_rational_equalities':True,
      'legacy_exact_rational_spectral_cells_embed_in_interval_theorem':True,
      'binary32_sqrt_target_libm_correspondence_closed':False,
      'binary32_pow_target_libm_correspondence_closed':False,
      'finite_tuner_candidate_interval_recurrence_composed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
