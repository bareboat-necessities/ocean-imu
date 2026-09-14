"""Binary32 production topology for the cached SpectralMSE q_eff^(1/14).

Shipping initializes/refreshed ``rs_qeff_pow_`` as

    qeff = 2.0f * r_a
    rs_qeff_pow_ = pow(qeff, 1.0f / 14.0f)

with ``r_a`` the configured reduced acceleration-noise density.  This cache is
used by every SpectralMSE target, so it may not remain a free theorem constant.

The ordinary multiply/divide operations are materialized exactly in binary32.
The cached pow result is retained as an explicit positive binary32 witness, and
its error interval is measured against the exact fourteenth root of the SAME
rounded qeff argument using exact rational bisection.  Target-libm correctness
remains open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_real_enclosure as ROOT

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
TWO=B.rn32(2); ONE=B.rn32(1); FOURTEEN=B.rn32(14); EXPONENT=B.div(ONE,FOURTEEN)
QUALIFICATION='OU3_ALT_QEFF_CACHE_BINARY32_V1'


def _root14(x,bits=ROOT.DEFAULT_BITS):
    return ROOT._root_enclosure(F(x),14,bits)


@dataclass(frozen=True)
class CacheWitness:
    r_a:F
    qeff_argument:F
    exponent:F
    result:F
    true_lo:F
    true_hi:F
    error_lo:F
    error_hi:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        for n in ('r_a','qeff_argument','exponent','result','true_lo','true_hi','error_lo','error_hi'):
            object.__setattr__(self,n,F(getattr(self,n)))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong qeff cache qualification')
        if not B.is_binary32(self.r_a) or self.r_a<=0: raise ValueError('qeff r_a must be positive binary32')
        if self.qeff_argument!=B.mul(TWO,self.r_a): raise ValueError('qeff pow argument detached from binary32 2*r_a')
        if self.exponent!=EXPONENT: raise ValueError('qeff pow exponent detached from binary32 1/14')
        if not B.is_binary32(self.result) or self.result<=0: raise ValueError('cached qeff pow result must be positive binary32')
        if not 0<self.true_lo<=self.true_hi: raise ValueError('invalid exact qeff root enclosure')
        if not self.true_lo**14<=self.qeff_argument<=self.true_hi**14:
            raise ValueError('qeff exact root enclosure detached from same argument')
        if self.error_lo!=self.result-self.true_hi or self.error_hi!=self.result-self.true_lo:
            raise ValueError('qeff cache error interval detached from witness')


def produce(*,r_a,pow_result,bits:int=ROOT.DEFAULT_BITS):
    ra=F(r_a); result=F(pow_result)
    if not B.is_binary32(ra) or ra<=0: raise ValueError('qeff r_a must be positive binary32')
    if not B.is_binary32(result) or result<=0: raise ValueError('cached qeff pow result must be positive binary32')
    arg=B.mul(TWO,ra); lo,hi=_root14(arg,bits)
    return CacheWitness(ra,arg,EXPONENT,result,lo,hi,result-hi,result-lo)


def default_r_a_binary32():
    """Literal source-order default: 0.0148f * 0.0148f * (1.0f/200.0f)."""
    noise=B.rn32(F(148,10000)); dt=B.div(B.rn32(1),B.rn32(200))
    return B.mul(B.mul(noise,noise),dt)


def qualify_candidate_cache(candidate_qeff_pow,cache:CacheWitness):
    if not isinstance(cache,CacheWitness): raise TypeError('CacheWitness required')
    q=F(candidate_qeff_pow)
    if not B.is_binary32(q): raise ValueError('candidate qeff cache must be actual binary32')
    if q!=cache.result: raise ValueError('candidate qeff cache detached from source-produced cached pow result')
    return cache


def _source_shape_matches():
    s=SOURCE.read_text()
    return all(x in s for x in (
      'constexpr float R_S_ACCEL_NOISE_DENSITY_DEFAULT = 0.0148f * 0.0148f * FREQ_SMOOTHER_DT;',
      'rs_qeff_pow_ = std::pow(2.0f * r_a, 1.0f / 14.0f);',
      'float rs_qeff_pow_            =',
      'std::pow(2.0f * R_S_ACCEL_NOISE_DENSITY_DEFAULT, 1.0f / 14.0f);'))


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_qeff_cache_source_shape_matches':_source_shape_matches(),
      'default_r_a_source_order_binary32_graph_materialized':True,
      'qeff_argument_binary32_2_times_r_a_materialized':True,
      'compiled_binary32_one_over_fourteen_exponent_materialized':True,
      'cached_pow_result_bound_to_same_machine_argument':True,
      'cached_pow_error_interval_against_exact_fourteenth_root_exposed':True,
      'candidate_qeff_cache_can_be_qualified_against_produced_cache':True,
      'target_libm_qeff_pow_correspondence_closed':False,
      'source_uniform_qeff_cache_supply_bound_closed':False,
      'binary32_SpectralMSE_target_correspondence_closed':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
