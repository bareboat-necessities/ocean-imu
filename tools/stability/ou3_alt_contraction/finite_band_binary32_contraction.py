"""One-step binary32 AdaptiveWaveBandPass relation under FP contraction choices.

The shipping source contains several scalar ``a*b + c*d`` expressions.  Unlike
the tuner EMA, a single separate-vs-FMA bit is not sufficient: each eligible
multiply-add may be contracted independently.  This module therefore returns
the finite set of binary32 values produced by all local no-reassociation
contraction choices for one already-qualified coefficient pair ``q_low,q_high``.

It covers the signal state and the exact unit-white covariance recurrence
``p00,p01,p11``.  Coefficient production (corner clamps and expf), persistence
through many samples, and the target compiler's no-reassociation/contraction
contract remain separate obligations.  No one outcome is selected by theorem
assumption.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B

SOURCE=Path(__file__).resolve().parents[3]/'src/tuner/AdaptiveWaveBandPass.h'
ZERO=B.rn32(0); ONE=B.rn32(1); TWO=B.rn32(2)
QUALIFICATION='OU3_ALT_BAND_BINARY32_CONTRACTION_V1'


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must be actual binary32')
    return q


def _uniq(values): return tuple(sorted(set(F(v) for v in values)))


def _sum_products(a,b,c,d):
    """All no-reassociation outcomes of a*b + c*d with optional one FMA."""
    ab=B.mul(a,b); cd=B.mul(c,d)
    return _uniq((B.add(ab,cd),B.fma(a,b,cd),B.fma(c,d,ab)))


def _add_product(acc,c,d):
    """All outcomes of acc + c*d when the new product may contract."""
    cd=B.mul(c,d)
    return _uniq((B.add(acc,cd),B.fma(c,d,acc)))


@dataclass(frozen=True)
class State:
    lowpass_low:F=ZERO
    band:F=ZERO
    p00:F=ZERO
    p01:F=ZERO
    p11:F=ZERO
    ready:bool=False
    def __post_init__(self):
        for n in ('lowpass_low','band','p00','p01','p11'):
            object.__setattr__(self,n,_q(getattr(self,n),n))
        if not isinstance(self.ready,bool): raise TypeError('literal band ready flag required')
        if self.p00<0 or self.p11<0: raise ValueError('band covariance diagonal must be nonnegative')


@dataclass(frozen=True)
class Envelope:
    predecessor:State
    x:F
    q_low:F
    q_high:F
    alpha_low:F
    alpha_high:F
    high_passed:F
    lowpass_values:tuple
    band_values:tuple
    p00_values:tuple
    p01_values:tuple
    p11_values:tuple
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.predecessor,State): raise TypeError('binary32 band predecessor required')
        for n in ('x','q_low','q_high','alpha_low','alpha_high','high_passed'):
            object.__setattr__(self,n,_q(getattr(self,n),n))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong band contraction qualification')
        for n in ('lowpass_values','band_values','p00_values','p01_values','p11_values'):
            vals=_uniq(getattr(self,n)); object.__setattr__(self,n,vals)
            if not vals or not all(B.is_binary32(v) for v in vals): raise ValueError('band contraction outcome set must contain binary32 values')
        if min(self.p00_values)<0 or min(self.p11_values)<0: raise ValueError('post-clamp covariance diagonal cannot be negative')

    def accepts(self,state:State):
        if not isinstance(state,State): return False
        return (state.ready and state.lowpass_low in self.lowpass_values and state.band in self.band_values
                and state.p00 in self.p00_values and state.p01 in self.p01_values and state.p11 in self.p11_values)


def step(predecessor:State,*,x,q_low,q_high):
    """All local contraction outcomes after valid corners/exp have produced q's."""
    if not isinstance(predecessor,State): raise TypeError('binary32 band State required')
    x=_q(x,'band input'); ql=_q(q_low,'q_low'); qh=_q(q_high,'q_high')
    if not 0<ql<=1 or not 0<qh<=1: raise ValueError('band decay coefficients must lie in (0,1]')
    al=B.sub(ONE,ql); ah=B.sub(ONE,qh)
    low_values=_sum_products(ql,predecessor.lowpass_low,al,x)
    hp=B.mul(ql,B.sub(x,predecessor.lowpass_low))
    band_values=_sum_products(qh,predecessor.band,ah,hp)

    a00=ql
    a10=-B.mul(ah,ql)
    a11=qh
    b0=al
    b1=B.mul(ah,ql)

    # p00 = a00*a00*p00 + b0*b0.  The first a00*a00 multiply cannot
    # contract (no add); the final product may contract with either add term.
    c00=B.mul(a00,a00)
    p00_raw=_sum_products(c00,predecessor.p00,b0,b0)
    p00_values=_uniq(max(ZERO,v) for v in p00_raw)

    # p01 = a00*(a10*p00 + a11*p01) + b0*b1.
    inner=_sum_products(a10,predecessor.p00,a11,predecessor.p01)
    p01=[]
    for v in inner:
        p01.extend(_sum_products(a00,v,b0,b1))
    p01_values=_uniq(p01)

    # p11 = a10*a10*p00 + 2*a10*a11*p01 + a11*a11*p11 + b1*b1.
    # Preserve the source's left-associated multiplicative coefficients, then
    # enumerate optional contraction of each product into the addition chain.
    c1=B.mul(a10,a10)
    c2=B.mul(B.mul(TWO,a10),a11)
    c3=B.mul(a11,a11)
    first=_sum_products(c1,predecessor.p00,c2,predecessor.p01)
    second=[]
    for v in first: second.extend(_add_product(v,c3,predecessor.p11))
    third=[]
    for v in _uniq(second): third.extend(_add_product(v,b1,b1))
    p11_values=_uniq(max(ZERO,v) for v in third)

    return Envelope(predecessor,x,ql,qh,al,ah,hp,low_values,band_values,p00_values,p01_values,p11_values)


def choose(envelope:Envelope,*,lowpass_low,band,p00,p01,p11):
    """Bind an actual one-step machine successor to the finite contraction set."""
    if not isinstance(envelope,Envelope): raise TypeError('band contraction Envelope required')
    s=State(lowpass_low,band,p00,p01,p11,True)
    if not envelope.accepts(s): raise ValueError('machine band successor outside all legal local contraction outcomes')
    return s


def _source_shape_matches():
    s=SOURCE.read_text()
    needles=('lowpass_low_ = q_low * low_prev + alpha_low * x;',
      'const float high_passed = q_low * (x - low_prev);',
      'band_ = q_high * band_prev + alpha_high * high_passed;',
      'const float p00_new = a00 * a00 * p00_ + b0 * b0;',
      'const float p01_new = a00 * (a10 * p00_ + a11 * p01_) + b0 * b1;',
      'const float p11_new = a10 * a10 * p00_',
      '+ 2.0f * a10 * a11 * p01_', '+ a11 * a11 * p11_', '+ b1 * b1;')
    return all(x in s for x in needles)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_adaptive_band_signal_and_covariance_source_shape_matches':_source_shape_matches(),
      'one_step_all_local_no_reassociation_contraction_outcomes_enumerated':True,
      'no_single_global_FMA_bit_assumed_for_band_polynomial':True,
      'actual_machine_successor_can_be_bound_to_finite_contraction_set':True,
      'band_corner_and_exp_binary32_coefficient_production_closed':False,
      'target_compiler_no_reassociation_contract_qualified':False,
      'persistent_band_machine_history_composed':False,
      'band_noise_floor_sqrt_machine_relation_composed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
