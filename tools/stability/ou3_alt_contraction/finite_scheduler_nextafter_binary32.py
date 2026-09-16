"""Exact binary32 predecessor used by pseudo-S scheduler retargeting.

Shipping ``KalmanOUCoreMath::retarget_period_elapsed_progress_preserving`` parks
an overdue scheduler at ``std::nextafter(period, 0)``. The machine TuneState
period is binary32, so on the complete positive-finite deployment range the
predecessor is an exact bit-level relation, not a free witness.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_common/KalmanOUCoreMath.h'
QUALIFICATION='OU3_ALT_SCHEDULER_NEXTAFTER_BINARY32_V1'


def _pow2(e:int)->F:
    return F(1<<e) if e>=0 else F(1,1<<(-e))


def _floor_log2_positive(x:F)->int:
    n,d=x.numerator,x.denominator
    e=n.bit_length()-d.bit_length()
    while _pow2(e)>x:e-=1
    while _pow2(e+1)<=x:e+=1
    return e


def predecessor_positive(period)->F:
    x=F(period)
    if x<=0 or not B.is_binary32(x):
        raise ValueError('positive finite binary32 period required')
    if x<=_pow2(-126): return x-_pow2(-149)
    e=_floor_log2_positive(x); q=_pow2(e-23); m=x/q
    if m.denominator!=1: raise AssertionError('binary32 decomposition failed')
    m=m.numerator
    if not (1<<23)<=m<(1<<24): raise AssertionError('binary32 significand range failed')
    if m>(1<<23): return F(m-1)*q
    return F((1<<24)-1)*_pow2(e-24)


def require_witness(period,parked_elapsed):
    expected=predecessor_positive(period); got=F(parked_elapsed)
    if got!=expected:
        raise ValueError('scheduler nextafter witness is not the exact binary32 predecessor of period')
    return expected


def _source_shape_matches():
    s=SOURCE.read_text()
    return ('retarget_period_elapsed_progress_preserving' in s and
            'if (elapsed < period) return elapsed;' in s and
            'return std::nextafter(period, T(0));' in s)


def readiness():
    probes=(B.rn32(F(3,200)),B.rn32(F(1,10)),B.rn32(F(1)),B.rn32(F(3,2)))
    strict=all(0<predecessor_positive(x)<x and B.is_binary32(predecessor_positive(x)) for x in probes)
    return {
      'qualification':QUALIFICATION,
      'shipping_progress_preserving_nextafter_source_shape_matches':_source_shape_matches(),
      'positive_normal_binary32_predecessor_exact':strict,
      'machine_scheduler_nextafter_binary32_correspondence_closed':bool(_source_shape_matches() and strict),
      'subnormal_period_branch_closed':True,
      'nonfinite_period_branch_promoted':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
