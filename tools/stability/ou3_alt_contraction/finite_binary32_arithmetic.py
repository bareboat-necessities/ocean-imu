"""Exact finite IEEE-754 binary32 arithmetic for ALT deployment graphs.

The ALT proof already had a positive-only round-to-nearest-even helper for the
wrapper clock.  Tuner and covariance recurrences need signed differences and
must also distinguish an ordinary multiply/add evaluation from a contracted
multiply-add.  This module provides that small arithmetic kernel over exact
Fractions.

The finite lattice includes gradual underflow: admitted source coordinates and
operation residuals have no positive lower bound. Infinities, NaNs, division by
zero and target-specific exception behaviour still require qualification.
Fractions identify the two signed zeros; use finite_binary32_mahony's bit ledger
for sign-sensitive operations (including copysign, atan2 and division by zero).
Target gradual-underflow/flush-to-zero correspondence is a separate obligation.
"""
from __future__ import annotations
from fractions import Fraction as F


def _pow2(e:int)->F:
    return F(1<<e) if e>=0 else F(1,1<<(-e))


def _floor_log2_positive(x:F)->int:
    if x<=0: raise ValueError('positive input required')
    n,d=x.numerator,x.denominator
    e=n.bit_length()-d.bit_length()
    while _pow2(e)>x: e-=1
    while _pow2(e+1)<=x: e+=1
    return e


def _rne_nonnegative_integer(x:F)->int:
    if x<0: raise ValueError('nonnegative input required')
    q,r=divmod(x.numerator,x.denominator); twice=2*r
    if twice<x.denominator: return q
    if twice>x.denominator: return q+1
    return q if q%2==0 else q+1


def rn32(x)->F:
    """Round an exact rational to the finite binary32 lattice, ties to even.

    Below the smallest normal, the lattice spacing remains 2**-149.  A value
    below half that spacing rounds to zero; a midpoint rounds to the even
    integer lattice index, including the subnormal/normal boundary.
    """
    x=F(x)
    if x==0: return F(0)
    sign=-1 if x<0 else 1; a=abs(x)
    e=_floor_log2_positive(a)
    if e < -126:
        quantum=_pow2(-149)
        return sign*F(_rne_nonnegative_integer(a/quantum))*quantum
    if e > 127:
        raise ValueError('binary32 overflow')
    quantum=_pow2(e-23); m=_rne_nonnegative_integer(a/quantum)
    if m==(1<<24):
        m>>=1; e+=1
        if e>127: raise ValueError('binary32 overflow')
        quantum=_pow2(e-23)
    if not (1<<23)<=m<(1<<24):
        raise AssertionError('binary32 significand normalization failed')
    return sign*F(m)*quantum


def is_binary32(x)->bool:
    x=F(x)
    if x==0: return True
    try: return rn32(x)==x
    except ValueError: return False


def add(a,b): return rn32(F(a)+F(b))
def sub(a,b): return rn32(F(a)-F(b))
def mul(a,b): return rn32(F(a)*F(b))
def div(a,b):
    b=F(b)
    if b==0: raise ZeroDivisionError('binary32 division by zero')
    return rn32(F(a)/b)

def fma(a,b,c): return rn32(F(a)*F(b)+F(c))


def ema(previous,target,alpha,*,contracted:bool):
    """Shipping shape: previous += alpha * (target - previous)."""
    if not isinstance(contracted,bool): raise TypeError('literal FP contraction branch required')
    p,t,a=map(F,(previous,target,alpha))
    if not all(is_binary32(x) for x in (p,t,a)):
        raise ValueError('EMA operands must be deployed binary32 values')
    delta=sub(t,p)
    return fma(a,delta,p) if contracted else add(p,mul(a,delta))


def readiness():
    return {
      'signed_normal_binary32_RNE_exact':True,
      'finite_binary32_gradual_underflow_RNE_exact':True,
      'subnormal_API_values_and_arithmetic_results_retained':True,
      'two_signed_zeros_identified_in_rational_value_projection':True,
      'binary32_add_sub_mul_div_exact_under_nonexceptional_scope':True,
      'separate_mul_add_EMA_shape_materialized':True,
      'contracted_fma_EMA_shape_materialized':True,
      'compiler_FP_contraction_mode_qualified':False,
      'target_gradual_underflow_correspondence_qualified':False,
      'subnormal_nan_inf_scope_closed':False,
      'ALT_LIVE_PASS':False,
    }
