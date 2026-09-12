"""Rigorous real enclosures for shipping exp/expm1 prediction roots.

Shipping evaluates exp and expm1 separately.  The theorem must therefore not
identify their rounded results bit-for-bit, but neither result may be a free
coefficient.  For 0 <= x <= 1 the alternating series gives

    1-x <= exp(-x) <= 1-x+x^2/2,
       -x <= expm1(-x) <= -x+x^2/2.

The current OU shipping argument satisfies x=h/tau <= 0.25 on the canonical
5 ms word because tau is clamped at >=0.02 s.  The active BA arguments are far
smaller (h/5000 and 2h/5000).  These real inequalities therefore bind the two
separate shipping transcendental outputs to their SAME source-owned arguments.

This is deliberately not a binary32/libm correspondence theorem.  The final
finite-precision proof must enclose the implementation's correctly/incorrectly
rounded library results inside a deployment residual; no equality expm1=exp-1
is assumed here.
"""
from __future__ import annotations
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M

ONE=F(1)


def enclosure(x):
    x=M.rational(x)
    if x < 0 or x > 1:
        raise ValueError('exp source argument outside certified [0,1] enclosure')
    lo=ONE-x
    hi=lo+x*x/F(2)
    return lo,hi,lo-ONE,hi-ONE


def validate_exp_expm1(x, exp_value, expm1_value):
    lo,hi,mlo,mhi=enclosure(x)
    e=M.rational(exp_value); m=M.rational(expm1_value)
    if not lo <= e <= hi:
        raise ValueError('exp root detached from source-owned argument')
    if not mlo <= m <= mhi:
        raise ValueError('expm1 root detached from source-owned argument')
    return e,m


def validate_ou(decay):
    # Avoid importing OU here to keep this primitive acyclic; structural fields
    # are checked by finite_ou_runtime_primitives before this theorem layer.
    x=M.rational(decay.h)/M.rational(decay.tau)
    validate_exp_expm1(x,decay.alpha,decay.em1)
    return decay


def validate_bias(decay, h):
    if not decay.active:
        if decay.phi_b != 1 or decay.em1_2 != 0:
            raise ValueError('held BA branch must consume no transcendental decay')
        return decay
    h=M.rational(h); tau=M.rational(decay.tau_b)
    # Mean uses exp(-h/tau_b); Q uses expm1(-2h/tau_b), so these are distinct
    # arguments and are checked separately rather than algebraically equated.
    elo,ehi,_,_=enclosure(h/tau)
    phi=M.rational(decay.phi_b)
    if not elo <= phi <= ehi:
        raise ValueError('BA exp root detached from source-owned argument')
    _,_,mlo,mhi=enclosure(2*h/tau)
    em=M.rational(decay.em1_2)
    if not mlo <= em <= mhi:
        raise ValueError('BA expm1 root detached from source-owned argument')
    return decay


def readiness():
    return {
      'OU_exp_root_real_enclosed_at_same_h_over_tau':True,
      'OU_expm1_root_real_enclosed_at_same_h_over_tau':True,
      'BA_exp_and_expm1_distinct_arguments_real_enclosed':True,
      'exp_expm1_bit_identity_assumed':False,
      'binary32_libm_correspondence_closed':False,
      'deployment_roundoff_supply_attached':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
