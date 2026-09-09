#!/usr/bin/env python3
"""Dependency-reduced small-x scaled OU covariance for P4 tube consumers.

BASE computes qbar(x) and then divides entries by powers of x.  On the small-x
Maclaurin branch this is algebraically exact but interval-dependency-expensive:
terms such as x^7/x^6 are evaluated as independent interval occurrences.  For
very narrow x~4.17e-4 cells that can hide the strict SPD lower bound and trigger
an exponential subdivision tree.

Here the division is performed SYMBOLICALLY on the existing Maclaurin terms.
The resulting polynomial is exactly the same truncation with each exponent
shifted by the similarity-scaling power.  No source domain or OU model changes.
Production callers can temporarily install this equivalent evaluator into the
BASE proof module while building the tube.
"""
from __future__ import annotations
import contextlib
import ou3_brmm_riccati_tube as BASE

I=BASE.I; poly=BASE.poly

def _p(x,*terms):return poly(x,tuple(terms))

def small_scaled(x):
    qvv=_p(x,(1,2/3),(2,-1/2),(3,7/30),(4,-1/12),(5,31/1260),(6,-1/160),(7,127/90720))
    qvp=_p(x,(1,1/4),(2,-1/6),(3,5/72),(4,-1/45),(5,17/2880),(6,-41/30240))
    qva=_p(x,(1,1),(2,-1),(3,7/12),(4,-1/4),(5,31/360),(6,-1/40),(7,127/20160),(8,-17/12096))
    qpp=_p(x,(1,1/10),(2,-1/18),(3,5/252),(4,-1/180),(5,17/12960))
    qpa=_p(x,(1,1/3),(2,-1/3),(3,11/60),(4,-13/180),(5,19/840),(6,-1/168),(7,247/181440))
    qaa=_p(x,(1,2),(2,-2),(3,4/3),(4,-2/3),(5,4/15),(6,-4/45),(7,8/315),(8,-2/315),(9,4/2835))
    qvS=_p(x,(1,1/15),(2,-1/24),(3,41/2520),(4,-7/1440),(5,109/90720))
    qpS=_p(x,(1,1/36),(2,-1/72),(3,13/2880),(4,-1/864))
    qSS=_p(x,(1,1/126),(2,-1/288),(3,13/12960))
    qSa=_p(x,(1,1/12),(2,-1/12),(3,2/45),(4,-1/60),(5,11/2240),(6,-73/60480))
    return [[qvv,qvp,qvS,qva],[qvp,qpp,qpS,qpa],[qvS,qpS,qSS,qSa],[qva,qpa,qSa,qaa]]

def step_scaled_q(x):
    if x.hi < BASE.BRANCH_X:return small_scaled(x)
    return BASE._ORIGINAL_STEP_SCALED_Q_FOR_P4(x) if hasattr(BASE,'_ORIGINAL_STEP_SCALED_Q_FOR_P4') else BASE.step_scaled_q(x)

@contextlib.contextmanager
def installed(max_depth=14):
    old_step=BASE.step_scaled_q;old_depth=BASE.MAX_X_SPLIT_DEPTH
    # Keep a stable fallback for non-small cells without recursive self-call.
    BASE._ORIGINAL_STEP_SCALED_Q_FOR_P4=old_step
    BASE.step_scaled_q=step_scaled_q;BASE.MAX_X_SPLIT_DEPTH=int(max_depth)
    try:yield BASE
    finally:
        BASE.step_scaled_q=old_step;BASE.MAX_X_SPLIT_DEPTH=old_depth
        try:delattr(BASE,'_ORIGINAL_STEP_SCALED_Q_FOR_P4')
        except AttributeError:pass

def build_base(*args,**kwargs):
    with installed():return BASE.build(*args,**kwargs)
