#!/usr/bin/env python3
"""Validated positive sqrt/log enclosures for OU-III proof work.

The P4 joint estimator needs the shipping WavePeriodEstimator moment-ratio map,
which contains sqrt and log.  This module keeps those operations out of ordinary
libm proof arithmetic.

For sqrt, a binary64 seed is accepted only after exact Fraction comparisons and
is moved with nextafter until it brackets the exact rational square root.
For log, positive binary64 endpoints are represented exactly as Fractions,
range-reduced by a power of two into [1,2), and evaluated with

    log(x) = 2 * sum_{j>=0} y^(2j+1)/(2j+1),  y=(x-1)/(x+1).

On [1,2), |y|<=1/3.  The omitted positive tail is bounded geometrically.  ln(2)
is enclosed by the same exact-rational series.  All final float endpoints are
directed by exact Fraction comparison where applicable and nextafter slack.
"""
from __future__ import annotations

from fractions import Fraction
import math

from ou3_interval import Interval, down, up


def _down_fraction(q: Fraction) -> float:
    f=float(q)
    if not math.isfinite(f): raise OverflowError("fraction does not fit binary64")
    if Fraction.from_float(f)>q: f=math.nextafter(f,-math.inf)
    return f


def _up_fraction(q: Fraction) -> float:
    f=float(q)
    if not math.isfinite(f): raise OverflowError("fraction does not fit binary64")
    if Fraction.from_float(f)<q: f=math.nextafter(f,math.inf)
    return f


def _sqrt_fraction_bounds(q: Fraction) -> tuple[float,float]:
    if q < 0: raise ValueError("sqrt requires nonnegative input")
    if q == 0: return 0.0,0.0
    seed=math.sqrt(float(q))
    lo=seed
    while Fraction.from_float(lo)*Fraction.from_float(lo) > q:
        lo=math.nextafter(lo,-math.inf)
    while True:
        nxt=math.nextafter(lo,math.inf)
        if Fraction.from_float(nxt)*Fraction.from_float(nxt) <= q: lo=nxt
        else: break
    hi=lo
    if Fraction.from_float(hi)*Fraction.from_float(hi) < q:
        hi=math.nextafter(hi,math.inf)
    return lo,hi


def sqrt_point(x: float) -> Interval:
    x=float(x)
    if not math.isfinite(x) or x<0: raise ValueError("sqrt point requires finite x>=0")
    lo,hi=_sqrt_fraction_bounds(Fraction.from_float(x))
    return Interval(lo,hi)


def sqrt_interval(x: Interval) -> Interval:
    if x.lo<0 or not math.isfinite(x.hi): raise ValueError("sqrt interval must be finite/nonnegative")
    return Interval(sqrt_point(x.lo).lo,sqrt_point(x.hi).hi)


def _log_unit_fraction(q: Fraction, terms: int=48) -> tuple[Fraction,Fraction]:
    if not (Fraction(1,1) <= q < Fraction(2,1)):
        raise ValueError("unit log argument must lie in [1,2)")
    y=(q-1)/(q+1)
    if y==0: return Fraction(0),Fraction(0)
    y2=y*y
    term=y
    s=Fraction(0)
    for j in range(terms):
        s += term/Fraction(2*j+1,1)
        term *= y2
    # Every omitted term has the sign of y.  Here y>=0 after reduction.
    # 1/(2j+1)<=1 and sum y^(2j+1) is geometric.
    tail = term/(Fraction(1,1)-y2)
    return 2*s, 2*(s+tail)


def _ln2_bounds(terms:int=64)->tuple[Fraction,Fraction]:
    # ln2 = 2 atanh(1/3)
    y=Fraction(1,3); y2=y*y; term=y; s=Fraction(0)
    for j in range(terms):
        s += term/Fraction(2*j+1,1)
        term *= y2
    tail=term/(Fraction(1,1)-y2)
    return 2*s,2*(s+tail)


def _log_fraction_bounds(q: Fraction) -> tuple[Fraction,Fraction]:
    if q<=0: raise ValueError("log requires positive input")
    # Exact binary range reduction q = m*2^k, m in [1,2).
    k=q.numerator.bit_length()-q.denominator.bit_length()
    if k>=0: m=q/Fraction(1<<k,1)
    else: m=q*Fraction(1<<(-k),1)
    while m<Fraction(1,1):
        m*=2; k-=1
    while m>=Fraction(2,1):
        m/=2; k+=1
    ml,mh=_log_unit_fraction(m)
    l2l,l2h=_ln2_bounds()
    if k>=0:
        return ml+k*l2l, mh+k*l2h
    return ml+k*l2h, mh+k*l2l


def log_point(x: float) -> Interval:
    x=float(x)
    if not math.isfinite(x) or not x>0: raise ValueError("log point requires finite x>0")
    lo,hi=_log_fraction_bounds(Fraction.from_float(x))
    return Interval(_down_fraction(lo),_up_fraction(hi))


def log_interval(x: Interval) -> Interval:
    if x.lo<=0 or not math.isfinite(x.hi): raise ValueError("log interval must be finite/positive")
    return Interval(log_point(x.lo).lo,log_point(x.hi).hi)


def self_test() -> list[str]:
    failures=[]
    for x in (0.0,1e-12,0.5,1.0,2.0,10.0,180.0):
        r=sqrt_point(x)
        if not r.contains(math.sqrt(x)): failures.append(f"sqrt containment {x}")
    for x in (0.05,0.5,1.0,2.0,6.0,180.0):
        r=log_point(x)
        if not r.contains(math.log(x)): failures.append(f"log containment {x}")
    return failures


if __name__=="__main__":
    f=self_test(); print({"failures":f}); raise SystemExit(bool(f))
