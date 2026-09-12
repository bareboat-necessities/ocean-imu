"""Exact IEEE-binary32 wrapper-clock recurrence on the canonical 5 ms grid.

The outer ``SeaStateFusion_OU_III`` wrapper stores ``t_`` as ``float`` and
updates it with ``t_ += dt``.  The finite ALT graph previously treated that
clock as exact real time.  This module closes the arithmetic recurrence for the
canonical 200 Hz source grid through the latest possible startup handoff
(150 s) plus one 600-transition contraction word (3 s).

This is an exact integer/rational IEEE-754 round-to-nearest-even calculation;
it does not call host floating point.  It proves the binary32 values produced
under the stated arithmetic profile, and quantifies their difference from the
ideal grid.  It deliberately does NOT claim indefinite lifetime: at much larger
magnitudes a 5 ms increment eventually loses resolution, and arbitrary dt / a
different rounding mode remain separate deployment obligations.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

DT_REAL=F(1,200)
STARTUP_TIMEOUT_STEPS=150*200
STARTUP_MIN_STEPS=8*200
WORD_STEPS=600
MAX_STEPS=STARTUP_TIMEOUT_STEPS+WORD_STEPS
SOURCE=Path(__file__).resolve().parents[2]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
QUALIFICATION='OU3_ALT_WRAPPER_CLOCK_BINARY32_CANONICAL_5MS_V1'


def _pow2(e:int)->F:
    return F(1<<e) if e>=0 else F(1,1<<(-e))


def _floor_log2(x:F)->int:
    if x<=0: raise ValueError('positive input required')
    n,d=x.numerator,x.denominator
    e=n.bit_length()-d.bit_length()
    while _pow2(e)>x: e-=1
    while _pow2(e+1)<=x: e+=1
    return e


def _rne_integer(x:F)->int:
    """Round nonnegative rational to nearest integer, ties to even."""
    if x<0: raise ValueError('nonnegative input required')
    q,r=divmod(x.numerator,x.denominator)
    twice=2*r
    if twice<x.denominator: return q
    if twice>x.denominator: return q+1
    return q if q%2==0 else q+1


def binary32_positive(x)->F:
    """Exact RN-even binary32 value for the positive normal range used here."""
    x=F(x)
    if x==0:return F(0)
    if x<0: raise ValueError('clock arithmetic is nonnegative')
    e=_floor_log2(x)
    if e < -126 or e > 127:
        raise ValueError('outside normal finite binary32 range used by this proof')
    quantum=_pow2(e-23)
    m=_rne_integer(x/quantum)
    if m==(1<<24):
        m>>=1; e+=1; quantum=_pow2(e-23)
    if not (1<<23)<=m<(1<<24):
        raise AssertionError('binary32 significand normalization failed')
    return F(m)*quantum


DT_FLOAT=binary32_positive(DT_REAL)

@dataclass(frozen=True)
class Report:
    dt_float:F
    startup_timeout_clock:F
    latest_word_end_clock:F
    max_absolute_grid_error:F
    max_absolute_grid_error_step:int
    max_three_second_elapsed_error:F
    max_three_second_elapsed_error_start_step:int
    all_updates_strictly_advance:bool
    shipping_source_shape_matches:bool


def _source_shape_matches()->bool:
    s=SOURCE.read_text()
    return (s.count('float t_ = 0.0f;')==1 and
            s.count('t_ += dt;')>=1 and
            s.count('float proxy_startup_timeout_sec = 150.0f;')==1)


def build()->Report:
    t=F(0); times=[t]
    max_err=F(0); max_k=0; strict=True
    for k in range(1,MAX_STEPS+1):
        nxt=binary32_positive(t+DT_FLOAT)
        if not nxt>t: strict=False
        t=nxt; times.append(t)
        err=abs(t-F(k,200))
        if err>max_err: max_err,max_k=err,k

    max_elapsed=F(0); max_start=STARTUP_MIN_STEPS
    for k in range(STARTUP_MIN_STEPS,STARTUP_TIMEOUT_STEPS+1):
        err=abs((times[k+WORD_STEPS]-times[k])-F(3))
        if err>max_elapsed: max_elapsed,max_start=err,k

    return Report(
        DT_FLOAT,
        times[STARTUP_TIMEOUT_STEPS],
        times[MAX_STEPS],
        max_err,max_k,
        max_elapsed,max_start,
        strict,
        _source_shape_matches())


def readiness():
    r=build()
    return {
      'qualification':QUALIFICATION,
      'canonical_real_dt_s':DT_REAL,
      'canonical_binary32_dt_s':r.dt_float,
      'shipping_outer_clock_is_float_and_incremented_by_dt':r.shipping_source_shape_matches,
      'all_binary32_clock_updates_strictly_advance_through_timeout_plus_one_word':r.all_updates_strictly_advance,
      'startup_timeout_clock_binary32':r.startup_timeout_clock,
      'latest_word_end_clock_binary32':r.latest_word_end_clock,
      'max_absolute_clock_vs_ideal_grid_error_s':r.max_absolute_grid_error,
      'max_absolute_clock_vs_ideal_grid_error_step':r.max_absolute_grid_error_step,
      'max_elapsed_error_over_any_3s_word_starting_between_8s_and_150s_s':r.max_three_second_elapsed_error,
      'max_elapsed_error_word_start_step':r.max_three_second_elapsed_error_start_step,
      'canonical_5ms_wrapper_clock_prefix_binary32_closed': bool(
          r.shipping_source_shape_matches and r.all_updates_strictly_advance),
      'arbitrary_dt_wrapper_clock_closed':False,
      'indefinite_wrapper_clock_lifetime_closed':False,
      'magnetic_counter_lifetime_closed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
