"""Exact binary64 deployment clock relation for periodic a_w covariance sync.

The OU-III implementation stores its inner ``time_`` and
``last_aw_cov_sync_sec_`` as ``double`` but receives ``dt`` and
``adapt_every_secs_`` as binary32 values. On the canonical ALT 200 Hz source
prefix this module evaluates those operations exactly under IEEE round-to-
nearest-even and proves the default periodic predicate keeps the same branch as
the exact 5 ms graph through the 150 s startup horizon plus one 600-edge word.

This is a finite-prefix arithmetic certificate, not an indefinite clock theorem
and not a replacement for target compiler/library qualification.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B32
from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as WRAPPER

DT_REAL=F(1,200)
DT_FLOAT=B32.rn32(DT_REAL)
ADAPT_REAL=F(1,10)
ADAPT_FLOAT=B32.rn32(ADAPT_REAL)
STARTUP_TIMEOUT_STEPS=WRAPPER.STARTUP_TIMEOUT_STEPS
WORD_STEPS=WRAPPER.WORD_STEPS
MAX_STEPS=WRAPPER.MAX_STEPS
SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
QUALIFICATION='OU3_ALT_AW_SYNC_CLOCK_BINARY64_CANONICAL_5MS_V1'


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
    q,r=divmod(x.numerator,x.denominator); twice=2*r
    if twice<x.denominator:return q
    if twice>x.denominator:return q+1
    return q if q%2==0 else q+1


def rn64(x)->F:
    """Round exact rational to nearest-even normal binary64 or zero."""
    x=F(x)
    if x==0:return F(0)
    sign=-1 if x<0 else 1; a=abs(x)
    e=_floor_log2_positive(a)
    if e < -1022 or e > 1023:
        raise ValueError('outside normal finite binary64 range used by this proof')
    quantum=_pow2(e-52); m=_rne_nonnegative_integer(a/quantum)
    if m==(1<<53):
        m>>=1; e+=1
        if e>1023: raise ValueError('binary64 overflow')
        quantum=_pow2(e-52)
    if not (1<<52)<=m<(1<<53):
        raise AssertionError('binary64 significand normalization failed')
    return sign*F(m)*quantum


def add64(a,b): return rn64(F(a)+F(b))
def sub64(a,b): return rn64(F(a)-F(b))


@lru_cache(maxsize=1)
def _times():
    out=[F(0)]; t=F(0)
    for _ in range(MAX_STEPS):
        t=add64(t,DT_FLOAT); out.append(t)
    return tuple(out)


def canonical_step_index(real_time)->int:
    t=F(real_time); q=t/DT_REAL
    if q.denominator!=1:
        raise ValueError('inner-clock lookup requires canonical 5 ms exact graph time')
    k=q.numerator
    if not 0<=k<=MAX_STEPS:
        raise ValueError('inner-clock lookup outside startup plus one-word prefix')
    return k


def clock_at_step(step:int)->F:
    if not isinstance(step,int) or isinstance(step,bool) or not 0<=step<=MAX_STEPS:
        raise ValueError('clock step outside certified prefix')
    return _times()[step]


def elapsed_steps(now_step:int,last_step:int)->F:
    if not (0<=last_step<=now_step<=MAX_STEPS):
        raise ValueError('invalid clock step interval')
    return sub64(clock_at_step(now_step),clock_at_step(last_step))


def deployed_due(now_step:int,last_step:int)->bool:
    return elapsed_steps(now_step,last_step)>ADAPT_FLOAT


def exact_due(now_step:int,last_step:int)->bool:
    return F(now_step-last_step,200)>ADAPT_REAL


def qualify_exact_tick(time,last_sync_time,adapt_every=ADAPT_REAL):
    """Bind one exact finite tick to the canonical deployed binary64 predicate."""
    if F(adapt_every)!=ADAPT_REAL:
        raise ValueError('current binary64 certificate is for shipping default 0.1 s cadence')
    now=canonical_step_index(time); last=canonical_step_index(last_sync_time)
    if last>now: raise ValueError('last sync cannot follow current time')
    e=exact_due(now,last); d=deployed_due(now,last)
    if e!=d:
        raise RuntimeError('binary64 aw-sync predicate diverges from exact finite branch')
    return {
      'now_step':now,'last_step':last,
      'exact_time':F(time),'exact_last_sync_time':F(last_sync_time),
      'deployed_time':clock_at_step(now),'deployed_last_sync_time':clock_at_step(last),
      'deployed_elapsed':elapsed_steps(now,last),'deployed_cadence':ADAPT_FLOAT,
      'exact_due':e,'deployed_due':d,
    }


def _source_shape_matches():
    s=SOURCE.read_text()
    return all(x in s for x in (
      'constexpr float ADAPT_EVERY_SECS           = 0.1f;',
      'time_ += dt;',
      'if (time_ - last_aw_cov_sync_sec_ <= adapt_every_secs_) return;',
      'last_aw_cov_sync_sec_ = time_;',
      'double last_aw_cov_sync_sec_ = 0.0;',
      'double time_;','time_(0.0)',
      'float adapt_every_secs_       = ADAPT_EVERY_SECS;'))


@dataclass(frozen=True)
class Report:
    max_20_step_elapsed:F
    min_21_step_elapsed:F
    max_20_start:int
    min_21_start:int
    strict_clock_advance:bool
    predicate_partition_closed:bool
    source_shape_matches:bool


def build()->Report:
    ts=_times(); strict=all(ts[k+1]>ts[k] for k in range(MAX_STEPS))
    mx20=F(-1); mi21=None; kmx=kmi=0
    for k in range(MAX_STEPS-21+1):
        e20=sub64(ts[k+20],ts[k])
        if e20>mx20: mx20,kmx=e20,k
        e21=sub64(ts[k+21],ts[k])
        if mi21 is None or e21<mi21: mi21,kmi=e21,k
    closed=bool(mx20<=ADAPT_FLOAT and mi21>ADAPT_FLOAT)
    return Report(mx20,mi21,kmx,kmi,strict,closed,_source_shape_matches())


def readiness():
    r=build()
    return {
      'qualification':QUALIFICATION,
      'canonical_dt_float':DT_FLOAT,
      'default_adapt_every_float':ADAPT_FLOAT,
      'shipping_inner_time_and_aw_sync_source_shape_matches':r.source_shape_matches,
      'binary64_clock_strictly_advances_through_startup_plus_word':r.strict_clock_advance,
      'max_deployed_elapsed_over_20_samples':r.max_20_step_elapsed,
      'min_deployed_elapsed_over_21_samples':r.min_21_step_elapsed,
      'default_due_partition_is_exactly_gap_ge_21_samples':r.predicate_partition_closed,
      'canonical_aw_sync_binary64_predicate_closed': bool(r.source_shape_matches and r.strict_clock_advance and r.predicate_partition_closed),
      'arbitrary_dt_inner_clock_closed':False,
      'mutable_adapt_every_setter_ancestry_closed':False,
      'indefinite_inner_clock_lifetime_closed':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
    }
