"""Binary32 coefficient producer for shipping AdaptiveWaveBandPass.

This materializes the scalar graph before the contraction-sensitive state step:
frequency/corner clamps, rounded ``2*pi*corner*dt`` arguments, the two distinct
``std::exp`` results, and the literal source construction
``alpha=1-exp; q=1-alpha``.  Exp witnesses are tied to tight exact-real
``exp(-x)`` enclosures through binary32 RNE cells; platform-libm correspondence
is deliberately still open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_band_variance_runtime as R
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_binary32 as E

SOURCE=Path(__file__).resolve().parents[3]/'src/tuner/AdaptiveWaveBandPass.h'
ZERO=B.rn32(0); ONE=B.rn32(1); TWO=B.rn32(2)
PI_F=B.rn32(F(314159265358979323846,10**20)); TWO_PI=B.mul(TWO,PI_F)
NYQUIST_FACTOR=B.rn32(F(45,100)); SPACING=B.rn32(F(105,100))
QUALIFICATION='OU3_ALT_BAND_COEFFICIENTS_BINARY32_V1'


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must be actual binary32')
    return q


@dataclass(frozen=True)
class Coefficients:
    active:bool
    f_ref:F
    dt:F
    upper:F
    low:F|None
    high:F|None
    x_low:F|None
    x_high:F|None
    exp_low:F|None
    exp_high:F|None
    alpha_low:F|None
    alpha_high:F|None
    q_low:F|None
    q_high:F|None
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.active,bool): raise TypeError('literal band active branch required')
        for n in ('f_ref','dt','upper'):
            object.__setattr__(self,n,_q(getattr(self,n),n))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong band coefficient qualification')
        optionals=('low','high','x_low','x_high','exp_low','exp_high','alpha_low','alpha_high','q_low','q_high')
        if self.active:
            for n in optionals:
                if getattr(self,n) is None: raise ValueError('active band coefficient branch missing '+n)
                object.__setattr__(self,n,_q(getattr(self,n),n))
            if not self.high>self.low: raise ValueError('active band corners not ordered')
            if not 0<self.q_low<=1 or not 0<self.q_high<=1: raise ValueError('invalid band q coefficients')
        elif any(getattr(self,n) is not None for n in optionals):
            raise ValueError('inactive band coefficient branch consumes no exp/corner witnesses')


def produce(cfg:R.BandConfig,*,f_ref,dt,exp_low=None,exp_high=None):
    if not isinstance(cfg,R.BandConfig): raise TypeError('BandConfig required')
    f=_q(f_ref,'band reference frequency'); h=_q(dt,'band dt')
    if f<=0 or h<=0: raise ValueError('positive band reference frequency/dt required')
    low_ratio=B.rn32(cfg.low_ratio); high_ratio=B.rn32(cfg.high_ratio)
    min_hz=B.rn32(cfg.min_hz); max_hz=B.rn32(cfg.max_hz)
    guard=B.div(NYQUIST_FACTOR,h); upper=min(max_hz,guard)
    if not upper>min_hz:
        if exp_low is not None or exp_high is not None: raise ValueError('inactive upper-limit branch consumes no exp witnesses')
        return Coefficients(False,f,h,upper,None,None,None,None,None,None,None,None,None,None)
    low=max(min_hz,B.mul(low_ratio,f)); low=min(low,B.div(upper,SPACING))
    high=min(upper,B.mul(high_ratio,f)); high=max(high,B.mul(low,SPACING)); high=min(high,upper)
    if not high>low:
        if exp_low is not None or exp_high is not None: raise ValueError('inactive corner-order branch consumes no exp witnesses')
        return Coefficients(False,f,h,upper,None,None,None,None,None,None,None,None,None,None)
    if exp_low is None or exp_high is None: raise ValueError('active band coefficient branch requires two distinct exp witnesses')
    xl=B.mul(B.mul(TWO_PI,low),h); xh=B.mul(B.mul(TWO_PI,high),h)
    el=_q(exp_low,'low-corner exp result'); eh=_q(exp_high,'high-corner exp result')
    for x,e,name in ((xl,el,'low'),(xh,eh,'high')):
        lo,hi=E.exp_minus_enclosure(x)
        if not E._interval_hits_rne_cell(lo,hi,e):
            raise ValueError(f'{name}-corner exp witness detached from SAME rounded argument RNE cell')
    al=B.sub(ONE,el); ah=B.sub(ONE,eh)
    ql=B.sub(ONE,al); qh=B.sub(ONE,ah)
    return Coefficients(True,f,h,upper,low,high,xl,xh,el,eh,al,ah,ql,qh)


def _source_shape_matches():
    s=SOURCE.read_text()
    needles=('const float nyquist_guard_hz = 0.45f / dt;',
      'float low_hz = std::max(min_hz_, low_ratio_ * f_ref_hz);',
      'low_hz = std::min(low_hz, upper_limit / 1.05f);',
      'float high_hz = std::min(upper_limit, high_ratio_ * f_ref_hz);',
      'high_hz = std::max(high_hz, low_hz * 1.05f);',
      'const float alpha_low = 1.0f - std::exp(-(2.0f * 3.14159265358979323846f) * low_hz * dt);',
      'const float alpha_high = 1.0f - std::exp(-(2.0f * 3.14159265358979323846f) * high_hz * dt);',
      'const float q_low = 1.0f - alpha_low;','const float q_high = 1.0f - alpha_high;')
    return all(x in s for x in needles)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_band_corner_and_decay_source_shape_matches':_source_shape_matches(),
      'corner_clamp_and_nyquist_guard_binary32_graph_materialized':True,
      'two_pi_corner_dt_arguments_binary32_source_order_materialized':True,
      'low_and_high_exp_calls_retained_distinct':True,
      'exp_results_bound_to_same_arguments_by_tight_real_enclosure_and_RNE_cell':True,
      'source_two_subtraction_alpha_q_rounding_materialized':True,
      'target_exp_libm_correspondence_closed':False,
      'reference_frequency_machine_production_closed':False,
      'persistent_band_machine_history_composed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
