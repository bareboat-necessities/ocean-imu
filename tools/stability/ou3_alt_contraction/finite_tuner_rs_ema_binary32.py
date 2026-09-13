"""Literal binary32 update of shipping tune_.RS_applied.

Shipping executes

    tune_.RS_applied += alpha_RS * (RS_t - tune_.RS_applied);

where all operands are float.  This module materializes the source-order
binary32 subtraction, multiplication, and addition with the shared RNE kernel.
``alpha_RS`` is an explicit same-event binary32 witness; its production from
adaptiveSmoothingHorizonSec and exp remains a separate deployment obligation.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
QUALIFICATION='OU3_ALT_RS_EMA_BINARY32_V1'


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must be actual binary32')
    return q


@dataclass(frozen=True)
class Step:
    previous:F
    target:F
    alpha:F
    delta:F
    increment:F
    next:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        for n in ('previous','target','alpha','delta','increment','next'):
            object.__setattr__(self,n,F(getattr(self,n)))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong RS EMA qualification')
        if not 0<=self.alpha<=1: raise ValueError('alpha_RS outside [0,1]')
        if self.delta!=B.sub(self.target,self.previous): raise ValueError('RS EMA delta detached from source-order subtraction')
        if self.increment!=B.mul(self.alpha,self.delta): raise ValueError('RS EMA increment detached from source-order multiply')
        if self.next!=B.add(self.previous,self.increment): raise ValueError('RS EMA next detached from source-order addition')


def step(previous,target,alpha):
    p=_q(previous,'previous RS'); t=_q(target,'target RS'); a=_q(alpha,'alpha_RS')
    if not 0<=a<=1: raise ValueError('alpha_RS outside [0,1]')
    d=B.sub(t,p); inc=B.mul(a,d); nxt=B.add(p,inc)
    return Step(p,t,a,d,inc,nxt)


def _source_shape_matches():
    s=SOURCE.read_text()
    return ('const float alpha_RS = 1.0f - std::exp(-dt / RS_sec);' in s and
            'tune_.RS_applied    += alpha_RS * (RS_t    - tune_.RS_applied);' in s)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_RS_EMA_source_shape_matches':_source_shape_matches(),
      'RS_target_minus_previous_binary32_subtraction_materialized':True,
      'alpha_times_delta_binary32_multiply_materialized':True,
      'previous_plus_increment_binary32_add_materialized':True,
      'alpha_RS_target_libm_production_closed':False,
      'source_uniform_alpha_RS_supply_bound_closed':False,
      'machine_RS_EMA_to_exact_interval_join_closed':False,
      'RS_commit_correspondence_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
