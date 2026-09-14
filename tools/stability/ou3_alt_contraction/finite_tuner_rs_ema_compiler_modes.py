"""Dual compiler-mode binary32 graph for shipping ``tune_.RS_applied`` EMA.

The C++ source expression

    tune_.RS_applied += alpha_RS * (RS_t - tune_.RS_applied);

has a required rounded subtraction for the parenthesized delta, followed by a
multiply/add pair that may be evaluated as two operations or contracted to an
FMA depending on the deployment toolchain.  No repository fact currently locks
that contraction mode.  This module therefore carries BOTH legal histories from
the same predecessor, target and source-owned alpha.

It does not claim which history shipping selects.  A later deployment fact may
choose one; until then theorem-facing ledgers must retain both.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_binary32 as A

SOURCE=Path(__file__).resolve().parents[3]/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
QUALIFICATION='OU3_ALT_RS_EMA_COMPILER_MODES_V1'


def _q(x,name):
    q=F(x)
    if not B.is_binary32(q): raise ValueError(f'{name} must be actual binary32')
    return q


@dataclass(frozen=True)
class Step:
    previous:F
    target:F
    alpha_source:A.Step
    delta:F
    next_separate:F
    next_fma:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        for n in ('previous','target','delta','next_separate','next_fma'):
            object.__setattr__(self,n,F(getattr(self,n)))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong RS compiler-mode qualification')
        if not isinstance(self.alpha_source,A.Step): raise TypeError('source-owned alpha step required')
        if not all(B.is_binary32(getattr(self,n)) for n in ('previous','target','delta','next_separate','next_fma')):
            raise ValueError('RS compiler-mode graph stores binary32 values only')
        if self.delta!=B.sub(self.target,self.previous):
            raise ValueError('RS delta detached from required rounded subtraction')
        a=self.alpha_source.alpha
        if self.next_separate!=B.add(self.previous,B.mul(a,self.delta)):
            raise ValueError('separate RS EMA result detached')
        if self.next_fma!=B.fma(a,self.delta,self.previous):
            raise ValueError('FMA RS EMA result detached')


def step(previous,target,alpha_source:A.Step):
    if not isinstance(alpha_source,A.Step): raise TypeError('source-owned alpha step required')
    p=_q(previous,'previous RS'); t=_q(target,'target RS')
    d=B.sub(t,p); a=alpha_source.alpha
    sep=B.add(p,B.mul(a,d)); fused=B.fma(a,d,p)
    return Step(p,t,alpha_source,d,sep,fused)


def committed(step_result:Step,*,contracted:bool):
    if not isinstance(step_result,Step): raise TypeError('RS compiler-mode step required')
    if not isinstance(contracted,bool): raise TypeError('literal compiler contraction mode required')
    return step_result.next_fma if contracted else step_result.next_separate


def _source_shape_matches():
    s=SOURCE.read_text()
    return ('tune_.RS_applied    += alpha_RS * (RS_t    - tune_.RS_applied);' in s)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'shipping_RS_EMA_source_shape_matches':_source_shape_matches(),
      'required_target_minus_previous_subtraction_materialized':True,
      'separate_multiply_add_history_materialized':True,
      'contracted_FMA_history_materialized':True,
      'both_histories_share_same_predecessor_target_and_source_owned_alpha':True,
      'shipping_compiler_FP_contraction_mode_qualified':False,
      'source_uniform_RS_roundoff_supply_bound_closed':False,
      'complete_word_RS_compiler_history_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
