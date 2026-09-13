"""One-step machine-versus-exact-real R_S EMA residual recurrence.

Inputs are deliberately separated:

* an exact-real predecessor interval for smoothed R_S;
* the actual binary32 predecessor value;
* a SAME-candidate SpectralMSE machine/real target join;
* exact-real alpha and an independently supplied binary32 alpha_RS.

The exact interval is mapped by the affine EMA with exact alpha.  The machine
state is advanced by the literal binary32 subtract/multiply/add graph.  Their
difference is retained as an exact rational interval.  No equality between the
two alpha values is assumed; ``alpha_supply`` is exposed explicitly.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate_interval as I
from tools.stability.ou3_alt_contraction import finite_tuner_rs_ema_binary32 as M
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_machine_real_join as T

QUALIFICATION='OU3_ALT_RS_EMA_MACHINE_REAL_JOIN_V1'


@dataclass(frozen=True)
class Join:
    previous_exact:I.IntervalTuneState
    previous_machine:F
    target:T.Join
    exact_alpha:F
    machine:M.Step
    exact_next_lo:F
    exact_next_hi:F
    previous_residual_lo:F
    previous_residual_hi:F
    alpha_supply:F
    next_residual_lo:F
    next_residual_hi:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        for n in ('previous_machine','exact_alpha','exact_next_lo','exact_next_hi',
                  'previous_residual_lo','previous_residual_hi','alpha_supply',
                  'next_residual_lo','next_residual_hi'):
            object.__setattr__(self,n,F(getattr(self,n)))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong RS EMA machine-real qualification')
        if not isinstance(self.previous_exact,I.IntervalTuneState) or not isinstance(self.target,T.Join) or not isinstance(self.machine,M.Step):
            raise TypeError('exact predecessor, target join and machine EMA step required')
        if self.machine.previous!=self.previous_machine or self.machine.target!=self.target.machine_target_RS:
            raise ValueError('machine EMA detached from predecessor or SAME target join')
        if not 0<=self.exact_alpha<=1: raise ValueError('exact alpha outside [0,1]')
        if self.previous_residual_lo!=self.previous_machine-self.previous_exact.RS_hi or self.previous_residual_hi!=self.previous_machine-self.previous_exact.RS_lo:
            raise ValueError('predecessor residual detached')
        if self.alpha_supply!=self.machine.alpha-self.exact_alpha: raise ValueError('alpha supply detached')
        if self.next_residual_lo!=self.machine.next-self.exact_next_hi or self.next_residual_hi!=self.machine.next-self.exact_next_lo:
            raise ValueError('next residual detached')


def join(previous_exact:I.IntervalTuneState,previous_machine,target:T.Join,*,exact_alpha,machine_alpha):
    if not isinstance(previous_exact,I.IntervalTuneState) or not isinstance(target,T.Join):
        raise TypeError('interval predecessor and machine-real target join required')
    pm=F(previous_machine); ae=F(exact_alpha); ma=F(machine_alpha)
    if not B.is_binary32(pm): raise ValueError('machine predecessor R_S must be actual binary32')
    if not 0<=ae<=1: raise ValueError('exact alpha outside [0,1]')
    # Exact affine interval image; all weights are nonnegative.
    elo=(1-ae)*previous_exact.RS_lo+ae*target.exact.target_RS_lo
    ehi=(1-ae)*previous_exact.RS_hi+ae*target.exact.target_RS_hi
    machine=M.step(pm,target.machine_target_RS,ma)
    return Join(previous_exact,pm,target,ae,machine,elo,ehi,
                pm-previous_exact.RS_hi,pm-previous_exact.RS_lo,
                machine.alpha-ae,machine.next-ehi,machine.next-elo)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'exact_RS_interval_EMA_image_materialized':True,
      'binary32_RS_EMA_step_composed_with_same_target_join':True,
      'machine_minus_exact_predecessor_RS_residual_carried':True,
      'machine_alpha_minus_exact_alpha_supply_carried':True,
      'machine_minus_exact_next_RS_residual_interval_exposed':True,
      'alpha_RS_target_libm_correspondence_closed':False,
      'source_uniform_target_and_alpha_supply_bounds_closed':False,
      'binary32_RS_commit_correspondence_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
