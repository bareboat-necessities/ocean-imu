"""Machine-versus-exact-real R_S EMA residual recurrence.

The exact interval is advanced with an exact-real alpha.  Deployment arithmetic
uses the SAME SpectralMSE target and source-owned binary32 alpha.  Because the
C++ multiply/add may be either separate or contracted, the theorem-facing
``join_compiler_modes`` carries both legal binary32 successors and therefore
both next residual intervals.  No compiler mode or libm correctness is guessed.

Older single-mode helpers remain for component regressions only.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_candidate_interval as I
from tools.stability.ou3_alt_contraction import finite_tuner_rs_alpha_binary32 as A
from tools.stability.ou3_alt_contraction import finite_tuner_rs_ema_binary32 as M
from tools.stability.ou3_alt_contraction import finite_tuner_rs_ema_compiler_modes as CM
from tools.stability.ou3_alt_contraction import finite_tuner_spectral_machine_real_join as T

QUALIFICATION='OU3_ALT_RS_EMA_MACHINE_REAL_JOIN_V3'


def _exact_image(previous_exact,target,ae):
    elo=(1-ae)*previous_exact.RS_lo+ae*target.exact.target_RS_lo
    ehi=(1-ae)*previous_exact.RS_hi+ae*target.exact.target_RS_hi
    return elo,ehi


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
    alpha_source:A.Step|None=None
    qualification:str=QUALIFICATION
    def __post_init__(self):
        for n in ('previous_machine','exact_alpha','exact_next_lo','exact_next_hi',
                  'previous_residual_lo','previous_residual_hi','alpha_supply',
                  'next_residual_lo','next_residual_hi'):
            object.__setattr__(self,n,F(getattr(self,n)))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong RS EMA machine-real qualification')
        if not isinstance(self.previous_exact,I.IntervalTuneState) or not isinstance(self.target,T.Join) or not isinstance(self.machine,M.Step):
            raise TypeError('exact predecessor, target join and machine EMA step required')
        if self.alpha_source is not None and not isinstance(self.alpha_source,A.Step):
            raise TypeError('alpha source must be deployment RS alpha step')
        if self.machine.previous!=self.previous_machine or self.machine.target!=self.target.machine_target_RS:
            raise ValueError('machine EMA detached from predecessor or SAME target join')
        if self.alpha_source is not None:
            if self.alpha_source.alpha!=self.machine.alpha: raise ValueError('machine EMA alpha detached from source-owned alpha step')
            if self.alpha_source.tau_target!=self.target.target.tau_target: raise ValueError('RS alpha horizon detached from SAME candidate tau target')
        if not 0<=self.exact_alpha<=1: raise ValueError('exact alpha outside [0,1]')
        if self.previous_residual_lo!=self.previous_machine-self.previous_exact.RS_hi or self.previous_residual_hi!=self.previous_machine-self.previous_exact.RS_lo:
            raise ValueError('predecessor residual detached')
        if self.alpha_supply!=self.machine.alpha-self.exact_alpha: raise ValueError('alpha supply detached')
        if self.next_residual_lo!=self.machine.next-self.exact_next_hi or self.next_residual_hi!=self.machine.next-self.exact_next_lo:
            raise ValueError('next residual detached')


@dataclass(frozen=True)
class CompilerModeJoin:
    previous_exact:I.IntervalTuneState
    previous_machine:F
    target:T.Join
    exact_alpha:F
    alpha_source:A.Step
    machine:CM.Step
    exact_next_lo:F
    exact_next_hi:F
    previous_residual_lo:F
    previous_residual_hi:F
    alpha_supply:F
    separate_residual_lo:F
    separate_residual_hi:F
    fma_residual_lo:F
    fma_residual_hi:F
    qualification:str=QUALIFICATION
    def __post_init__(self):
        for n in ('previous_machine','exact_alpha','exact_next_lo','exact_next_hi','previous_residual_lo',
                  'previous_residual_hi','alpha_supply','separate_residual_lo','separate_residual_hi',
                  'fma_residual_lo','fma_residual_hi'):
            object.__setattr__(self,n,F(getattr(self,n)))
        if self.qualification!=QUALIFICATION: raise ValueError('wrong RS compiler-mode join qualification')
        if not isinstance(self.previous_exact,I.IntervalTuneState) or not isinstance(self.target,T.Join):
            raise TypeError('exact predecessor and target join required')
        if not isinstance(self.alpha_source,A.Step) or not isinstance(self.machine,CM.Step):
            raise TypeError('source-owned alpha and dual compiler-mode step required')
        if self.alpha_source.tau_target!=self.target.target.tau_target:
            raise ValueError('RS alpha horizon detached from SAME candidate tau target')
        if self.machine.previous!=self.previous_machine or self.machine.target!=self.target.machine_target_RS or self.machine.alpha_source!=self.alpha_source:
            raise ValueError('compiler-mode machine step detached from same predecessor/target/alpha')
        if self.previous_residual_lo!=self.previous_machine-self.previous_exact.RS_hi or self.previous_residual_hi!=self.previous_machine-self.previous_exact.RS_lo:
            raise ValueError('predecessor residual detached')
        if self.alpha_supply!=self.alpha_source.alpha-self.exact_alpha: raise ValueError('alpha supply detached')
        if self.separate_residual_lo!=self.machine.next_separate-self.exact_next_hi or self.separate_residual_hi!=self.machine.next_separate-self.exact_next_lo:
            raise ValueError('separate residual detached')
        if self.fma_residual_lo!=self.machine.next_fma-self.exact_next_hi or self.fma_residual_hi!=self.machine.next_fma-self.exact_next_lo:
            raise ValueError('FMA residual detached')


def _join(previous_exact,previous_machine,target,exact_alpha,machine_alpha,alpha_source):
    if not isinstance(previous_exact,I.IntervalTuneState) or not isinstance(target,T.Join):
        raise TypeError('interval predecessor and machine-real target join required')
    pm=F(previous_machine); ae=F(exact_alpha); ma=F(machine_alpha)
    if not B.is_binary32(pm): raise ValueError('machine predecessor R_S must be actual binary32')
    if not 0<=ae<=1: raise ValueError('exact alpha outside [0,1]')
    elo,ehi=_exact_image(previous_exact,target,ae)
    machine=M.step(pm,target.machine_target_RS,ma)
    return Join(previous_exact,pm,target,ae,machine,elo,ehi,
                pm-previous_exact.RS_hi,pm-previous_exact.RS_lo,
                machine.alpha-ae,machine.next-ehi,machine.next-elo,alpha_source)


def join(previous_exact:I.IntervalTuneState,previous_machine,target:T.Join,*,exact_alpha,machine_alpha):
    return _join(previous_exact,previous_machine,target,exact_alpha,machine_alpha,None)


def join_from_alpha_step(previous_exact:I.IntervalTuneState,previous_machine,target:T.Join,*,exact_alpha,alpha_step:A.Step):
    if not isinstance(alpha_step,A.Step): raise TypeError('source-owned RS alpha step required')
    if alpha_step.tau_target!=target.target.tau_target: raise ValueError('RS alpha step detached from SAME candidate tau target')
    return _join(previous_exact,previous_machine,target,exact_alpha,alpha_step.alpha,alpha_step)


def join_compiler_modes(previous_exact:I.IntervalTuneState,previous_machine,target:T.Join,*,exact_alpha,alpha_step:A.Step):
    """Theorem-facing join retaining both legal compiler evaluation histories."""
    if not isinstance(previous_exact,I.IntervalTuneState) or not isinstance(target,T.Join) or not isinstance(alpha_step,A.Step):
        raise TypeError('interval predecessor, target join and source-owned alpha step required')
    pm=F(previous_machine); ae=F(exact_alpha)
    if not B.is_binary32(pm): raise ValueError('machine predecessor R_S must be actual binary32')
    if not 0<=ae<=1: raise ValueError('exact alpha outside [0,1]')
    if alpha_step.tau_target!=target.target.tau_target: raise ValueError('RS alpha step detached from SAME candidate tau target')
    elo,ehi=_exact_image(previous_exact,target,ae)
    machine=CM.step(pm,target.machine_target_RS,alpha_step)
    return CompilerModeJoin(previous_exact,pm,target,ae,alpha_step,machine,elo,ehi,
        pm-previous_exact.RS_hi,pm-previous_exact.RS_lo,alpha_step.alpha-ae,
        machine.next_separate-ehi,machine.next_separate-elo,
        machine.next_fma-ehi,machine.next_fma-elo)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'exact_RS_interval_EMA_image_materialized':True,
      'binary32_RS_EMA_step_composed_with_same_target_join':True,
      'source_owned_binary32_RS_alpha_can_feed_persistent_recurrence':True,
      'same_candidate_tau_owns_RS_alpha_horizon':True,
      'separate_and_FMA_RS_histories_both_joined_to_same_exact_interval':True,
      'compiler_mode_not_guessed_by_theorem_join':True,
      'machine_minus_exact_predecessor_RS_residual_carried':True,
      'machine_alpha_minus_exact_alpha_supply_carried':True,
      'machine_minus_exact_next_RS_residual_interval_exposed':True,
      'alpha_RS_target_libm_correspondence_closed':False,
      'shipping_compiler_FP_contraction_mode_qualified':False,
      'source_uniform_target_and_alpha_supply_bounds_closed':False,
      'binary32_RS_commit_correspondence_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
