"""Coherent dual-compiler machine TuneState product for ALT.

The shipping tuner owns one struct ``TuneState{tau_applied,sigma_applied,
RS_applied}``.  Separate scalar ledgers are insufficient unless they are proven
to describe the SAME sample history.  This product composes the existing tau,
sigma and R_S deployment ledgers and enforces, for each global compiler mode:

* sigma uses the qualified common-alpha object whose TauStep is exactly the tau
  ledger step for that mode;
* the SpectralMSE R_S candidate uses the same mode's tau_target and sigma_target;
* sigma and R_S use the same deployment configuration;
* all three persistent ledgers advance one update together;
* one common pending bit is attached only after the whole TuneState successor is
  formed.

Thus a caller cannot splice tau from one machine history, sigma from another,
and R_S from a third while still entering the theorem product.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as TAU
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_deployment_ledger as SIG
from tools.stability.ou3_alt_contraction import finite_tuner_rs_deployment_ledger as RS

QUALIFICATION='OU3_ALT_MACHINE_TUNESTATE_PRODUCT_V1'


@dataclass(frozen=True)
class State:
    tau:TAU.State
    sigma:SIG.State
    rs:RS.State
    pending:bool=False
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.tau,TAU.State) or not isinstance(self.sigma,SIG.State) or not isinstance(self.rs,RS.State):
            raise TypeError('tau, sigma and R_S deployment ledger states required')
        if not isinstance(self.pending,bool): raise TypeError('literal common pending bit required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine TuneState qualification')
        if not (self.tau.updates==self.sigma.updates==self.rs.updates):
            raise ValueError('machine TuneState scalar ledgers have different update counts')


@dataclass(frozen=True)
class Result:
    state:State
    tau:TAU.StepResult
    sigma:SIG.Result
    rs:RS.Result
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.tau,TAU.StepResult) or not isinstance(self.sigma,SIG.Result) or not isinstance(self.rs,RS.Result):
            raise TypeError('machine TuneState result components required')
        if self.state.tau!=self.tau.state or self.state.sigma!=self.sigma.state or self.state.rs!=self.rs.state:
            raise ValueError('machine TuneState successor detached from scalar ledgers')
        _check_mode(self.tau.separate_step,self.sigma.separate,self.rs.separate,'separate')
        _check_mode(self.tau.fma_step,self.sigma.fma,self.rs.fma,'fma')


def _check_mode(tau_step,sigma_track,rs_track,mode):
    if sigma_track.mode!=mode or rs_track.mode!=mode: raise ValueError('machine TuneState compiler modes crossed')
    if sigma_track.alpha.step!=tau_step:
        raise ValueError('sigma common-alpha source detached from same-mode tau step')
    if rs_track.target.target.tau_target!=tau_step.tau_target:
        raise ValueError('R_S SpectralMSE tau target detached from same-mode tau target')
    if rs_track.target.target.sigma_target!=sigma_track.target.sigma_target:
        raise ValueError('R_S SpectralMSE sigma target detached from same-mode sigma target')
    if rs_track.target.cfg!=sigma_track.target.cfg:
        raise ValueError('sigma and R_S detached from same deployment config')
    if rs_track.alpha.tau_target!=tau_step.tau_target:
        raise ValueError('R_S alpha horizon detached from same-mode tau target')


def initial():
    return State(TAU.initial(),SIG.initial(),RS.initial(),False)


def compose_after_sample(previous:State,*,tau:TAU.StepResult,sigma:SIG.Result,rs:RS.Result,pending_after):
    """Join already-executed same-sample scalar recurrences into one TuneState."""
    if not isinstance(previous,State): raise TypeError('machine TuneState predecessor required')
    if not isinstance(tau,TAU.StepResult) or not isinstance(sigma,SIG.Result) or not isinstance(rs,RS.Result):
        raise TypeError('tau, sigma and R_S sample results required')
    if not isinstance(pending_after,bool): raise TypeError('literal pending-after bit required')
    if tau.separate_step.previous!=previous.tau.separate or tau.fma_step.previous!=previous.tau.fma:
        raise ValueError('tau result detached from machine TuneState predecessor')
    if sigma.separate.previous!=previous.sigma.separate or sigma.fma.previous!=previous.sigma.fma:
        raise ValueError('sigma result detached from machine TuneState predecessor')
    if rs.separate.ema.previous!=previous.rs.separate or rs.fma.ema.previous!=previous.rs.fma:
        raise ValueError('R_S result detached from machine TuneState predecessor')
    expected=previous.tau.updates+1
    if not (tau.state.updates==sigma.state.updates==rs.state.updates==expected):
        raise ValueError('machine TuneState ledgers did not advance together')
    _check_mode(tau.separate_step,sigma.separate,rs.separate,'separate')
    _check_mode(tau.fma_step,sigma.fma,rs.fma,'fma')
    return Result(State(tau.state,sigma.state,rs.state,pending_after),tau,sigma,rs)


def hold(previous:State):
    """Cold/no-adapt sample: the complete machine TuneState is identity."""
    if not isinstance(previous,State): raise TypeError('machine TuneState State required')
    return previous


def readiness():
    return {
      'shipping_tau_sigma_RS_persistent_machine_states_composed':True,
      'all_three_scalar_ledgers_advance_with_one_common_update_index':True,
      'sigma_common_alpha_is_same_mode_tau_step':True,
      'SpectralMSE_RS_uses_same_mode_tau_and_sigma_targets':True,
      'sigma_and_RS_share_same_deployment_config_per_mode':True,
      'one_common_pending_bit_attached_after_whole_TuneState_successor':True,
      'cold_or_nonadapting_complete_TuneState_identity_materialized':True,
      'next_boundary_common_machine_commit_attached':False,
      'startup_frontend_machine_TuneState_product_attached':False,
      'Live_600_step_machine_TuneState_product_attached':False,
      'all_target_libm_correspondence_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
