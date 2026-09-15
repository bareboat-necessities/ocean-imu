"""Persistent actual-machine AdaptiveWaveBandPass ledger for ALT.

The compiler's local contraction choices are not selected by the theorem.
Instead, one concrete machine state is carried and every active successor must
belong to ``finite_band_binary32_contraction``'s complete local outcome set for
the same source-owned coefficient object and input.  This prevents exponential
history enumeration while retaining every locally legal contraction choice.

Inactive corner branches are literal identity and consume no successor witness.
The target compiler/libm correspondence that says the deployed execution really
belongs to these per-step relations remains fail-closed.
"""
from __future__ import annotations

from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as HORIZON
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_band_binary32_contraction as C
from tools.stability.ou3_alt_contraction import finite_band_coefficients_binary32 as Q

QUALIFICATION='OU3_ALT_BAND_MACHINE_LEDGER_V1'
MAX_SAMPLES=HORIZON.MAX_STEPS


@dataclass(frozen=True)
class State:
    machine:C.State=C.State()
    samples:int=0
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.machine,C.State): raise TypeError('actual machine band State required')
        if not isinstance(self.samples,int) or not 0<=self.samples<=MAX_SAMPLES: raise ValueError('band machine sample count outside bounded startup+word horizon')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong band machine ledger qualification')


@dataclass(frozen=True)
class StepResult:
    before:State
    state:State
    coefficients:Q.Coefficients
    envelope:C.Envelope|None
    def __post_init__(self):
        if not isinstance(self.before,State) or not isinstance(self.state,State) or not isinstance(self.coefficients,Q.Coefficients):
            raise TypeError('band machine ledger step objects required')
        if self.state.samples!=self.before.samples+1: raise ValueError('band machine ledger sample count did not advance once')
        if self.coefficients.active:
            if not isinstance(self.envelope,C.Envelope): raise TypeError('active band step requires contraction envelope')
            if self.envelope.predecessor!=self.before.machine: raise ValueError('band contraction envelope detached from ledger predecessor')
            if (self.envelope.q_low,self.envelope.q_high)!=(self.coefficients.q_low,self.coefficients.q_high):
                raise ValueError('band contraction envelope detached from same coefficient producer')
            if not self.envelope.accepts(self.state.machine): raise ValueError('carried machine successor outside certified contraction envelope')
        else:
            if self.envelope is not None: raise ValueError('inactive band step consumes no contraction envelope')
            if self.state.machine!=self.before.machine: raise ValueError('inactive adaptive-band branch must preserve machine state exactly')


def initial(): return State()


def step(state:State,coefficients:Q.Coefficients,*,x,successor:C.State|None=None):
    if not isinstance(state,State) or not isinstance(coefficients,Q.Coefficients): raise TypeError('band machine State and coefficient relation required')
    if state.samples>=MAX_SAMPLES: raise ValueError('band machine ledger exceeded bounded startup+word horizon')
    if coefficients.active:
        if successor is None: raise ValueError('active machine band step requires actual successor witness')
        env=C.step(state.machine,x=x,q_low=coefficients.q_low,q_high=coefficients.q_high)
        if not isinstance(successor,C.State) or not env.accepts(successor):
            raise ValueError('actual machine band successor outside same-step contraction set')
        nxt=State(successor,state.samples+1)
        return StepResult(state,nxt,coefficients,env)
    if successor is not None: raise ValueError('inactive machine band step consumes no successor witness')
    nxt=State(state.machine,state.samples+1)
    return StepResult(state,nxt,coefficients,None)


def readiness():
    return {
      'qualification':QUALIFICATION,
      'actual_machine_band_state_persists_without_contraction_history_branch_explosion':True,
      'each_active_successor_must_belong_to_same_step_finite_contraction_set':True,
      'inactive_band_branch_is_literal_machine_identity':True,
      'bounded_startup_plus_600_sample_horizon_enforced':True,
      'band_coefficient_source_relation_consumed':True,
      'target_compiler_execution_membership_in_local_contraction_relation_closed':False,
      'target_exp_libm_correspondence_closed':False,
      'band_noise_floor_machine_consumer_attached':False,
      'startup_frontend_machine_history_attached':False,
      'Live_600_step_machine_history_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
