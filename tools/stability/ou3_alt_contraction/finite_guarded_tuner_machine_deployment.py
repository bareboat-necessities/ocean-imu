"""Guarded startup frontend with coherent WPE + whole machine TuneState.

This layer strengthens ``finite_guarded_tuner_wpe_tau_deployment`` without
reimplementing its shipping-order WPE/tau graph.  It carries the sigma and R_S
machine ledgers alongside that already-rooted tau history and joins all three
only when the exact frontend actually produces a tuner candidate.

Machine sigma/R_S targets need not equal the exact-real candidate targets.
Their differences are retained explicitly as deployment supplies.  What is
forbidden is ancestry splicing: tau, sigma and R_S must advance from the same
whole predecessor, with one common update index and the exact candidate's
pending bit.

Cold/noncandidate samples preserve sigma/R_S and pending while the lower WPE
history may advance.  Startup reachability and upstream binary32/libm
correspondence remain open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_guarded_tuner_wpe_tau_deployment as BASE
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as MACHINE
from tools.stability.ou3_alt_contraction import finite_tuner_sigma_deployment_ledger as SIGMA
from tools.stability.ou3_alt_contraction import finite_tuner_rs_deployment_ledger as RS


@dataclass(frozen=True)
class ModeSupply:
    sigma_target:F
    RS_target:F
    def __post_init__(self):
        object.__setattr__(self,'sigma_target',F(self.sigma_target))
        object.__setattr__(self,'RS_target',F(self.RS_target))


@dataclass(frozen=True)
class State:
    base:BASE.State
    machine:MACHINE.State
    def __post_init__(self):
        if not isinstance(self.base,BASE.State) or not isinstance(self.machine,MACHINE.State):
            raise TypeError('coherent startup WPE/tau base and machine TuneState required')
        if self.machine.tau != self.base.tau:
            raise ValueError('whole machine TuneState tau detached from startup WPE/tau history')
        if not (self.machine.sigma.updates==self.machine.rs.updates==self.machine.tau.updates):
            raise ValueError('startup machine scalar update counts differ')


@dataclass(frozen=True)
class Result:
    state:State
    base:BASE.Result
    machine_result:MACHINE.Result|None
    separate_supply:ModeSupply|None
    fma_supply:ModeSupply|None
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.base,BASE.Result):
            raise TypeError('startup machine result requires strengthened state and lower result')
        if self.base.tau_step is None:
            if self.machine_result is not None or self.separate_supply is not None or self.fma_supply is not None:
                raise ValueError('noncandidate startup sample cannot carry machine candidate result')
        else:
            if not isinstance(self.machine_result,MACHINE.Result):
                raise TypeError('candidate startup sample requires coherent machine TuneState result')
            if self.machine_result.tau != self.base.tau_step:
                raise ValueError('whole machine result detached from lower same-sample tau result')
            if self.state.machine != self.machine_result.state:
                raise ValueError('startup machine successor detached from coherent result')


def initial(frontend):
    base=BASE.initial(frontend)
    machine=MACHINE.initial()
    if machine.tau != base.tau:
        raise AssertionError('shipping machine tau seed disagrees with guarded startup tau seed')
    return State(base,machine)


def step(state:State,raw,*,sigma_result:SIGMA.Result|None=None,rs_result:RS.Result|None=None,**kwargs):
    if not isinstance(state,State): raise TypeError('guarded startup machine State required')
    lower=BASE.step(state.base,raw,**kwargs)
    cand=lower.frontend.tuner.candidate
    if lower.tau_step is None:
        if sigma_result is not None or rs_result is not None:
            raise ValueError('Cold/noncandidate startup branch consumes no sigma/R_S machine results')
        machine=MACHINE.State(lower.state.tau,state.machine.sigma,state.machine.rs,state.machine.pending)
        return Result(State(lower.state,machine),lower,None,None,None)
    if cand is None:
        raise AssertionError('tau machine update exists without exact frontend candidate')
    if not isinstance(sigma_result,SIGMA.Result) or not isinstance(rs_result,RS.Result):
        raise TypeError('candidate startup sample requires same-sample sigma and R_S machine results')
    composed=MACHINE.compose_after_sample(state.machine,tau=lower.tau_step,
                                          sigma=sigma_result,rs=rs_result,
                                          pending_after=cand.pending_after)
    ss=ModeSupply(sigma_result.separate.target.sigma_target-F(cand.sigma_target),
                  rs_result.separate.target.target.RS_target-F(cand.RS_target))
    fs=ModeSupply(sigma_result.fma.target.sigma_target-F(cand.sigma_target),
                  rs_result.fma.target.target.RS_target-F(cand.RS_target))
    return Result(State(lower.state,composed.state),lower,composed,ss,fs)


def readiness():
    b=BASE.readiness(); m=MACHINE.readiness()
    return {
      'guarded_startup_WPE_tau_product_consumed':True,
      'construction_roots_whole_tau_sigma_RS_machine_TuneState':True,
      'machine_tau_is_exact_same_history_as_startup_WPE_tau_product':True,
      'Cold_noncandidate_samples_preserve_sigma_RS_and_pending':True,
      'candidate_samples_require_tau_sigma_RS_to_advance_from_one_predecessor':True,
      'exact_candidate_pending_bit_controls_whole_machine_successor':True,
      'machine_minus_exact_sigma_and_RS_target_supplies_retained_per_mode':True,
      'startup_global_compiler_track_coherence_includes_WPE_tau_sigma_RS': bool(
          b['startup_global_compiler_track_coherence_includes_WPE_log_and_tau_states'] and
          m['shipping_tau_sigma_RS_persistent_machine_states_composed']),
      'upstream_sigma_binary32_source_production_closed':False,
      'upstream_RS_binary32_source_production_closed':False,
      'all_target_libm_correspondence_closed':False,
      'every_admitted_startup_history_reaches_TunerReady_with_this_product':False,
      'goLive_carries_this_exact_WPE_machine_TuneState_product':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,
    }
