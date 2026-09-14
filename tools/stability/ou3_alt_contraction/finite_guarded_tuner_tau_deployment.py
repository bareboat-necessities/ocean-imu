"""Guard-persistent startup/frontend prefix with binary32 tau deployment ledger.

This wraps ``finite_guarded_tuner_prefix`` rather than reimplementing startup
logic. The lower prefix owns the literal guard -> private Mahony -> band/stats
-> candidate -> WPE order. This layer carries the machine tau recurrence from
construction.

In addition to event-by-event operand identity, State enforces the induction
relation between the exact-real frontend TuneState and each deployed compiler
track. Both start at the same shipping 1.1f seed. If the current ledger has seen
n adapting updates, each machine tau must lie within n*2^-20 s of the exact
frontend tau. This follows from the local deployment residual bound and the
nonexpansive exact EMA error recurrence (multiplier 1-alpha in [0,1]). A random
machine tau can therefore no longer be paired with an unrelated frontend state
at goLive.

WPE->StoredFrequency and target-libm exp provenance remain explicit upstream
obligations. This module closes composition/invariant shape, not those machine
facts or startup reachability.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_guarded_tuner_prefix as FRONT
from tools.stability.ou3_alt_contraction import finite_tuner_candidate as CAND
from tools.stability.ou3_alt_contraction import finite_tuner_frequency_binary32 as FREQ
from tools.stability.ou3_alt_contraction import finite_tuner_tau_deployment_ledger as LEDGER
from tools.stability.ou3_alt_contraction import finite_tuner_tau_roundoff_bound as ROUND
from tools.stability.ou3_alt_contraction import finite_shipping_tau_target_binary32 as TARGET


@dataclass(frozen=True)
class State:
    frontend:FRONT.State
    tau:LEDGER.State
    def __post_init__(self):
        if not isinstance(self.frontend,FRONT.State) or not isinstance(self.tau,LEDGER.State):
            raise TypeError('guarded frontend and tau deployment ledger required')
        exact=F(self.frontend.tuner.tune.tau_applied)
        if exact<=0: raise ValueError('exact frontend tau must stay positive')
        allowance=self.tau.updates*ROUND.UNIFORM_RESIDUAL_MAX
        if abs(self.tau.separate-exact)>allowance:
            raise ValueError('separate machine tau detached from exact frontend recurrence envelope')
        if abs(self.tau.fma-exact)>allowance:
            raise ValueError('FMA machine tau detached from exact frontend recurrence envelope')


@dataclass(frozen=True)
class Result:
    state:State
    frontend:FRONT.Result
    tau_step:LEDGER.StepResult|None


def initial(frontend:FRONT.State):
    # State's zero-update invariant requires the exact frontend seed to be the
    # same compiled 1.1f value as the deployed ledger.
    return State(frontend,LEDGER.initial())


def step(state:State,raw,*,dt,stored_frequency:FREQ.StoredFrequency|None=None,
         tau_exp_decay=None,**kwargs):
    if not isinstance(state,State): raise TypeError('frontend tau deployment State required')
    cfg=kwargs.get('candidate_cfg')
    out=FRONT.step(state.frontend,raw,dt=dt,**kwargs)
    cand=out.tuner.candidate
    if cand is None:
        if stored_frequency is not None or tau_exp_decay is not None:
            raise ValueError('Cold/noncandidate branch consumes no tau deployment witnesses')
        return Result(State(out.state,state.tau),out,None)
    if not isinstance(cfg,CAND.CandidateConfig):
        raise TypeError('post-Cold candidate requires persistent CandidateConfig')
    LEDGER.qualify_config(cfg)
    if not isinstance(stored_frequency,FREQ.StoredFrequency) or tau_exp_decay is None:
        raise TypeError('post-Cold candidate requires stored frequency and binary32 tau exp witness')
    ema=kwargs.get('ema')
    if not isinstance(ema,CAND.EmaWitness):
        raise TypeError('post-Cold candidate requires exact tuner EmaWitness')
    if F(ema.decay_tau_sigma)!=F(tau_exp_decay):
        raise ValueError('binary32 tau exp witness detached from exact candidate EMA decay')
    target=TARGET.evaluate(stored_frequency.stored_hz)
    if F(cand.frequency)!=F(stored_frequency.stored_hz):
        raise ValueError('startup candidate frequency detached from stored binary32 tuner frequency')
    if F(cand.tau_target)!=F(target.exact_target):
        raise ValueError('startup candidate tau target detached from source-locked machine-operand target')
    tstep=LEDGER.step(state.tau,stored_frequency,cfg,dt=dt,exp_decay=tau_exp_decay)
    # Construction of the successor State checks the accumulated exact-vs-
    # machine envelope using the incremented ledger count.
    return Result(State(out.state,tstep.state),out,tstep)


def readiness():
    return {
      'guard_Mahony_band_stats_candidate_WPE_order_reused_without_reimplementation':True,
      'shipping_tau_ledger_begins_at_same_1p1f_exact_frontend_seed':True,
      'Cold_noncandidate_branch_preserves_tau_ledger':True,
      'postCold_candidate_advances_both_global_compiler_tau_tracks':True,
      'startup_candidate_frequency_target_decay_share_tau_ledger_operands':True,
      'startup_each_tau_update_retains_both_roundoff_certificates':True,
      'exact_frontend_vs_each_machine_tau_track_error_envelope_inductively_checked':True,
      'arbitrary_machine_tau_frontend_splice_forbidden':True,
      'upstream_WPE_to_StoredFrequency_binary32_correspondence_closed':False,
      'tuner_exp_libm_binary32_correspondence_closed':False,
      'admitted_startup_source_history_carries_this_frontend_product':False,
      'goLive_carries_resulting_tau_ledger_into_Live_product':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,
      'ALT_LIVE_PASS':False,
    }
