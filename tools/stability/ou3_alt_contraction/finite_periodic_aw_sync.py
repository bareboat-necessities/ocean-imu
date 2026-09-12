"""Finite deployed periodic a_w covariance synchronization recurrence.

The shipping default is ``periodic_aw_cov_sync_=true``,
``congruent_aw_cov_sync_=false`` and the MEKF default non-legacy synchronization.
On that path ``periodic_aw_cov_sync_tick_`` does not rewrite posterior covariance
at the end of the current sample.  When

    time - last_aw_cov_sync_sec > adapt_every_secs

it calls ``synchronize_aw_covariance_to_stationary()``, which stores the current
stationary Sigma_aw as a target and sets ``aw_covariance_floor_pending_=true``.
The positive-part floor is consumed inside the NEXT prediction.

Alternative legacy block replacement and congruent immediate synchronization are
real shipping configuration branches but are deliberately not substituted by the
deployed queued-floor relation here; they fail closed until separately composed.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P


def R(x): return P.rational(x)

@dataclass(frozen=True)
class State:
    pending: bool=False
    last_sync_time: F=F(0)
    def __post_init__(self):
        if not isinstance(self.pending,bool): raise TypeError('literal aw-floor pending bit required')
        t=R(self.last_sync_time)
        if t<0: raise ValueError('nonnegative aw-sync clock required')
        object.__setattr__(self,'last_sync_time',t)

@dataclass(frozen=True)
class Result:
    state: State
    requested_now: bool


def tick(state:State,*,time,adapt_every,live,enabled=True,congruent=False,legacy=False):
    if not isinstance(state,State): raise TypeError('periodic aw-sync state required')
    if not all(isinstance(x,bool) for x in (live,enabled,congruent,legacy)):
        raise TypeError('literal aw-sync policy branches required')
    time,adapt_every=R(time),R(adapt_every)
    if time<state.last_sync_time or adapt_every<0: raise ValueError('monotone time and nonnegative cadence required')
    if congruent or legacy:
        raise NotImplementedError('non-default immediate aw synchronization branch is not represented by queued-floor lemma')
    if not enabled or not live or time-state.last_sync_time<=adapt_every:
        return Result(state,False)
    return Result(State(True,time),True)


def consume_at_prediction(state:State):
    """Prediction consumes/clears the queued request regardless of eigensolver outcome."""
    if not isinstance(state,State): raise TypeError('periodic aw-sync state required')
    return State(False,state.last_sync_time) if state.pending else state


def readiness():
    return {
      'deployed_periodic_queue_predicate_materialized':True,
      'queued_floor_request_persists_to_next_prediction':True,
      'prediction_consumes_pending_request':True,
      'queued_floor_target_is_current_active_Sigma_aw':True,
      'legacy_immediate_replacement_branch_attached':False,
      'congruent_immediate_sync_branch_attached':False,
      'clock_binary64_roundoff_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
