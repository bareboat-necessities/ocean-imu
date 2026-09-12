"""Finite deployed periodic a_w covariance synchronization recurrence.

The shipping default is ``periodic_aw_cov_sync_=true``,
``congruent_aw_cov_sync_=false`` and non-legacy MEKF synchronization.  On that
path the end-of-sample tick snapshots the THEN-current stationary covariance
into ``aw_covariance_floor_target_`` and sets a pending bit.  The target is
therefore historical state: a tuner commit at the next IMU boundary may change
``Sigma_aw_stat`` but must NOT change the already queued floor target.

The positive-part floor consumes that snapshotted target inside the next
prediction.  Legacy block replacement and congruent immediate synchronization
are real configurable branches but fail closed here until separately composed.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M


def R(x): return P.rational(x)

@dataclass(frozen=True)
class State:
    pending: bool=False
    last_sync_time: F=F(0)
    target: tuple|None=None
    def __post_init__(self):
        if not isinstance(self.pending,bool): raise TypeError('literal aw-floor pending bit required')
        t=R(self.last_sync_time)
        if t<0: raise ValueError('nonnegative aw-sync clock required')
        object.__setattr__(self,'last_sync_time',t)
        if self.target is not None:
            A=M.mat(self.target,3,3)
            if A != M.transpose(A): raise ValueError('symmetric queued aw-floor target required')
            object.__setattr__(self,'target',tuple(map(tuple,A)))
        if self.pending and self.target is None:
            raise ValueError('pending aw-floor request must retain its snapshotted target')

@dataclass(frozen=True)
class Result:
    state: State
    requested_now: bool


def tick(state:State,*,time,adapt_every,live,active_sigma=None,
         enabled=True,congruent=False,legacy=False):
    if not isinstance(state,State): raise TypeError('periodic aw-sync state required')
    if not all(isinstance(x,bool) for x in (live,enabled,congruent,legacy)):
        raise TypeError('literal aw-sync policy branches required')
    time,adapt_every=R(time),R(adapt_every)
    if time<state.last_sync_time or adapt_every<0: raise ValueError('monotone time and nonnegative cadence required')
    if congruent or legacy:
        raise NotImplementedError('non-default immediate aw synchronization branch is not represented by queued-floor lemma')
    due=enabled and live and time-state.last_sync_time>adapt_every
    if not due:
        if active_sigma is not None:
            raise ValueError('not-due aw-sync branch consumes no stationary-covariance operand')
        return Result(state,False)
    if active_sigma is None: raise ValueError('due aw-sync branch requires current active stationary covariance')
    A=M.mat(active_sigma,3,3)
    if A != M.transpose(A): raise ValueError('symmetric current active stationary covariance required')
    return Result(State(True,time,tuple(map(tuple,A))),True)


def floor_target(state:State):
    if not isinstance(state,State): raise TypeError('periodic aw-sync state required')
    return state.target if state.pending else None


def consume_at_prediction(state:State):
    """Prediction clears the request regardless of floor eigensolver outcome."""
    if not isinstance(state,State): raise TypeError('periodic aw-sync state required')
    return State(False,state.last_sync_time,None) if state.pending else state


def readiness():
    return {
      'deployed_periodic_queue_predicate_materialized':True,
      'queued_floor_request_persists_to_next_prediction':True,
      'queued_floor_target_snapshotted_at_request':True,
      'next_boundary_tuner_commit_cannot_rewrite_queued_target':True,
      'prediction_consumes_pending_request':True,
      'legacy_immediate_replacement_branch_attached':False,
      'congruent_immediate_sync_branch_attached':False,
      'clock_binary64_roundoff_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
