"""Next-sample pending boundary for dual-compiler SpectralMSE R_S histories.

Shipping updates ``tune_.RS_applied`` on physical sample k and only sets
``online_tune_apply_pending_``.  At the beginning of IMU sample k+1,
``apply_pending_online_tune_`` consumes the ALREADY-CARRIED TuneState and, when
Live, ``apply_RS_tune_`` writes the active S pseudo-measurement covariance.

This module keeps that one-sample staging exact for both persistent compiler
histories.  It never consumes a same-sample candidate result directly: callers
must first advance ``finite_tuner_rs_deployment_ledger.State`` and carry it to
the next boundary.

MAG/HOLD events do not run this IMU boundary and therefore preserve both stored
R_S tracks and any pending bit by identity.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_tuner_rs_commit_binary32 as C
from tools.stability.ou3_alt_contraction import finite_tuner_rs_deployment_ledger as L

QUALIFICATION='OU3_ALT_RS_BOUNDARY_DEPLOYMENT_V1'


@dataclass(frozen=True)
class State:
    ledger:L.State
    pending:bool=False
    def __post_init__(self):
        if not isinstance(self.ledger,L.State): raise TypeError('persistent RS ledger required')
        if not isinstance(self.pending,bool): raise TypeError('literal pending bit required')


@dataclass(frozen=True)
class Boundary:
    state:State
    separate_commit:C.Commit|None
    fma_commit:C.Commit|None
    consumed:bool
    def __post_init__(self):
        if not isinstance(self.state,State) or not isinstance(self.consumed,bool):
            raise TypeError('RS boundary successor and consumed flag required')
        if self.consumed:
            if not isinstance(self.separate_commit,C.Commit) or not isinstance(self.fma_commit,C.Commit):
                raise TypeError('consumed Live RS boundary requires both compiler commits')
            if self.state.pending: raise ValueError('consumed RS pending bit must clear')
            if self.separate_commit.stored_RS!=self.state.ledger.separate or self.fma_commit.stored_RS!=self.state.ledger.fma:
                raise ValueError('RS commits detached from carried ledger snapshot')
        elif self.separate_commit is not None or self.fma_commit is not None:
            raise ValueError('nonconsumed RS boundary cannot carry commits')


def begin(ledger:L.State|None=None,*,pending=False):
    return State(L.initial() if ledger is None else ledger,pending)


def after_sample(state:State,result:L.Result,*,pending_after):
    """Carry sample-k candidate successor; do NOT commit it on sample k."""
    if not isinstance(state,State) or not isinstance(result,L.Result):
        raise TypeError('RS boundary state and ledger result required')
    if result.separate.ema.previous!=state.ledger.separate or result.fma.ema.previous!=state.ledger.fma:
        raise ValueError('RS sample result detached from boundary predecessor')
    if not isinstance(pending_after,bool): raise TypeError('literal pending-after bit required')
    return State(result.state,pending_after)


def imu_boundary(state:State,*,live,min_RS,max_RS,rs_scale=1,x_factor=C.DEFAULT_X_FACTOR,y_factor=C.DEFAULT_Y_FACTOR):
    """Apply pending candidate only at the following IMU boundary."""
    if not isinstance(state,State): raise TypeError('RS boundary State required')
    if not isinstance(live,bool): raise TypeError('literal Live branch required')
    if not state.pending or not live:
        # Pre-Live pending is consumed by the OU transaction but does not write
        # R_S.  This RS-only ledger therefore records no covariance command;
        # the common full tuner boundary owns the pending-bit transaction.
        return Boundary(State(state.ledger,False if state.pending else state.pending),None,None,False)
    sc=C.commit(state.ledger.separate,min_RS=min_RS,max_RS=max_RS,rs_scale=rs_scale,x_factor=x_factor,y_factor=y_factor)
    fc=C.commit(state.ledger.fma,min_RS=min_RS,max_RS=max_RS,rs_scale=rs_scale,x_factor=x_factor,y_factor=y_factor)
    return Boundary(State(state.ledger,False),sc,fc,True)


def mag_or_hold(state:State):
    if not isinstance(state,State): raise TypeError('RS boundary State required')
    return state


def readiness():
    return {
      'sample_k_RS_candidate_is_carried_without_same_sample_commit':True,
      'pending_RS_is_consumed_only_at_following_IMU_boundary':True,
      'following_Live_boundary_commits_carried_separate_and_FMA_snapshots':True,
      'MAG_and_HOLD_are_literal_RS_ledger_and_pending_identity':True,
      'preLive_RS_boundary_performs_no_RS_covariance_write':True,
      'full_common_tuner_pending_transaction_attached':False,
      'active_MEKF_RS_Eigen_correspondence_closed':False,
      'startup_frontend_RS_machine_history_attached':False,
      'Live_600_step_RS_machine_history_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
