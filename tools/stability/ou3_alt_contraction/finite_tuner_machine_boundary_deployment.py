"""Common pending/commit boundary for the dual-compiler machine TuneState.

Shipping owns one ``TuneState{tau_applied,sigma_applied,RS_applied}`` and one
``online_tune_apply_pending_`` flag.  The machine-state product already proves
that tau, sigma and R_S are advanced on one coherent sample history.  This
module closes the next step in that transaction: a pending snapshot is consumed
ONLY at the following IMU boundary, tau/sigma are committed from that exact
carried TuneState, and Live additionally applies R_S from the same snapshot.

The two compiler histories remain global alternatives; no per-sample compiler
switch is introduced.  MAG/HOLD are exact identity on this machine transaction.
The final Eigen ``set_RS_noise`` execution correspondence remains open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_binary32_arithmetic as B
from tools.stability.ou3_alt_contraction import finite_tuner_machine_tunestate_product as T
from tools.stability.ou3_alt_contraction import finite_tuner_rs_commit_binary32 as RSC

QUALIFICATION='OU3_ALT_MACHINE_TUNESTATE_BOUNDARY_V1'


@dataclass(frozen=True)
class Snapshot:
    mode:str
    tau:F
    sigma:F
    RS:F
    rs_commit:RSC.Commit|None
    def __post_init__(self):
        if self.mode not in ('separate','fma'): raise ValueError('invalid compiler mode')
        for n in ('tau','sigma','RS'):
            object.__setattr__(self,n,F(getattr(self,n)))
            if not B.is_binary32(getattr(self,n)):
                raise ValueError('machine TuneState boundary stores binary32 scalars only')
        if self.rs_commit is not None:
            if not isinstance(self.rs_commit,RSC.Commit): raise TypeError('RS commit object required')
            if self.rs_commit.stored_RS!=self.RS:
                raise ValueError('R_S commit detached from same machine TuneState snapshot')


@dataclass(frozen=True)
class Boundary:
    state:T.State
    separate:Snapshot|None
    fma:Snapshot|None
    consumed:bool
    live:bool
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.state,T.State): raise TypeError('machine TuneState successor required')
        if not isinstance(self.consumed,bool) or not isinstance(self.live,bool):
            raise TypeError('literal boundary branch flags required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine boundary qualification')
        if self.consumed:
            if not isinstance(self.separate,Snapshot) or not isinstance(self.fma,Snapshot):
                raise TypeError('consumed boundary requires both compiler snapshots')
            if self.state.pending: raise ValueError('consumed common pending bit must clear')
            if self.separate.mode!='separate' or self.fma.mode!='fma':
                raise ValueError('compiler histories crossed at boundary')
            if self.separate.tau!=self.state.tau.separate or self.fma.tau!=self.state.tau.fma:
                raise ValueError('tau commit detached from carried machine TuneState')
            if self.separate.sigma!=self.state.sigma.separate or self.fma.sigma!=self.state.sigma.fma:
                raise ValueError('sigma commit detached from carried machine TuneState')
            if self.separate.RS!=self.state.rs.separate or self.fma.RS!=self.state.rs.fma:
                raise ValueError('R_S commit detached from carried machine TuneState')
            if self.live and (self.separate.rs_commit is None or self.fma.rs_commit is None):
                raise ValueError('Live common boundary requires R_S application for both compiler tracks')
            if not self.live and (self.separate.rs_commit is not None or self.fma.rs_commit is not None):
                raise ValueError('pre-Live common boundary must not write R_S covariance')
        elif self.separate is not None or self.fma is not None:
            raise ValueError('nonconsumed boundary cannot expose commit snapshots')


def _snapshot(state:T.State,mode,*,live,min_RS,max_RS,rs_scale,x_factor,y_factor):
    if mode=='separate':
        tau,sigma,rs=state.tau.separate,state.sigma.separate,state.rs.separate
    elif mode=='fma':
        tau,sigma,rs=state.tau.fma,state.sigma.fma,state.rs.fma
    else: raise ValueError('invalid compiler mode')
    rc=None
    if live:
        rc=RSC.commit(rs,min_RS=min_RS,max_RS=max_RS,rs_scale=rs_scale,
                      x_factor=x_factor,y_factor=y_factor)
    return Snapshot(mode,tau,sigma,rs,rc)


def imu_boundary(state:T.State,*,live,min_RS,max_RS,rs_scale=1,
                 x_factor=RSC.DEFAULT_X_FACTOR,y_factor=RSC.DEFAULT_Y_FACTOR):
    """Consume the carried common pending bit at the following IMU boundary."""
    if not isinstance(state,T.State): raise TypeError('machine TuneState State required')
    if not isinstance(live,bool): raise TypeError('literal Live branch required')
    if not state.pending:
        return Boundary(state,None,None,False,live)
    nxt=T.State(state.tau,state.sigma,state.rs,False)
    s=_snapshot(nxt,'separate',live=live,min_RS=min_RS,max_RS=max_RS,
                rs_scale=rs_scale,x_factor=x_factor,y_factor=y_factor)
    f=_snapshot(nxt,'fma',live=live,min_RS=min_RS,max_RS=max_RS,
                rs_scale=rs_scale,x_factor=x_factor,y_factor=y_factor)
    return Boundary(nxt,s,f,True,live)


def mag_or_hold(state:T.State):
    if not isinstance(state,T.State): raise TypeError('machine TuneState State required')
    return state


def readiness():
    return {
      'full_tau_sigma_RS_machine_TuneState_pending_transaction_composed':True,
      'pending_consumed_only_at_following_IMU_boundary':True,
      'tau_and_sigma_commit_from_same_carried_machine_snapshot':True,
      'Live_RS_commit_uses_same_carried_machine_snapshot':True,
      'preLive_pending_clears_without_RS_covariance_write':True,
      'separate_and_FMA_histories_remain_global_and_coherent':True,
      'MAG_and_HOLD_are_literal_machine_TuneState_identity':True,
      'Eigen_set_RS_noise_execution_correspondence_closed':False,
      'startup_frontend_machine_TuneState_product_attached':False,
      'Live_600_step_machine_TuneState_product_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
