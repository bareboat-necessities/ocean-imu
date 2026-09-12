"""Compatibility bridge from staged tuner commit to finite runtime events.

Shipping ``apply_ou_tune_`` installs diagonal Sigma_aw, tau and a tau-derived
pseudo-S period.  Changing that period does not reset the periodic scheduler:
``set_pseudo_update_period_s`` preserves elapsed credit when it is below the new
period and otherwise parks it at binary ``nextafter(period,0)`` so the next
valid sample services the owed update.  This module keeps that scheduler
retarget on the same active-parameter ancestry.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_tuner_commit as T
from tools.stability.ou3_alt_contraction import finite_prediction_runtime as PR
from tools.stability.ou3_alt_contraction import finite_post_prediction as POST


@dataclass(frozen=True)
class ActiveParameters:
    tau: object
    Sigma_aw: tuple
    pseudo_period: object
    R_S: tuple | None
    def __post_init__(self):
        sig=M.mat(self.Sigma_aw,3,3)
        if sig != M.transpose(sig): raise ValueError('symmetric active Sigma_aw required')
        object.__setattr__(self,'Sigma_aw',tuple(map(tuple,sig)))
        object.__setattr__(self,'tau',P.rational(self.tau))
        object.__setattr__(self,'pseudo_period',P.rational(self.pseudo_period))
        if self.tau<=0 or self.pseudo_period<=0: raise ValueError('positive active tau/period required')

    @classmethod
    def from_commit(cls,result:T.CommitResult):
        if not isinstance(result,T.CommitResult): raise TypeError('completed tuner commit required')
        return cls(result.tau,result.Sigma_aw,result.pseudo_period,result.R_S)

    def require_prediction(self,*,ou,qaxis:PR.QAxisBranch):
        if ou.tau != self.tau: raise ValueError('prediction tau detached from committed TuneState')
        if qaxis.correlated: raise ValueError('set_aw_stationary_std forces shipping independent-axis Q branch')
        if qaxis.sigma_aw != self.Sigma_aw: raise ValueError('prediction Sigma_aw detached from committed stationary covariance')
        return True

    def require_scheduler(self,scheduler):
        if scheduler.period != self.pseudo_period: raise ValueError('scheduler period detached from committed tau/cadence')
        return True

    def require_RS(self,R_S):
        if self.R_S is None: raise ValueError('no Live applied R_S in this committed parameter state')
        if tuple(map(tuple,M.mat(R_S,3,3))) != self.R_S: raise ValueError('measurement R_S detached from same committed TuneState')
        return True


@dataclass(frozen=True)
class NextafterParkWitness:
    """Deployment witness for std::nextafter(new_period,0) on overdue retarget."""
    parked_elapsed: object
    def __post_init__(self):
        x=P.rational(self.parked_elapsed)
        if x<0: raise ValueError('nonnegative parked elapsed required')
        object.__setattr__(self,'parked_elapsed',x)


def retarget_scheduler(active:ActiveParameters,scheduler:POST.Scheduler,*,park:NextafterParkWitness|None=None):
    """Literal control graph of set_pseudo_update_period_s scheduler retarget.

    The exact bit-pattern relation of ``nextafter`` is deliberately not claimed
    here.  On the overdue branch a witness must provide the deployed predecessor
    of the new period; finite-precision closure must later prove that witness.
    """
    if not isinstance(active,ActiveParameters) or not isinstance(scheduler,POST.Scheduler):
        raise TypeError('active parameters and scheduler required')
    newp=active.pseudo_period
    if scheduler.elapsed < newp:
        if park is not None: raise ValueError('credit-preserving retarget consumes no nextafter witness')
        return POST.Scheduler(newp,scheduler.elapsed,scheduler.tolerance)
    if park is None: raise ValueError('overdue period retarget requires nextafter witness')
    x=park.parked_elapsed
    if not 0<=x<newp: raise ValueError('nextafter witness must park immediately below positive period')
    return POST.Scheduler(newp,x,scheduler.tolerance)


def readiness():
    return {
      'committed_tau_to_prediction_tau':True,
      'committed_stationary_covariance_to_Qaxis':True,
      'tuner_commit_forces_independent_Qaxis_branch':True,
      'committed_period_to_scheduler':True,
      'period_change_preserves_subdeadline_elapsed_credit':True,
      'overdue_period_change_nextafter_branch_materialized':True,
      'committed_live_RS_to_S_measurement':True,
      'TuneState_frontend_history_attached':False,
      'nextafter_binary32_ancestry_attached':False,
      'commit_binary32_roundoff_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
