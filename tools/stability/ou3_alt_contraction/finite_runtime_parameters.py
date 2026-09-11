"""Compatibility bridge from staged tuner commit to finite runtime events.

Shipping ``apply_ou_tune_`` calls ``set_aw_stationary_std``, which installs a
diagonal stationary covariance and sets ``aw_process_correlated_=false``.  This
bridge prevents the finite prediction word from combining that commit with an
impossible correlated-Qaxis branch or a different tau/Sigma/period/R_S.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_tuner_commit as T
from tools.stability.ou3_alt_contraction import finite_prediction_runtime as PR


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


def readiness():
    return {
      'committed_tau_to_prediction_tau':True,
      'committed_stationary_covariance_to_Qaxis':True,
      'tuner_commit_forces_independent_Qaxis_branch':True,
      'committed_period_to_scheduler':True,
      'committed_live_RS_to_S_measurement':True,
      'TuneState_frontend_history_attached':False,
      'commit_binary32_roundoff_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
