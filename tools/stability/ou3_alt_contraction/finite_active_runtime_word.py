"""Bind one committed ActiveParameters state to shipping runtime events.

The tuner-boundary proof produces ``ActiveParameters``.  This layer makes that
object a mandatory ancestor of the ordinary prediction and S-service prefix:
OU tau and stationary Sigma_aw must match prediction, the Qaxis branch must be
the independent branch installed by ``set_aw_stationary_std``, scheduler period
must match committed pseudo cadence, and a due Live S event must use the same
committed R_S.

This is ancestry/composition only.  Runtime transcendental, PSD/eigensolver,
LDLT and sensor/source disturbance bounds remain open.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_prediction_runtime as PR
from tools.stability.ou3_alt_contraction import finite_post_prediction as POST
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MR


@dataclass(frozen=True)
class Predicted:
    active: ACTIVE.ActiveParameters
    state: object


def prediction_from_active(active:ACTIVE.ActiveParameters,state,segment,sample,*,
                           angular,Qbase,ou,bias,qaxis:PR.QAxisBranch,
                           use_exact_attitude_Q=True,
                           attitude_first_ldlt_success=True,
                           attitude_second_ldlt_success=None):
    if not isinstance(active,ACTIVE.ActiveParameters): raise TypeError('committed ActiveParameters required')
    active.require_prediction(ou=ou,qaxis=qaxis)
    out=PR.prediction_from_raw(state,segment,sample,angular=angular,Qbase=Qbase,
        ou=ou,bias=bias,qaxis=qaxis,use_exact_attitude_Q=use_exact_attitude_Q,
        attitude_first_ldlt_success=attitude_first_ldlt_success,
        attitude_second_ldlt_success=attitude_second_ldlt_success)
    return Predicted(active,out)


def post_prediction_from_active(predicted:Predicted,*,h,pending_aw_floor,aw_floor_target,
                                scheduler:POST.Scheduler,floor_solver_success=None,
                                floor_eigenvectors=None,floor_eigenvalues=None):
    if not isinstance(predicted,Predicted): raise TypeError('active-rooted prediction required')
    predicted.active.require_scheduler(scheduler)
    return POST.post_prediction_prefix(predicted.state,h=h,pending_aw_floor=pending_aw_floor,
        aw_floor_target=aw_floor_target,scheduler=scheduler,
        floor_solver_success=floor_solver_success,floor_eigenvectors=floor_eigenvectors,
        floor_eigenvalues=floor_eigenvalues)


def service_S_from_active(active:ACTIVE.ActiveParameters,prefix:POST.PostPrediction,*,
                          ldlt:MR.SafeLDLT|None=None,alpha=1,radius=None):
    if not isinstance(active,ACTIVE.ActiveParameters) or not isinstance(prefix,POST.PostPrediction):
        raise TypeError('active parameters and post-prediction prefix required')
    if not prefix.S_service_due:
        if ldlt is not None: raise ValueError('not-due active S branch consumes no LDLT witness')
        return POST.ServicedPostPrediction(prefix.state,prefix,None)
    if active.R_S is None: raise ValueError('due S service has no Live R_S in committed active state')
    active.require_RS(active.R_S)
    kwargs={'R_S':active.R_S,'ldlt':ldlt,'alpha':alpha}
    if radius is not None: kwargs['radius']=radius
    return POST.service_S_if_due(prefix,**kwargs)


def readiness():
    return {
      'committed_active_tau_sigma_to_raw_prediction':True,
      'committed_active_period_to_scheduler':True,
      'committed_active_live_RS_to_due_S_service':True,
      'same_history_tuner_parameters_attached_at_event_interface':True,
      'independent_Qaxis_branch_enforced_after_tuner_commit':True,
      'runtime_transcendental_and_solver_finite_precision_attached':False,
      'sensor_residual_source_bounds_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
