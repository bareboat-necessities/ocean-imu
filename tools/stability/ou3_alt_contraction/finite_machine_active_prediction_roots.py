"""Source-owned prediction roots driven by carried machine ActiveParameters.

The exact admitted Live word already derives attitude, OU, Q-axis and BA roots
from one source-owned physical segment. Deployment differs only in which
applied tuner parameters feed the OU/Q-axis coefficient paths. This module
reuses the same source/packet ancestry and the existing primitive validators,
but substitutes one carried global compiler-mode ``ActiveParameters`` object for
``tau`` and ``Sigma_aw``.

A subtle ordering point matters: a pending tuner boundary is consumed BEFORE
prediction on the same IMU sample. Therefore the exact ActiveParameters used by
the executed prediction may differ from the ActiveParameters stored in the
pre-boundary source-word state. ``build`` accepts that executed exact active
state explicitly; omitting it retains the historical no-boundary/component
behavior. The exact-vs-machine join is always formed against the selected
executed active state, not silently against a stale pre-boundary shadow.

This is not a parallel physical history and it does not execute a second source
transition. It is a coefficient relation on the SAME event. Scheduler and R_S
measurement effects are separate obligations.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as WORD
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_attitude_runtime as ATT
from tools.stability.ou3_alt_contraction import finite_source_bound_attitude_trig as TRIG
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_ou_runtime_primitives as OU
from tools.stability.ou3_alt_contraction import finite_prediction_runtime as PRED
from tools.stability.ou3_alt_contraction import finite_qaxis_binary32_branch as QB
from tools.stability.ou3_alt_contraction import finite_source_bound_prediction_word as BASE
from tools.stability.ou3_alt_contraction import finite_runtime_parameters as ACTIVE
from tools.stability.ou3_alt_contraction import finite_active_parameter_machine_real_join as JOIN

QUALIFICATION='OU3_ALT_MACHINE_ACTIVE_PREDICTION_ROOTS_V2'


@dataclass(frozen=True)
class Roots:
    active:ACTIVE.ActiveParameters
    active_join:JOIN.Join
    angular:ATT.AngularRuntime
    ou:OU.OUDecay
    qaxis:PRED.QAxisBranch
    bias:OU.BiasDecay
    qualification:str=QUALIFICATION
    def __post_init__(self):
        if not isinstance(self.active,ACTIVE.ActiveParameters) or not isinstance(self.active_join,JOIN.Join):
            raise TypeError('machine active parameters and exact/machine join required')
        if self.active_join.machine!=self.active:
            raise ValueError('machine prediction roots detached from active-parameter join')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine prediction-root qualification')
        self.active.require_prediction(ou=self.ou,qaxis=self.qaxis)


def build(state:WORD.State,physical:SOURCE.QualifiedPhysicalSegment,
          raw:SENSOR.RawImuSample,machine_active:ACTIVE.ActiveParameters,*,mode,
          exact_active:ACTIVE.ActiveParameters|None=None,
          ou_alpha,ou_em1,bias_phi=None,bias_em1_2=None,
          angular_full=None,angular_half=None,
          qaxis_marginal_exp=None,qaxis_final_exp=None,
          qaxis_marginal_psd,qaxis_final_psd):
    if not isinstance(state,WORD.State) or not isinstance(physical,SOURCE.QualifiedPhysicalSegment):
        raise TypeError('source-owning Live state and checked physical segment required')
    if not isinstance(raw,SENSOR.RawImuSample) or not isinstance(machine_active,ACTIVE.ActiveParameters):
        raise TypeError('same-source raw packet and carried machine ActiveParameters required')
    if physical.root!=state.source.root or physical.witness.ordinal!=state.source.next_ordinal:
        raise ValueError('machine prediction roots detached from next carried source transition')
    core=state.live.live.live.mekf; segment=physical.segment
    if segment.before!=core.reference or raw.physical!=segment.before:
        raise ValueError('machine prediction roots detached from current physical predecessor')

    if exact_active is None:
        exact_active=state.live.live.live.active
    if not isinstance(exact_active,ACTIVE.ActiveParameters):
        raise TypeError('executed exact ActiveParameters required')
    aj=JOIN.join(exact_active,machine_active,mode)
    omega_hat=raw.required_bias_corrected_relation(core.z[3:6])
    angular=ATT.AngularRuntime(tuple(omega_hat),segment.h,full=angular_full,half=angular_half)
    TRIG.validate(angular)
    ou=OU.OUDecay(segment.h,machine_active.tau,ou_alpha,em1=ou_em1)
    EXP.validate_ou(ou)
    coefficient_branch=QB.branch(machine_active.tau,segment.h)
    covariance_exp=BASE._qaxis_exp_pair(coefficient_branch,qaxis_marginal_exp,qaxis_final_exp)
    qaxis=PRED.QAxisBranch(False,machine_active.Sigma_aw,
                           tuple(qaxis_marginal_psd),tuple(qaxis_final_psd),
                           BASE.SHIPPING_FLOAT_EPSILON,coefficient_branch,covariance_exp)
    machine_active.require_prediction(ou=ou,qaxis=qaxis)
    bias=BASE._bias_root(core,h=segment.h,bias_phi=bias_phi,bias_em1_2=bias_em1_2)
    return Roots(machine_active,aj,angular,ou,qaxis,bias)


def readiness():
    return {
      'same_source_physical_segment_and_raw_packet_as_exact_prediction_required':True,
      'machine_applied_tau_drives_OU_decay_and_Qaxis_branch_selection':True,
      'machine_applied_stationary_Sigma_drives_Qaxis_covariance':True,
      'attitude_and_BA_roots_reuse_same_source_owned_event':True,
      'small_general_Qaxis_branch_is_recomputed_from_machine_tau_not_shadow_tau':True,
      'exact_vs_machine_active_parameter_join_retained_with_roots':True,
      'executed_post_boundary_exact_active_can_be_bound_explicitly':True,
      'machine_prediction_root_relation_attached':True,
      'machine_root_effect_injected_into_joint24_event_relation':False,
      'machine_pseudo_period_scheduler_effect_attached':False,
      'machine_RS_measurement_effect_attached':False,
      'OU_Qaxis_target_libm_and_Eigen_correspondence_closed':False,
      'source_uniform_machine_coefficient_supply_bound_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
