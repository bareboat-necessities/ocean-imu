"""Bind prediction coefficient roots to the same source-owning Live word.

``finite_source_bound_live_word`` owns the next physical transition, raw IMU
packet, active TuneState-derived parameters and static runtime config. This
layer removes independently injectable coefficient roots while retaining the
actual shipping transcendental topology:

* attitude angular rate is reconstructed from the exact raw packet and current
  gyro-bias error;
* OU h/tau is reconstructed from the exact source segment and carried active
  tau; the distinct exp(-x) and expm1(-x) results remain arithmetic witnesses
  of that same argument until binary32/libm correspondence is enclosed;
* the independent-axis Q branch is reconstructed from carried Sigma_aw and its
  regularization epsilon is fixed by Kalman3D_Wave_OU_III<float> to 2^-23;
* residual accelerometer-bias OU ownership is the shipping core's fixed
  tau_bacc_=5000 s and Q_bacc_=2.5e-7 I. H18 has phi_b=1 exactly. A21 retains
  separate exp(-h/tau_b) and expm1(-2h/tau_b) arithmetic results.

No equality between separately rounded exp/expm1 calls is silently assumed.
Their common mathematical arguments are source-owned here; their deployment
error relation remains fail-closed. This is not a storage certificate.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as WORD
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_attitude_runtime as ATT
from tools.stability.ou3_alt_contraction import finite_ou_runtime_primitives as OU
from tools.stability.ou3_alt_contraction import finite_prediction_runtime as PRED

SHIPPING_BA_TAU = F(5000)
SHIPPING_BA_Q = tuple(tuple(F(1,4_000_000) if i == j else F(0)
                            for j in range(3)) for i in range(3))
SHIPPING_FLOAT_EPSILON = F(1, 1 << 23)


@dataclass(frozen=True)
class Roots:
    angular: ATT.AngularRuntime
    ou: OU.OUDecay
    qaxis: PRED.QAxisBranch
    bias: OU.BiasDecay


def _bias_root(core, *, bias_phi=None, bias_em1_2=None):
    """Literal shipping BA prediction topology with no fake exp identity."""
    active = core.mode == 'A'
    if active:
        if bias_phi is None or bias_em1_2 is None:
            raise ValueError('active A21 bias prediction requires exp and expm1 arithmetic witnesses')
        phi = bias_phi
    else:
        if bias_phi is not None or bias_em1_2 is not None:
            raise ValueError('held H18 bias prediction consumes no exp/expm1 arithmetic witness')
        phi = F(1)
    return OU.BiasDecay(active, SHIPPING_BA_TAU, phi, SHIPPING_BA_Q,
                        em1_2=bias_em1_2)


def build(state: WORD.State, physical: SOURCE.QualifiedPhysicalSegment,
          raw: SENSOR.RawImuSample, *, ou_alpha, ou_em1,
          bias_phi=None, bias_em1_2=None,
          angular_full: ATT.TrigWitness | None = None,
          angular_half: ATT.TrigWitness | None = None,
          qaxis_marginal_psd, qaxis_final_psd):
    """Construct prediction roots from one exact source predecessor."""
    if not isinstance(state, WORD.State):
        raise TypeError('source-owning Live state required')
    if not isinstance(physical, SOURCE.QualifiedPhysicalSegment):
        raise TypeError('checked next physical source transition required')
    if not isinstance(raw, SENSOR.RawImuSample):
        raise TypeError('same-source raw IMU packet required')
    if physical.root != state.source.root:
        raise ValueError('prediction source transition detached from carried root')
    if physical.witness.ordinal != state.source.next_ordinal:
        raise ValueError('prediction roots require exactly the next source ordinal')
    core = state.live.live.live.mekf
    segment = physical.segment
    if segment.before != core.reference or raw.physical != segment.before:
        raise ValueError('prediction roots detached from current physical predecessor')

    omega_hat = raw.required_bias_corrected_relation(core.z[3:6])
    angular = ATT.AngularRuntime(tuple(omega_hat), segment.h,
                                 full=angular_full, half=angular_half)

    active = state.live.live.live.active
    ou = OU.OUDecay(segment.h, active.tau, ou_alpha, em1=ou_em1)

    qaxis = PRED.QAxisBranch(False, active.Sigma_aw,
                             tuple(qaxis_marginal_psd), tuple(qaxis_final_psd),
                             SHIPPING_FLOAT_EPSILON)
    active.require_prediction(ou=ou, qaxis=qaxis)

    bias = _bias_root(core, bias_phi=bias_phi, bias_em1_2=bias_em1_2)
    return Roots(angular, ou, qaxis, bias)


def imu_step(state: WORD.State, *, witness: SOURCE.StepWitness,
             segment, raw: SENSOR.RawImuSample, packet_id: str,
             ou_alpha, ou_em1, bias_phi=None, bias_em1_2=None,
             angular_full=None, angular_half=None,
             qaxis_marginal_psd, qaxis_final_psd,
             **dynamic):
    """Execute one source-owning IMU edge with source-bound prediction roots."""
    forbidden = {'angular','ou','qaxis','bias','machine_epsilon'} & set(dynamic)
    if forbidden:
        raise TypeError('source-bound prediction roots cannot be overridden '+repr(sorted(forbidden)))
    if not isinstance(state, WORD.State):
        raise TypeError('source-owning Live state required')
    if not isinstance(witness, SOURCE.StepWitness):
        raise TypeError('checked source step witness required')
    if witness.ordinal != state.source.next_ordinal:
        raise ValueError('prediction roots require exactly the next source ordinal')
    physical = SOURCE.QualifiedPhysicalSegment(state.source.root, witness, segment)
    roots = build(state, physical, raw, ou_alpha=ou_alpha,ou_em1=ou_em1,
                  bias_phi=bias_phi,bias_em1_2=bias_em1_2,
                  angular_full=angular_full, angular_half=angular_half,
                  qaxis_marginal_psd=qaxis_marginal_psd,
                  qaxis_final_psd=qaxis_final_psd)
    return WORD.imu_step(state, witness=witness, segment=segment, raw=raw,
                         packet_id=packet_id, angular=roots.angular,
                         ou=roots.ou, qaxis=roots.qaxis, bias=roots.bias, **dynamic)


def readiness():
    lower = WORD.readiness()
    return {
      'source_owned_next_transition_required': True,
      'attitude_omega_reconstructed_from_same_raw_packet_and_current_bias_error': True,
      'OU_h_tau_argument_from_same_source_and_active_TuneState': True,
      'OU_exp_and_expm1_shipping_results_retained_separately': True,
      'Qaxis_Sigma_aw_from_same_carried_active_TuneState': True,
      'shipping_independent_Qaxis_branch_fixed_by_active_parameter_contract': True,
      'Qaxis_machine_epsilon_bound_to_shipping_binary32': True,
      'prediction_angular_OU_Qaxis_roots_bound_to_same_source_continuation': True,
      'shipping_BA_tau_Q_defaults_bound_at_prediction_entry': True,
      'shipping_BA_hold_active_branch_derived_from_current_MEKF_mode': True,
      'held_H18_BA_phi_exactly_one_without_transcendental_witness': True,
      'active_A21_BA_exp_and_expm1_results_retained_separately': True,
      'caller_cannot_override_prediction_roots': True,
      'sample_zero_origin_bridge_preserved': lower['sample_zero_startup_to_checked_outer_endpoint_bridge_closed'],
      'OU_exp_expm1_binary32_relation_closed': False,
      'BA_exp_expm1_binary32_relation_closed': False,
      'attitude_trig_arithmetic_deployment_closed': False,
      'Qaxis_PSD_regularization_deployment_closed': False,
      'accelerometer_bias_prediction_root_bound_here': True,
      'all_estimator_coefficients_bound_to_same_source_continuation': False,
      'quantitative_sensor_residual_ISS_envelope_attached': False,
      'source_uniform_complete_600_step_word_qualified': False,
      'storage_search_allowed': False,
      'ALT_LIVE_PASS': False,
      'ALT_STARTUP_PASS': False,
      'ALT_END_TO_END_PASS': False,
    }
