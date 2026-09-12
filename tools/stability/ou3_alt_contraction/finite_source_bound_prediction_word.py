"""Bind shipping IMU coefficient roots to the same source-owning Live word.

``finite_source_bound_live_word`` owns the next physical transition, raw IMU
packet, active TuneState-derived parameters and static runtime config. This
layer removes independently injectable coefficient roots while retaining the
actual shipping transcendental topology:

* attitude angular rate comes from the exact raw packet/current gyro-bias error;
  nonzero-rate sin/cos witnesses must also lie in rigorous rational enclosures
  at that SAME source-owned rotation angle;
* OU h/tau comes from the source segment and carried TuneState; distinct
  exp(-x) and expm1(-x) results are separately constrained by rigorous real
  enclosures at that SAME h/tau argument;
* Qaxis Sigma_aw comes from TuneState and regularization epsilon is fixed by
  the shipping ``Kalman3D_Wave_OU_III<float>`` instantiation to 2^-23. Its
  polynomial/general formula decision is derived from exact binary32 rounding
  of those same source-owned h/tau operands, rather than from exact-real h/tau;
* residual BA tau/Q are the current shipping core defaults. H18 consumes no BA
  transcendental; A21 retains distinct exp and expm1 results and constrains each
  to its literal h/tau_b or 2h/tau_b source argument;
* accelerometer temperature-model coefficient is the shipping core default
  k_a=(0.002,0.002,0.002) at tempC_ref=35 C. The caller supplies only the
  actual temperature input; it cannot replace k_a or lever-arm conditioning.

Binary32/libm coefficient correspondence, tuner-commit storage correspondence,
temperature-history admissibility, PSD solver outcomes and quantitative
roundoff remain open. This is not storage or a stability certificate.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as WORD
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_attitude_runtime as ATT
from tools.stability.ou3_alt_contraction import finite_source_bound_attitude_trig as TRIG
from tools.stability.ou3_alt_contraction import finite_source_bound_exp_enclosure as EXP
from tools.stability.ou3_alt_contraction import finite_ou_runtime_primitives as OU
from tools.stability.ou3_alt_contraction import finite_prediction_runtime as PRED
from tools.stability.ou3_alt_contraction import finite_qaxis_binary32_branch as QB

SHIPPING_BA_TAU = F(5000)
SHIPPING_BA_Q = tuple(tuple(F(1,4_000_000) if i == j else F(0)
                            for j in range(3)) for i in range(3))
SHIPPING_FLOAT_EPSILON = F(1, 1 << 23)
SHIPPING_TEMP_REF_C = F(35)
SHIPPING_KA = (F(1,500), F(1,500), F(1,500))
ZERO3 = (F(0),F(0),F(0))


@dataclass(frozen=True)
class Roots:
    angular: ATT.AngularRuntime
    ou: OU.OUDecay
    qaxis: PRED.QAxisBranch
    bias: OU.BiasDecay


def _bias_root(core, *, h=None, bias_phi=None, bias_em1_2=None):
    active = core.mode == 'A'
    if active:
        if bias_phi is None or bias_em1_2 is None or h is None:
            raise ValueError('active A21 bias prediction requires h plus exp and expm1 arithmetic witnesses')
        phi = bias_phi
    else:
        if bias_phi is not None or bias_em1_2 is not None:
            raise ValueError('held H18 bias prediction consumes no exp/expm1 arithmetic witness')
        phi = F(1)
    out=OU.BiasDecay(active, SHIPPING_BA_TAU, phi, SHIPPING_BA_Q,
                     em1_2=bias_em1_2)
    if active:
        EXP.validate_bias(out,h)
    return out


def _accel_conditioning(temperature_c):
    temp=SENSOR.R(temperature_c)
    return SENSOR.AccelConditioning(temp-SHIPPING_TEMP_REF_C, SHIPPING_KA, ZERO3)


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
    TRIG.validate(angular)
    active = state.live.live.live.active
    ou = OU.OUDecay(segment.h, active.tau, ou_alpha, em1=ou_em1)
    EXP.validate_ou(ou)
    coefficient_branch=QB.branch(active.tau,segment.h)
    qaxis = PRED.QAxisBranch(False, active.Sigma_aw,
                             tuple(qaxis_marginal_psd), tuple(qaxis_final_psd),
                             SHIPPING_FLOAT_EPSILON, coefficient_branch)
    active.require_prediction(ou=ou, qaxis=qaxis)
    bias = _bias_root(core,h=segment.h,bias_phi=bias_phi,bias_em1_2=bias_em1_2)
    return Roots(angular, ou, qaxis, bias)


def imu_step(state: WORD.State, *, witness: SOURCE.StepWitness,
             segment, raw: SENSOR.RawImuSample, packet_id: str,
             temperature_c, ou_alpha, ou_em1,
             bias_phi=None, bias_em1_2=None,
             angular_full=None, angular_half=None,
             qaxis_marginal_psd, qaxis_final_psd,
             **dynamic):
    """Execute one source-owning IMU edge with source-bound model roots."""
    forbidden = {'angular','ou','qaxis','bias','machine_epsilon','accel_conditioning'} & set(dynamic)
    if forbidden:
        raise TypeError('source-bound IMU roots cannot be overridden '+repr(sorted(forbidden)))
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
    conditioning=_accel_conditioning(temperature_c)
    return WORD.imu_step(state, witness=witness, segment=segment, raw=raw,
                         packet_id=packet_id, angular=roots.angular,
                         ou=roots.ou, qaxis=roots.qaxis, bias=roots.bias,
                         accel_conditioning=conditioning, **dynamic)


def readiness():
    lower = WORD.readiness(); trig=TRIG.readiness(); exp=EXP.readiness(); qb=QB.readiness()
    return {
      'source_owned_next_transition_required': True,
      'attitude_omega_reconstructed_from_same_raw_packet_and_current_bias_error': True,
      'attitude_trig_full_half_bound_to_same_source_owned_rotation_angle': trig['trig_full_half_angles_derived_from_same_angular_rate_and_step'],
      'detached_attitude_unit_circle_points_rejected': trig['detached_unit_circle_points_rejected'],
      'OU_h_tau_argument_from_same_source_and_active_TuneState': True,
      'OU_exp_and_expm1_shipping_results_retained_separately': True,
      'OU_exp_expm1_real_enclosed_at_same_source_argument': bool(
          exp['OU_exp_root_real_enclosed_at_same_h_over_tau'] and
          exp['OU_expm1_root_real_enclosed_at_same_h_over_tau']),
      'Qaxis_Sigma_aw_from_same_carried_active_TuneState': True,
      'shipping_independent_Qaxis_branch_fixed_by_active_parameter_contract': True,
      'Qaxis_machine_epsilon_bound_to_shipping_binary32': True,
      'Qaxis_formula_branch_derived_from_binary32_rounding_of_same_source_tau_h': bool(
          qb['shipping_nested_and_final_Qaxis_share_literal_tau_h_branch_shape'] and
          qb['small_general_comparison_binary32_attached']),
      'Qaxis_source_tau_commit_binary32_correspondence_closed':False,
      'prediction_angular_OU_Qaxis_roots_bound_to_same_source_continuation': True,
      'shipping_BA_tau_Q_defaults_bound_at_prediction_entry': True,
      'shipping_BA_hold_active_branch_derived_from_current_MEKF_mode': True,
      'held_H18_BA_phi_exactly_one_without_transcendental_witness': True,
      'active_A21_BA_exp_and_expm1_results_retained_separately': True,
      'active_A21_BA_exp_expm1_real_enclosed_at_literal_arguments':exp['BA_exp_and_expm1_distinct_arguments_real_enclosed'],
      'accelerometer_temperature_coefficient_bound_to_shipping_default': True,
      'zero_lever_accelerometer_scope_enforced_at_source_bound_entry': True,
      'temperature_is_explicit_per_sample_input_not_free_model_coefficient': True,
      'caller_cannot_override_prediction_or_accel_model_roots': True,
      'sample_zero_origin_bridge_preserved': lower['sample_zero_startup_to_checked_outer_endpoint_bridge_closed'],
      'attitude_one_radian_local_angle_guard_retained_for_every_prefix': False,
      'temperature_history_admissibility_attached': False,
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
