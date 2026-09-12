"""Bind prediction coefficient roots to the same source-owning Live word.

``finite_source_bound_live_word`` already owns the next physical transition,
raw IMU packet, active TuneState-derived parameters and static runtime config.
This layer removes independently injectable prediction roots:

* attitude angular rate is reconstructed from the exact raw packet and current
  gyro-bias error;
* OU h/tau is reconstructed from the exact source segment and carried active
  tau, leaving only the exp(-h/tau) arithmetic witness;
* the shipping independent-axis Q branch is reconstructed from the carried
  active Sigma_aw, leaving only PSD/regularization branch witnesses;
* residual accelerometer-bias OU ownership is the current shipping core's fixed
  tau_bacc_=5000 s and Q_bacc_=2.5e-7 I.  H18/held therefore has phi_b=1
  exactly; A21 consumes only the still-open exp(-h/tau_bacc) arithmetic value.

The transcendental values and LDLT/PSD outcomes are still conditional runtime
witnesses. Frontend transcendental arithmetic, solver finite precision and
quantitative sensor forcing remain open. Therefore this is not a storage
certificate and does not promote any ALT theorem gate.
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
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M

# Current shipping Kalman3D_Wave_OU_III defaults.  The wrapper does not call the
# public BA-model setters.  These constants are therefore part of the present
# shipping implementation being proved, not a theorem tuning surface.
SHIPPING_BA_TAU = F(5000)
SHIPPING_BA_Q = tuple(tuple(F(1,4_000_000) if i == j else F(0)
                            for j in range(3)) for i in range(3))


@dataclass(frozen=True)
class Roots:
    angular: ATT.AngularRuntime
    ou: OU.OUDecay
    qaxis: PRED.QAxisBranch
    bias: OU.BiasDecay


def _bias_root(core, *, bias_phi=None):
    """Literal current shipping BA prediction ownership.

    The H18 branch has no exponential evaluation and consumes no witness.
    A21 consumes the deployed exp value, but h and tau are not caller-owned.
    Its libm/binary32 correspondence remains a separate fail-closed obligation.
    """
    active = core.mode == 'A'
    if active:
        if bias_phi is None:
            raise ValueError('active A21 bias prediction requires exp arithmetic witness')
        phi = bias_phi
    else:
        if bias_phi is not None:
            raise ValueError('held H18 bias prediction consumes no exp arithmetic witness')
        phi = F(1)
    return OU.BiasDecay(active, SHIPPING_BA_TAU, phi, SHIPPING_BA_Q)


def build(state: WORD.State, physical: SOURCE.QualifiedPhysicalSegment,
          raw: SENSOR.RawImuSample, *, ou_alpha, bias_phi=None,
          angular_full: ATT.TrigWitness | None = None,
          angular_half: ATT.TrigWitness | None = None,
          qaxis_marginal_psd, qaxis_final_psd, machine_epsilon):
    """Construct prediction roots from one exact source predecessor.

    No physical, TuneState or BA-model coefficient may be supplied by the
    caller. Remaining arguments are arithmetic/branch witnesses whose
    deployment correspondence is deliberately left open.
    """
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

    # Exact shipping relation after de-heel and subtraction of the estimated
    # gyro bias: omega_hat = omega_sample + e_bg + n_g.
    omega_hat = raw.required_bias_corrected_relation(core.z[3:6])
    angular = ATT.AngularRuntime(tuple(omega_hat), segment.h,
                                 full=angular_full, half=angular_half)

    # h and tau are owned by the source step and current committed active state.
    # alpha remains the finite exp arithmetic witness until deployment closure.
    active = state.live.live.live.active
    ou = OU.OUDecay(segment.h, active.tau, ou_alpha)

    # set_aw_stationary_std selects the shipping independent-axis branch.  The
    # stationary covariance cannot be replaced per event: it is the carried
    # active Sigma_aw from the same TuneState ancestry.
    qaxis = PRED.QAxisBranch(False, active.Sigma_aw,
                             tuple(qaxis_marginal_psd), tuple(qaxis_final_psd),
                             machine_epsilon)
    active.require_prediction(ou=ou, qaxis=qaxis)

    # BA tau and driving density are core-owned shipping defaults on the current
    # wrapper path; only the A21 exponential arithmetic value remains open.
    bias = _bias_root(core, bias_phi=bias_phi)
    return Roots(angular, ou, qaxis, bias)


def imu_step(state: WORD.State, *, witness: SOURCE.StepWitness,
             segment, raw: SENSOR.RawImuSample, packet_id: str,
             ou_alpha, bias_phi=None, angular_full=None, angular_half=None,
             qaxis_marginal_psd, qaxis_final_psd, machine_epsilon,
             **dynamic):
    """Execute one source-owning IMU edge with source-bound prediction roots."""
    forbidden = {'angular','ou','qaxis','bias'} & set(dynamic)
    if forbidden:
        raise TypeError('source-bound prediction roots cannot be overridden '+repr(sorted(forbidden)))
    if not isinstance(state, WORD.State):
        raise TypeError('source-owning Live state required')
    if not isinstance(witness, SOURCE.StepWitness):
        raise TypeError('checked source step witness required')
    # Ordinal ownership is a theorem-product invariant. Reject a skipped or
    # duplicated transition before constructing any lower-level physical cell,
    # so a malformed caller cannot select a different validation path first.
    if witness.ordinal != state.source.next_ordinal:
        raise ValueError('prediction roots require exactly the next source ordinal')
    physical = SOURCE.QualifiedPhysicalSegment(state.source.root, witness, segment)
    roots = build(state, physical, raw, ou_alpha=ou_alpha,bias_phi=bias_phi,
                  angular_full=angular_full, angular_half=angular_half,
                  qaxis_marginal_psd=qaxis_marginal_psd,
                  qaxis_final_psd=qaxis_final_psd,
                  machine_epsilon=machine_epsilon)
    return WORD.imu_step(state, witness=witness, segment=segment, raw=raw,
                         packet_id=packet_id, angular=roots.angular,
                         ou=roots.ou, qaxis=roots.qaxis, bias=roots.bias, **dynamic)


def readiness():
    lower = WORD.readiness()
    return {
      'source_owned_next_transition_required': True,
      'attitude_omega_reconstructed_from_same_raw_packet_and_current_bias_error': True,
      'OU_h_from_same_source_transition': True,
      'OU_tau_from_same_carried_active_TuneState': True,
      'Qaxis_Sigma_aw_from_same_carried_active_TuneState': True,
      'shipping_independent_Qaxis_branch_fixed_by_active_parameter_contract': True,
      'prediction_angular_OU_Qaxis_roots_bound_to_same_source_continuation': True,
      'shipping_BA_tau_Q_defaults_bound_at_prediction_entry': True,
      'shipping_BA_hold_active_branch_derived_from_current_MEKF_mode': True,
      'held_H18_BA_phi_exactly_one_without_exp_witness': True,
      'active_A21_BA_phi_only_remaining_exp_arithmetic_witness': True,
      'caller_cannot_override_prediction_roots': True,
      'sample_zero_origin_bridge_preserved': lower['sample_zero_startup_to_checked_outer_endpoint_bridge_closed'],
      'OU_exp_arithmetic_deployment_closed': False,
      'BA_exp_arithmetic_deployment_closed': False,
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
