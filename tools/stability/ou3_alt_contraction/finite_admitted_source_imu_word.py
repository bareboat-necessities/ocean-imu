"""Strongest current theorem-facing ALT IMU edge.

This composes, on one event, the pieces that were previously parallel:

* one quantified admitted COMPLETE-BRMM history;
* one quantified admitted BIAS0/1/2 history;
* the same kth 5 ms restriction of both histories;
* one source-qualified raw IMU packet;
* prediction/model roots reconstructed by ``finite_source_bound_prediction_word``;
* the exact post-guard IMU/thermal forcing consumed by the executed prefix.

Thus the event direction is

  admitted histories -> kth restriction -> raw packet -> shipping coefficients
                     -> shipping IMU successor + retained forcing.

Dynamic libm/Eigen/roundoff witnesses are still conditional. Sensor residuals
are retained as ISS supply but are not yet given a theorem-level admissibility
class. This module does not authorize storage or promote an ALT gate.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as ADLIVE
from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as ABRMM
from tools.stability.ou3_alt_contraction import finite_admitted_bias_history as ABIAS
from tools.stability.ou3_alt_contraction import finite_source_bound_prediction_word as PRED
from tools.stability.ou3_alt_contraction import finite_source_bound_imu_forcing as FORCE
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR


@dataclass(frozen=True)
class Result:
    state: ADLIVE.State
    event: object
    forcing: FORCE.ImuForcing


def _forcing(pred_result, temperature_c):
    prefix=FORCE._executed_prefix(pred_result)
    guarded=prefix.guarded
    if not isinstance(guarded,SENSOR.GuardedImuSample):
        raise TypeError('executed admitted source word lost guarded same-packet sample')
    conditioning=PRED._accel_conditioning(temperature_c)
    obs=SENSOR.finite_accel_core_observation(guarded,conditioning)
    thermal=tuple(conditioning.k_a_hat_internal[i]*conditioning.temperature_delta
                  for i in range(3))
    return FORCE.ImuForcing(
        guarded.raw.gyro_residual_internal,
        guarded.effective_accel_residual_internal,
        thermal,obs.nu_acc,conditioning.temperature_delta)


def imu_step(state:ADLIVE.State, *, restricted:ABRMM.RestrictedSegment,
             bias_restricted:ABIAS.RestrictedBiasStep,
             witness:SOURCE.StepWitness, raw:SENSOR.RawImuSample, packet_id:str,
             temperature_c, **dynamic):
    if not isinstance(state,ADLIVE.State):
        raise TypeError('joint admitted BRMM/BIAS Live state required')
    if not isinstance(restricted,ABRMM.RestrictedSegment) or not isinstance(bias_restricted,ABIAS.RestrictedBiasStep):
        raise TypeError('BRMM and BIAS kth restrictions required')
    if restricted.history != state.admitted_history:
        raise ValueError('IMU restriction detached from carried admitted BRMM history')
    if bias_restricted.history != state.bias_history:
        raise ValueError('IMU restriction detached from carried admitted BIAS history')
    if restricted.ordinal != bias_restricted.ordinal or restricted.ordinal != witness.ordinal:
        raise ValueError('BRMM/BIAS/source restriction ordinals differ')
    if restricted.segment != bias_restricted.segment:
        raise ValueError('BRMM and BIAS restrictions must share exact physical segment')

    # Consume the universal restriction theorem before constructing any model
    # coefficient root or executing the shipping event.
    ABRMM.qualify_step(state.live_word.source.root,restricted,witness)

    pred=PRED.imu_step(state.live_word,witness=witness,segment=restricted.segment,
                       raw=raw,packet_id=packet_id,temperature_c=temperature_c,
                       **dynamic)
    forcing=_forcing(pred,temperature_c)
    nxt=ADLIVE.State(pred.state,state.admitted_history,state.bias_history)
    return Result(nxt,pred.event,forcing)


def readiness():
    source=ADLIVE.readiness(); pred=PRED.readiness(); force=FORCE.readiness()
    return {
      'admitted_BRMM_and_BIAS_histories_carried_jointly': bool(
          source['universal_admitted_COMPLETE_BRMM_history_carried_in_Live_product']
          and source['admitted_BIAS_history_carried_in_same_Live_product']),
      'same_kth_BRMM_BIAS_restriction_consumed_before_shipping_execution':True,
      'same_restriction_drives_source_bound_prediction_coefficients':True,
      'prediction_roots_caller_override_forbidden':pred['caller_cannot_override_prediction_or_accel_model_roots'],
      'same_executed_packet_forcing_retained':force['forcing_supply_derived_from_executed_source_owned_IMU_event'],
      'physical_moments_sensor_forcing_and_model_roots_share_one_event':True,
      'sensor_residual_admissibility_attached':False,
      'startup_sample_zero_admitted_history_equality_closed':False,
      'all_runtime_arithmetic_deployment_correspondence_closed':False,
      'complete_600_step_shipping_word_composed':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
