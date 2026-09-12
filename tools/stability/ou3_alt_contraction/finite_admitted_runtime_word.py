"""Master admitted-history/runtime edge for the OU-III ALT proof.

This layer closes a structural bypass between two stronger proof products:

* ``finite_admitted_source_live_word`` carries one quantified admitted
  COMPLETE-BRMM history and one quantified admitted BIAS0/1/2 history, requiring
  each IMU event to consume the same kth physical restriction;
* ``finite_source_bound_imu_forcing`` binds the shipping prediction/model roots
  to that source-owned step and derives the actual same-event IMU ISS supply;
* ``finite_source_bound_mag_forcing`` derives the correlated magnetic ISS supply
  from the same source-owned magnetic event.

An IMU edge through this module therefore cannot choose between "admitted
source" and "strong runtime" paths: it must satisfy both.  This still does NOT
close startup sample-zero equality, deployment exp/expm1/trig/Eigen/LDLT
correspondence, finite clocks/counters, or a complete 600-step hybrid word.
Storage remains forbidden.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as ADMITTED
from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as BRMM
from tools.stability.ou3_alt_contraction import finite_admitted_bias_history as BIAS
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_source_bound_imu_forcing as IMU
from tools.stability.ou3_alt_contraction import finite_source_bound_mag_forcing as MAG


@dataclass(frozen=True)
class ImuResult:
    state: ADMITTED.State
    event: object
    forcing: IMU.ImuForcing
    runtime_word: object


@dataclass(frozen=True)
class MagResult:
    state: ADMITTED.State
    event: object
    forcing: MAG.MagForcing | None
    runtime_word: object


def _qualify_step(state:ADMITTED.State, *, restricted:BRMM.RestrictedSegment,
                  bias_restricted:BIAS.RestrictedBiasStep,
                  witness:SOURCE.StepWitness):
    """Require one kth restriction of both carried admitted histories."""
    if not isinstance(state,ADMITTED.State):
        raise TypeError('admitted BRMM/BIAS Live state required')
    if not isinstance(restricted,BRMM.RestrictedSegment):
        raise TypeError('kth admitted COMPLETE-BRMM restriction required')
    if not isinstance(bias_restricted,BIAS.RestrictedBiasStep):
        raise TypeError('kth admitted BIAS restriction required')
    if restricted.history != state.admitted_history:
        raise ValueError('IMU restriction detached from carried admitted history')
    if bias_restricted.history != state.bias_history:
        raise ValueError('BIAS restriction detached from carried admitted BIAS history')
    if restricted.ordinal != bias_restricted.ordinal or restricted.ordinal != witness.ordinal:
        raise ValueError('BRMM/BIAS/source restriction ordinals differ')
    if restricted.segment != bias_restricted.segment:
        raise ValueError('BRMM and BIAS restrictions must describe the same physical segment')
    # This applies the universal COMPLETE-BRMM restriction theorem and finite
    # necessary constraints.  RestrictedBiasStep has already checked the
    # selected admitted BIAS recurrence on the exact same segment.
    return BRMM.qualify_step(state.live_word.source.root,restricted,witness)


def imu_step(state:ADMITTED.State, *, restricted:BRMM.RestrictedSegment,
             bias_restricted:BIAS.RestrictedBiasStep,
             witness:SOURCE.StepWitness, raw:SENSOR.RawImuSample, packet_id:str,
             **runtime):
    """Execute one admitted-source IMU step through the strongest runtime edge."""
    qualified=_qualify_step(state,restricted=restricted,
                            bias_restricted=bias_restricted,witness=witness)
    if qualified.segment != restricted.segment:
        raise AssertionError('qualified admitted segment identity lost')
    out=IMU.imu_step(state.live_word,witness=witness,segment=restricted.segment,
                     raw=raw,packet_id=packet_id,**runtime)
    nxt=ADMITTED.State(out.state,state.admitted_history,state.bias_history)
    return ImuResult(nxt,out.word.event,out.forcing,out)


def mag_step(state:ADMITTED.State, **kwargs):
    """Execute magnetic event while preserving both admitted histories."""
    if not isinstance(state,ADMITTED.State):
        raise TypeError('admitted BRMM/BIAS Live state required')
    out=MAG.mag_step(state.live_word,**kwargs)
    nxt=ADMITTED.State(out.state,state.admitted_history,state.bias_history)
    return MagResult(nxt,out.word.event,out.forcing,out)


def set_hold(state:ADMITTED.State, *, hold):
    """External hold event does not consume either admitted physical history."""
    if not isinstance(state,ADMITTED.State):
        raise TypeError('admitted BRMM/BIAS Live state required')
    out=ADMITTED.set_hold(state,hold=hold)
    return out


def readiness():
    admitted=ADMITTED.readiness(); imu=IMU.readiness(); mag=MAG.readiness()
    return {
      'admitted_COMPLETE_BRMM_history_and_strong_runtime_joined':True,
      'admitted_BIAS_history_and_strong_runtime_joined':True,
      'each_IMU_event_requires_same_admitted_BRMM_BIAS_and_source_ordinal':True,
      'same_admitted_segment_drives_source_bound_prediction_and_measurement':True,
      'prediction_model_roots_cannot_bypass_admitted_history_edge':True,
      'same_event_IMU_ISS_supply_attached_to_admitted_step':imu['forcing_supply_derived_from_executed_source_owned_IMU_event'],
      'same_event_magnetic_ISS_supply_attached_to_admitted_product':mag['magnetic_effective_residual_derived_from_same_qualified_event'],
      'BIAS_generating_history_attached':admitted['BIAS_generating_history_attached'],
      'finite_tokens_used_as_source_membership_oracle':False,
      'sensor_or_temperature_amplitude_bound_invented':False,
      'startup_sample_zero_equal_to_admitted_history_restriction_proved':False,
      'bounded_input_history_qualified':False,
      'one_radian_attitude_guard_closed_for_every_admitted_prefix':False,
      'deployment_exp_expm1_trig_Eigen_LDLT_closed':False,
      'deployment_roundoff_supply_attached':False,
      'complete_600_step_shipping_word_composed_from_restrictions':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
