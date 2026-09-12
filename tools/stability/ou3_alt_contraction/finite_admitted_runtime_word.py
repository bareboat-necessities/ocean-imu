"""Master admitted-history/runtime edge for the OU-III ALT proof.

This theorem-facing product requires, simultaneously:
* one quantified admitted COMPLETE-BRMM history;
* one quantified admitted BIAS0/1/2 history;
* the canonical sample-zero restriction of that admitted physical history,
  equal to the fresh startup joint24 physical Reference;
* source-bound shipping prediction/model roots;
* same-event IMU and magnetic ISS forcing coordinates.

Live magnetic successors now use the source-qualified dual-clock edge: exact
binary32 outer-wrapper time for outer magnetic control/calibration and physical
inner-MEKF time for the Kalman call. The startup history entering this master is
still produced by the older startup interleave, so end-to-end dual-clock startup
ancestry remains fail-closed. Storage remains forbidden.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as ADMITTED
from tools.stability.ou3_alt_contraction import finite_admitted_startup_origin as ORIGIN
from tools.stability.ou3_alt_contraction import finite_admitted_brmm_restriction as BRMM
from tools.stability.ou3_alt_contraction import finite_admitted_bias_history as BIAS
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_source_bound_imu_forcing as IMU
from tools.stability.ou3_alt_contraction import finite_source_bound_mag_dual_clock as MAG
from tools.stability.ou3_alt_contraction import finite_source_bound_attitude_trig as TRIG
from tools.stability.ou3_alt_contraction import finite_mag_call_schedule as MAGSCHEDULE
from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as CLOCK
from tools.stability.ou3_alt_contraction import finite_magnetic_wrapper_clock as MAGCLOCK


@dataclass(frozen=True)
class State:
    admitted: ADMITTED.State
    origin: ORIGIN.RestrictedOrigin
    def __post_init__(self):
        if not isinstance(self.admitted,ADMITTED.State) or not isinstance(self.origin,ORIGIN.RestrictedOrigin):
            raise TypeError('admitted Live product and canonical admitted origin required')
        if self.origin.history != self.admitted.admitted_history:
            raise ValueError('admitted runtime origin detached from carried admitted history')
        root=self.admitted.live_word.source.root
        if root.history_id != self.origin.history.history_id or root.live_origin != self.origin.endpoint.live_origin:
            raise ValueError('finite source root detached from admitted runtime origin')
        BRMM.qualify_origin(root,self.origin)
        if self.admitted.live_word.source.steps:
            if self.admitted.live_word.source.steps[0].segment.before != self.origin.endpoint:
                raise ValueError('runtime source chain did not start at admitted sample zero')
        else:
            ref=self.admitted.live_word.live.live.live.mekf.reference
            if ref != self.origin.endpoint:
                raise ValueError('fresh runtime Reference is not admitted sample zero')


@dataclass(frozen=True)
class ImuResult:
    state: State
    event: object
    forcing: IMU.ImuForcing
    runtime_word: object


@dataclass(frozen=True)
class MagResult:
    state: State
    event: object
    forcing: object
    runtime_word: object


def bind_startup(admitted:ADMITTED.State, origin:ORIGIN.RestrictedOrigin):
    ORIGIN.bind(admitted,origin)
    return State(admitted,origin)


def _qualify_step(state:State, *, restricted:BRMM.RestrictedSegment,
                  bias_restricted:BIAS.RestrictedBiasStep,
                  witness:SOURCE.StepWitness):
    if not isinstance(state,State): raise TypeError('origin-bound admitted runtime state required')
    admitted=state.admitted
    if not isinstance(restricted,BRMM.RestrictedSegment): raise TypeError('kth admitted COMPLETE-BRMM restriction required')
    if not isinstance(bias_restricted,BIAS.RestrictedBiasStep): raise TypeError('kth admitted BIAS restriction required')
    if restricted.history != admitted.admitted_history:
        raise ValueError('IMU restriction detached from carried admitted history')
    if bias_restricted.history != admitted.bias_history:
        raise ValueError('BIAS restriction detached from carried admitted BIAS history')
    if restricted.ordinal != bias_restricted.ordinal or restricted.ordinal != witness.ordinal:
        raise ValueError('BRMM/BIAS/source restriction ordinals differ')
    if restricted.segment != bias_restricted.segment:
        raise ValueError('BRMM and BIAS restrictions must describe the same physical segment')
    if restricted.ordinal == 1 and restricted.segment.before != state.origin.endpoint:
        raise ValueError('first admitted restriction does not start at admitted sample zero')
    return BRMM.qualify_step(admitted.live_word.source.root,restricted,witness)


def imu_step(state:State, *, restricted:BRMM.RestrictedSegment,
             bias_restricted:BIAS.RestrictedBiasStep,
             witness:SOURCE.StepWitness, raw:SENSOR.RawImuSample, packet_id:str,
             **runtime):
    qualified=_qualify_step(state,restricted=restricted,bias_restricted=bias_restricted,witness=witness)
    if qualified.segment != restricted.segment: raise AssertionError('qualified admitted segment identity lost')
    admitted=state.admitted
    out=IMU.imu_step(admitted.live_word,witness=witness,segment=restricted.segment,
                     raw=raw,packet_id=packet_id,**runtime)
    nxt=ADMITTED.State(out.state,admitted.admitted_history,admitted.bias_history)
    return ImuResult(State(nxt,state.origin),out.word.event,out.forcing,out)


def mag_step(state:State, **kwargs):
    if not isinstance(state,State): raise TypeError('origin-bound admitted runtime state required')
    admitted=state.admitted
    out=MAG.mag_step(admitted.live_word,**kwargs)
    nxt=ADMITTED.State(out.state,admitted.admitted_history,admitted.bias_history)
    return MagResult(State(nxt,state.origin),out.word.event,out.forcing,out)


def set_hold(state:State, *, hold):
    if not isinstance(state,State): raise TypeError('origin-bound admitted runtime state required')
    out=ADMITTED.set_hold(state.admitted,hold=hold)
    return State(out.state,state.origin)


def readiness():
    admitted=ADMITTED.readiness(); origin=ORIGIN.readiness(); imu=IMU.readiness(); mag=MAG.readiness()
    trig=TRIG.readiness(); counter=MAGSCHEDULE.counter_lifetime(); clock=CLOCK.readiness(); magclock=MAGCLOCK.readiness()
    return {
      'admitted_COMPLETE_BRMM_history_and_strong_runtime_joined':True,
      'admitted_BIAS_history_and_strong_runtime_joined':True,
      'canonical_admitted_sample_zero_origin_persists_in_runtime_product':True,
      'first_restriction_forced_to_start_at_admitted_sample_zero':True,
      'startup_sample_zero_equal_to_admitted_history_restriction_proved':origin['startup_sample_zero_equal_to_admitted_history_restriction_proved'],
      'each_IMU_event_requires_same_admitted_BRMM_BIAS_and_source_ordinal':True,
      'same_admitted_segment_drives_source_bound_prediction_and_measurement':True,
      'prediction_model_roots_cannot_bypass_admitted_history_edge':True,
      'same_event_IMU_ISS_supply_attached_to_admitted_step':imu['forcing_supply_derived_from_executed_source_owned_IMU_event'],
      'same_event_magnetic_ISS_supply_attached_to_admitted_product':mag['same_event_correlated_magnetic_ISS_forcing_retained'],
      'BIAS_generating_history_attached':admitted['BIAS_generating_history_attached'],
      'finite_tokens_used_as_source_membership_oracle':False,
      'sensor_or_temperature_amplitude_bound_invented':False,
      'startup_reachability_for_every_admitted_history_proved':False,
      'bounded_input_history_qualified':False,
      'global_exact_real_attitude_trig_ancestry_closed':trig['global_finite_rational_angle_enclosure_available'],
      'one_radian_attitude_guard_closed_for_every_admitted_prefix':not trig['one_radian_local_angle_guard_required'],
      'deployment_exp_expm1_trig_Eigen_LDLT_closed':False,
      'canonical_5ms_wrapper_clock_prefix_binary32_closed':clock['canonical_5ms_wrapper_clock_prefix_binary32_closed'],
      'wrapper_clock_arbitrary_dt_closed':clock['arbitrary_dt_wrapper_clock_closed'],
      'wrapper_clock_indefinite_lifetime_closed':clock['indefinite_wrapper_clock_lifetime_closed'],
      'outer_magnetic_wrapper_clock_operands_materialized':magclock['exact_binary32_wrapper_timestamp_derived_from_canonical_physical_endpoint'],
      'dual_clock_live_magnetic_edge_composed':mag['outer_binary32_and_inner_physical_magnetic_clocks_composed_on_Live_edge'],
      'dual_clock_startup_history_feeds_master':mag['startup_dual_clock_history_feeds_this_master_edge'],
      'dual_clock_magnetic_word_composed':bool(
          mag['outer_binary32_and_inner_physical_magnetic_clocks_composed_on_Live_edge']
          and mag['startup_dual_clock_history_feeds_this_master_edge']),
      'deployment_roundoff_supply_attached':False,
      'mag_schedule_supplies_uniform_call_count_upper':counter.uniform_call_count_upper is not None,
      'shipping_signed_mag_counter_lifetime_closed':counter.no_signed_overflow_proved,
      'complete_600_step_shipping_word_composed_from_restrictions':False,
      'successive_600_step_words_tiled_without_restarting_Live_origin':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
