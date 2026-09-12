"""Source-owning finite Live product for the ALT proof.

The lower ``finite_live_interleave`` layer composes exact shipping event
relations but accepts a source-qualified segment/endpoint supplied by its caller.
That is still too weak for the universal 600-step theorem: a caller could present
individually qualified transitions without proving that they are successive
members of one carried O^601_BRMM / BIAS history.

This module makes the source continuation part of the product state itself.
Every theorem-facing IMU event appends exactly the next source transition before
executing the shipping event.  Thus source-cell ancestry, primitive ancestry,
physical endpoint equality, BIAS-family token, sensor-residual histories, full
MEKF/frontend/calibration state and the next physical endpoint advance together.
Magnetic and external-hold events do not advance the physical source and are
forced to remain at the current carried endpoint.

The initial Live endpoint is inherited from the already-qualified startup
handoff.  This module does not newly prove startup membership in O^601_BRMM; that
startup-to-source-root bridge remains explicit and fail-closed.  Runtime
transcendentals/solver roundoff and quantitative sensor ISS bounds also remain
open.  Therefore the finite-storage master and every ALT theorem gate remain
false.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_live_interleave as LIVE
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR


@dataclass(frozen=True)
class State:
    live: LIVE.State
    source: SOURCE.Continuation
    sensor_root: SOURCE.SensorDisturbanceRoot
    def __post_init__(self):
        if not isinstance(self.live,LIVE.State) or not isinstance(self.source,SOURCE.Continuation):
            raise TypeError('interleaved Live state and persistent source continuation required')
        if not isinstance(self.sensor_root,SOURCE.SensorDisturbanceRoot):
            raise TypeError('persistent sensor disturbance root required')
        if self.sensor_root.source_root != self.source.root:
            raise ValueError('sensor histories detached from carried source continuation')
        core=self.live.live.live.mekf
        root=self.source.root
        if core.reference.history_id != root.history_id:
            raise ValueError('carried source history detached from current MEKF physical reference')
        if core.reference.live_origin != root.live_origin:
            raise ValueError('carried source root restarted one-time Live origin')
        if self.source.steps:
            if self.source.steps[-1].segment.after != core.reference:
                raise ValueError('Live product state detached from last carried source endpoint')
        else:
            # Sample zero comes from the startup handoff; membership of that
            # handoff in the universal outer source remains a separate bridge.
            if core.reference.time != root.live_origin:
                raise ValueError('empty source continuation must sit at fresh Live origin')


@dataclass(frozen=True)
class Result:
    state: State
    event: object


def from_live(live:LIVE.State, root:SOURCE.SourceRoot,
              sensor_root:SOURCE.SensorDisturbanceRoot):
    if not isinstance(live,LIVE.State) or not isinstance(root,SOURCE.SourceRoot):
        raise TypeError('Live product and source root required')
    if sensor_root.source_root != root:
        raise ValueError('sensor histories detached from source root')
    return State(live,SOURCE.begin(root),sensor_root)


def imu_step(state:State, *, witness:SOURCE.StepWitness,
             segment:PHYS.PhysicalSegment, raw:SENSOR.RawImuSample,
             packet_id:str, **kwargs):
    """Append the next admitted physical transition and execute one IMU event."""
    if not isinstance(state,State): raise TypeError('source-owning Live state required')
    if witness.ordinal != state.source.next_ordinal:
        raise ValueError('IMU event must consume exactly the next source ordinal')
    nxt_source=SOURCE.append(state.source,witness=witness,segment=segment)
    qualified=nxt_source.steps[-1]
    packet=SOURCE.qualify_raw_imu(qualified,state.sensor_root,raw,packet_id)
    event=LIVE.imu_step_source_qualified(state.live,packet,**kwargs)
    return Result(State(event.state,nxt_source,state.sensor_root),event)


def mag_step(state:State, **kwargs):
    """Execute an async magnetic event at the current carried source endpoint.

    Before the first post-handoff IMU transition there is not yet a transition-
    derived source endpoint in this module.  That sample-zero/startup-source
    bridge remains explicitly open instead of being manufactured from equality
    of time/history labels.
    """
    if not isinstance(state,State): raise TypeError('source-owning Live state required')
    if not state.source.steps:
        raise NotImplementedError('sample-zero startup-to-COMPLETE-BRMM magnetic endpoint bridge remains open')
    endpoint=SOURCE.endpoint(state.source.steps[-1],'after')
    event=LIVE.mag_step_source_qualified(state.live,endpoint,**kwargs)
    return Result(State(event.state,state.source,state.sensor_root),event)


def set_hold(state:State, *, hold):
    if not isinstance(state,State): raise TypeError('source-owning Live state required')
    event=LIVE.set_hold(state.live,hold=hold)
    return Result(State(event.state,state.source,state.sensor_root),event)


def readiness():
    lower=LIVE.readiness(); src=SOURCE.readiness()
    return {
      'source_continuation_is_part_of_theorem_product_state':True,
      'every_IMU_event_appends_exactly_next_source_ordinal':True,
      'source_cell_and_physical_primitive_chains_advance_with_filter_state':True,
      'raw_IMU_sensor_histories_advance_with_same_source_transition':True,
      'magnetic_and_hold_events_preserve_current_source_endpoint':True,
      'qualified_post_first_IMU_magnetic_entry_available':lower['source_qualified_async_magnetic_endpoint_entry_available'],
      'one_bias_parameter_token_carried_over_word':src['one_bias_family_parameter_token_over_word_required'],
      'sample_zero_startup_to_COMPLETE_BRMM_endpoint_bridge_closed':False,
      'quantitative_sensor_residual_ISS_envelope_attached':False,
      'finite_estimator_coefficients_bound_to_same_source_continuation':False,
      'all_runtime_arithmetic_deployment_correspondence_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
