"""Source-qualified dual-clock magnetic edge for the ALT master word.

This is the theorem-facing Live magnetometer successor after the outer/inner
clock split. It preserves the admitted/source-checked physical endpoint, uses
shipping's exact binary32 outer wrapper clock for magnetic gating/refinement/
continuous calibration, retains physical inner-MEKF time for the actual Kalman
measurement call, advances the separate deployment call-schedule clock in
physical time, and exposes the same correlated magnetic ISS forcing.

Startup history still has to be rooted through the dual-clock startup composer
before the complete master can claim an end-to-end dual-clock magnetic word.
"""
from __future__ import annotations
from dataclasses import dataclass, replace

from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as WORD
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_live_interleave as LIVE
from tools.stability.ou3_alt_contraction import finite_live_magnetic_dual_clock as DUAL
from tools.stability.ou3_alt_contraction import finite_mag_call_schedule as SCHEDULE
from tools.stability.ou3_alt_contraction import finite_source_bound_mag_forcing as FORCE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR


@dataclass(frozen=True)
class Result:
    word: WORD.Result
    forcing: FORCE.MagForcing | None

    @property
    def state(self): return self.word.state


def _endpoint(state:WORD.State):
    if state.source.steps:
        return SOURCE.endpoint(state.source.steps[-1],'after')
    return SOURCE.origin_endpoint(state.source.root,state.live.live.live.mekf.reference)


def _assert_endpoint(state:WORD.State, endpoint):
    endpoint_types=(SOURCE.QualifiedPhysicalEndpoint,SOURCE.QualifiedPhysicalOrigin)
    if not isinstance(endpoint,endpoint_types):
        raise TypeError('source-checked magnetic physical endpoint required')
    inter=state.live; core=inter.live.live.mekf; root=endpoint.root
    if root.history_id != core.reference.history_id or root.history_id != inter.magnetic.memory.history_id:
        raise ValueError('qualified magnetic endpoint detached from persistent physical history')
    if root.live_origin != core.reference.live_origin:
        raise ValueError('qualified magnetic endpoint restarted Live origin')
    if endpoint.endpoint != core.reference:
        raise ValueError('asynchronous magnetic call not attached to current admitted endpoint')
    return True


def _forcing(event):
    if event.qualification is None:
        if event.effective_residual is not None:
            raise AssertionError('gated magnetic call unexpectedly produced effective residual')
        return None
    source=event.qualification.sample
    active=event.state.active
    memory=event.state.memory
    ref=event.filter.reference
    field_difference=tuple(source.model.world_field[i]-active.model.world_reference[i]
                           for i in range(3))
    rotated=tuple(SENSOR.q_rotate(ref.q_world_to_body,field_difference))
    correction=tuple(-x for x in memory.applied.total_bias)
    return FORCE.MagForcing(rotated,source.model.hard_iron_body,correction,
                            source.residual_body,event.effective_residual)


def mag_step(state:WORD.State, **kwargs):
    if not isinstance(state,WORD.State):
        raise TypeError('source-owning Live word required')
    endpoint=_endpoint(state); _assert_endpoint(state,endpoint)
    inter=state.live
    event=DUAL.live_call(inter.magnetic,inter.live.live.mekf,
                         inter.live.live.tuner.vertical,**kwargs)

    # MAG-CALL-SCHEDULE-v1 is a deployment cadence condition in physical time;
    # it is deliberately distinct from the binary32 outer-wrapper clock used
    # inside DUAL.live_call.
    clock=inter.clock
    if event.measurement is not None and event.measurement.wrapper_attempted:
        clock=SCHEDULE.record_call(clock,inter.schedule,time=event.filter.reference.time)

    live=replace(inter.live,live=replace(inter.live.live,mekf=event.filter))
    nxt_inter=LIVE.State(live,event.state,clock,inter.schedule)
    nxt=WORD.State(nxt_inter,state.source,state.sensor_root,state.bias_history_id,state.runtime)
    wrapped=WORD.Result(nxt,LIVE.Result(nxt_inter,event))
    return Result(wrapped,_forcing(event))


def readiness():
    d=DUAL.readiness(); f=FORCE.readiness()
    return {
      'source_checked_async_endpoint_required_before_dual_clock_magnetic_edge':True,
      'outer_binary32_and_inner_physical_magnetic_clocks_composed_on_Live_edge':d['dual_clock_magnetic_word_composed'],
      'physical_MAG_CALL_SCHEDULE_clock_kept_distinct_from_wrapper_binary32_clock':True,
      'same_event_correlated_magnetic_ISS_forcing_retained':f['magnetic_effective_residual_derived_from_same_qualified_event'],
      'source_and_Live_origin_persist_across_dual_clock_magnetic_edge':True,
      'startup_dual_clock_history_feeds_this_master_edge':False,
      'binary32_exp_solver_roundoff_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
