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
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_source_bound_live_word as WORD
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_live_interleave as LIVE
from tools.stability.ou3_alt_contraction import finite_live_magnetic_dual_clock as DUAL
from tools.stability.ou3_alt_contraction import finite_source_bound_mag_forcing as FORCE


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


def mag_step(state:WORD.State, **kwargs):
    if not isinstance(state,WORD.State):
        raise TypeError('source-owning Live word required')
    endpoint=_endpoint(state); _assert_endpoint(state,endpoint)
    # Consume the shared composer, including ungauged waiting, actual north
    # acquisition, service-clock initialization and saturated counter projection.
    out=LIVE.mag_step(state.live,**kwargs)
    nxt=WORD.State(out.state,state.source,state.sensor_root,state.bias_history_id,state.runtime)
    return Result(WORD.Result(nxt,out),FORCE.from_event(out.event))


def readiness():
    d=DUAL.readiness(); f=FORCE.readiness()
    return {
      'source_checked_async_endpoint_required_before_dual_clock_magnetic_edge':True,
      'outer_binary32_and_inner_physical_magnetic_clocks_composed_on_Live_edge':d['dual_clock_magnetic_word_composed'],
      'physical_MAG_CALL_SCHEDULE_clock_kept_distinct_from_wrapper_binary32_clock':True,
      'same_event_correlated_magnetic_ISS_forcing_retained':f['magnetic_effective_residual_derived_from_same_qualified_event'],
      'source_and_Live_origin_persist_across_dual_clock_magnetic_edge':True,
      'ungauged_north_and_saturated_counter_use_shared_event_composer':True,
      'startup_dual_clock_history_feeds_this_master_edge':False,
      'binary32_exp_solver_roundoff_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
