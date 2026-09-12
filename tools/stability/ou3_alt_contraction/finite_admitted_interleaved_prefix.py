"""Finite interleaving induction over the strongest admitted-source ALT edges.

This is the structural prefix composer needed before a source-uniform 600-step
word can exist. It folds the strongest current theorem-facing events without
restarting state:

* IMU: admitted BRMM+BIAS kth restriction -> source-bound coefficients ->
  shipping successor + retained IMU forcing;
* MAG: same current admitted endpoint (or explicit t_L origin at sample zero) ->
  shipping magnetic successor + retained correlated magnetic forcing;
* HOLD: literal external hold-control event, preserving both admitted histories.

The composer proves by construction that only IMU events advance the source
ordinal and that every successor is the predecessor of the next event. It does
NOT quantify/close all deployment arithmetic witnesses or sensor admissibility,
so a 600-transition trace is representable but not yet a certified universal
shipping word. No storage search is authorized.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as LIVE
from tools.stability.ou3_alt_contraction import finite_admitted_source_imu_word as IMU
from tools.stability.ou3_alt_contraction import finite_admitted_source_mag_word as MAG
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE


@dataclass(frozen=True)
class EventRecord:
    kind: str
    source_steps_before: int
    source_steps_after: int
    def __post_init__(self):
        if self.kind not in ('imu','mag','hold'):
            raise ValueError('unknown interleaved event kind')
        expected=self.source_steps_before+(1 if self.kind=='imu' else 0)
        if self.source_steps_after != expected:
            raise ValueError('event/source ordinal accounting mismatch')


@dataclass(frozen=True)
class State:
    live: LIVE.State
    events: tuple[EventRecord,...]=()
    def __post_init__(self):
        if not isinstance(self.live,LIVE.State):
            raise TypeError('joint admitted Live state required')
        object.__setattr__(self,'events',tuple(self.events))
        n=len(self.live.live_word.source.steps)
        if n>SOURCE.TRANSITIONS:
            raise ValueError('more than canonical 600 source transitions')
        if self.events:
            if self.events[-1].source_steps_after != n:
                raise ValueError('event ledger detached from carried source continuation')
            if sum(e.kind=='imu' for e in self.events) != n:
                raise ValueError('IMU event count detached from source transition count')

    @property
    def imu_steps(self): return len(self.live.live_word.source.steps)
    @property
    def complete_source_horizon(self): return self.imu_steps==SOURCE.TRANSITIONS


def begin(live:LIVE.State): return State(live,())


def imu_step(state:State, **kwargs):
    if not isinstance(state,State): raise TypeError('interleaved admitted prefix required')
    before=state.imu_steps
    out=IMU.imu_step(state.live,**kwargs)
    record=EventRecord('imu',before,len(out.state.live_word.source.steps))
    return State(out.state,state.events+(record,)),out


def mag_step(state:State, **kwargs):
    if not isinstance(state,State): raise TypeError('interleaved admitted prefix required')
    before=state.imu_steps
    out=MAG.mag_step(state.live,**kwargs)
    record=EventRecord('mag',before,len(out.state.live_word.source.steps))
    return State(out.state,state.events+(record,)),out


def set_hold(state:State, *, hold):
    if not isinstance(state,State): raise TypeError('interleaved admitted prefix required')
    before=state.imu_steps
    out=LIVE.set_hold(state.live,hold=hold)
    record=EventRecord('hold',before,len(out.state.live_word.source.steps))
    return State(out.state,state.events+(record,)),out


def readiness():
    i=IMU.readiness(); m=MAG.readiness()
    return {
      'finite_successor_induction_over_IMU_MAG_HOLD_closed':True,
      'only_IMU_advances_admitted_source_ordinal':True,
      'IMU_successor_uses_strong_admitted_source_edge':i['physical_moments_sensor_forcing_and_model_roots_share_one_event'],
      'MAG_successor_retains_same_event_correlated_forcing':m['same_event_correlated_magnetic_forcing_retained'],
      'sample_zero_MAG_requires_explicit_admitted_origin':m['fresh_sample_zero_magnetic_edge_requires_explicit_admitted_tL_origin'],
      'both_admitted_histories_persist_through_all_three_event_types':True,
      'canonical_600_transition_source_horizon_representable':True,
      'all_event_arithmetic_witnesses_source_uniformly_qualified':False,
      'sensor_disturbance_admissibility_attached':False,
      'startup_reachability_to_admitted_fresh_state_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
