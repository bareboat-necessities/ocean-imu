"""Interleaved admitted ALT prefix carrying one bounded IMU ISS history.

This is the strongest current structural prefix.  It reuses the admitted
BRMM/BIAS/origin IMU-MAG-HOLD induction and adds one persistent bounded IMU
forcing history.  IMU consumes the exact kth forcing through
``finite_admitted_disturbance_imu_word``; asynchronous MAG and HOLD preserve the
same disturbance-history object and consume no IMU disturbance ordinal.

This remains finite structural composition.  Deployment arithmetic, counter
lifetime, startup reachability and a storage inequality are still open.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_interleaved_prefix as BASE
from tools.stability.ou3_alt_contraction import finite_admitted_disturbance_imu_word as IMU
from tools.stability.ou3_alt_contraction import finite_admitted_imu_disturbance as DIST


@dataclass(frozen=True)
class State:
    prefix: BASE.State
    disturbance: DIST.BoundedHistory

    def __post_init__(self):
        if not isinstance(self.prefix,BASE.State) or not isinstance(self.disturbance,DIST.BoundedHistory):
            raise TypeError('interleaved prefix and bounded IMU disturbance history required')
        if self.disturbance.sensor_root != self.prefix.live.live_word.sensor_root:
            raise ValueError('interleaved ISS history detached from carried sensor residual histories')

    @property
    def imu_steps(self): return self.prefix.imu_steps
    @property
    def complete_source_horizon(self): return self.prefix.complete_source_horizon


def begin(prefix:BASE.State, disturbance:DIST.BoundedHistory):
    return State(prefix,disturbance)


def imu_step(state:State, **kwargs):
    if not isinstance(state,State): raise TypeError('ISS interleaved prefix required')
    before=state.prefix.imu_steps
    product=IMU.State(state.prefix.live,state.disturbance)
    out=IMU.imu_step(product,**kwargs)
    record=BASE.EventRecord('imu',before,len(out.state.live.live_word.source.steps))
    nxt_prefix=BASE.State(out.state.live,state.prefix.origin,state.prefix.events+(record,))
    return State(nxt_prefix,out.state.disturbance),out


def mag_step(state:State, **kwargs):
    if not isinstance(state,State): raise TypeError('ISS interleaved prefix required')
    nxt,event=BASE.mag_step(state.prefix,**kwargs)
    return State(nxt,state.disturbance),event


def set_hold(state:State, *, hold):
    if not isinstance(state,State): raise TypeError('ISS interleaved prefix required')
    nxt,event=BASE.set_hold(state.prefix,hold=hold)
    return State(nxt,state.disturbance),event


def readiness():
    b=BASE.readiness(); d=IMU.readiness()
    return {
      'admitted_IMU_MAG_HOLD_successor_induction_consumed':b['finite_successor_induction_over_IMU_MAG_HOLD_closed'],
      'one_bounded_IMU_ISS_history_persists_through_entire_prefix':True,
      'only_IMU_consumes_kth_bounded_forcing_restriction':True,
      'MAG_and_HOLD_preserve_bounded_IMU_history_without_consumption':True,
      'Racc_not_reinterpreted_as_pathwise_bound':d['Racc_covariance_not_used_as_pathwise_bound'],
      'canonical_600_transition_source_horizon_representable':b['canonical_600_transition_source_horizon_representable'],
      'bounded_input_history_structurally_attached':True,
      'all_event_arithmetic_witnesses_source_uniformly_qualified':False,
      'magnetic_counter_lifetime_closed':False,
      'startup_reachability_to_admitted_fresh_state_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
