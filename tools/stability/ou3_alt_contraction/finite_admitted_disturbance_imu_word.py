"""Strongest admitted-source IMU edge with theorem-level bounded ISS history.

This layer carries one bounded disturbance history as theorem state and binds the
exact forcing produced by ``finite_admitted_source_imu_word`` at ordinal k to
that same history.  The event itself remains the already-composed shipping edge;
this wrapper adds no filter approximation and no numeric noise envelope.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_source_imu_word as IMU
from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as LIVE
from tools.stability.ou3_alt_contraction import finite_admitted_imu_disturbance as DIST


@dataclass(frozen=True)
class State:
    live: LIVE.State
    disturbance: DIST.BoundedHistory

    def __post_init__(self):
        if not isinstance(self.live,LIVE.State) or not isinstance(self.disturbance,DIST.BoundedHistory):
            raise TypeError('admitted Live state and bounded IMU disturbance history required')
        if self.disturbance.sensor_root != self.live.live_word.sensor_root:
            raise ValueError('bounded IMU history detached from carried sensor residual histories')


@dataclass(frozen=True)
class Result:
    state: State
    event: object
    forcing: object
    restriction: DIST.RestrictedForcing


def imu_step(state:State, **kwargs):
    if not isinstance(state,State):
        raise TypeError('admitted Live + bounded disturbance product required')
    out=IMU.imu_step(state.live,**kwargs)
    ordinal=out.state.live_word.source.steps[-1].witness.ordinal
    restricted=DIST.bind(state.disturbance,
                         sensor_root=out.state.live_word.sensor_root,
                         ordinal=ordinal,forcing=out.forcing)
    return Result(State(out.state,state.disturbance),out.event,out.forcing,restricted)


def readiness():
    i=IMU.readiness(); d=DIST.readiness()
    return {
      'strongest_admitted_source_IMU_edge_consumed':True,
      'bounded_IMU_ISS_history_carried_as_theorem_state':d['arbitrary_bounded_IMU_ISS_history_quantifier_available'],
      'same_persistent_sensor_root_required':True,
      'same_kth_executed_forcing_charged_to_history':True,
      'Racc_covariance_not_used_as_pathwise_bound':d['Racc_covariance_not_used_as_pathwise_noise_bound'],
      'same_BRMM_BIAS_event_and_forcing_retained':i['physical_moments_sensor_forcing_and_model_roots_share_one_event'],
      'deployment_roundoff_supply_attached':False,
      'complete_600_step_shipping_word_composed':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
