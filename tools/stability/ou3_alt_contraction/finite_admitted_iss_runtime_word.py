"""Master admitted runtime word with one bounded IMU ISS history.

This composes ``finite_admitted_runtime_word`` with the theorem-level bounded
IMU forcing history.  It is the strongest current runtime product: source
admission, BIAS admission, sample-zero origin, shipping event algebra and one
persistent bounded disturbance history coexist in the same state.

The bound W remains a symbolic theorem parameter.  Racc is not a pathwise cap.
MAG/HOLD do not consume IMU forcing ordinals.  Deployment arithmetic, startup
reachability, magnetic counter lifetime and storage remain open.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_admitted_runtime_word as BASE
from tools.stability.ou3_alt_contraction import finite_admitted_imu_disturbance as DIST


@dataclass(frozen=True)
class State:
    runtime: BASE.State
    disturbance: DIST.BoundedHistory

    def __post_init__(self):
        if not isinstance(self.runtime,BASE.State) or not isinstance(self.disturbance,DIST.BoundedHistory):
            raise TypeError('admitted runtime state and bounded IMU disturbance history required')
        if self.disturbance.sensor_root != self.runtime.admitted.live_word.sensor_root:
            raise ValueError('master ISS history detached from carried sensor residual histories')


@dataclass(frozen=True)
class ImuResult:
    state: State
    event: object
    forcing: object
    restriction: DIST.RestrictedForcing
    runtime_word: object


def bind_startup(runtime:BASE.State, disturbance:DIST.BoundedHistory):
    return State(runtime,disturbance)


def imu_step(state:State, **kwargs):
    if not isinstance(state,State): raise TypeError('master ISS runtime state required')
    out=BASE.imu_step(state.runtime,**kwargs)
    ordinal=out.state.admitted.live_word.source.steps[-1].witness.ordinal
    restricted=DIST.bind(state.disturbance,
                         sensor_root=out.state.admitted.live_word.sensor_root,
                         ordinal=ordinal,forcing=out.forcing)
    return ImuResult(State(out.state,state.disturbance),out.event,out.forcing,restricted,out.runtime_word)


def mag_step(state:State, **kwargs):
    if not isinstance(state,State): raise TypeError('master ISS runtime state required')
    out=BASE.mag_step(state.runtime,**kwargs)
    return State(out.state,state.disturbance),out


def set_hold(state:State, *, hold):
    if not isinstance(state,State): raise TypeError('master ISS runtime state required')
    out=BASE.set_hold(state.runtime,hold=hold)
    return State(out,state.disturbance)


def readiness():
    b=BASE.readiness(); d=DIST.readiness()
    return {
      **b,
      'bounded_input_history_qualified':True,
      'bounded_IMU_ISS_history_is_part_of_master_runtime_state':True,
      'same_kth_executed_IMU_supply_charged_to_symbolic_W':True,
      'Racc_covariance_not_used_as_pathwise_noise_bound':d['Racc_covariance_not_used_as_pathwise_noise_bound'],
      'finite_horizon_probability_used_to_prune_disturbances':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
