"""Theorem-level bounded IMU disturbance history for ALT ISS.

Configured Racc is a covariance and MUST NOT be reinterpreted as a pathwise
sensor-noise cap.  The stability theorem instead quantifies over an arbitrary
bounded disturbance history.  This module makes that ISS quantifier explicit:

* one persistent pair of gyro/accelerometer residual-history identities is
  inherited from ``SensorDisturbanceRoot``;
* one persistent temperature/model history identity is carried beside them;
* ``supply_norm_upper`` is a symbolic theorem parameter, not a tuned numeric
  sensor envelope;
* each kth event is charged by the exact forcing vector produced by the
  executed shipping event after the vibration guard and temperature model.

Thus a later storage theorem may state an ultimate bound as a function of W,
where ``||w_k|| <= W``.  No probability, sigma multiple, covariance-consistency
assumption, replay, or independent disturbance box is introduced here.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_source_continuation as SOURCE
from tools.stability.ou3_alt_contraction import finite_source_bound_imu_forcing as FORCE

QUALIFICATION='OU3_ALT_BOUNDED_IMU_ISS_HISTORY_V1'
TRANSITIONS=SOURCE.TRANSITIONS


def _norm2(v):
    q=tuple(M.vec(v,len(v)))
    return sum((x*x for x in q),F(0))


@dataclass(frozen=True)
class BoundedHistory:
    """One arbitrary bounded IMU/model forcing history under the theorem quantifier."""
    sensor_root: SOURCE.SensorDisturbanceRoot
    temperature_model_history_id: str
    supply_norm_upper: F
    qualification: str=QUALIFICATION

    def __post_init__(self):
        if not isinstance(self.sensor_root,SOURCE.SensorDisturbanceRoot):
            raise TypeError('persistent sensor disturbance root required')
        if not isinstance(self.temperature_model_history_id,str) or not self.temperature_model_history_id:
            raise ValueError('persistent temperature/model history identity required')
        b=M.rational(self.supply_norm_upper)
        if b < 0:
            raise ValueError('nonnegative symbolic ISS forcing bound required')
        if self.qualification != QUALIFICATION:
            raise ValueError('wrong bounded IMU disturbance theorem qualification')
        object.__setattr__(self,'supply_norm_upper',b)


@dataclass(frozen=True)
class RestrictedForcing:
    """The kth actual executed forcing value of one bounded theorem history."""
    history: BoundedHistory
    ordinal: int
    forcing: FORCE.ImuForcing

    def __post_init__(self):
        if not isinstance(self.history,BoundedHistory) or not isinstance(self.forcing,FORCE.ImuForcing):
            raise TypeError('bounded history and exact executed IMU forcing required')
        if not isinstance(self.ordinal,int) or isinstance(self.ordinal,bool) or not 1 <= self.ordinal <= TRANSITIONS:
            raise ValueError('IMU disturbance ordinal must be one of the 600 source transitions')
        w=self.forcing.supply_vector
        if _norm2(w) > self.history.supply_norm_upper**2:
            raise ValueError('executed IMU ISS forcing exceeds carried theorem bound')


def bind(history:BoundedHistory, *, sensor_root:SOURCE.SensorDisturbanceRoot,
         ordinal:int, forcing:FORCE.ImuForcing):
    """Bind the exact executed event forcing to the same persistent source IDs."""
    if not isinstance(history,BoundedHistory):
        raise TypeError('bounded IMU disturbance history required')
    if history.sensor_root != sensor_root:
        raise ValueError('bounded IMU history detached from carried sensor residual histories')
    return RestrictedForcing(history,ordinal,forcing)


def readiness():
    return {
      'Racc_covariance_not_used_as_pathwise_noise_bound':True,
      'arbitrary_bounded_IMU_ISS_history_quantifier_available':True,
      'symbolic_supply_bound_not_fixed_sigma_multiple':True,
      'persistent_gyro_and_accel_history_ids_carried':True,
      'persistent_temperature_model_history_id_carried':True,
      'kth_bound_applies_to_exact_executed_post_guard_thermal_forcing':True,
      'finite_horizon_probability_used_to_prune_disturbances':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
    }
