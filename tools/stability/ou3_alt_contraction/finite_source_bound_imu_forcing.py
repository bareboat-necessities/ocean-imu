"""Retain actual IMU disturbance coordinates for the future ALT ISS supply.

The source-bound prediction word owns the physical transition, raw packet and
shipping accelerometer temperature model.  A storage proof must not replace
sensor/model disturbances by independent boxes or silently absorb them into
state error.  This layer therefore derives the forcing record from the SAME
executed IMU event:

* ``gyro_residual`` is the raw packet's persistent gyro-source residual;
* ``guarded_accel_residual`` is the exact descendant after the shipping
  vibration guard, so guard dynamics are not reclassified as free noise;
* ``thermal_model`` is k_a*(tempC-35 C) using the source-locked shipping k_a;
* ``nu_acc`` is the exact accelerometer measurement-model forcing consumed by
  the finite core, equal to guarded_accel_residual-thermal_model.

No amplitude bound is imposed.  These are ISS supply coordinates: a later
storage inequality may charge their actual magnitude.  This module does not
claim bounded-input admissibility, finite precision, contraction or storage.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_source_bound_prediction_word as WORD
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M


@dataclass(frozen=True)
class ImuForcing:
    gyro_residual_internal: tuple
    guarded_accel_residual_internal: tuple
    thermal_model_internal: tuple
    nu_acc_internal: tuple
    temperature_delta_c: object

    def __post_init__(self):
        for name in ('gyro_residual_internal','guarded_accel_residual_internal',
                     'thermal_model_internal','nu_acc_internal'):
            object.__setattr__(self,name,tuple(M.vec(getattr(self,name),3)))
        object.__setattr__(self,'temperature_delta_c',M.rational(self.temperature_delta_c))
        expected=tuple(self.guarded_accel_residual_internal[i]-self.thermal_model_internal[i]
                       for i in range(3))
        if self.nu_acc_internal != expected:
            raise ValueError('accelerometer forcing detached from guarded residual/thermal model')

    @property
    def supply_vector(self):
        """Canonical real forcing coordinates, without an invented weighting."""
        return (self.gyro_residual_internal + self.nu_acc_internal
                + self.thermal_model_internal)


@dataclass(frozen=True)
class Result:
    word: object
    forcing: ImuForcing

    @property
    def state(self):
        return self.word.state


def _executed_prefix(word_result):
    """Reach the literal finite_live_imu_prefix.Result through composed layers."""
    try:
        return word_result.event.event.live
    except AttributeError as exc:
        raise TypeError('source-bound result missing executed Live IMU prefix') from exc


def imu_step(*args, temperature_c, **kwargs):
    """Execute source-bound IMU word and expose its actual disturbance supply."""
    out=WORD.imu_step(*args,temperature_c=temperature_c,**kwargs)
    prefix=_executed_prefix(out)
    guarded=prefix.guarded
    if not isinstance(guarded,SENSOR.GuardedImuSample):
        raise TypeError('executed source-bound word lost guarded same-packet sample')
    conditioning=WORD._accel_conditioning(temperature_c)
    obs=SENSOR.finite_accel_core_observation(guarded,conditioning)
    thermal=tuple(conditioning.k_a_hat_internal[i]*conditioning.temperature_delta
                  for i in range(3))
    forcing=ImuForcing(
        guarded.raw.gyro_residual_internal,
        guarded.effective_accel_residual_internal,
        thermal,
        obs.nu_acc,
        conditioning.temperature_delta)
    return Result(out,forcing)


def readiness():
    lower=WORD.readiness()
    return {
      'same_packet_gyro_residual_retained_as_ISS_supply':True,
      'post_guard_accel_residual_derived_not_free':True,
      'temperature_model_term_retained_separately':True,
      'accelerometer_nu_equals_guarded_residual_minus_thermal_term':True,
      'forcing_supply_derived_from_executed_source_owned_IMU_event':True,
      'physical_source_moments_remain_separate_from_sensor_supply':True,
      'sensor_or_temperature_amplitude_bound_invented':False,
      'bounded_input_history_qualified':False,
      'deployment_roundoff_supply_attached':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
      'lower_model_roots_source_bound': bool(
          lower['prediction_angular_OU_Qaxis_roots_bound_to_same_source_continuation']
          and lower['accelerometer_temperature_coefficient_bound_to_shipping_default']),
    }
