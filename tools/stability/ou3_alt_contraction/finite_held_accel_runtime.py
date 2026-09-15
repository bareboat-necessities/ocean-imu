"""Exact held-accelerometer bridge across one physical prediction segment.

SeaStateFusionFilter_OU_III calls ``time_update(gyro,dt)`` and then reuses the
same already-conditioned ``acc_in`` in ``measurement_update_acc_only``.  The raw
sensor packet is therefore rooted at the physical predecessor while the finite
measurement graph is rooted at the post-prediction physical endpoint.

Do not erase that one-step hold mismatch.  For a guarded packet from endpoint k
and target endpoint k+1 this module proves

  y = f_{k+1} + beta_{k+1} + nu_hold

with

  nu_hold = n_a,k - k_hat*dT + Delta_guard
            + (f_k + beta_k - f_{k+1} - beta_{k+1}),

where ``Delta_guard = D_h acc_in - D_h acc_raw``.  Every term is a descendant
of the same PhysicalSegment/raw packet/guard result.  Bounding this forcing is a
later COMPLETE-BRMM/source task; it is not declared small here.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR


@dataclass(frozen=True)
class HeldObservation:
    observed: tuple
    nu_acc: tuple
    guard_delta: tuple
    physical_hold_forcing: tuple
    def __post_init__(self):
        for n in ('observed','nu_acc','guard_delta','physical_hold_forcing'):
            object.__setattr__(self,n,tuple(M.vec(getattr(self,n),3)))


def _specific_force(ref,gravity_world):
    inertial=[ref.acceleration[i]-gravity_world[i] for i in range(3)]
    return tuple(SENSOR.q_rotate(ref.q_world_to_body,inertial))


def observation_from_conditioned(raw:SENSOR.RawImuSample,conditioned_accel_body,
                               segment:PHYS.PhysicalSegment,
                               conditioning:SENSOR.AccelConditioning):
    """Held observation from an explicitly attached conditioned body sample.

    Used by the binary32 guard relation. ``conditioned_accel_body`` is not a
    free disturbance: the caller owns the same-packet machine guard recurrence.
    """
    if not isinstance(raw,SENSOR.RawImuSample) or not isinstance(segment,PHYS.PhysicalSegment):
        raise TypeError('raw sample and same-history PhysicalSegment required')
    if not isinstance(conditioning,SENSOR.AccelConditioning): raise TypeError('AccelConditioning required')
    if conditioning.lever_internal != (0,0,0):
        raise ValueError('current finite measurement theorem is the declared zero-lever branch')
    if raw.physical != segment.before:
        raise ValueError('held accelerometer packet is not rooted at this segment predecessor')
    conditioned=tuple(M.vec(conditioned_accel_body,3))

    before,after=segment.before,segment.after
    for name in ('history_id','bias_root','bias_family'):
        if hasattr(before,name) or hasattr(after,name):
            if not (hasattr(before,name) and hasattr(after,name) and getattr(before,name)==getattr(after,name)):
                raise ValueError(name+' restarts across held accelerometer segment')

    t=conditioning.temperature_delta; k=conditioning.k_a_hat_internal
    modeled=tuple(k[i]*t for i in range(3))
    internal=SENSOR.mv3(raw.deheel_body_to_internal,conditioned)
    observed=tuple(internal[i]-modeled[i] for i in range(3))
    guard_delta=tuple(internal[i]-raw.internal_accel[i] for i in range(3))
    fb=_specific_force(before,raw.gravity_world); fa=_specific_force(after,raw.gravity_world)
    hold=tuple(fb[i]+before.beta[i]-fa[i]-after.beta[i] for i in range(3))
    nu=tuple(raw.accel_residual_internal[i]-modeled[i]+guard_delta[i]+hold[i] for i in range(3))
    expected=tuple(fa[i]+after.beta[i]+nu[i] for i in range(3))
    if observed != expected:
        raise AssertionError('held accelerometer same-history decomposition lost')
    return HeldObservation(observed,nu,guard_delta,hold)


def observation(sample:SENSOR.GuardedImuSample,segment:PHYS.PhysicalSegment,
                conditioning:SENSOR.AccelConditioning):
    if not isinstance(sample,SENSOR.GuardedImuSample):
        raise TypeError('guarded sample required')
    return observation_from_conditioned(sample.raw,sample.conditioned_accel_body,segment,conditioning)


def readiness():
    return {
      'same_raw_guarded_sample_reused_after_prediction':True,
      'pre_to_post_physical_specific_force_hold_term_retained':True,
      'pre_to_post_true_bias_change_retained':True,
      'guard_delta_retained_separately_from_raw_sensor_residual':True,
      'explicit_conditioned_sample_adapter_available_for_machine_guard':True,
      'temperature_model_retained_in_post_prediction_residual':True,
      'held_accel_COMPLETE_BRMM_bound_attached':False,
      'guard_and_temperature_binary32_ancestry_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
