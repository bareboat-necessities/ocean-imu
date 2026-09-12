"""Literal safe-LDLT control graph around finite accepted measurements.

The shipping accelerometer path is now represented at two levels.  Component
entries retain same-endpoint raw/guarded identities.  The full IMU-sample entry
uses one guarded predecessor packet after prediction and consumes the exact
``finite_held_accel_runtime`` observation relative to the segment successor,
plus the exact pre-measurement ``finite_racc_runtime`` covariance.

This preserves shipping's predict-then-reuse-held-acceleration ordering without
pretending the held sample was measured at the post-prediction endpoint.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_prediction_graph as P
from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_racc_runtime as RACC
from tools.stability.ou3_alt_contraction import finite_held_accel_runtime as HELD
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS

BUMP_SCALE = F(1,10**6)


@dataclass(frozen=True)
class SafeLDLT:
    first_success: bool
    second_success: bool | None
    noise_scale: F
    machine_epsilon: F
    def __post_init__(self):
        if not isinstance(self.first_success,bool): raise TypeError('literal first LDLT outcome required')
        ns,eps=P.rational(self.noise_scale),P.rational(self.machine_epsilon)
        if ns < 0 or eps <= 0: raise ValueError('nonnegative noise scale and positive epsilon required')
        if self.first_success:
            if self.second_success is not None: raise ValueError('successful first factorization has no retry outcome')
        elif not isinstance(self.second_success,bool): raise TypeError('failed first factorization requires literal retry outcome')
        object.__setattr__(self,'noise_scale',ns); object.__setattr__(self,'machine_epsilon',eps)
    @property
    def bump(self): return max(self.machine_epsilon,BUMP_SCALE*(self.noise_scale+1))
    @property
    def accepted(self): return self.first_success or bool(self.second_success)
    @property
    def innovation_shift(self): return F(0) if self.first_success else self.bump


@dataclass(frozen=True)
class MeasurementRuntimeResult:
    state: CORE.State
    accepted: bool
    branch: SafeLDLT
    accepted_graph: CORE.Accepted | None


def measurement(state, kind, *, ldlt:SafeLDLT, **kwargs):
    if not isinstance(ldlt,SafeLDLT): raise TypeError('safe-LDLT runtime branch required')
    if not ldlt.accepted: return MeasurementRuntimeResult(state,False,ldlt,None)
    accepted=CORE.measurement(state,kind,innovation_shift=ldlt.innovation_shift,**kwargs)
    return MeasurementRuntimeResult(accepted.state,True,ldlt,accepted)


def _gravity(sample,gravity):
    if gravity is None:
        if sample.gravity_world[0] or sample.gravity_world[1]:
            raise ValueError('finite accelerometer core currently uses NED scalar gravity on z')
        return sample.gravity_world[2]
    g=P.rational(gravity)
    if sample.gravity_world != (0,0,g): raise ValueError('core gravity detached from sensor physical model')
    return g


def _accel_event(state,sample,conditioning,*,ldlt,R,gravity=None,guarded=False,**kwargs):
    if sample.physical != state.reference: raise ValueError('accelerometer packet detached from SAME physical reference')
    if guarded: SENSOR.assert_guarded_acc_measurement_input(sample,sample.conditioned_accel_body)
    else: SENSOR.assert_acc_measurement_input(sample,sample.raw_accel_body)
    packet=SENSOR.finite_accel_core_observation(sample,conditioning)
    g=_gravity(sample,gravity)
    if 'observed' in kwargs or 'kind' in kwargs: raise TypeError('accelerometer observation/kind are owned by sensor bridge')
    return measurement(state,'accelerometer',ldlt=ldlt,R=R,observed=packet.observed,gravity=g,**kwargs)


def accelerometer_from_raw(state,sample:SENSOR.RawImuSample,conditioning:SENSOR.AccelConditioning,*,
                           ldlt:SafeLDLT,R,gravity=None,**kwargs):
    if not isinstance(sample,SENSOR.RawImuSample): raise TypeError('RawImuSample required')
    return _accel_event(state,sample,conditioning,ldlt=ldlt,R=R,gravity=gravity,guarded=False,**kwargs)


def accelerometer_from_guarded(state,sample:SENSOR.GuardedImuSample,conditioning:SENSOR.AccelConditioning,*,
                               ldlt:SafeLDLT,R,gravity=None,**kwargs):
    if not isinstance(sample,SENSOR.GuardedImuSample): raise TypeError('GuardedImuSample required')
    return _accel_event(state,sample,conditioning,ldlt=ldlt,R=R,gravity=gravity,guarded=True,**kwargs)


def accelerometer_from_guarded_racc(state,sample:SENSOR.GuardedImuSample,
                                    conditioning:SENSOR.AccelConditioning,
                                    racc:RACC.Result,*,ldlt:SafeLDLT,gravity=None,**kwargs):
    if not isinstance(sample,SENSOR.GuardedImuSample): raise TypeError('GuardedImuSample required')
    if not isinstance(racc,RACC.Result): raise TypeError('finite Racc runtime result required')
    if 'R' in kwargs: raise TypeError('Racc covariance is owned by finite Racc runtime result')
    return _accel_event(state,sample,conditioning,ldlt=ldlt,R=racc.covariance,
                        gravity=gravity,guarded=True,**kwargs)


def accelerometer_from_held_guarded_racc(state,segment:PHYS.PhysicalSegment,
                                         sample:SENSOR.GuardedImuSample,
                                         conditioning:SENSOR.AccelConditioning,
                                         racc:RACC.Result,*,ldlt:SafeLDLT,
                                         gravity=None,**kwargs):
    """Full shipping post-prediction accelerometer event.

    ``sample`` is rooted at ``segment.before`` and must be the same guarded packet
    already used by Mahony. ``state`` must be rooted at ``segment.after``.  The
    observation stays held; HELD derives its exact residual relative to the new
    physical endpoint instead of changing the sample's provenance.
    """
    if not isinstance(segment,PHYS.PhysicalSegment) or not isinstance(sample,SENSOR.GuardedImuSample):
        raise TypeError('PhysicalSegment and GuardedImuSample required')
    if not isinstance(racc,RACC.Result): raise TypeError('finite Racc runtime result required')
    if sample.physical != segment.before or state.reference != segment.after:
        raise ValueError('held accelerometer must span the exact prediction predecessor/successor')
    SENSOR.assert_guarded_acc_measurement_input(sample,sample.conditioned_accel_body)
    packet=HELD.observation(sample,segment,conditioning)
    g=_gravity(sample,gravity)
    if 'R' in kwargs or 'observed' in kwargs or 'kind' in kwargs:
        raise TypeError('held accelerometer observation and Racc are owned by shipping runtime ancestry')
    return measurement(state,'accelerometer',ldlt=ldlt,R=racc.covariance,
                       observed=packet.observed,gravity=g,**kwargs)


def readiness():
    return {
      'first_LDLT_success_branch':True,
      'single_diagonal_bump_retry_branch':True,
      'double_LDLT_failure_rejection_branch':True,
      'rejected_measurement_preserves_state_covariance':True,
      'same_retry_shift_used_by_gain_and_Joseph':True,
      'accelerometer_observation_from_same_raw_packet':True,
      'accelerometer_deheel_and_temperature_removal_attached':True,
      'accelerometer_observation_from_same_guarded_packet':True,
      'accelerometer_guard_deheel_and_temperature_removal_attached':True,
      'guard_effective_residual_is_derived_not_free':True,
      'shipping_accelerometer_Racc_from_same_runtime_result':True,
      'held_predecessor_sample_to_postprediction_measurement_attached':True,
      'held_physical_evolution_forcing_retained':True,
      'temperature_and_k_a_runtime_ancestry_attached':False,
      'guard_exp_sqrt_binary32_ancestry_attached':False,
      'Racc_hypot_sqrt_binary32_ancestry_attached':False,
      'nominal_Racc_stage_ancestry_attached':False,
      'held_accel_COMPLETE_BRMM_bound_attached':False,
      'noise_scale_frobenius_runtime_source_attached':False,
      'Eigen_LDLT_outcomes_finite_precision_attached':False,
      'machine_epsilon_deployment_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
