"""Finite raw-IMU/source provenance for the ALT shipping word.

No numerical sensor-noise bound is invented here. Shipping receives gyro and
accelerometer vectors in physical body B. The raw packet satisfies

  D_h gyro_raw = omega_sample_B' + b_g_true + n_g,
  D_h acc_raw  = R_true(a_true-g) + beta + n_a_internal.

Shipping does not, however, feed ``acc_raw`` directly to private Mahony or the
MEKF. ``AccelVibrationGuard::step`` first produces one conditioned body-frame
sample ``acc_in`` and every accelerometer consumer sees that same descendant.
This module therefore keeps the physical raw packet immutable and represents the
guarded descendant separately.  Its effective post-guard residual is an exact
identity, not a new free disturbance:

  n_a_guard = D_h acc_in - R_true(a_true-g) - beta.

For the zero-lever theorem branch, the finite accelerometer core then removes
its modeled temperature-bias term from that guarded sample:

  y_core = D_h acc_in - k_a_hat * (T-T_ref),
  nu_acc = n_a_guard - k_a_hat * (T-T_ref).

Thus physical sensor residual, deterministic vibration conditioning, and the
measurement-model residual remain distinguishable on one same-history graph.
Their numerical source bounds, guard exp/sqrt binary32 ancestry, temperature
parameters and deployment conversion roundoff remain open rather than invented.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as VERT
from tools.stability.ou3_alt_contraction import finite_accel_guard_runtime as GUARD


def R(x): return M.rational(x)


def q_rotate(q, v):
    q=M.vec(q,4); v=M.vec(v,3); n2=M.dot(q,q)
    if not n2: raise ValueError('zero physical quaternion')
    w,x,y,z=q
    return [
        ((w*w+x*x-y*y-z*z)*v[0] + 2*(x*y-w*z)*v[1] + 2*(x*z+w*y)*v[2])/n2,
        (2*(x*y+w*z)*v[0] + (w*w-x*x+y*y-z*z)*v[1] + 2*(y*z-w*x)*v[2])/n2,
        (2*(x*z-w*y)*v[0] + 2*(y*z+w*x)*v[1] + (w*w-x*x-y*y+z*z)*v[2])/n2,
    ]


def mv3(A,v): return tuple(M.mv(M.mat(A,3,3),M.vec(v,3)))

IDENTITY3=((F(1),F(0),F(0)),(F(0),F(1),F(0)),(F(0),F(0),F(1)))


@dataclass(frozen=True)
class RawImuSample:
    physical: PHYS.PhysicalKinematics
    omega_sample_internal: tuple
    gyro_residual_internal: tuple
    accel_residual_internal: tuple
    raw_gyro_body: tuple
    raw_accel_body: tuple
    gravity_world: tuple=(F(0),F(0),F(980665,100000))
    deheel_body_to_internal: tuple=IDENTITY3

    def __post_init__(self):
        if not isinstance(self.physical,PHYS.PhysicalKinematics): raise TypeError('same physical kinematics object required')
        for name in ('omega_sample_internal','gyro_residual_internal','accel_residual_internal','raw_gyro_body','raw_accel_body','gravity_world'):
            object.__setattr__(self,name,tuple(M.vec(getattr(self,name),3)))
        D=tuple(map(tuple,M.mat(self.deheel_body_to_internal,3,3))); object.__setattr__(self,'deheel_body_to_internal',D)
        if M.mm([list(r) for r in D],M.transpose([list(r) for r in D])) != M.eye(3): raise ValueError('de-heel map must be orthogonal')
        gyro_internal=mv3(D,self.raw_gyro_body)
        expected_g=tuple(self.omega_sample_internal[i]+self.physical.gyro_bias[i]+self.gyro_residual_internal[i] for i in range(3))
        if gyro_internal != expected_g: raise ValueError('de-heeled gyro detached from physical rate/true bias/source residual')
        accel_internal=mv3(D,self.raw_accel_body)
        inertial=[self.physical.acceleration[i]-self.gravity_world[i] for i in range(3)]
        f_internal=q_rotate(self.physical.q_world_to_body,inertial)
        expected_a=tuple(f_internal[i]+self.physical.beta[i]+self.accel_residual_internal[i] for i in range(3))
        if accel_internal != expected_a: raise ValueError('de-heeled accelerometer detached from physical specific force/true bias/source residual')

    @property
    def internal_gyro(self): return mv3(self.deheel_body_to_internal,self.raw_gyro_body)
    @property
    def internal_accel(self): return mv3(self.deheel_body_to_internal,self.raw_accel_body)

    def bias_corrected_internal_gyro(self,e_bg):
        e=M.vec(e_bg,3); return tuple(self.internal_gyro[i]-(self.physical.gyro_bias[i]-e[i]) for i in range(3))
    def required_bias_corrected_relation(self,e_bg):
        e=M.vec(e_bg,3); return tuple(self.omega_sample_internal[i]+e[i]+self.gyro_residual_internal[i] for i in range(3))


@dataclass(frozen=True)
class GuardedImuSample:
    """Exact descendant of one raw packet after AccelVibrationGuard::step."""
    raw: RawImuSample
    guard: GUARD.Result
    def __post_init__(self):
        if not isinstance(self.raw,RawImuSample) or not isinstance(self.guard,GUARD.Result):
            raise TypeError('raw packet and finite guard result required')

    @property
    def physical(self): return self.raw.physical
    @property
    def raw_gyro_body(self): return self.raw.raw_gyro_body
    @property
    def raw_accel_body(self): return self.raw.raw_accel_body
    @property
    def conditioned_accel_body(self): return self.guard.output
    @property
    def gravity_world(self): return self.raw.gravity_world
    @property
    def deheel_body_to_internal(self): return self.raw.deheel_body_to_internal
    @property
    def internal_gyro(self): return self.raw.internal_gyro
    @property
    def internal_accel(self): return mv3(self.deheel_body_to_internal,self.conditioned_accel_body)
    @property
    def effective_accel_residual_internal(self):
        inertial=[self.physical.acceleration[i]-self.gravity_world[i] for i in range(3)]
        physical=q_rotate(self.physical.q_world_to_body,inertial)
        return tuple(self.internal_accel[i]-physical[i]-self.physical.beta[i] for i in range(3))


def guarded_sample(raw:RawImuSample,guard_state:GUARD.State,guard_cfg:GUARD.Config,*,dt,
                   decay:GUARD.DecayWitness|None=None,rms:GUARD.RmsWitness|None=None):
    if not isinstance(raw,RawImuSample): raise TypeError('RawImuSample required')
    result=GUARD.step(guard_state,guard_cfg,raw.raw_accel_body,dt,decay=decay,rms=rms)
    out=GuardedImuSample(raw,result)
    # Re-establish the conditioned physical identity explicitly at construction.
    effective=out.effective_accel_residual_internal
    inertial=[raw.physical.acceleration[i]-raw.gravity_world[i] for i in range(3)]
    physical=q_rotate(raw.physical.q_world_to_body,inertial)
    expected=tuple(physical[i]+raw.physical.beta[i]+effective[i] for i in range(3))
    if out.internal_accel != expected: raise AssertionError('guarded accelerometer physical identity lost')
    return out


@dataclass(frozen=True)
class AccelConditioning:
    """Shipping accelerometer model inputs at one measurement event."""
    temperature_delta: F
    k_a_hat_internal: tuple
    lever_internal: tuple=(F(0),F(0),F(0))
    def __post_init__(self):
        object.__setattr__(self,'temperature_delta',R(self.temperature_delta))
        object.__setattr__(self,'k_a_hat_internal',tuple(M.vec(self.k_a_hat_internal,3)))
        object.__setattr__(self,'lever_internal',tuple(M.vec(self.lever_internal,3)))


@dataclass(frozen=True)
class AccelCoreObservation:
    observed: tuple
    nu_acc: tuple
    def __post_init__(self):
        object.__setattr__(self,'observed',tuple(M.vec(self.observed,3)))
        object.__setattr__(self,'nu_acc',tuple(M.vec(self.nu_acc,3)))


def assert_prediction_gyro(sample:RawImuSample,state,gyro_body_raw):
    graw=tuple(M.vec(gyro_body_raw,3))
    if graw != sample.raw_gyro_body: raise ValueError('MEKF prediction gyro detached from SAME raw IMU packet')
    if tuple(state.reference.gyro_bias) != sample.physical.gyro_bias: raise ValueError('prediction physical predecessor differs from raw sensor predecessor')
    if sample.bias_corrected_internal_gyro(state.z[3:6]) != sample.required_bias_corrected_relation(state.z[3:6]): raise AssertionError('omega_hat identity lost after de-heel')
    return sample.internal_gyro


def vertical_step_from_raw(vertical_state:VERT.State,vertical_cfg:VERT.Config,sample:RawImuSample,**kwargs):
    """Legacy unguarded identity branch; full shipping composition uses guarded."""
    if 'gyro' in kwargs or 'acc' in kwargs: raise TypeError('private vertical IMU inputs are owned by RawImuSample')
    return VERT.step(vertical_state,vertical_cfg,gyro=sample.raw_gyro_body,acc=sample.raw_accel_body,**kwargs)


def vertical_step_from_guarded(vertical_state:VERT.State,vertical_cfg:VERT.Config,sample:GuardedImuSample,**kwargs):
    if not isinstance(sample,GuardedImuSample): raise TypeError('GuardedImuSample required')
    if 'gyro' in kwargs or 'acc' in kwargs: raise TypeError('private vertical IMU inputs are owned by guarded sample')
    return VERT.step(vertical_state,vertical_cfg,gyro=sample.raw_gyro_body,acc=sample.conditioned_accel_body,**kwargs)


def assert_acc_measurement_input(sample:RawImuSample,acc_input_raw):
    """Legacy raw identity assertion retained for existing lower-level tests."""
    a=tuple(M.vec(acc_input_raw,3))
    if a != sample.raw_accel_body: raise ValueError('MEKF accelerometer input detached from SAME raw IMU packet')
    return sample.internal_accel


def assert_guarded_acc_measurement_input(sample:GuardedImuSample,acc_input):
    a=tuple(M.vec(acc_input,3))
    if a != sample.conditioned_accel_body: raise ValueError('MEKF accelerometer input detached from SAME guarded IMU descendant')
    return sample.internal_accel


def finite_accel_core_observation(sample:RawImuSample|GuardedImuSample,conditioning:AccelConditioning):
    """Exact zero-lever bridge to finite_core accelerometer convention."""
    if not isinstance(conditioning,AccelConditioning): raise TypeError('AccelConditioning required')
    if conditioning.lever_internal != (0,0,0):
        raise ValueError('current finite measurement theorem is the declared zero-lever branch')
    if not isinstance(sample,(RawImuSample,GuardedImuSample)): raise TypeError('raw or guarded IMU sample required')
    t=conditioning.temperature_delta; k=conditioning.k_a_hat_internal
    modeled=tuple(k[i]*t for i in range(3))
    observed=tuple(sample.internal_accel[i]-modeled[i] for i in range(3))
    residual=(sample.effective_accel_residual_internal if isinstance(sample,GuardedImuSample)
              else sample.accel_residual_internal)
    nu=tuple(residual[i]-modeled[i] for i in range(3))
    inertial=[sample.physical.acceleration[i]-sample.gravity_world[i] for i in range(3)]
    physical=q_rotate(sample.physical.q_world_to_body,inertial)
    expected=tuple(physical[i]+sample.physical.beta[i]+nu[i] for i in range(3))
    if observed != expected: raise AssertionError('conditioned finite accelerometer observation identity lost')
    return AccelCoreObservation(observed,nu)


def readiness():
    return {
      'raw_body_to_internal_deheel_map_materialized':True,
      'internal_gyro_physical_bias_residual_identity':True,
      'internal_accel_physical_specific_force_bias_residual_identity':True,
      'guarded_accel_is_exact_descendant_of_same_raw_packet':True,
      'guarded_effective_residual_not_free_source':True,
      'same_raw_gyro_private_vertical_and_prediction_API':True,
      'same_guarded_accel_private_vertical_and_measurement_API':True,
      'omega_hat_equals_omega_sample_plus_e_bg_plus_n_g_after_deheel':True,
      'guarded_accel_to_temperature_removed_core_observation_bridge':True,
      'nu_acc_from_same_guarded_residual_and_temperature_model':True,
      'zero_lever_theorem_branch_enforced':True,
      'deheel_sincos_binary32_ancestry_attached':False,
      'guard_exp_sqrt_binary32_ancestry_attached':False,
      'temperature_and_k_a_runtime_ancestry_attached':False,
      'sensor_residual_source_bounds_attached':False,
      'binary32_sensor_conversion_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
