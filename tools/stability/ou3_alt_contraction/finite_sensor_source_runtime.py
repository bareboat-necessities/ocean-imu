"""Finite raw-IMU/source provenance for the ALT shipping word.

No numerical sensor-noise bound is invented here. Shipping receives gyro and
accelerometer vectors in physical body B. The OU-III MEKF de-heels them into
virtual body B', while the private VerticalAccelComplementary observer consumes
the original raw B-frame pair.

  D_h gyro_raw = omega_sample_B' + b_g_true + n_g,
  D_h acc_raw  = R_true(a_true-g) + beta + n_a_internal.

For the zero-lever theorem branch, shipping's finite accelerometer core removes
its modeled temperature-bias term before the already-proved measurement graph:

  y_core = D_h acc_raw - k_a_hat * (T-T_ref),
  nu_acc = n_a_internal - k_a_hat * (T-T_ref).

Thus ``nu_acc`` is now an exact descendant of the same raw sample and the same
temperature/model parameters; its source bound and those parameters' runtime
ancestry remain open rather than being invented here.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as VERT


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
    if 'gyro' in kwargs or 'acc' in kwargs: raise TypeError('private vertical IMU inputs are owned by RawImuSample')
    return VERT.step(vertical_state,vertical_cfg,gyro=sample.raw_gyro_body,acc=sample.raw_accel_body,**kwargs)


def assert_acc_measurement_input(sample:RawImuSample,acc_input_raw):
    a=tuple(M.vec(acc_input_raw,3))
    if a != sample.raw_accel_body: raise ValueError('MEKF accelerometer input detached from SAME raw IMU packet')
    return sample.internal_accel


def finite_accel_core_observation(sample:RawImuSample,conditioning:AccelConditioning):
    """Exact zero-lever bridge to finite_core accelerometer convention."""
    if not isinstance(conditioning,AccelConditioning): raise TypeError('AccelConditioning required')
    if conditioning.lever_internal != (0,0,0):
        raise ValueError('current finite measurement theorem is the declared zero-lever branch')
    t=conditioning.temperature_delta; k=conditioning.k_a_hat_internal
    modeled=tuple(k[i]*t for i in range(3))
    observed=tuple(sample.internal_accel[i]-modeled[i] for i in range(3))
    nu=tuple(sample.accel_residual_internal[i]-modeled[i] for i in range(3))
    inertial=[sample.physical.acceleration[i]-sample.gravity_world[i] for i in range(3)]
    physical=q_rotate(sample.physical.q_world_to_body,inertial)
    expected=tuple(physical[i]+sample.physical.beta[i]+nu[i] for i in range(3))
    if observed != expected: raise AssertionError('raw->temperature-removed finite accelerometer observation identity lost')
    return AccelCoreObservation(observed,nu)


def readiness():
    return {
      'raw_body_to_internal_deheel_map_materialized':True,
      'internal_gyro_physical_bias_residual_identity':True,
      'internal_accel_physical_specific_force_bias_residual_identity':True,
      'same_raw_gyro_private_vertical_and_prediction_API':True,
      'same_raw_accel_private_vertical_and_measurement_API':True,
      'omega_hat_equals_omega_sample_plus_e_bg_plus_n_g_after_deheel':True,
      'raw_accel_to_temperature_removed_core_observation_bridge':True,
      'nu_acc_from_same_raw_residual_and_temperature_model':True,
      'zero_lever_theorem_branch_enforced':True,
      'deheel_sincos_binary32_ancestry_attached':False,
      'temperature_and_k_a_runtime_ancestry_attached':False,
      'sensor_residual_source_bounds_attached':False,
      'binary32_sensor_conversion_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
