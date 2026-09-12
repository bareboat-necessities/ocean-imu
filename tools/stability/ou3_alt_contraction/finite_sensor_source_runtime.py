"""Finite raw-IMU/source provenance for the ALT shipping word.

This module does NOT invent numerical sensor-noise bounds.  It binds one raw IMU
sample to the physical kinematics already carried by ``PhysicalKinematics`` and
to explicit source residuals:

  gyro_raw = omega_sample_body + b_g_true + n_g,
  acc_raw  = R_true (a_true - g_world) + beta + n_a_raw.

The same raw gyro is then the only gyro allowed at the private Mahony and MEKF
prediction entries.  The same raw accelerometer is the only accelerometer
allowed at the private Mahony and MEKF measurement input boundary.

``n_a_raw`` is deliberately NOT equated to the finite measurement theorem's
``nu_acc``.  Shipping's measurement model adds the estimated temperature term
and other conditioning/model terms; their mismatch belongs in ``nu_acc`` and
requires a separate exact bridge.  Source admission/bounds and binary32 sensor
conversion remain open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as VERT


def R(x): return M.rational(x)


def q_rotate(q, v):
    """Projective world->body quaternion rotation, exact in real arithmetic."""
    q=M.vec(q,4); v=M.vec(v,3)
    n2=M.dot(q,q)
    if not n2: raise ValueError('zero physical quaternion')
    w,x,y,z=q
    # R(q)v with division by ||q||^2 so PhysicalKinematics may be projective.
    return [
        ((w*w+x*x-y*y-z*z)*v[0] + 2*(x*y-w*z)*v[1] + 2*(x*z+w*y)*v[2])/n2,
        (2*(x*y+w*z)*v[0] + (w*w-x*x+y*y-z*z)*v[1] + 2*(y*z-w*x)*v[2])/n2,
        (2*(x*z-w*y)*v[0] + 2*(y*z+w*x)*v[1] + (w*w-x*x-y*y+z*z)*v[2])/n2,
    ]


@dataclass(frozen=True)
class RawImuSample:
    physical: PHYS.PhysicalKinematics
    omega_sample_body: tuple
    gyro_residual: tuple
    accel_residual_raw: tuple
    raw_gyro: tuple
    raw_accel: tuple
    gravity_world: tuple=(F(0),F(0),F(980665,100000))

    def __post_init__(self):
        if not isinstance(self.physical,PHYS.PhysicalKinematics):
            raise TypeError('same physical kinematics object required')
        for name in ('omega_sample_body','gyro_residual','accel_residual_raw','raw_gyro','raw_accel','gravity_world'):
            object.__setattr__(self,name,tuple(M.vec(getattr(self,name),3)))
        omega=self.omega_sample_body; bg=self.physical.gyro_bias; ng=self.gyro_residual
        expected_g=tuple(omega[i]+bg[i]+ng[i] for i in range(3))
        if self.raw_gyro != expected_g:
            raise ValueError('raw gyro detached from physical rate/true bias/source residual')
        inertial=[self.physical.acceleration[i]-self.gravity_world[i] for i in range(3)]
        fbody=q_rotate(self.physical.q_world_to_body,inertial)
        expected_a=tuple(fbody[i]+self.physical.beta[i]+self.accel_residual_raw[i] for i in range(3))
        if self.raw_accel != expected_a:
            raise ValueError('raw accelerometer detached from physical specific force/true bias/source residual')

    def bias_corrected_gyro(self, e_bg):
        """Actual shipping gyro after subtracting b_g_hat=b_g_true-e_bg."""
        e=M.vec(e_bg,3)
        return tuple(self.raw_gyro[i]-(self.physical.gyro_bias[i]-e[i]) for i in range(3))

    def required_bias_corrected_relation(self, e_bg):
        e=M.vec(e_bg,3)
        return tuple(self.omega_sample_body[i]+e[i]+self.gyro_residual[i] for i in range(3))


def assert_prediction_gyro(sample:RawImuSample,state,gyro_body):
    """Reject a prediction fed by any gyro other than this same raw sample."""
    g=tuple(M.vec(gyro_body,3))
    if g != sample.raw_gyro: raise ValueError('MEKF prediction gyro detached from SAME raw IMU sample')
    if tuple(state.reference.gyro_bias) != sample.physical.gyro_bias:
        raise ValueError('prediction physical predecessor differs from raw sensor predecessor')
    if sample.bias_corrected_gyro(state.z[3:6]) != sample.required_bias_corrected_relation(state.z[3:6]):
        raise AssertionError('omega_hat=omega_sample+e_bg+n_g identity lost')
    return g


def vertical_step_from_raw(vertical_state:VERT.State,vertical_cfg:VERT.Config,sample:RawImuSample,**kwargs):
    """Private Mahony consumes exactly this raw gyro/accelerometer pair."""
    if 'gyro' in kwargs or 'acc' in kwargs:
        raise TypeError('private vertical IMU inputs are owned by RawImuSample')
    return VERT.step(vertical_state,vertical_cfg,gyro=sample.raw_gyro,acc=sample.raw_accel,**kwargs)


def assert_acc_measurement_input(sample:RawImuSample,acc_input):
    """Wrapper passes the same acc_in to private vertical and MEKF acc update."""
    a=tuple(M.vec(acc_input,3))
    if a != sample.raw_accel: raise ValueError('MEKF accelerometer input detached from SAME raw IMU sample')
    return a


def readiness():
    return {
      'raw_gyro_physical_bias_residual_identity':True,
      'raw_accel_physical_specific_force_bias_residual_identity':True,
      'same_raw_gyro_private_vertical_and_prediction':True,
      'same_raw_accel_private_vertical_and_measurement_input':True,
      'omega_hat_equals_omega_sample_plus_e_bg_plus_n_g':True,
      'raw_accel_residual_to_nu_acc_temperature_conditioning_bridge':False,
      'sensor_residual_source_bounds_attached':False,
      'binary32_sensor_conversion_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
