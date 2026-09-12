"""Finite raw-IMU/source provenance for the ALT shipping word.

No numerical sensor-noise bound is invented here.  Shipping receives gyro and
accelerometer vectors in the physical sensor/body frame B.  The OU-III MEKF
rotates each through its steady-wind de-heel map D_h into the internal virtual
body frame B' before attitude/accelerometer modeling, while the private
VerticalAccelComplementary observer consumes the original raw B-frame sample.

This module therefore binds one raw packet by

  D_h gyro_raw = omega_sample_B' + b_g_true + n_g,
  D_h acc_raw  = R_true(a_true-g) + beta + n_a_internal.

The same raw packet feeds private Mahony and the wrapper MEKF APIs; the same
D_h-transformed packet feeds the finite MEKF source identities.  The de-heel
matrix's runtime sin/cos ancestry and all sensor residual bounds remain open.

``n_a_internal`` is deliberately NOT identified with the finite measurement
``nu_acc`` yet: the MEKF predicted force contains estimated temperature-bias and
possible model terms, so their exact difference must be retained in a later
conditioning bridge.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as VERT


def R(x): return M.rational(x)


def q_rotate(q, v):
    """Projective world->B' quaternion rotation, exact in real arithmetic."""
    q=M.vec(q,4); v=M.vec(v,3)
    n2=M.dot(q,q)
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
        if not isinstance(self.physical,PHYS.PhysicalKinematics):
            raise TypeError('same physical kinematics object required')
        for name in ('omega_sample_internal','gyro_residual_internal','accel_residual_internal',
                     'raw_gyro_body','raw_accel_body','gravity_world'):
            object.__setattr__(self,name,tuple(M.vec(getattr(self,name),3)))
        D=tuple(map(tuple,M.mat(self.deheel_body_to_internal,3,3)))
        object.__setattr__(self,'deheel_body_to_internal',D)
        # Shipping de-heel is a rotation.  Keep orthogonality as a hard real
        # relation; its sin/cos finite-precision ancestry remains separate.
        if M.mm([list(r) for r in D],M.transpose([list(r) for r in D])) != M.eye(3):
            raise ValueError('de-heel map must be an orthogonal shipping frame rotation')

        gyro_internal=mv3(D,self.raw_gyro_body)
        omega=self.omega_sample_internal; bg=self.physical.gyro_bias; ng=self.gyro_residual_internal
        expected_g=tuple(omega[i]+bg[i]+ng[i] for i in range(3))
        if gyro_internal != expected_g:
            raise ValueError('de-heeled gyro detached from physical rate/true bias/source residual')

        accel_internal=mv3(D,self.raw_accel_body)
        inertial=[self.physical.acceleration[i]-self.gravity_world[i] for i in range(3)]
        f_internal=q_rotate(self.physical.q_world_to_body,inertial)
        expected_a=tuple(f_internal[i]+self.physical.beta[i]+self.accel_residual_internal[i] for i in range(3))
        if accel_internal != expected_a:
            raise ValueError('de-heeled accelerometer detached from physical specific force/true bias/source residual')

    @property
    def internal_gyro(self): return mv3(self.deheel_body_to_internal,self.raw_gyro_body)
    @property
    def internal_accel(self): return mv3(self.deheel_body_to_internal,self.raw_accel_body)

    def bias_corrected_internal_gyro(self, e_bg):
        e=M.vec(e_bg,3)
        return tuple(self.internal_gyro[i]-(self.physical.gyro_bias[i]-e[i]) for i in range(3))

    def required_bias_corrected_relation(self, e_bg):
        e=M.vec(e_bg,3)
        return tuple(self.omega_sample_internal[i]+e[i]+self.gyro_residual_internal[i] for i in range(3))


def assert_prediction_gyro(sample:RawImuSample,state,gyro_body_raw):
    """Wrapper packet must be this raw B-frame gyro; return shipping B' gyro."""
    graw=tuple(M.vec(gyro_body_raw,3))
    if graw != sample.raw_gyro_body: raise ValueError('MEKF prediction gyro detached from SAME raw IMU packet')
    if tuple(state.reference.gyro_bias) != sample.physical.gyro_bias:
        raise ValueError('prediction physical predecessor differs from raw sensor predecessor')
    if sample.bias_corrected_internal_gyro(state.z[3:6]) != sample.required_bias_corrected_relation(state.z[3:6]):
        raise AssertionError('omega_hat=omega_sample+e_bg+n_g identity lost after de-heel')
    return sample.internal_gyro


def vertical_step_from_raw(vertical_state:VERT.State,vertical_cfg:VERT.Config,sample:RawImuSample,**kwargs):
    """Private Mahony consumes the original raw B-frame pair, before de-heel."""
    if 'gyro' in kwargs or 'acc' in kwargs:
        raise TypeError('private vertical IMU inputs are owned by RawImuSample')
    return VERT.step(vertical_state,vertical_cfg,gyro=sample.raw_gyro_body,acc=sample.raw_accel_body,**kwargs)


def assert_acc_measurement_input(sample:RawImuSample,acc_input_raw):
    """Wrapper passes the same raw acc_in to private vertical and MEKF API."""
    a=tuple(M.vec(acc_input_raw,3))
    if a != sample.raw_accel_body: raise ValueError('MEKF accelerometer input detached from SAME raw IMU packet')
    return sample.internal_accel


def readiness():
    return {
      'raw_body_to_internal_deheel_map_materialized':True,
      'internal_gyro_physical_bias_residual_identity':True,
      'internal_accel_physical_specific_force_bias_residual_identity':True,
      'same_raw_gyro_private_vertical_and_prediction_API':True,
      'same_raw_accel_private_vertical_and_measurement_API':True,
      'omega_hat_equals_omega_sample_plus_e_bg_plus_n_g_after_deheel':True,
      'deheel_sincos_binary32_ancestry_attached':False,
      'raw_accel_residual_to_nu_acc_temperature_conditioning_bridge':False,
      'sensor_residual_source_bounds_attached':False,
      'binary32_sensor_conversion_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
