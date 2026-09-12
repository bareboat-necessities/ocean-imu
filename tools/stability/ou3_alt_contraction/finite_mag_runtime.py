"""Finite asynchronous magnetic source / Rmag / measurement relation for ALT.

One external updateMag packet is rooted at the current physical endpoint and
satisfies the literal zero-hard-iron theorem branch

    D_h m_raw = R_true B_world + n_m .

The same de-heeled packet is the MEKF observation.  Rmag is generated only from
the configured per-axis magnetic standard deviations.  Shipping rejects a
non-finite or <=1e-6-norm packet before SafeLDLT; in this exact-rational layer
non-finite deployment values remain a binary32/source obligation, while the
norm threshold branch is represented exactly.

This is event algebra, not a COMPLETE-BRMM magnetic-noise bound and not a proof
that the configured world reference came from an admitted startup history.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MR
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS


def R(x): return M.rational(x)

@dataclass(frozen=True)
class Model:
    world_reference:tuple
    sigma_internal:tuple
    def __post_init__(self):
        b=tuple(M.vec(self.world_reference,3)); s=tuple(M.vec(self.sigma_internal,3))
        if any(x<0 for x in s): raise ValueError('nonnegative magnetic standard deviations required')
        object.__setattr__(self,'world_reference',b); object.__setattr__(self,'sigma_internal',s)
    @property
    def covariance(self):
        return tuple(tuple(self.sigma_internal[i]*self.sigma_internal[i] if i==j else F(0)
                           for j in range(3)) for i in range(3))

@dataclass(frozen=True)
class Sample:
    physical:PHYS.PhysicalKinematics
    raw_body:tuple
    residual_internal:tuple
    model:Model
    deheel_body_to_internal:tuple=SENSOR.IDENTITY3
    def __post_init__(self):
        if not isinstance(self.physical,PHYS.PhysicalKinematics) or not isinstance(self.model,Model):
            raise TypeError('physical endpoint and magnetic model required')
        raw=tuple(M.vec(self.raw_body,3)); n=tuple(M.vec(self.residual_internal,3))
        D=tuple(map(tuple,M.mat(self.deheel_body_to_internal,3,3)))
        if M.mm([list(r) for r in D],M.transpose([list(r) for r in D])) != M.eye(3):
            raise ValueError('mag de-heel map must be orthogonal')
        object.__setattr__(self,'raw_body',raw); object.__setattr__(self,'residual_internal',n)
        object.__setattr__(self,'deheel_body_to_internal',D)
        obs=SENSOR.mv3(D,raw)
        physical=tuple(SENSOR.q_rotate(self.physical.q_world_to_body,self.model.world_reference))
        expected=tuple(physical[i]+n[i] for i in range(3))
        if obs != expected: raise ValueError('mag packet detached from physical field/source residual')
    @property
    def observed_internal(self): return SENSOR.mv3(self.deheel_body_to_internal,self.raw_body)
    @property
    def norm2(self): return M.dot(self.observed_internal,self.observed_internal)

@dataclass(frozen=True)
class Result:
    state:CORE.State
    attempted_measurement:bool
    accepted:bool
    measurement:MR.MeasurementRuntimeResult|None
    sample:Sample


def update(state:CORE.State,sample:Sample,*,ldlt:MR.SafeLDLT|None=None,alpha=1,radius=F(2,5)):
    if not isinstance(state,CORE.State) or not isinstance(sample,Sample):
        raise TypeError('finite filter state and magnetic sample required')
    if state.reference != sample.physical:
        raise ValueError('mag packet detached from current physical endpoint')
    # C++ checks mag_norm > 1e-6. Avoid sqrt while preserving the same branch.
    if sample.norm2 <= F(1,10**12):
        if ldlt is not None: raise ValueError('sanity-rejected mag packet consumes no LDLT branch')
        return Result(state,False,False,None,sample)
    if not isinstance(ldlt,MR.SafeLDLT):
        raise TypeError('valid mag packet requires safe-LDLT branch')
    m=MR.measurement(state,'magnetometer',ldlt=ldlt,R=sample.model.covariance,
                     observed=sample.observed_internal,
                     magnetic_reference=sample.model.world_reference,alpha=alpha,radius=radius)
    return Result(m.state,True,m.accepted,m,sample)


def readiness():
    return {
      'mag_raw_body_to_internal_deheel_materialized':True,
      'mag_physical_world_reference_residual_identity':True,
      'same_mag_packet_supplies_observation_and_residual':True,
      'Rmag_diagonal_from_configured_sigma_only':True,
      'mag_norm_gt_1e6_threshold_inverse_square_materialized':True,
      'safe_LDLT_accept_retry_reject_reused':True,
      'mag_world_reference_startup_ancestry_attached':False,
      'mag_sigma_constructor_runtime_ancestry_attached':False,
      'mag_noise_source_bound_attached':False,
      'deheel_and_norm_binary32_attached':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
