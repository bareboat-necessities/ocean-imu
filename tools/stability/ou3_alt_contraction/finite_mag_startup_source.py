"""Physical raw startup magnetometer source relation for ALT.

The outer OU-III wrapper feeds *uncorrected* ``mag_body_ned`` to the provisional
startup MagAutoTuner.  For one physical endpoint this module retains

    m_raw_B = R_true(W->B) B_true_W + b_HI_B + n_m_B.

The world field and body-fixed hard-iron offset belong to one persistent model
root; they are not reselected per magnetic call.  The source residual is kept
explicit and unbounded here.  Refinement/Live hard-iron subtraction is a later
relation and must not be back-propagated into this startup packet.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_mag_startup_prefix as PREFIX


def R(x): return M.rational(x)

@dataclass(frozen=True)
class PhysicalEndpoint:
    time:F
    q_world_to_body:tuple
    history_id:str
    def __post_init__(self):
        t=R(self.time); q=tuple(M.vec(self.q_world_to_body,4))
        if t<0 or M.dot(q,q)==0: raise ValueError('valid startup physical endpoint required')
        if not isinstance(self.history_id,str) or not self.history_id:
            raise ValueError('persistent startup physical history id required')
        object.__setattr__(self,'time',t); object.__setattr__(self,'q_world_to_body',q)

@dataclass(frozen=True)
class Model:
    world_field:tuple
    hard_iron_body:tuple=(F(0),F(0),F(0))
    model_root:str='startup-mag-model'
    def __post_init__(self):
        object.__setattr__(self,'world_field',tuple(M.vec(self.world_field,3)))
        object.__setattr__(self,'hard_iron_body',tuple(M.vec(self.hard_iron_body,3)))
        if not isinstance(self.model_root,str) or not self.model_root:
            raise ValueError('persistent magnetic model root required')

@dataclass(frozen=True)
class Sample:
    physical:PhysicalEndpoint
    model:Model
    residual_body:tuple
    raw_body:tuple
    packet_id:str
    def __post_init__(self):
        if not isinstance(self.physical,PhysicalEndpoint) or not isinstance(self.model,Model):
            raise TypeError('startup physical endpoint and magnetic model required')
        n=tuple(M.vec(self.residual_body,3)); raw=tuple(M.vec(self.raw_body,3))
        predicted=tuple(SENSOR.q_rotate(self.physical.q_world_to_body,self.model.world_field))
        expected=tuple(predicted[i]+self.model.hard_iron_body[i]+n[i] for i in range(3))
        if raw!=expected:
            raise ValueError('startup raw magnetic packet detached from true field/hard-iron/source residual')
        if not isinstance(self.packet_id,str) or not self.packet_id:
            raise ValueError('startup magnetic packet id required')
        object.__setattr__(self,'residual_body',n); object.__setattr__(self,'raw_body',raw)
    @property
    def predicted_body_field(self):
        return tuple(SENSOR.q_rotate(self.physical.q_world_to_body,self.model.world_field))
    def packet(self):
        return PREFIX.Packet(self.physical.time,self.raw_body,self.packet_id)


def make_sample(physical:PhysicalEndpoint,model:Model,residual_body,packet_id):
    n=tuple(M.vec(residual_body,3))
    p=tuple(SENSOR.q_rotate(physical.q_world_to_body,model.world_field))
    raw=tuple(p[i]+model.hard_iron_body[i]+n[i] for i in range(3))
    return Sample(physical,model,n,raw,packet_id)


def assert_same_model_root(samples):
    samples=tuple(samples)
    if not samples: raise ValueError('nonempty startup magnetic source sequence required')
    root=samples[0].model.model_root
    field=samples[0].model.world_field; hard=samples[0].model.hard_iron_body
    history=samples[0].physical.history_id
    for s in samples:
        if not isinstance(s,Sample): raise TypeError('startup magnetic samples required')
        if s.model.model_root!=root or s.model.world_field!=field or s.model.hard_iron_body!=hard:
            raise ValueError('startup magnetic field/hard-iron model restarted within one history')
        if s.physical.history_id!=history:
            raise ValueError('startup magnetic physical history restarted')
    return True


def readiness():
    return {
      'startup_raw_mag_equals_true_rotated_field_plus_hard_iron_plus_residual':True,
      'startup_stream_is_uncorrected_before_MagAutoTuner':True,
      'world_field_and_body_hard_iron_have_persistent_model_root':True,
      'startup_physical_history_id_persists_across_mag_samples':True,
      'source_residual_not_absorbed_into_learned_reference':True,
      'startup_mag_noise_bound_attached':False,
      'world_field_physical_bound_attached':False,
      'hard_iron_physical_bound_attached':False,
      'sensor_conversion_binary32_attached':False,
      'complete_word_finite_identity':False,
      'ALT_STARTUP_PASS':False,
    }
