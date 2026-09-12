"""Startup MagAutoTuner-ready transition for the ALT finite word.

When the provisional startup tuner becomes ready before Live, shipping:
  1. writes its gauge-fixed world reference through setMagWorldRef_;
  2. obtains the yaw gauge of the SAME accepted mean;
  3. stores pending_yaw_abs = wrapPi(-yaw_gauge) for later MEKF handoff.

The finite-real graph represents the pending absolute yaw by its unit Z-rotation
quaternion rather than by a free angle scalar.  It is tied exactly to the same
accepted mean direction used by getYawGaugeCorrectionRad().  Native
atan2/wrapPi/AngleAxis binary32 correspondence remains open.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_mag_tuner_default as TUNER
from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as TILT
from tools.stability.ou3_alt_contraction import finite_mag_runtime as MAG
from tools.stability.ou3_alt_contraction import finite_mag_reference_runtime as REF
from tools.stability.ou3_alt_contraction import finite_mag_reference_events as EVENTS


@dataclass(frozen=True)
class Config:
    sigma_internal:tuple
    producer_root:str
    def __post_init__(self):
        s=tuple(M.vec(self.sigma_internal,3))
        if any(x<0 for x in s): raise ValueError('nonnegative constructor magnetic sigma required')
        if not isinstance(self.producer_root,str) or not self.producer_root:
            raise ValueError('persistent magnetic reference producer root required')
        object.__setattr__(self,'sigma_internal',s)

@dataclass(frozen=True)
class PendingYaw:
    """Absolute yaw Z quaternion equivalent to wrapPi(-atan2(mean.y,mean.x))."""
    q_abs:tuple
    mean_xy:tuple
    def __post_init__(self):
        q=tuple(M.vec(self.q_abs,4)); xy=tuple(M.vec(self.mean_xy,2))
        if M.dot(q,q)!=1 or q[1]!=0 or q[2]!=0:
            raise ValueError('unit absolute-yaw Z quaternion required')
        object.__setattr__(self,'q_abs',q); object.__setattr__(self,'mean_xy',xy)

@dataclass(frozen=True)
class Result:
    active_reference:REF.State
    pending_yaw:PendingYaw
    tuner_state:TUNER.State


def transition(tuner:TUNER.State,cfg:Config,*,yaw_half:TILT.YawHalfWitness):
    if not isinstance(tuner,TUNER.State) or not isinstance(cfg,Config):
        raise TypeError('ready startup tuner and reference config required')
    if not tuner.ready or tuner.mean is None or tuner.world_reference is None:
        raise ValueError('startup magnetic reference write requires ready MagAutoTuner state')
    mx,my,_=tuner.mean
    if not isinstance(yaw_half,TILT.YawHalfWitness) or yaw_half.c!=mx or yaw_half.s!=my:
        raise ValueError('startup yaw-gauge witness detached from same accepted magnetic mean')
    # q(+yaw_gauge)=(ch,0,0,sh), so absolute -yaw uses the conjugate.
    q_abs=(yaw_half.cos_half,F(0),F(0),-yaw_half.sin_half)
    pending=PendingYaw(q_abs,(mx,my))
    model=MAG.Model(tuner.world_reference,cfg.sigma_internal)
    witness=EVENTS.WriteWitness('startup_provisional',model,cfg.producer_root)
    active=EVENTS.first_write(witness)
    return Result(active,pending,tuner)


def readiness():
    return {
      'startup_ready_reference_from_same_MagAutoTuner_state':True,
      'startup_reference_write_uses_generation_zero_setMagWorldRef_event':True,
      'constructor_sigma_Rmag_retained_at_first_reference_write':True,
      'pending_absolute_yaw_tied_to_same_accepted_mean':True,
      'pending_yaw_not_applied_to_MEKF_before_handoff':True,
      'atan2_wrapPi_AngleAxis_binary32_attached':False,
      'magnetic_source_envelope_declared':False,
      'complete_word_finite_identity':False,
      'ALT_STARTUP_PASS':False,
    }
