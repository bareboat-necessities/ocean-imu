"""Post-eligibility startup magnetometer edge for the ALT finite word.

Once the outer startup gravity/quality gate has admitted an updateMag call,
shipping uses the persistent startup Mahony quaternion as the attitude source,
computes

    dt_mag = t-last_mag_t  if last_mag_t is finite and t>last_mag_t
             cfg.mag_sample_dt otherwise,

stores ``last_mag_sample_t = t`` unconditionally, yaw-strips the SAME Mahony
quaternion, and calls the default MagAutoTuner with the same raw body magnetic
packet.  The Mahony state itself is not advanced by this asynchronous event.

This module intentionally starts *after* the gravity-alignment/eligibility gate.
That gate and the physical raw magnetic source relation remain named open
obligations; neither is represented by a free acceptance boolean here.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_measurement_graph as M
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as VERT
from tools.stability.ou3_alt_contraction import finite_mag_tuner_default as TUNER
from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as TILT
from tools.stability.ou3_alt_contraction import finite_mag_gauge_fix as GAUGE


def R(x): return M.rational(x)


@dataclass(frozen=True)
class Packet:
    wrapper_time:F
    raw_body:tuple
    packet_id:str
    def __post_init__(self):
        t=R(self.wrapper_time)
        if t<0: raise ValueError('nonnegative startup wrapper time required')
        if not isinstance(self.packet_id,str) or not self.packet_id:
            raise ValueError('persistent startup magnetic packet id required')
        object.__setattr__(self,'wrapper_time',t)
        object.__setattr__(self,'raw_body',tuple(M.vec(self.raw_body,3)))


@dataclass(frozen=True)
class State:
    proxy:VERT.State
    tuner:TUNER.State=TUNER.State()
    last_mag_time:F|None=None
    def __post_init__(self):
        if not isinstance(self.proxy,VERT.State) or not isinstance(self.tuner,TUNER.State):
            raise TypeError('persistent Mahony proxy and magnetic tuner states required')
        if self.last_mag_time is not None:
            t=R(self.last_mag_time)
            if t<0: raise ValueError('nonnegative previous magnetic wrapper time required')
            object.__setattr__(self,'last_mag_time',t)


@dataclass(frozen=True)
class Result:
    state:State
    tuner_step:TUNER.StepResult
    dt_mag:F
    packet:Packet


def eligible_update(state:State,cfg:TUNER.Config,packet:Packet,*,sample_dt,
                    boat_q_norm:TILT.SqrtWitness,
                    yaw_half:TILT.YawHalfWitness|None,
                    mag_norm:TUNER.SqrtWitness|None=None,
                    mean_norm:TUNER.SqrtWitness|None=None,
                    horizontal_sqrt:GAUGE.HorizontalSqrt|None=None):
    if not isinstance(state,State) or not isinstance(cfg,TUNER.Config) or not isinstance(packet,Packet):
        raise TypeError('startup mag state/config/packet required')
    fallback=R(sample_dt)
    if fallback<0: raise ValueError('nonnegative configured mag sample dt required')
    t=packet.wrapper_time
    dt_mag=(t-state.last_mag_time
            if state.last_mag_time is not None and t>state.last_mag_time
            else fallback)
    # Wrapper stores last_mag_sample_t_ before calling addSampleWithTiltQuatDt,
    # hence this clock mutation survives both tuner acceptance and rejection.
    out=TUNER.step_from_boat_quaternion(
        state.tuner,cfg,q_boat_bw=state.proxy.q,mag_body=packet.raw_body,dt=dt_mag,
        boat_q_norm=boat_q_norm,yaw_half=yaw_half,mag_norm=mag_norm,
        mean_norm=mean_norm,horizontal_sqrt=horizontal_sqrt)
    return Result(State(state.proxy,out.state,t),out,dt_mag,packet)


def readiness():
    return {
      'startup_mag_reads_persistent_Mahony_quaternion':True,
      'startup_mag_event_does_not_advance_Mahony_state':True,
      'last_mag_sample_time_persisted_across_events':True,
      'strict_time_advance_else_sample_dt_fallback_materialized':True,
      'clock_update_survives_tuner_rejection':True,
      'same_raw_body_packet_reaches_default_tuner':True,
      'gravity_alignment_eligibility_gate_attached':False,
      'startup_updateMag_delay_and_have_last_imu_gate_attached':False,
      'startup_raw_mag_physical_source_relation_attached':False,
      'startup_mag_call_schedule_source_attached':False,
      'deployment_finite_precision_closed':False,
      'complete_word_finite_identity':False,
      'ALT_STARTUP_PASS':False,
    }
