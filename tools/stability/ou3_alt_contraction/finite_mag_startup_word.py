"""Compose startup magnetic gravity admission with the finite tuner prefix.

The wrapper has two distinct event types:
  * IMU samples advance the persistent gravity-alignment certificate;
  * asynchronous updateMag calls inspect that certificate and, only when the
    literal wrapper admission gate opens, advance the magnetic tuner/clock.

This product prevents proof code from invoking the post-eligibility tuner edge
without an actual admission result from the same persistent gate state.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_mag_gravity_gate as GATE
from tools.stability.ou3_alt_contraction import finite_mag_startup_prefix as MAG
from tools.stability.ou3_alt_contraction import finite_mag_tuner_default as TUNER
from tools.stability.ou3_alt_contraction import finite_mag_tilt_frame as TILT
from tools.stability.ou3_alt_contraction import finite_mag_gauge_fix as GAUGE


@dataclass(frozen=True)
class State:
    gate:GATE.State
    mag:MAG.State
    def __post_init__(self):
        if not isinstance(self.gate,GATE.State) or not isinstance(self.mag,MAG.State):
            raise TypeError('startup gravity-gate and magnetic-prefix states required')

@dataclass(frozen=True)
class ImuResult:
    state:State
    gate:GATE.ImuResult

@dataclass(frozen=True)
class MagResult:
    state:State
    admission:GATE.AdmissionResult
    magnetic:MAG.Result|None


def imu_gate_step(state:State,cfg:GATE.Config,*,acc_body,gyro_body,dt,
                  lpf_exp:GATE.LpfExpWitness|None,
                  lpf_norm:GATE.SqrtWitness|None=None,
                  horizontal_norm:GATE.SqrtWitness|None=None,
                  gyro_norm:GATE.SqrtWitness|None=None):
    """Advance only the gravity certificate from the persistent proxy attitude."""
    out=GATE.imu_step(state.gate,cfg,q_proxy_bw=state.mag.proxy.q,
                      acc_body=acc_body,gyro_body=gyro_body,dt=dt,
                      lpf_exp=lpf_exp,lpf_norm=lpf_norm,
                      horizontal_norm=horizontal_norm,gyro_norm=gyro_norm)
    return ImuResult(State(out.state,state.mag),out)


def update_mag_call(state:State,gate_cfg:GATE.Config,tuner_cfg:TUNER.Config,
                    packet:MAG.Packet,*,begun,have_last_imu,mag_ref_set=False,
                    sample_dt,
                    boat_q_norm:TILT.SqrtWitness|None=None,
                    yaw_half:TILT.YawHalfWitness|None=None,
                    mag_norm:TUNER.SqrtWitness|None=None,
                    mean_norm:TUNER.SqrtWitness|None=None,
                    horizontal_sqrt:GAUGE.HorizontalSqrt|None=None):
    """One startup wrapper updateMag call; tuner operands are consumed iff admitted."""
    admission=GATE.startup_mag_admission(
        state.gate,gate_cfg,wrapper_time=packet.wrapper_time,begun=begun,
        have_last_imu=have_last_imu,mag_ref_set=mag_ref_set)
    if not admission.admitted:
        if any(x is not None for x in (boat_q_norm,yaw_half,mag_norm,mean_norm,horizontal_sqrt)):
            raise ValueError('nonadmitted startup updateMag consumes no tuner arithmetic witnesses')
        return MagResult(State(admission.state,state.mag),admission,None)
    if not isinstance(boat_q_norm,TILT.SqrtWitness):
        raise TypeError('admitted startup updateMag requires proxy-quaternion norm witness')
    out=MAG.eligible_update(
        state.mag,tuner_cfg,packet,sample_dt=sample_dt,
        boat_q_norm=boat_q_norm,yaw_half=yaw_half,mag_norm=mag_norm,
        mean_norm=mean_norm,horizontal_sqrt=horizontal_sqrt)
    return MagResult(State(admission.state,out.state),admission,out)


def readiness():
    return {
      'IMU_gravity_gate_and_async_mag_events_are_separate':True,
      'gravity_gate_reads_same_persistent_proxy_as_mag_tuner':True,
      'nonadmitted_updateMag_cannot_advance_tuner_or_clock':True,
      'admitted_updateMag_composes_literal_wrapper_gate_to_tuner':True,
      'gravity_gate_and_tuner_states_persist_jointly':True,
      'startup_raw_mag_physical_source_relation_attached':False,
      'startup_mag_call_schedule_source_attached':False,
      'gravity_gate_binary32_closed':False,
      'complete_word_finite_identity':False,
      'ALT_STARTUP_PASS':False,
    }
