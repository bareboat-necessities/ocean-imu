"""Compose one asynchronous updateMag event with magnetic source and BA gate.

This event occurs between IMU samples and advances no physical primitive.  The
wrapper clock must equal the current finite physical endpoint time.  Once the
with-mag/delay gate opens, the exact magnetic packet is passed through the
shipping sanity/SafeLDLT measurement relation and that successor is then passed
to the wrapper count/unlock recurrence.  Crucially, wrapper counting does not
read magnetic acceptance.

The world magnetic reference and Rmag are persistent active state. updateMag
cannot choose either per event. How startup/refinement/continuous correction
produce a new active generation remains a separate, fail-closed obligation.
"""
from __future__ import annotations
from dataclasses import dataclass

from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_measurement_runtime as MR
from tools.stability.ou3_alt_contraction import finite_mag_runtime as MAG
from tools.stability.ou3_alt_contraction import finite_mag_bias_gate as GATE
from tools.stability.ou3_alt_contraction import finite_mag_reference_runtime as MAGREF

@dataclass(frozen=True)
class State:
    filter:CORE.State
    control:GATE.State
    magnetic_active:MAGREF.State
    def __post_init__(self):
        if not isinstance(self.filter,CORE.State) or not isinstance(self.control,GATE.State):
            raise TypeError('finite filter and magnetometer control states required')
        if not isinstance(self.magnetic_active,MAGREF.State):
            raise TypeError('persistent active magnetic model required')

@dataclass(frozen=True)
class Result:
    state:State
    gate:GATE.GateResult
    magnetic:MAG.Result|None
    wrapper_attempted:bool
    measurement_accepted:bool


def update_mag_call(state:State,cfg:GATE.Config,*,time,live,
                    sample:MAG.Sample|None=None,ldlt:MR.SafeLDLT|None=None,
                    alpha=1,radius=None):
    if not isinstance(state,State) or not isinstance(cfg,GATE.Config):
        raise TypeError('async-mag state and config required')
    if state.filter.reference.time != GATE.R(time):
        raise ValueError('async wrapper clock detached from current physical endpoint')
    attempted=cfg.with_mag and GATE.R(time)>=cfg.mag_delay
    if not attempted:
        if sample is not None or ldlt is not None or radius is not None:
            raise ValueError('delay/disabled updateMag branch consumes no magnetic proof operands')
        gate=GATE.update_mag_call(state.control,state.filter,cfg,time=time,live=live)
        return Result(State(gate.filter_state,gate.state,state.magnetic_active),gate,None,False,False)
    if not isinstance(sample,MAG.Sample):
        raise TypeError('attempted updateMag requires same-endpoint magnetic packet')
    if sample.physical != state.filter.reference:
        raise ValueError('async magnetic packet detached from filter endpoint')
    MAGREF.require_same(state.magnetic_active,sample)
    mkw={'ldlt':ldlt,'alpha':alpha}
    if radius is not None: mkw['radius']=radius
    magnetic=MAG.update(state.filter,sample,**mkw)
    gate=GATE.update_mag_call(state.control,state.filter,cfg,time=time,live=live,
                              measurement_state=magnetic.state)
    return Result(State(gate.filter_state,gate.state,state.magnetic_active),gate,magnetic,True,magnetic.accepted)


def readiness():
    return {
      'async_updateMag_no_physical_time_advance':True,
      'wrapper_clock_equals_current_physical_endpoint':True,
      'delay_gate_precedes_magnetic_packet_consumption':True,
      'same_packet_to_mag_measurement_then_bias_gate':True,
      'sanity_and_LDLT_rejections_still_increment_wrapper_count':True,
      'H18_to_A21_release_independent_of_mag_acceptance':True,
      'persistent_active_world_reference_and_Rmag_bound_to_event':True,
      'per_event_free_world_reference_forbidden':True,
      'per_event_free_Rmag_forbidden':True,
      'async_mag_call_schedule_source_attached':False,
      'mag_world_reference_startup_ancestry_attached':False,
      'mag_reference_refinement_and_continuous_update_attached':False,
      'mag_noise_source_bound_attached':False,
      'deployment_finite_precision_closed':False,
      'complete_word_finite_identity':False,
      'ALT_LIVE_PASS':False,
    }
