"""Exact paired prediction displacement from machine tuner coefficients.

This is the finite relation needed before machine R_S/S-service can be composed
soundly. From one physical segment/raw IMU packet, execute paired predictions:

* the source-owned exact-shadow prediction roots; and
* one global-compiler machine-active prediction-root family.

The legacy local comparison starts both from one predecessor and requires
identical attitude/BA roots. The persistent continuation instead supplies a
distinct full machine predecessor; its source-owned roots must carry that same
state and derive its own corrected gyro and H18/A21 bias branch. This retains
prior mean/covariance defects and permits different nominal attitudes/charts.
Both finite states and their chart indices remain in the relation: subtracting
their coordinate vectors does not identify different atlas charts physically.

The output retains all 24 joint state coordinates and every entry of the full
21x21 covariance displacement.  It is an exact finite relation, not a norm box
and not a storage certificate.  Post-prediction floor/scheduler/S service and
later accelerometer correction are deliberately not executed here.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction as F

from tools.stability.ou3_alt_contraction import finite_core as CORE
from tools.stability.ou3_alt_contraction import finite_prediction_runtime as PRED
from tools.stability.ou3_alt_contraction import finite_source_bound_prediction_word as EXACTROOT
from tools.stability.ou3_alt_contraction import finite_machine_active_prediction_roots as MACHROOT
from tools.stability.ou3_alt_contraction import finite_sensor_source_runtime as SENSOR
from tools.stability.ou3_alt_contraction import finite_physical_prediction as PHYS

QUALIFICATION='OU3_ALT_MACHINE_PREDICTION_DISPLACEMENT_V1'


def _subvec(a,b): return tuple(F(x)-F(y) for x,y in zip(a,b))
def _submat(a,b): return tuple(tuple(F(x)-F(y) for x,y in zip(ra,rb)) for ra,rb in zip(a,b))


@dataclass(frozen=True)
class Supply:
    z:tuple
    covariance:tuple
    def __post_init__(self):
        z=tuple(F(x) for x in self.z)
        p=tuple(tuple(F(x) for x in row) for row in self.covariance)
        if len(z)!=24 or len(p)!=21 or any(len(row)!=21 for row in p):
            raise ValueError('full joint24 and 21x21 prediction displacement required')
        object.__setattr__(self,'z',z); object.__setattr__(self,'covariance',p)


@dataclass(frozen=True)
class Relation:
    predecessor:CORE.State
    exact:CORE.State
    machine:CORE.State
    exact_roots:EXACTROOT.Roots
    machine_roots:MACHROOT.Roots
    supply:Supply
    qualification:str=QUALIFICATION
    machine_predecessor:CORE.State|None=None
    def __post_init__(self):
        if not isinstance(self.predecessor,CORE.State) or not isinstance(self.exact,CORE.State) or not isinstance(self.machine,CORE.State):
            raise TypeError('predecessor and both paired prediction states required')
        if not isinstance(self.exact_roots,EXACTROOT.Roots) or not isinstance(self.machine_roots,MACHROOT.Roots):
            raise TypeError('exact and machine prediction roots required')
        if self.qualification!=QUALIFICATION: raise ValueError('wrong machine prediction-displacement qualification')
        machine_pre=self.predecessor if self.machine_predecessor is None else self.machine_predecessor
        if not isinstance(machine_pre,CORE.State): raise TypeError('full machine predecessor required')
        if machine_pre.reference!=self.predecessor.reference or self.exact.reference!=self.machine.reference:
            raise ValueError('paired predictions must share physical predecessor and successor')
        if self.machine_roots.predecessor!=machine_pre:
            raise ValueError('machine prediction roots detached from carried predecessor')
        if (machine_pre==self.predecessor and self.machine_roots.machine_packet is None
                and self.exact.q_hat!=self.machine.q_hat):
            raise ValueError('tuner-only prediction displacement unexpectedly changed nominal attitude')
        if self.supply.z!=_subvec(self.machine.z,self.exact.z):
            raise ValueError('joint24 prediction supply detached from paired successors')
        if self.supply.covariance!=_submat(self.machine.covariance,self.exact.covariance):
            raise ValueError('covariance prediction supply detached from paired successors')


def compare(predecessor:CORE.State,segment:PHYS.PhysicalSegment,raw:SENSOR.RawImuSample,
            exact_roots:EXACTROOT.Roots,machine_roots:MACHROOT.Roots,*,Qbase,
            machine_predecessor:CORE.State|None=None,
            machine_attitude_solver=None,
            use_exact_attitude_Q=True,attitude_first_ldlt_success=True,
            attitude_second_ldlt_success=None):
    if not isinstance(predecessor,CORE.State) or not isinstance(segment,PHYS.PhysicalSegment) or not isinstance(raw,SENSOR.RawImuSample):
        raise TypeError('one predecessor, physical segment and raw packet required')
    if not isinstance(exact_roots,EXACTROOT.Roots) or not isinstance(machine_roots,MACHROOT.Roots):
        raise TypeError('exact and machine source-owned prediction roots required')
    machine_pre=predecessor if machine_predecessor is None else machine_predecessor
    if not isinstance(machine_pre,CORE.State) or machine_pre.reference!=predecessor.reference:
        raise ValueError('machine predecessor detached from same physical source')
    if machine_roots.predecessor!=machine_pre:
        raise ValueError('machine prediction roots detached from carried predecessor')
    packet=machine_roots.machine_packet
    if packet is not None and packet.raw!=raw:
        raise ValueError('machine prediction API projection detached from same raw packet')
    if machine_pre==predecessor and packet is None:
        if exact_roots.angular!=machine_roots.angular:
            raise ValueError('tuner prediction supply cannot absorb attitude-root discrepancy')
        if exact_roots.bias!=machine_roots.bias:
            raise ValueError('tuner prediction supply cannot absorb BA-root discrepancy')
    if machine_roots.active_join.exact.require_prediction(ou=exact_roots.ou,qaxis=exact_roots.qaxis) is not True:
        raise AssertionError('exact active/root consistency failed')

    common=dict(Qbase=Qbase,use_exact_attitude_Q=use_exact_attitude_Q,
                attitude_first_ldlt_success=attitude_first_ldlt_success,
                attitude_second_ldlt_success=attitude_second_ldlt_success)
    machine_common=dict(common)
    if machine_attitude_solver is not None:
        allowed={'attitude_first_ldlt_success','attitude_second_ldlt_success'}
        if not isinstance(machine_attitude_solver,dict) or set(machine_attitude_solver)-allowed:
            raise TypeError('only same-mode attitude LDLT outcomes may differ')
        machine_common.update(machine_attitude_solver)
    exact=PRED.prediction_from_raw(predecessor,segment,raw,angular=exact_roots.angular,
                                  ou=exact_roots.ou,bias=exact_roots.bias,qaxis=exact_roots.qaxis,**common)
    if packet is None:
        machine=PRED.prediction_from_raw(machine_pre,segment,raw,angular=machine_roots.angular,
                                        ou=machine_roots.ou,bias=machine_roots.bias,qaxis=machine_roots.qaxis,**machine_common)
    else:
        gyro,_=MACHROOT.projected_gyro(packet,machine_pre)
        machine=PRED.prediction(machine_pre,segment,gyro_body=gyro,angular=machine_roots.angular,
                                ou=machine_roots.ou,bias=machine_roots.bias,qaxis=machine_roots.qaxis,**machine_common)
    return Relation(predecessor,exact,machine,exact_roots,machine_roots,
                    Supply(_subvec(machine.z,exact.z),_submat(machine.covariance,exact.covariance)),
                    machine_predecessor=machine_pre)


def readiness():
    return {
      'one_predecessor_segment_packet_shared_by_exact_and_machine_prediction':True,
      'attitude_and_BA_roots_forbidden_from_hiding_inside_tuner_supply':True,
      'machine_tau_effect_on_all_joint24_prediction_coordinates_retained':True,
      'machine_tau_Sigma_effect_on_full_21x21_covariance_retained':True,
      'nominal_attitude_and_physical_successor_identity_checked':True,
      'full_prediction_coefficient_displacement_relation_attached':True,
      'distinct_compiler_CORE_predecessor_and_source_owned_roots_supported':True,
      'preceding_machine_covariance_and_bias_error_not_reset_to_shadow':True,
      'post_prediction_floor_scheduler_and_S_service_reexecuted_from_machine_state':False,
      'accelerometer_measurement_reexecuted_from_machine_state':False,
      'source_uniform_machine_prediction_supply_bound_closed':False,
      'machine_root_effect_injected_into_complete_admitted_event_relation':False,
      'all_target_libm_and_Eigen_correspondence_closed':False,
      'source_uniform_complete_600_step_word_qualified':False,
      'storage_search_allowed':False,
      'ALT_LIVE_PASS':False,
    }
