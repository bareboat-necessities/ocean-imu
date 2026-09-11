#!/usr/bin/env python3
"""Physical-reference relations required by the ALT same-trajectory theorem.

This module does not create a new sea model. It consumes the exact shared
COMPLETE-BRMM prediction-forcing/source-sector interfaces and supplies finite
reference terms that a differential-only word omits:

* S=0 uses r_S=e_S-S_phys, not e_S;
* bias prediction uses one physical beta history;
* H18 accelerometer residual retains held accelerometer-bias error.
"""
from __future__ import annotations
from typing import Sequence
from ou3_interval import Interval
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_brmm_physical_prediction_forcing as PRED_FORCE
import ou3_p4_brmm_physical_acceleration_witness_sector as PRED_SECTOR
import ou3_p4_S_origin_literal_prefix_invariant as S_ORIGIN
import ou3_p4_bias0_family as BIAS0
import ou3_p4_bias1_family as BIAS1
import ou3_p4_bias2_family as BIAS2

QUALIFICATION="OU3_ALT_PHYSICAL_REFERENCE_FORCING_V3"

def _vec3(x: Sequence[Interval], name: str):
    if len(x)!=3 or any(not isinstance(v,Interval) for v in x): raise ValueError(f"{name} must be Interval[3]")
    return tuple(x)
def _shape(A):
    r=len(A);c=len(A[0]) if r else 0
    if any(len(row)!=c for row in A): raise ValueError("ragged matrix")
    return r,c
def _matvec(A,x):
    r,c=_shape(A)
    if c!=len(x): raise ValueError("matrix/vector dimension mismatch")
    out=[]
    for i in range(r):
        v=Interval.point(0)
        for j in range(c): v=v+A[i][j]*x[j]
        out.append(v)
    return out

def validate_primitive_reference(primitive: KERNEL.WavePrimitivePayload) -> None:
    if not isinstance(primitive,KERNEL.WavePrimitivePayload): raise TypeError("typed physical primitive required")
    if not primitive.generator_id or not primitive.primitive_in_id or not primitive.primitive_out_id: raise ValueError("physical primitive ancestry missing")
    if not primitive.live_origin_id: raise ValueError("one-time Live origin missing")
    _vec3(primitive.centered_S,"centered_S")

def s_zero_residual(error_S,primitive):
    validate_primitive_reference(primitive);e=_vec3(error_S,"e_S");s=_vec3(primitive.centered_S,"centered_S")
    return tuple(e[i]-s[i] for i in range(3))
def s_zero_residual_operator():
    z=Interval.point(0);one=Interval.point(1);neg=Interval.point(-1)
    return [[one if i==j else z for j in range(3)]+[neg if i==j else z for j in range(3)] for i in range(3)]
def bias_error_prediction(error_before,beta_before,beta_after,phi_hat):
    e=_vec3(error_before,"e_ba");b0=_vec3(beta_before,"beta_before");b1=_vec3(beta_after,"beta_after")
    if not isinstance(phi_hat,Interval): raise TypeError("phi_hat must be Interval")
    return tuple(phi_hat*e[i]+b1[i]-phi_hat*b0[i] for i in range(3))
def bias_prediction_operator(phi_hat):
    if not isinstance(phi_hat,Interval): raise TypeError("phi_hat must be Interval")
    z=Interval.point(0);one=Interval.point(1)
    return [[phi_hat if i==j else z for j in range(3)]+[-phi_hat if i==j else z for j in range(3)]+[one if i==j else z for j in range(3)] for i in range(3)]
def h18_accelerometer_residual(base_residual,held_bias_error):
    r=_vec3(base_residual,"base residual");e=_vec3(held_bias_error,"held e_ba")
    return tuple(r[i]+e[i] for i in range(3))

def _upstream_contracts():
    expected={"prediction_forcing":"OU3_P4_BRMM_SAME_HISTORY_PHYSICAL_PREDICTION_FORCING_V1","prediction_sector":"OU3_P4_BRMM_PHYSICAL_ACCELERATION_WITNESS_SECTOR_V1","S_origin":"OU3_P4_SHARED_S_ORIGIN_LITERAL_PREFIX_INVARIANT_V1"}
    actual={"prediction_forcing":getattr(PRED_FORCE,"QUALIFICATION",None),"prediction_sector":getattr(PRED_SECTOR,"QUALIFICATION",None),"S_origin":getattr(S_ORIGIN,"QUALIFICATION",None)}
    if actual!=expected: raise RuntimeError(f"shared physical contract changed: {actual!r}")
    for name,module in (("BIAS0",BIAS0),("BIAS1",BIAS1),("BIAS2",BIAS2)):
        if not callable(getattr(module,"build",None)) or not callable(getattr(module,"validate",None)): raise RuntimeError(name+" definition/validator interface missing")
    actual["bias_interfaces"]="BIAS0/BIAS1/BIAS2 build+validate";return actual

def build():
    upstream=_upstream_contracts()
    return {"qualification":QUALIFICATION,"canonical_source":"COMPLETE_BRMM_NORMAL_LIVE_WORD","upstream_contracts":upstream,
      "shared_exact_BRMM_vs_OU_prediction_forcing_required":True,"shared_exact_BRMM_vs_OU_prediction_forcing_consumed":True,
      "shared_joint_15D_acceleration_witness_sector_required":True,"shared_joint_15D_acceleration_witness_sector_consumed":True,
      "shared_one_time_S_origin_prefix_invariant_required":True,"shared_one_time_S_origin_prefix_invariant_consumed":True,
      "all_BIAS0_BIAS1_BIAS2_definitions_required_separately":True,"S_zero_residual_identity":"r_S=e_S-S_phys",
      "bias_error_identity":"e_ba+=phi_hat*e_ba+beta+-phi_hat*beta","H18_accelerometer_bias_identity":"r_acc_H18=r_acc_18+e_ba_held",
      "physical_reference_algebra_closed":True,"homogeneous_OU_error_prediction_substituted_for_physical_truth":False,"S_phys_replaced_by_zero":False,
      "legacy_300m_s_error_ball_used":False,"independent_bias_driver_used":False,"replay_or_unreachable_perturbation_used":False,
      "attached_to_every_literal_source_cell":False,"complete_word_master_closed":False,"ALT_LIVE_PASS":False,
      "next_obligation":"attach finite S and bias reference terms to every estimator-owned physical H18/A21 event while retaining upstream correlated q15 prediction forcing; then assemble the literal physical word before storage search"}
def validate(d):
    f=[]
    if d.get("qualification")!=QUALIFICATION:f.append("qualification mismatch")
    if d.get("canonical_source")!="COMPLETE_BRMM_NORMAL_LIVE_WORD":f.append("canonical source changed")
    for k in ("shared_exact_BRMM_vs_OU_prediction_forcing_required","shared_joint_15D_acceleration_witness_sector_required","shared_one_time_S_origin_prefix_invariant_required","all_BIAS0_BIAS1_BIAS2_definitions_required_separately","physical_reference_algebra_closed"):
        if d.get(k) is not True:f.append(k+" not true")
    for k in ("homogeneous_OU_error_prediction_substituted_for_physical_truth","S_phys_replaced_by_zero","legacy_300m_s_error_ball_used","independent_bias_driver_used","replay_or_unreachable_perturbation_used","attached_to_every_literal_source_cell","complete_word_master_closed","ALT_LIVE_PASS"):
        if d.get(k) is not False:f.append(k+" not false")
    return f
