#!/usr/bin/env python3
"""Theorem-facing same-history COMPLETE-BRMM source-cover contract.

The complete frontend predecessor family and same-history coefficient relation
are now closed.  A theorem event cell still must carry one joint tuple
(state,P,H/R geometry,true bias,active tuner/scheduler state); no detached
P/H/R/K or f/sigma/tau/T_S/R_S coordinate is permitted.  Literal event/radial
coverage and endpoint/every-prefix storage remain fail-closed here.
"""
from __future__ import annotations
import argparse,math,json
from dataclasses import dataclass
from pathlib import Path
import ou3_brmm_physical_wave_source as WAVE
from typing import Sequence
from ou3_interval import Interval
import ou3_brmm_complete_source as COMPLETE
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_brmm_full_normal_live_word as WORD
import ou3_brmm_finite_window_primitive_qualification as PRIMITIVE
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_bias1_family as BIAS1
import ou3_p4_complete_brmm_universal_target_relation as TARGET_REL
import ou3_p4_joint_brmm_frontend_transition as BRMM_JOINT

REPO=Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
SCHEMA=10
QUALIFICATION='OU3_P4_COMPLETE_BRMM_SAME_HISTORY_SOURCE_COVER_CONTRACT_V10'
EVENT_KINDS=('prediction','aw_floor','S_zero','accelerometer','magnetometer','H_to_A')

def _shape(A):
    r=len(A);c=len(A[0]) if r else 0
    if any(len(row)!=c for row in A):raise ValueError('ragged matrix')
    return r,c

def _vec(x,n,name):
    if len(x)!=n or any(not isinstance(v,Interval) for v in x):raise ValueError(f'{name} must be Interval[{n}]')

def _mat(A,r,c,name):
    if _shape(A)!=(r,c) or any(not isinstance(v,Interval) for row in A for v in row):raise ValueError(f'{name} must be Interval[{r}x{c}]')

@dataclass(frozen=True)
class SourceCoverCell:
    source_token:str;predecessor_token:str|None;mode:str;sample_index:int;event_ordinal:int;kind:str
    state:Sequence[Interval];P:Sequence[Sequence[Interval]];dt_s:Interval;tau_applied_s:Interval;sigma_aw_mps2:Interval;pseudo_elapsed_s:Interval
    R:Sequence[Sequence[Interval]]|None=None;R_provenance:str|None=None;f_hat:Sequence[Interval]|None=None
    R_hat:Sequence[Sequence[Interval]]|None=None;m_body:Sequence[Interval]|None=None;true_bias:Sequence[Interval]|None=None
    bias_projection_limit:float|None=None;radial_scale:Interval|None=None;estimator_source_token:str|None=None
    estimator_predecessor_token:str|None=None;estimator_generated_coefficients:bool=False
    wave_primitive:KERNEL.WavePrimitivePayload|None=None

def _diag_cov(std):
    if len(std)!=3 or any(not isinstance(x,Interval) or x.lo<=0 for x in std):raise ValueError('positive applied R_S std required')
    z=Interval.point(0)
    return [[std[i].square() if i==j else z for j in range(3)] for i in range(3)]

def estimator_owns_event_token(c):
    e=c.estimator_source_token
    return bool(isinstance(e,str) and e and (c.source_token==e or c.source_token.startswith(e+':e')))

def source_cell_from_joint_image(image:BRMM_JOINT.Image,*,mode,sample_index,event_ordinal,kind,state,P,dt_s,pseudo_elapsed_s,radial_scale,
    event_source_token=None,event_predecessor_token=None,R=None,R_provenance=None,f_hat=None,R_hat=None,m_body=None,true_bias=None,bias_projection_limit=None,wave_primitive=None):
    if not isinstance(image,BRMM_JOINT.Image):raise TypeError('joint estimator image required')
    active=image.active_schedule_for_current_riccati
    if active.tau.lo<=0 or active.sigma.lo<=0 or active.rs_base.lo<=0:raise ValueError('positive active schedule required')
    src=event_source_token or image.source_token; pred=image.predecessor_token if event_predecessor_token is None else event_predecessor_token
    if src!=image.source_token and not src.startswith(image.source_token+':e'):raise ValueError('detached event token')
    if kind=='S_zero':
        if R is not None or R_provenance is not None:raise ValueError('S R must come from active schedule')
        R=_diag_cov(image.actual_rs_std_xyz_for_current_riccati);R_provenance=EVENTS.ACTUAL_RS_PROVENANCE
    c=SourceCoverCell(src,pred,mode,sample_index,event_ordinal,kind,state,P,dt_s,active.tau,active.sigma,pseudo_elapsed_s,R,R_provenance,f_hat,R_hat,m_body,true_bias,bias_projection_limit,radial_scale,image.source_token,image.predecessor_token,True,wave_primitive)
    f=validate_cell(c,require_estimator_provenance=True)
    if f:raise ValueError('invalid joint source cell: '+repr(f))
    return c

def validate_cell(c,*,require_estimator_provenance=False):
    f=[]
    if not isinstance(c.source_token,str) or not c.source_token:f.append('missing source token')
    if c.predecessor_token is not None and not isinstance(c.predecessor_token,str):f.append('invalid predecessor')
    if c.mode not in ('H','A'):return f+['mode must be H/A']
    n=18 if c.mode=='H' else 21
    if c.kind not in EVENT_KINDS:f.append('unsupported event kind')
    if not isinstance(c.sample_index,int) or c.sample_index<0 or not isinstance(c.event_ordinal,int) or c.event_ordinal<0:f.append('invalid event index')
    try:_vec(c.state,n,'state');_mat(c.P,n,n,'P')
    except ValueError as e:f.append(str(e))
    for name,x in (('dt',c.dt_s),('tau',c.tau_applied_s),('sigma',c.sigma_aw_mps2),('pseudo_elapsed',c.pseudo_elapsed_s)):
        if not isinstance(x,Interval) or not(math.isfinite(x.lo) and math.isfinite(x.hi)):f.append(name+' nonfinite')
    if isinstance(c.dt_s,Interval) and c.dt_s.lo<=0:f.append('dt nonpositive')
    if isinstance(c.tau_applied_s,Interval) and c.tau_applied_s.lo<=0:f.append('tau nonpositive')
    if isinstance(c.sigma_aw_mps2,Interval) and c.sigma_aw_mps2.lo<=0:f.append('sigma nonpositive')
    if not isinstance(c.radial_scale,Interval):f.append('radial missing')
    elif c.radial_scale.lo<0 or c.radial_scale.hi>1:f.append('radial outside [0,1]')
    if c.wave_primitive is not None and not isinstance(c.wave_primitive,KERNEL.WavePrimitivePayload):f.append("invalid physical wave ancestry payload")
    if require_estimator_provenance:
        if not c.estimator_generated_coefficients:f.append('coefficients not estimator generated')
        if not estimator_owns_event_token(c):f.append('detached estimator token')
        if not isinstance(c.estimator_predecessor_token,str) or not c.estimator_predecessor_token:f.append('estimator predecessor missing')
    if c.kind in ('S_zero','accelerometer','magnetometer'):
        if c.R is None:f.append('Joseph R missing')
        else:
            try:_mat(c.R,3,3,'R')
            except ValueError as e:f.append(str(e))
    if c.kind=='S_zero' and c.R_provenance!=EVENTS.ACTUAL_RS_PROVENANCE:f.append('actual R_S provenance missing')
    if c.kind=='accelerometer':
        try:_vec(c.f_hat or (),3,'f_hat');_mat(c.R_hat or (),3,3,'R_hat')
        except ValueError as e:f.append(str(e))
    if c.kind=='magnetometer':
        try:_vec(c.m_body or (),3,'m_body')
        except ValueError as e:f.append(str(e))
    if c.mode=='A':
        try:_vec(c.true_bias or (),3,'true_bias')
        except ValueError as e:f.append(str(e))
        if not(c.bias_projection_limit is not None and math.isfinite(float(c.bias_projection_limit)) and float(c.bias_projection_limit)>0):f.append('projection radius missing')
    return list(dict.fromkeys(f))

def joseph_event_kwargs(c):
    f=validate_cell(c)
    if f:raise ValueError('invalid source cell: '+repr(f))
    if c.kind not in ('S_zero','accelerometer','magnetometer'):raise ValueError('not Joseph event')
    return {'mode':c.mode,'state':c.state,'P':c.P,'R':c.R,'kind':c.kind,'f_hat':c.f_hat,'R_hat':c.R_hat,'m_body':c.m_body,'R_provenance':c.R_provenance,'bias_true':c.true_bias,'bias_projection_limit':c.bias_projection_limit}

def build(domain_path=DEFAULT_DOMAIN):
    path=Path(domain_path).resolve();source=COMPLETE.build(path);literal=WORD.build(path);entry=ENTRY.build();bias1=BIAS1.build();target=TARGET_REL.build(path);joint=BRMM_JOINT.build();primitive=PRIMITIVE.build(path)
    bad={'source':COMPLETE.validate(source),'word':WORD.validate(literal),'entry':ENTRY.validate(entry),'BIAS1':BIAS1.validate(bias1),'target':TARGET_REL.validate(target),'joint':BRMM_JOINT.validate(joint),'primitive':PRIMITIVE.validate(primitive)};bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('source-cover prerequisites failed: '+repr(bad))
    realization=source['BRMM_dynamic_realization'];pm=primitive['uniform_physical_primitives'];mom=primitive['finite_window_moment_consequences']
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'cover_cell_API':'SourceCoverCell + source_cell_from_joint_image + validate_cell + joseph_event_kwargs',
      'theorem_cell_coefficients_must_come_from_joint_estimator_image':True,'estimator_ancestry_distinct_from_literal_event_ancestry':True,
      'event_tokens_must_be_owned_by_estimator_token':True,'multiple_literal_events_may_share_one_estimator_image':True,'free_tau_sigma_constructor_for_theorem_forbidden':True,
      'S_zero_R_derived_from_same_active_schedule':True,'same_history_source_token_required':True,'literal_predecessor_relation_required':True,
      'reachable_shipping_Riccati_P_required':True,'same_cell_state_P_geometry_R_tuner_scheduler_required':True,'actual_applied_RS_provenance_required':True,
      'A21_same_source_true_bias_and_projection_required':True,'all_radial_segments_zero_to_full_hard_entry_required':True,
      'every_valid_IMU_prediction_and_accelerometer_required':bool(literal['every_valid_imu_sample_requires_prediction'] and literal['every_valid_imu_sample_requires_accelerometer_Joseph']),
      'all_due_S_and_asynchronous_vector_events_required':bool(literal['S_scheduler_is_executed_not_replaced_by_selected_four'] and literal['magnetometer_is_asynchronous_external_event_family']),
      'complete_BRMM_single_history_required':bool(realization['single_history_required']),'same_history_drives_translation_rotation_frontend_tuner_geometry':bool(realization['same_realization_drives_translation_rotation_frontend_tuner_geometry']),
      'hard_entry_full_declared_scale_required':bool(entry['full_declared_scale_enforced']),'BRMM_uniform_finite_window_primitive_qualification_closed':True,
      'BRMM_uniform_V_m_mps':pm['V_m_norm_upper_mps'],'BRMM_uniform_P_m_m':pm['P_m_norm_upper_m'],'BRMM_independent_S_m_used':False,
      'BRMM_S_dot_equals_p_same_history_required':True,
      'physical_wave_generator_contract':WAVE.build(),
      'physical_wave_primitive_phi_and_Live_phi_same_history_required':True,
      'physical_generator_constraints_consumed_in_joint24_master':False,'BRMM_DeltaS_3s_norm_upper_m_s':mom['Delta_S_norm_upper_from_uniform_position_m_s'],
      'BIAS1_one_common_root_driver_history_required':bool(bias1['one_root_one_parameter_history_required']),
      'independent_P_H_R_K_boxes_forbidden':True,'independent_tuner_RS_schedule_forbidden':True,
      'independent_frequency_sigma_coordinates_forbidden':target['independent_frequency_sigma_coordinates_forbidden'],'same_signal_history_must_generate_period_and_sigma':target['same_signal_history_must_generate_period_and_sigma'],
      'wave_period_moments_must_be_retained':target['wave_period_moments_must_be_retained'],'period_scaled_sigma_band_must_be_retained':target['period_scaled_sigma_band_must_be_retained'],
      'sigma_variance_statistic_must_be_retained':target['sigma_variance_statistic_must_be_retained'],'period_frequency_reciprocal_relation_must_be_retained':target['period_frequency_reciprocal_relation_must_be_retained'],
      'coefficient_outer_postimage_cover_available':target['coefficient_outer_postimage_cover_available'],'coarse_frequency_sigma_rectangle_may_promote_P4':False,
      'coefficient_target_inclusion_closed':target['coefficient_target_inclusion_closed'],'joint_estimator_relation_materialized':target['joint_estimator_relation_materialized_here'],
      'joint_estimator_physical_BRMM_attachment_closed':target['joint_transition_physical_BRMM_attachment_closed'],'complete_BRMM_predecessor_family_covered':target['complete_BRMM_predecessor_family_covered'],
      'tau_TS_SpectralMSE_RS_same_cell_relation_available':target['tau_TS_RS_are_correlated_same_cell_images'],'raw_effective_sigma_split_preserved':target['raw_and_effective_sigma_coordinates_preserved'],
      'trajectory_replay_or_pinned_generator_may_establish_uniform_cover':False,'finite_harmonic_source_may_replace_BRMM':False,'marginal_Pbar_only_cover_forbidden':True,
      'point_trace_schema_can_regress_interface_only':True,'point_trace_can_promote_source_uniform_cover':False,
      'source_cover_transition_operator_materialized':False,'source_cover_all_BRMM_continuations_covered':False,'source_cover_every_radial_segment_covered':False,
      'source_cover_all_Joseph_cells_correlated':False,'SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED':False,'P4_promoted_here':False,
      'next_obligation':'propagate universal estimator images through reachable covariance/event/radial lineage, then close endpoint/every-prefix augmented LDLT and first-exit retention'}

def validate(d):
    f=WAVE.validate(d.get("physical_wave_generator_contract",{}))
    if d.get("physical_wave_primitive_phi_and_Live_phi_same_history_required") is not True:f.append("physical generator primitive ancestry lost")
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema mismatch')
    true=('theorem_cell_coefficients_must_come_from_joint_estimator_image','estimator_ancestry_distinct_from_literal_event_ancestry','event_tokens_must_be_owned_by_estimator_token','multiple_literal_events_may_share_one_estimator_image','free_tau_sigma_constructor_for_theorem_forbidden','S_zero_R_derived_from_same_active_schedule','same_history_source_token_required','literal_predecessor_relation_required','reachable_shipping_Riccati_P_required','same_cell_state_P_geometry_R_tuner_scheduler_required','actual_applied_RS_provenance_required','A21_same_source_true_bias_and_projection_required','all_radial_segments_zero_to_full_hard_entry_required','every_valid_IMU_prediction_and_accelerometer_required','all_due_S_and_asynchronous_vector_events_required','complete_BRMM_single_history_required','same_history_drives_translation_rotation_frontend_tuner_geometry','hard_entry_full_declared_scale_required','BRMM_uniform_finite_window_primitive_qualification_closed','BRMM_S_dot_equals_p_same_history_required','BIAS1_one_common_root_driver_history_required','independent_P_H_R_K_boxes_forbidden','independent_tuner_RS_schedule_forbidden','independent_frequency_sigma_coordinates_forbidden','same_signal_history_must_generate_period_and_sigma','wave_period_moments_must_be_retained','period_scaled_sigma_band_must_be_retained','sigma_variance_statistic_must_be_retained','period_frequency_reciprocal_relation_must_be_retained','coefficient_outer_postimage_cover_available','coefficient_target_inclusion_closed','joint_estimator_relation_materialized','joint_estimator_physical_BRMM_attachment_closed','complete_BRMM_predecessor_family_covered','tau_TS_SpectralMSE_RS_same_cell_relation_available','raw_effective_sigma_split_preserved','marginal_Pbar_only_cover_forbidden','point_trace_schema_can_regress_interface_only')
    for k in true:
        if d.get(k) is not True:f.append(k+' not true')
    false=('physical_generator_constraints_consumed_in_joint24_master','BRMM_independent_S_m_used','coarse_frequency_sigma_rectangle_may_promote_P4','trajectory_replay_or_pinned_generator_may_establish_uniform_cover','finite_harmonic_source_may_replace_BRMM','point_trace_can_promote_source_uniform_cover','source_cover_transition_operator_materialized','source_cover_all_BRMM_continuations_covered','source_cover_every_radial_segment_covered','source_cover_all_Joseph_cells_correlated','SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED','P4_promoted_here')
    for k in false:
        if d.get(k) is not False:f.append(k+' not false')
    if not(0<float(d.get('BRMM_uniform_V_m_mps',0))<=5):f.append('V_m invalid')
    if not(0<float(d.get('BRMM_uniform_P_m_m',0))<8):f.append('P_m invalid')
    if not(0<float(d.get('BRMM_DeltaS_3s_norm_upper_m_s',0))<23):f.append('DeltaS invalid')
    return list(dict.fromkeys(f))

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'contract_ready':not f,'physical_attachment':d['joint_estimator_physical_BRMM_attachment_closed'],'predecessor':d['complete_BRMM_predecessor_family_covered'],'coefficient':d['coefficient_target_inclusion_closed'],'uniform_cover':d['SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
