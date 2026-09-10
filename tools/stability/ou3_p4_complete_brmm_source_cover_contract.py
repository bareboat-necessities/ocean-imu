#!/usr/bin/env python3
"""Theorem-facing same-history COMPLETE-BRMM source-cover contract for OU-III P4.

One theorem cell carries one outward SAME-CELL tuple

    (state,P,H/geometry,R,true-bias,active tuner/scheduler state)

owned by one joint estimator image.  Frequency, sigma, tau, T_S and R_S are
functional images of that same physical history; no detached coefficient box is
permitted.  This contract also consumes the uniform finite-window BRMM primitive
qualification (uniform P_m,V_m and derived S_dot=p relation).

The coefficient target relation and its physical BRMM attachment are now closed.
The global predecessor/radial/literal-event cover is intentionally fail-closed
until every continuation and every prefix is propagated and certified.
"""
from __future__ import annotations
import argparse, math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from ou3_interval import Interval
import ou3_brmm_complete_source as COMPLETE
import ou3_brmm_full_normal_live_word as WORD
import ou3_brmm_finite_window_primitive_qualification as PRIMITIVE
import ou3_p4_complete_brmm_differential_events as EVENTS
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_bias1_family as BIAS1
import ou3_p4_complete_brmm_universal_target_relation as TARGET_REL
import ou3_p4_joint_brmm_frontend_transition as BRMM_JOINT

REPO=Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN=REPO/'tools'/'stability'/'ou3_proof_operating_domain.json'
SCHEMA=7
QUALIFICATION='OU3_P4_COMPLETE_BRMM_SAME_HISTORY_SOURCE_COVER_CONTRACT_V7'
EVENT_KINDS=('prediction','aw_floor','S_zero','accelerometer','magnetometer','H_to_A')

def _shape(A):
    r=len(A); c=len(A[0]) if r else 0
    if any(len(row)!=c for row in A): raise ValueError('ragged source-cover matrix')
    return r,c

def _interval_vector(x,n,name):
    if len(x)!=n or any(not isinstance(v,Interval) for v in x):
        raise ValueError(f'{name} must be an outward Interval vector of length {n}')

def _interval_matrix(A,r,c,name):
    if _shape(A)!=(r,c) or any(not isinstance(v,Interval) for row in A for v in row):
        raise ValueError(f'{name} must be an outward {r}x{c} Interval matrix')

@dataclass(frozen=True)
class SourceCoverCell:
    source_token:str
    predecessor_token:str|None
    mode:str
    sample_index:int
    event_ordinal:int
    kind:str
    state:Sequence[Interval]
    P:Sequence[Sequence[Interval]]
    dt_s:Interval
    tau_applied_s:Interval
    sigma_aw_mps2:Interval
    pseudo_elapsed_s:Interval
    R:Sequence[Sequence[Interval]]|None=None
    R_provenance:str|None=None
    f_hat:Sequence[Interval]|None=None
    R_hat:Sequence[Sequence[Interval]]|None=None
    m_body:Sequence[Interval]|None=None
    true_bias:Sequence[Interval]|None=None
    bias_projection_limit:float|None=None
    radial_scale:Interval|None=None
    estimator_source_token:str|None=None
    estimator_predecessor_token:str|None=None
    estimator_generated_coefficients:bool=False

def _diag_cov_from_std(std_xyz):
    if len(std_xyz)!=3 or any(not isinstance(x,Interval) or x.lo<=0 for x in std_xyz):
        raise ValueError('applied R_S std must be three positive outward intervals')
    z=Interval.point(0.0)
    return [[std_xyz[i].square() if i==j else z for j in range(3)] for i in range(3)]

def estimator_owns_event_token(cell:SourceCoverCell)->bool:
    est=cell.estimator_source_token
    return bool(isinstance(est,str) and est and (cell.source_token==est or cell.source_token.startswith(est+':e')))

def source_cell_from_joint_image(image:BRMM_JOINT.Image,*,mode:str,sample_index:int,event_ordinal:int,kind:str,
    state:Sequence[Interval],P:Sequence[Sequence[Interval]],dt_s:Interval,pseudo_elapsed_s:Interval,
    radial_scale:Interval,event_source_token:str|None=None,event_predecessor_token:str|None=None,
    R:Sequence[Sequence[Interval]]|None=None,R_provenance:str|None=None,f_hat:Sequence[Interval]|None=None,
    R_hat:Sequence[Sequence[Interval]]|None=None,m_body:Sequence[Interval]|None=None,
    true_bias:Sequence[Interval]|None=None,bias_projection_limit:float|None=None)->SourceCoverCell:
    if not isinstance(image,BRMM_JOINT.Image):
        raise TypeError('theorem source cell requires a BRMM-attached joint estimator image')
    active=image.active_schedule_for_current_riccati
    if active.tau.lo<=0 or active.sigma.lo<=0 or active.rs_base.lo<=0:
        raise ValueError('joint image lost positive current applied schedule')
    src=event_source_token or image.source_token
    pred=image.predecessor_token if event_predecessor_token is None else event_predecessor_token
    if src!=image.source_token and not src.startswith(image.source_token+':e'):
        raise ValueError('event token must remain owned by the joint estimator image')
    if kind=='S_zero':
        if R is not None or R_provenance is not None:
            raise ValueError('S=0 R must come only from the joint estimator active schedule')
        R=_diag_cov_from_std(image.actual_rs_std_xyz_for_current_riccati)
        R_provenance=EVENTS.ACTUAL_RS_PROVENANCE
    cell=SourceCoverCell(source_token=src,predecessor_token=pred,mode=mode,sample_index=sample_index,
        event_ordinal=event_ordinal,kind=kind,state=state,P=P,dt_s=dt_s,tau_applied_s=active.tau,
        sigma_aw_mps2=active.sigma,pseudo_elapsed_s=pseudo_elapsed_s,R=R,R_provenance=R_provenance,
        f_hat=f_hat,R_hat=R_hat,m_body=m_body,true_bias=true_bias,bias_projection_limit=bias_projection_limit,
        radial_scale=radial_scale,estimator_source_token=image.source_token,
        estimator_predecessor_token=image.predecessor_token,estimator_generated_coefficients=True)
    failures=validate_cell(cell,require_estimator_provenance=True)
    if failures: raise ValueError('joint-image source cell invalid: '+repr(failures))
    return cell

def validate_cell(cell:SourceCoverCell,*,require_estimator_provenance=False):
    f=[]
    if not isinstance(cell.source_token,str) or not cell.source_token:f.append('missing common source-history token')
    if cell.predecessor_token is not None and not isinstance(cell.predecessor_token,str):f.append('invalid predecessor token')
    if cell.mode not in ('H','A'):return f+['mode must be H/A']
    n=18 if cell.mode=='H' else 21
    if cell.kind not in EVENT_KINDS:f.append('unsupported event kind')
    if not isinstance(cell.sample_index,int) or cell.sample_index<0:f.append('invalid sample index')
    if not isinstance(cell.event_ordinal,int) or cell.event_ordinal<0:f.append('invalid event ordinal')
    try:_interval_vector(cell.state,n,'state');_interval_matrix(cell.P,n,n,'P')
    except ValueError as exc:f.append(str(exc))
    for name,x in (('dt_s',cell.dt_s),('tau_applied_s',cell.tau_applied_s),('sigma_aw_mps2',cell.sigma_aw_mps2),('pseudo_elapsed_s',cell.pseudo_elapsed_s)):
        if not isinstance(x,Interval) or not(math.isfinite(x.lo) and math.isfinite(x.hi)):f.append(name+' is not a finite outward interval')
    if isinstance(cell.dt_s,Interval) and cell.dt_s.lo<=0:f.append('dt_s lost positivity')
    if isinstance(cell.tau_applied_s,Interval) and cell.tau_applied_s.lo<=0:f.append('tau_applied_s lost positivity')
    if isinstance(cell.sigma_aw_mps2,Interval) and cell.sigma_aw_mps2.lo<=0:f.append('sigma_aw_mps2 lost positivity')
    if cell.radial_scale is None or not isinstance(cell.radial_scale,Interval):f.append('finite-error radial scale interval missing')
    elif cell.radial_scale.lo<0 or cell.radial_scale.hi>1:f.append('radial scale must lie inside [0,1]')
    if require_estimator_provenance:
        if not cell.estimator_generated_coefficients:f.append('theorem cell coefficients were not generated by joint estimator')
        if not estimator_owns_event_token(cell):f.append('event token is detached from joint estimator source token')
        if not isinstance(cell.estimator_predecessor_token,str) or not cell.estimator_predecessor_token:f.append('joint estimator predecessor missing')
    if cell.kind in ('S_zero','accelerometer','magnetometer'):
        if cell.R is None:f.append('Joseph event missing same-cell R')
        else:
            try:_interval_matrix(cell.R,3,3,'R')
            except ValueError as exc:f.append(str(exc))
    if cell.kind=='S_zero' and cell.R_provenance!=EVENTS.ACTUAL_RS_PROVENANCE:f.append('S=0 event lacks actual applied SpectralMSE R_S provenance')
    if cell.kind=='accelerometer':
        try:_interval_vector(cell.f_hat or (),3,'f_hat');_interval_matrix(cell.R_hat or (),3,3,'R_hat')
        except ValueError as exc:f.append(str(exc))
    if cell.kind=='magnetometer':
        try:_interval_vector(cell.m_body or (),3,'m_body')
        except ValueError as exc:f.append(str(exc))
    if cell.mode=='A':
        try:_interval_vector(cell.true_bias or (),3,'true_bias')
        except ValueError as exc:f.append(str(exc))
        if not(cell.bias_projection_limit is not None and math.isfinite(float(cell.bias_projection_limit)) and float(cell.bias_projection_limit)>0):
            f.append('A21 cell missing deployed bias projection radius')
    return list(dict.fromkeys(f))

def joseph_event_kwargs(cell):
    failures=validate_cell(cell)
    if failures: raise ValueError('invalid source-cover cell: '+repr(failures))
    if cell.kind not in ('S_zero','accelerometer','magnetometer'):raise ValueError('cell is not a Joseph event')
    return {'mode':cell.mode,'state':cell.state,'P':cell.P,'R':cell.R,'kind':cell.kind,'f_hat':cell.f_hat,
            'R_hat':cell.R_hat,'m_body':cell.m_body,'R_provenance':cell.R_provenance,
            'bias_true':cell.true_bias,'bias_projection_limit':cell.bias_projection_limit}

def build(domain_path:Path=DEFAULT_DOMAIN):
    path=Path(domain_path).resolve()
    source=COMPLETE.build(path); literal=WORD.build(path); entry=ENTRY.build(); bias1=BIAS1.build()
    target=TARGET_REL.build(path); joint=BRMM_JOINT.build(); primitive=PRIMITIVE.build(path)
    bad={'complete_source':COMPLETE.validate(source),'literal_word':WORD.validate(literal),'hard_entry':ENTRY.validate(entry),
         'BIAS1':BIAS1.validate(bias1),'target_relation':TARGET_REL.validate(target),'BRMM_joint':BRMM_JOINT.validate(joint),
         'primitive':PRIMITIVE.validate(primitive)}
    bad={k:v for k,v in bad.items() if v}
    if bad: raise RuntimeError('source-cover contract prerequisites failed: '+repr(bad))
    realization=source['BRMM_dynamic_realization']; pm=primitive['uniform_physical_primitives']; moments=primitive['finite_window_moment_consequences']
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'cover_cell_API':'SourceCoverCell + source_cell_from_joint_image + validate_cell + joseph_event_kwargs',
      'theorem_cell_coefficients_must_come_from_joint_estimator_image':True,'estimator_ancestry_distinct_from_literal_event_ancestry':True,
      'event_tokens_must_be_owned_by_estimator_token':True,'multiple_literal_events_may_share_one_estimator_image':True,
      'free_tau_sigma_constructor_for_theorem_forbidden':True,'S_zero_R_derived_from_same_active_schedule':True,
      'same_history_source_token_required':True,'literal_predecessor_relation_required':True,'reachable_shipping_Riccati_P_required':True,
      'same_cell_state_P_geometry_R_tuner_scheduler_required':True,'actual_applied_RS_provenance_required':True,
      'A21_same_source_true_bias_and_projection_required':True,'all_radial_segments_zero_to_full_hard_entry_required':True,
      'every_valid_IMU_prediction_and_accelerometer_required':bool(literal['every_valid_imu_sample_requires_prediction'] and literal['every_valid_imu_sample_requires_accelerometer_Joseph']),
      'all_due_S_and_asynchronous_vector_events_required':bool(literal['S_scheduler_is_executed_not_replaced_by_selected_four'] and literal['magnetometer_is_asynchronous_external_event_family']),
      'complete_BRMM_single_history_required':bool(realization['single_history_required']),
      'same_history_drives_translation_rotation_frontend_tuner_geometry':bool(realization['same_realization_drives_translation_rotation_frontend_tuner_geometry']),
      'hard_entry_full_declared_scale_required':bool(entry['full_declared_scale_enforced']),
      'BRMM_uniform_finite_window_primitive_qualification_closed':True,'BRMM_uniform_V_m_mps':pm['V_m_norm_upper_mps'],
      'BRMM_uniform_P_m_m':pm['P_m_norm_upper_m'],'BRMM_independent_S_m_used':False,
      'BRMM_S_dot_equals_p_same_history_required':True,'BRMM_DeltaS_3s_norm_upper_m_s':moments['Delta_S_norm_upper_from_uniform_position_m_s'],
      'BIAS1_one_common_root_driver_history_required':bool(bias1['one_root_one_parameter_history_required']),
      'independent_P_H_R_K_boxes_forbidden':True,'independent_tuner_RS_schedule_forbidden':True,
      'independent_frequency_sigma_coordinates_forbidden':bool(target['independent_frequency_sigma_coordinates_forbidden']),
      'same_signal_history_must_generate_period_and_sigma':bool(target['same_signal_history_must_generate_period_and_sigma']),
      'wave_period_moments_must_be_retained':bool(target['wave_period_moments_must_be_retained']),
      'period_scaled_sigma_band_must_be_retained':bool(target['period_scaled_sigma_band_must_be_retained']),
      'sigma_variance_statistic_must_be_retained':bool(target['sigma_variance_statistic_must_be_retained']),
      'period_frequency_reciprocal_relation_must_be_retained':bool(target['period_frequency_reciprocal_relation_must_be_retained']),
      'coarse_frequency_sigma_rectangle_may_promote_P4':False,'coefficient_target_inclusion_closed':bool(target['coefficient_target_inclusion_closed']),
      'joint_estimator_relation_materialized':bool(target['joint_estimator_relation_materialized_here']),
      'joint_estimator_physical_BRMM_attachment_closed':bool(target['joint_transition_physical_BRMM_attachment_closed']),
      'tau_TS_SpectralMSE_RS_same_cell_relation_available':bool(target['tau_TS_RS_are_correlated_same_cell_images']),
      'raw_effective_sigma_split_preserved':bool(target['raw_and_effective_sigma_coordinates_preserved']),
      'trajectory_replay_or_pinned_generator_may_establish_uniform_cover':False,'finite_harmonic_source_may_replace_BRMM':False,
      'marginal_Pbar_only_cover_forbidden':True,'point_trace_schema_can_regress_interface_only':True,'point_trace_can_promote_source_uniform_cover':False,
      'source_cover_transition_operator_materialized':False,'source_cover_all_BRMM_continuations_covered':False,
      'source_cover_every_radial_segment_covered':False,'source_cover_all_Joseph_cells_correlated':False,
      'SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED':False,'P4_promoted_here':False,
      'next_obligation':'propagate the qualified P_m/V_m/S_dot=p source through every estimator-owned predecessor/radial/event prefix; then close common endpoint/every-prefix augmented LDLT and first-exit retention'}

def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    true_keys=('theorem_cell_coefficients_must_come_from_joint_estimator_image','estimator_ancestry_distinct_from_literal_event_ancestry',
      'event_tokens_must_be_owned_by_estimator_token','multiple_literal_events_may_share_one_estimator_image','free_tau_sigma_constructor_for_theorem_forbidden',
      'S_zero_R_derived_from_same_active_schedule','same_history_source_token_required','literal_predecessor_relation_required','reachable_shipping_Riccati_P_required',
      'same_cell_state_P_geometry_R_tuner_scheduler_required','actual_applied_RS_provenance_required','A21_same_source_true_bias_and_projection_required',
      'all_radial_segments_zero_to_full_hard_entry_required','every_valid_IMU_prediction_and_accelerometer_required','all_due_S_and_asynchronous_vector_events_required',
      'complete_BRMM_single_history_required','same_history_drives_translation_rotation_frontend_tuner_geometry','hard_entry_full_declared_scale_required',
      'BRMM_uniform_finite_window_primitive_qualification_closed','BRMM_S_dot_equals_p_same_history_required','BIAS1_one_common_root_driver_history_required',
      'independent_P_H_R_K_boxes_forbidden','independent_tuner_RS_schedule_forbidden','independent_frequency_sigma_coordinates_forbidden',
      'same_signal_history_must_generate_period_and_sigma','wave_period_moments_must_be_retained','period_scaled_sigma_band_must_be_retained',
      'sigma_variance_statistic_must_be_retained','period_frequency_reciprocal_relation_must_be_retained','coefficient_target_inclusion_closed',
      'joint_estimator_relation_materialized','joint_estimator_physical_BRMM_attachment_closed','tau_TS_SpectralMSE_RS_same_cell_relation_available',
      'raw_effective_sigma_split_preserved','marginal_Pbar_only_cover_forbidden','point_trace_schema_can_regress_interface_only')
    for k in true_keys:
        if d.get(k) is not True:f.append(k+' not true')
    false_keys=('BRMM_independent_S_m_used','coarse_frequency_sigma_rectangle_may_promote_P4','trajectory_replay_or_pinned_generator_may_establish_uniform_cover',
      'finite_harmonic_source_may_replace_BRMM','point_trace_can_promote_source_uniform_cover','source_cover_transition_operator_materialized',
      'source_cover_all_BRMM_continuations_covered','source_cover_every_radial_segment_covered','source_cover_all_Joseph_cells_correlated',
      'SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED','P4_promoted_here')
    for k in false_keys:
        if d.get(k) is not False:f.append(k+' not false')
    if not(0<float(d.get('BRMM_uniform_V_m_mps',0))<=5.0):f.append('uniform V_m invalid')
    if not(0<float(d.get('BRMM_uniform_P_m_m',0))<8.0):f.append('uniform P_m invalid')
    if not(0<float(d.get('BRMM_DeltaS_3s_norm_upper_m_s',0))<23.0):f.append('correlated DeltaS bound invalid')
    return list(dict.fromkeys(f))

def main():
    import json
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d)
    d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'contract_ready':not f,'primitive_closed':d['BRMM_uniform_finite_window_primitive_qualification_closed'],
      'V_m':d['BRMM_uniform_V_m_mps'],'P_m':d['BRMM_uniform_P_m_m'],'DeltaS3':d['BRMM_DeltaS_3s_norm_upper_m_s'],
      'coefficient_target_inclusion':d['coefficient_target_inclusion_closed'],'physical_attachment':d['joint_estimator_physical_BRMM_attachment_closed'],
      'uniform_cover_closed':d['SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
