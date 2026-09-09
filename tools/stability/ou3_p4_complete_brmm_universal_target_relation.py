#!/usr/bin/env python3
"""Coarse target enclosure plus same-BRMM joint-estimator provenance for P4.

The clamp rectangle remains sanity evidence only.  The theorem-facing route is
now also materialized on one already-admitted BRMM/private-Mahony predecessor:
the same physical vertical sample drives stillness, period-scaled band/moments,
raw (f,tau,sigma,T_S,R_S), staging, and the current WPE update.  Complete source
inclusion remains false until that operator is propagated over every admitted
BRMM/stillness/radial predecessor cell.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_p4_complete_brmm_target_cell as TARGET
import ou3_p4_complete_brmm_adaptive_transition as ADAPT
import ou3_p4_joint_estimator_transition as JOINT
import ou3_p4_joint_brmm_frontend_transition as BRMM_JOINT

REPO=Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN=REPO/'tools'/'stability'/'ou3_proof_operating_domain.json'
SCHEMA=4
QUALIFICATION='OU3_P4_COMPLETE_BRMM_COARSE_TARGET_PLUS_JOINT_RELATION_V4'

def build(domain_path=DEFAULT_DOMAIN):
    dynamic=DYNAMIC.build(domain_path);df=DYNAMIC.validate(dynamic)
    target=TARGET.build(domain_path);tf=TARGET.validate(target)
    adapt=ADAPT.build(domain_path);af=ADAPT.validate(adapt)
    joint=JOINT.build();jf=JOINT.validate(joint)
    bj=BRMM_JOINT.build();bf=BRMM_JOINT.validate(bj)
    bad={k:v for k,v in (('dynamic',df),('target',tf),('adaptive',af),('joint',jf),('brmm_joint',bf)) if v}
    if bad:raise RuntimeError('target-relation prerequisites failed: '+repr(bad))
    parity=dynamic['source_parity'];inv=dynamic['dynamic_invariant']
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'finite_frontend_history_required':True,
      'all_shipping_tuning_frequency_outputs_enclosed':bool(parity['normal_live_measured_period_selector_present']),
      'all_shipping_raw_sigma_targets_enclosed':bool(parity['raw_sigma_target_uses_1e_minus6_variance_floor']),
      'quiet_continuation_enclosed':bool(dynamic['normal_live_contract']['quiet_case_admitted']),
      'frequency_raw_sigma_target_rectangle':{'frequency_hz':inv['tuning_frequency_hz'],'sigma_target_raw_mps2':inv['sigma_target_raw_mps2']},
      'tau_TS_RS_are_correlated_same_cell_images':bool(target['full_dynamic_target_rectangle_image_valid']),
      'SpectralMSE_RS_same_cell_map':bool(target['SpectralMSE_RS_is_same_tau_sigma_image']),
      'independent_RS_target_forbidden':bool(target['independent_RS_target_forbidden']),
      'active_schedule_connected_by_shipping_EMA':bool(adapt['EMA_candidate_transition_available']),
      'active_schedule_connected_by_staged_commit':bool(adapt['staged_commit_transition_available']),
      'pseudo_event_incidence_connected_by_scheduler':bool(adapt['scheduler_due_not_due_branch_preserved']),
      'raw_and_effective_sigma_coordinates_preserved':bool(adapt['raw_and_effective_sigma_coordinates_preserved']),
      'full_rectangle_is_coarse_sanity_enclosure_only':True,'target_pair_may_span_full_rectangle_in_theorem':False,
      'independent_frequency_sigma_coordinates_forbidden':True,'same_signal_history_must_generate_period_and_sigma':True,
      'wave_period_moments_must_be_retained':True,'period_scaled_sigma_band_must_be_retained':True,'sigma_variance_statistic_must_be_retained':True,'period_frequency_reciprocal_relation_must_be_retained':True,
      'joint_estimator_relation_materialized_here':bool(joint['joint_f_tau_sigma_TS_RS_image_materialized']),
      'joint_transition_rejects_independent_f_sigma':not bool(joint['independent_f_sigma_rectangle_accepted_by_transition']),
      'joint_transition_requires_split_on_dependency_loss':bool(joint['unresolved_wide_cells_require_source_split']),
      'same_BRMM_sample_joint_frontend_transition_materialized':bool(bj['same_private_Mahony_vertical_drives_stillness_tuner_and_WPE'] and bj['joint_raw_f_tau_sigma_TS_RS_tuple_emitted']),
      'same_BRMM_current_applied_schedule_retained':bool(bj['current_Riccati_active_schedule_precedes_current_measurement'] and bj['actual_applied_RS_xyz_retained_from_same_active_schedule']),
      'joint_transition_physical_BRMM_attachment_closed':False,
      'complete_BRMM_predecessor_family_covered':False,'coefficient_target_inclusion_closed':False,
      'physical_vector_geometry_transition_closed_here':False,'complete_source_cover_closed_here':False,'P4_promoted_here':False,
    }
def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('finite_frontend_history_required','all_shipping_tuning_frequency_outputs_enclosed','all_shipping_raw_sigma_targets_enclosed','quiet_continuation_enclosed','tau_TS_RS_are_correlated_same_cell_images','SpectralMSE_RS_same_cell_map','independent_RS_target_forbidden','active_schedule_connected_by_shipping_EMA','active_schedule_connected_by_staged_commit','pseudo_event_incidence_connected_by_scheduler','raw_and_effective_sigma_coordinates_preserved','full_rectangle_is_coarse_sanity_enclosure_only','independent_frequency_sigma_coordinates_forbidden','same_signal_history_must_generate_period_and_sigma','wave_period_moments_must_be_retained','period_scaled_sigma_band_must_be_retained','sigma_variance_statistic_must_be_retained','period_frequency_reciprocal_relation_must_be_retained','joint_estimator_relation_materialized_here','joint_transition_rejects_independent_f_sigma','joint_transition_requires_split_on_dependency_loss','same_BRMM_sample_joint_frontend_transition_materialized','same_BRMM_current_applied_schedule_retained'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('target_pair_may_span_full_rectangle_in_theorem','joint_transition_physical_BRMM_attachment_closed','complete_BRMM_predecessor_family_covered','coefficient_target_inclusion_closed','physical_vector_geometry_transition_closed_here','complete_source_cover_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    box=d.get('frequency_raw_sigma_target_rectangle',{})
    for k in ('frequency_hz','sigma_target_raw_mps2'):
        x=box.get(k,[])
        if len(x)!=2 or not 0<float(x[0])<=float(x[1]):f.append(k+' invalid')
    return f
def main():
    p=argparse.ArgumentParser();p.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build(a.domain);f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'joint':d['joint_estimator_relation_materialized_here'],'BRMM_transition':d['same_BRMM_sample_joint_frontend_transition_materialized'],'coefficient_closed':d['coefficient_target_inclusion_closed'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
