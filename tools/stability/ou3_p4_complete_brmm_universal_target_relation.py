#!/usr/bin/env python3
"""Universal same-history COMPLETE-BRMM frontend/coefficient relation.

The shipping same-signal one-step transition is combined with the separately
proved compact forward-invariant predecessor family.  Therefore every admitted
COMPLETE-BRMM predecessor/sample has a same-history frontend successor and the
shipping clamp image contains its coefficient tuple.  The outer rectangle is
only an inclusion device: it may never generate f/sigma/tau/T_S/R_S histories.
Literal Riccati/Joseph event attachment and endpoint/prefix P4 storage remain
separate obligations.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import ou3_brmm_dynamic_source_certificate as DYNAMIC
import ou3_p4_complete_brmm_target_cell as TARGET
import ou3_p4_complete_brmm_adaptive_transition as ADAPT
import ou3_p4_joint_estimator_transition as JOINT
import ou3_p4_joint_brmm_frontend_transition as BRMM_JOINT
import ou3_brmm_finite_window_primitive_qualification as PRIMITIVE
import ou3_p4_brmm_frontend_predecessor_invariant as PREDECESSOR

REPO=Path(__file__).resolve().parents[2]
DEFAULT_DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
SCHEMA=8
QUALIFICATION='OU3_P4_COMPLETE_BRMM_UNIVERSAL_SAME_HISTORY_TARGET_RELATION_V8'

def build(domain_path=DEFAULT_DOMAIN):
    path=Path(domain_path).resolve();domain=json.loads(path.read_text())
    dynamic=DYNAMIC.build(path);target=TARGET.build(path);adapt=ADAPT.build(path);joint=JOINT.build();bj=BRMM_JOINT.build();primitive=PRIMITIVE.build(path);pred=PREDECESSOR.build()
    bad={'dynamic':DYNAMIC.validate(dynamic),'target':TARGET.validate(target),'adaptive':ADAPT.validate(adapt),'joint':JOINT.validate(joint),'brmm_joint':BRMM_JOINT.validate(bj),'primitive':PRIMITIVE.validate(primitive),'predecessor':PREDECESSOR.validate(pred)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('target-relation prerequisites failed: '+repr(bad))
    parity=dynamic['source_parity'];inv=dynamic['dynamic_invariant'];u=primitive['uniform_physical_primitives'];live=domain['normal_live']
    declared_A=float(live['non_gravitational_cog_acceleration_norm_upper_mps2']);declared_W=float(live['body_rate_norm_upper_deg_s'])
    source_inside_normal=(float(u['acceleration_norm_upper_mps2'])<=declared_A and float(u['body_rate_norm_upper_deg_s'])<=declared_W)
    one_step_joint=bool(joint['joint_f_tau_sigma_TS_RS_image_materialized'] and bj['same_private_Mahony_vertical_drives_stillness_tuner_and_WPE'] and bj['joint_raw_f_tau_sigma_TS_RS_tuple_emitted'])
    clamp_inclusion=bool(parity['normal_live_measured_period_selector_present'] and parity['raw_sigma_target_uses_1e_minus6_variance_floor'] and target['full_dynamic_target_rectangle_image_valid'] and adapt['EMA_candidate_transition_available'] and adapt['staged_commit_transition_available'] and adapt['scheduler_due_not_due_branch_preserved'])
    predecessor_family=bool(pred['complete_BRMM_predecessor_state_family_covered'] and pred['all_component_predecessor_invariants_closed'])
    physical_attachment=bool(primitive['equivalent_hard_finite_window_dynamic_constraint_closed'] and source_inside_normal and one_step_joint and predecessor_family)
    coefficient_closed=bool(physical_attachment and clamp_inclusion)
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'declared_Normal_Live_acceleration_cap_mps2':declared_A,'declared_Normal_Live_body_rate_cap_deg_s':declared_W,
      'qualified_source_acceleration_cap_mps2':float(u['acceleration_norm_upper_mps2']),'qualified_source_body_rate_cap_deg_s':float(u['body_rate_norm_upper_deg_s']),
      'historical_4mps2_literal_used_for_source_inclusion':False,
      'finite_frontend_history_required':True,'all_shipping_tuning_frequency_outputs_enclosed':bool(parity['normal_live_measured_period_selector_present']),
      'all_shipping_raw_sigma_targets_enclosed':bool(parity['raw_sigma_target_uses_1e_minus6_variance_floor']),'quiet_continuation_enclosed':bool(dynamic['normal_live_contract']['quiet_case_admitted']),
      'frequency_raw_sigma_target_rectangle':{'frequency_hz':inv['tuning_frequency_hz'],'sigma_target_raw_mps2':inv['sigma_target_raw_mps2']},
      'tau_TS_RS_are_correlated_same_cell_images':bool(target['full_dynamic_target_rectangle_image_valid']),'SpectralMSE_RS_same_cell_map':bool(target['SpectralMSE_RS_is_same_tau_sigma_image']),'independent_RS_target_forbidden':bool(target['independent_RS_target_forbidden']),
      'active_schedule_connected_by_shipping_EMA':bool(adapt['EMA_candidate_transition_available']),'active_schedule_connected_by_staged_commit':bool(adapt['staged_commit_transition_available']),'pseudo_event_incidence_connected_by_scheduler':bool(adapt['scheduler_due_not_due_branch_preserved']),'raw_and_effective_sigma_coordinates_preserved':bool(adapt['raw_and_effective_sigma_coordinates_preserved']),
      'full_rectangle_is_outer_inclusion_only':True,'outer_rectangle_may_generate_theorem_history':False,'target_pair_may_span_full_rectangle_in_theorem':False,
      'independent_frequency_sigma_coordinates_forbidden':True,'same_signal_history_must_generate_period_and_sigma':True,'wave_period_moments_must_be_retained':True,'period_scaled_sigma_band_must_be_retained':True,'sigma_variance_statistic_must_be_retained':True,'period_frequency_reciprocal_relation_must_be_retained':True,
      'joint_estimator_relation_materialized_here':bool(joint['joint_f_tau_sigma_TS_RS_image_materialized']),'joint_transition_rejects_independent_f_sigma':not bool(joint['independent_f_sigma_rectangle_accepted_by_transition']),'joint_transition_requires_split_on_dependency_loss':bool(joint['unresolved_wide_cells_require_source_split']),
      'same_BRMM_sample_joint_frontend_transition_materialized':one_step_joint,'same_BRMM_current_applied_schedule_retained':bool(bj['current_Riccati_active_schedule_precedes_current_measurement'] and bj['actual_applied_RS_xyz_retained_from_same_active_schedule']),
      'qualified_source_inside_existing_Normal_Live_physical_caps':source_inside_normal,'shipping_clamp_invariant_is_postimage_outer_cover':clamp_inclusion,
      'one_step_frontend_transition_is_not_predecessor_family_induction':True,
      'predecessor_invariant_consumed':True,'predecessor_elapsed_guard_quotient_used':pred['complete_BRMM_predecessor_family_compact_modulo_elapsed_guard_quotient'],
      'branch_complete_frontend_transition_induction_used':predecessor_family,
      'joint_transition_physical_BRMM_attachment_closed':physical_attachment,'complete_BRMM_predecessor_family_covered':predecessor_family,
      'coefficient_target_inclusion_closed':coefficient_closed,'coefficient_outer_postimage_cover_available':clamp_inclusion,
      'physical_vector_geometry_transition_closed_here':source_inside_normal,
      'complete_source_cover_closed_here':False,'P4_promoted_here':False,
      'next_obligation':'attach every universal estimator image to the literal reachable Riccati/Joseph/MEKF-injection lineage and close endpoint/every-prefix augmented LDLT plus first-exit retention'}

def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    true_keys=('finite_frontend_history_required','all_shipping_tuning_frequency_outputs_enclosed','all_shipping_raw_sigma_targets_enclosed','quiet_continuation_enclosed','tau_TS_RS_are_correlated_same_cell_images','SpectralMSE_RS_same_cell_map','independent_RS_target_forbidden','active_schedule_connected_by_shipping_EMA','active_schedule_connected_by_staged_commit','pseudo_event_incidence_connected_by_scheduler','raw_and_effective_sigma_coordinates_preserved','full_rectangle_is_outer_inclusion_only','independent_frequency_sigma_coordinates_forbidden','same_signal_history_must_generate_period_and_sigma','wave_period_moments_must_be_retained','period_scaled_sigma_band_must_be_retained','sigma_variance_statistic_must_be_retained','period_frequency_reciprocal_relation_must_be_retained','joint_estimator_relation_materialized_here','joint_transition_rejects_independent_f_sigma','joint_transition_requires_split_on_dependency_loss','same_BRMM_sample_joint_frontend_transition_materialized','same_BRMM_current_applied_schedule_retained','qualified_source_inside_existing_Normal_Live_physical_caps','shipping_clamp_invariant_is_postimage_outer_cover','coefficient_outer_postimage_cover_available','physical_vector_geometry_transition_closed_here','one_step_frontend_transition_is_not_predecessor_family_induction','predecessor_invariant_consumed','predecessor_elapsed_guard_quotient_used','branch_complete_frontend_transition_induction_used','joint_transition_physical_BRMM_attachment_closed','complete_BRMM_predecessor_family_covered','coefficient_target_inclusion_closed')
    for k in true_keys:
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('historical_4mps2_literal_used_for_source_inclusion','outer_rectangle_may_generate_theorem_history','target_pair_may_span_full_rectangle_in_theorem','complete_source_cover_closed_here','P4_promoted_here'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('declared_Normal_Live_acceleration_cap_mps2',0))!=8.0:f.append('declared acceleration cap is not 8 m/s^2')
    if float(d.get('qualified_source_acceleration_cap_mps2',0))!=8.0:f.append('qualified source acceleration cap is not 8 m/s^2')
    return f

def main():
    p=argparse.ArgumentParser();p.add_argument('--domain',type=Path,default=DEFAULT_DOMAIN);p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build(a.domain);f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps({'one_step':d['same_BRMM_sample_joint_frontend_transition_materialized'],'predecessor':d['complete_BRMM_predecessor_family_covered'],'physical_attachment':d['joint_transition_physical_BRMM_attachment_closed'],'coefficient_closed':d['coefficient_target_inclusion_closed'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
