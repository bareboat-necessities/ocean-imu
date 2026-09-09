#!/usr/bin/env python3
"""Single fail-closed conditional mathematical P4 promotion gate.

Structural proof machinery may validate while closure booleans remain false.
P4 promotion requires every mathematical closure bit simultaneously; an open
obligation is reported as a blocker, not misclassified as malformed tooling.
Target-toolchain/device qualification remains separate from mathematical P4.

The P4 theorem family is intentionally narrow and fail-closed. A promotion is
valid only for the deployed COMPLETE BRMM source with all physical accelerometer
bias families BIAS0/BIAS1/BIAS2 and the literal signal-generated adaptive
coefficient history. In particular tau, sigma, T_S and R_S are not independent
proof coordinates: tau and sigma are estimated from the same admitted input
signal history, T_S is derived from tau, R_S is the deployed dependent image of
(tau,sigma,T_S), and candidate/active EMA plus staged commit/scheduler semantics
must be retained through each actual Riccati/Joseph event.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_bias1_family as BIAS1
import ou3_p4_p3_execution_admission as P3A
import ou3_p4_exact_chord_signed_master_bridge as BRIDGE
import ou3_p4_reset_domain_binding as RESETBIND
import ou3_p4_kalman_reset_binary32_iss as FP
import ou3_p4_complete_brmm_source_cover_contract as SOURCE

QUALIFICATION='OU3_P4_FINAL_FAIL_CLOSED_GATE_V6'
REQUIRED_BIAS_FAMILIES=('BIAS0','BIAS1','BIAS2')


def build():
    e=ENTRY.build();b=BIAS1.build();p=P3A.build();g=BRIDGE.build();rb=RESETBIND.build();fp=FP.build();src=SOURCE.build()
    bad={'entry':ENTRY.validate(e),'bias1':BIAS1.validate(b),'p3':P3A.validate(p),
         'bridge':BRIDGE.validate(g),'reset_binding':RESETBIND.validate(rb),'fp':FP.validate(fp),
         'source_contract':SOURCE.validate(src)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('final P4 prerequisite validation failed: '+repr(bad))

    admissions=bool(e['P4_may_use_as_regional_entry_hypothesis'] and b['BIAS1_SOURCE_ADMISSION_PASS'] and p['P3_SCOPED_EXECUTION_ADMISSION_PASS'])
    backbone=bool(g['universal_full_entry_finite_angle_differential_backbone_closed'])
    graph_ready=bool(
      g['exact_chord_parameterized_reset_affine_graph_ready_for_augmented_master']
      and g['parameterized_same_cell_reset_IQC_primitive_ready']
      and g['affine_hard_entry_correction_and_prefix_IQC_ready']
      and g['joint_ISS_master_ready_for_BIAS1_and_roundoff']
      and g['physical_BIAS1_projection_and_Joseph_prerequisites_ready'])

    # The existing branch has a genuine BIAS1 conditional graph. BIAS0 and
    # BIAS2 are deliberately not inferred from it. They require their own
    # source-uniform same-history certificates before P4 may promote.
    bias_family_closed={
      'BIAS0':False,
      'BIAS1':bool(admissions and graph_ready),
      'BIAS2':False,
    }
    all_bias_families_closed=all(bias_family_closed[x] for x in REQUIRED_BIAS_FAMILIES)

    adaptive_source_contract=bool(
      src['complete_BRMM_single_history_required']
      and src['same_history_drives_translation_rotation_frontend_tuner_geometry']
      and src['theorem_cell_coefficients_must_come_from_joint_estimator_image']
      and src['same_signal_history_must_generate_period_and_sigma']
      and src['period_frequency_reciprocal_relation_must_be_retained']
      and src['period_scaled_sigma_band_must_be_retained']
      and src['sigma_variance_statistic_must_be_retained']
      and src['tau_TS_SpectralMSE_RS_same_cell_relation_available']
      and src['raw_effective_sigma_split_preserved']
      and src['S_zero_R_derived_from_same_active_schedule']
      and src['actual_applied_RS_provenance_required']
      and src['independent_frequency_sigma_coordinates_forbidden']
      and src['independent_tuner_RS_schedule_forbidden'])

    # Materialization/coverage, not just API declarations, is required.
    adaptive_source_cover_closed=bool(
      adaptive_source_contract
      and src['joint_estimator_relation_materialized']
      and src['joint_estimator_physical_BRMM_attachment_closed']
      and src['SOURCE_UNIFORM_COMPLETE_BRMM_COVER_CLOSED'])

    correction=bool(rb['source_uniform_same_graph_correction_domain_closed'])
    reset_iqc=bool(rb['homogeneous_same_cell_reset_IQC_consumed'])
    reset_gain=bool(rb['reset_gain_uniform_over_zero_to_correction_ceiling'])
    endpoint=bool(g['source_uniform_exact_graph_endpoint_augmented_LDLT_closed'])
    prefixes=bool(g['source_uniform_exact_graph_every_prefix_augmented_LDLT_closed'])
    hard_prefix=bool(g['same_graph_every_prefix_hard_domain_retention_closed'])
    arithmetic=bool(fp['full_shipping_Kalman_reset_finite_precision_enclosure_closed_conditionally'] and fp['additive_ISS_channel_complete_for_conditional_P4'])
    platform_qualified=bool(fp['deployment_finite_precision_qualification_closed'])

    motion=bool(admissions and backbone and graph_ready and all_bias_families_closed and adaptive_source_cover_closed
                and correction and reset_iqc and reset_gain and endpoint and prefixes and hard_prefix and arithmetic)
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','declared_domain_shrunk':False,
      'filter_changed':False,'quality_gates_changed':False,'P3_delta':1e-18,
      'required_bias_families':list(REQUIRED_BIAS_FAMILIES),
      'bias_family_source_uniform_same_history_closed':bias_family_closed,
      'all_required_bias_families_closed':all_bias_families_closed,
      'BRMM_required_for_any_P4_proof':True,
      'signal_generated_tau_sigma_required':True,
      'tau_sigma_TS_to_RS_interdependence_required':True,
      'adaptation_schedule_candidate_active_commit_scheduler_required':True,
      'independent_tau_sigma_TS_RS_boxes_can_promote':False,
      'independent_frequency_sigma_coordinates_can_promote':False,
      'adaptive_same_signal_source_contract_ready':adaptive_source_contract,
      'adaptive_same_signal_source_uniform_cover_closed':adaptive_source_cover_closed,
      'hard_entry_set_admitted_without_covariance_membership':admissions,
      'universal_full_entry_finite_angle_differential_backbone_closed':backbone,
      'H18_finite_angle_worst_LDLT_pivot_lower':g['universal_H18_worst_LDLT_pivot_lower'],
      'A21_first_active_ba_margin_lower':g['universal_A21_first_active_ba_margin_lower'],
      'exact_chord_projection_BIAS1_same_history_graph_ready':graph_ready,
      'same_cell_Joseph_reset_domain_binding_consumed':True,
      'source_uniform_same_graph_correction_domain_closed':correction,
      'homogeneous_same_cell_reset_IQC_consumed':reset_iqc,
      'reset_gain_uniform_over_zero_to_correction_ceiling':reset_gain,
      'reset_binding_modes':rb['modes'],
      'source_uniform_endpoint_augmented_LDLT_closed':endpoint,
      'source_uniform_every_prefix_augmented_LDLT_closed':prefixes,
      'source_uniform_every_prefix_hard_domain_retention_closed':hard_prefix,
      'conditional_full_shipping_finite_precision_additive_ISS_closed':arithmetic,
      'target_toolchain_finite_precision_qualified':platform_qualified,
      'P4_DEPLOYMENT_PASS':bool(motion and platform_qualified),
      'point_capture_can_promote':False,'rowwise_coefficient_boxes_can_promote':False,'differential_backbone_alone_can_promote':False,
      'P4_MOTION_PASS':motion,'P4_PASS':motion,'P5_MAY_START':motion,
      'remaining_mathematical_P4_blockers':[x for x,ok in (
        ('source-uniform BIAS0 same-history physical-driver family',bias_family_closed['BIAS0']),
        ('source-uniform BIAS1 same-history physical-driver family',bias_family_closed['BIAS1']),
        ('source-uniform BIAS2 same-history physical-driver family',bias_family_closed['BIAS2']),
        ('COMPLETE BRMM same-signal estimator/source cover with signal-derived tau and sigma',adaptive_source_cover_closed),
        ('interdependent (tau,sigma,T_S)->R_S with literal adaptation/commit/scheduler history',adaptive_source_cover_closed),
        ('source-uniform same-cell Joseph correction/reset domain',correction and reset_iqc and reset_gain),
        ('source-uniform exact-graph endpoint augmented LDLT',endpoint),
        ('source-uniform exact-graph every-prefix augmented LDLT',prefixes),
        ('same exact graph every-prefix hard-domain retention',hard_prefix)) if not ok],
      'remaining_deployment_blockers':[x for x,ok in (
        ('target-toolchain qualification of the declared binary32 LDLT/libm forward postconditions',platform_qualified),) if not ok]}


def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    if d.get('canonical_source')!='COMPLETE_BRMM_NORMAL_LIVE_WORD':f.append('source changed')
    if float(d.get('P3_delta',0))!=1e-18:f.append('P3 delta changed')
    if tuple(d.get('required_bias_families',()))!=REQUIRED_BIAS_FAMILIES:f.append('required bias family set changed')
    if set(d.get('bias_family_source_uniform_same_history_closed',{}))!=set(REQUIRED_BIAS_FAMILIES):f.append('bias family closure map incomplete')
    for k in ('BRMM_required_for_any_P4_proof','signal_generated_tau_sigma_required',
              'tau_sigma_TS_to_RS_interdependence_required','adaptation_schedule_candidate_active_commit_scheduler_required',
              'adaptive_same_signal_source_contract_ready','hard_entry_set_admitted_without_covariance_membership',
              'universal_full_entry_finite_angle_differential_backbone_closed','exact_chord_projection_BIAS1_same_history_graph_ready',
              'same_cell_Joseph_reset_domain_binding_consumed','conditional_full_shipping_finite_precision_additive_ISS_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('declared_domain_shrunk','filter_changed','quality_gates_changed','point_capture_can_promote',
              'rowwise_coefficient_boxes_can_promote','differential_backbone_alone_can_promote',
              'independent_tau_sigma_TS_RS_boxes_can_promote','independent_frequency_sigma_coordinates_can_promote'):
        if d.get(k) is not False:f.append(k+' not false')
    for k in ('H18_finite_angle_worst_LDLT_pivot_lower','A21_first_active_ba_margin_lower'):
        if float(d.get(k,0))<=0:f.append(k+' not positive')

    required=all(bool(d[k]) for k in (
      'hard_entry_set_admitted_without_covariance_membership','universal_full_entry_finite_angle_differential_backbone_closed',
      'exact_chord_projection_BIAS1_same_history_graph_ready','all_required_bias_families_closed',
      'adaptive_same_signal_source_uniform_cover_closed','source_uniform_same_graph_correction_domain_closed',
      'homogeneous_same_cell_reset_IQC_consumed','reset_gain_uniform_over_zero_to_correction_ceiling',
      'source_uniform_endpoint_augmented_LDLT_closed','source_uniform_every_prefix_augmented_LDLT_closed',
      'source_uniform_every_prefix_hard_domain_retention_closed','conditional_full_shipping_finite_precision_additive_ISS_closed'))
    if bool(d.get('P4_MOTION_PASS'))!=required:f.append('P4 motion promotion is not exact conjunction')
    if bool(d.get('P4_PASS'))!=required or bool(d.get('P5_MAY_START'))!=required:f.append('P4/P5 promotion mismatch')
    if bool(d.get('P4_DEPLOYMENT_PASS'))!=(required and bool(d.get('target_toolchain_finite_precision_qualified'))):f.append('deployment promotion mismatch')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f;a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(json.dumps(d,sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
