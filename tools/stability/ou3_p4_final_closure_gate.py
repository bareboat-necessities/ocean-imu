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
import ou3_p4_bias0_family as BIAS0
import ou3_p4_bias1_family as BIAS1
import ou3_p4_bias2_family as BIAS2
import ou3_p4_bias_family_joint_iss_supply as BIASISS
import ou3_p4_p3_execution_admission as P3A
import ou3_p4_exact_chord_signed_master_bridge as BRIDGE
import ou3_p4_reset_domain_binding as RESETBIND
import ou3_p4_kalman_reset_binary32_iss as FP
import ou3_p4_complete_brmm_source_cover_contract as SOURCE
import ou3_brmm_infinite_continuation as INFINITE

QUALIFICATION='OU3_P4_FINAL_FAIL_CLOSED_GATE_V9'
REQUIRED_BIAS_FAMILIES=('BIAS0','BIAS1','BIAS2')


def build():
    e=ENTRY.build();b0=BIAS0.build();b=BIAS1.build();b2=BIAS2.build();bs=BIASISS.build()
    p=P3A.build();g=BRIDGE.build();rb=RESETBIND.build();fp=FP.build();src=SOURCE.build();infinite=INFINITE.build()
    bad={'entry':ENTRY.validate(e),'bias0':BIAS0.validate(b0),'bias1':BIAS1.validate(b),
         'bias2':BIAS2.validate(b2),'bias_family_supply':BIASISS.validate(bs),'p3':P3A.validate(p),
         'bridge':BRIDGE.validate(g),'reset_binding':RESETBIND.validate(rb),'fp':FP.validate(fp),
         'source_contract':SOURCE.validate(src),'indefinite_source':INFINITE.validate(infinite)}
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

    # Every family is admitted and lifted from its OWN physical-driver module.
    # Projection/Joseph prerequisites are family-parametric; positive BIAS2
    # separation is optional for this bounded-bias/practical-motion objective.
    family_admission={
      'BIAS0':bool(b0['BIAS0_SOURCE_ADMISSION_PASS']),
      'BIAS1':bool(b['BIAS1_SOURCE_ADMISSION_PASS']),
      'BIAS2':bool(b2['BIAS2_SOURCE_ADMISSION_PASS']),
    }
    family_supply_materialized={
      x:bool(bs['each_family_pushed_through_deployed_24state_event_lift']
             and bs['family_supply'][x]['event_lift']['same_w_column_shared_by_error_and_truth'])
      for x in REQUIRED_BIAS_FAMILIES}
    # The bridge's projection/Joseph prerequisites are family-parametric: the
    # radial projection map carries no driver term, the Joseph/reset gains come
    # from the reachable P/H/R cell, and the three families share one true-bias
    # envelope.  Each family therefore reaches the same-history graph on its own
    # certificate rather than by inheriting BIAS1's.
    family_graph_ready={x:bool(graph_ready and g['physical_projection_and_Joseph_prerequisites_ready_per_family'][x])
                        for x in REQUIRED_BIAS_FAMILIES}
    # BIAS2 admits phi_true=1.  The declared objective is bounded bias error
    # plus practical ISS of the other 18, and the bias half comes from the
    # closed radial projection sector, which needs no relaxation root and no
    # separation constant.  mu_sep would only sharpen motion gains, so it is an
    # optional sharpener here and is reported as such rather than as a blocker.
    bias2_separation_closed=bool(b2['source_uniform_separation_closed'])
    bias2_separation_required=bool(b2['separation_sector_required_for_bounded_bias_objective'])
    bias_family_closed={
      x:bool(admissions and family_admission[x] and family_supply_materialized[x] and family_graph_ready[x])
      for x in REQUIRED_BIAS_FAMILIES}
    if bias2_separation_required and not bias2_separation_closed:
        bias_family_closed['BIAS2']=False
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

    # A hypothetical Cartesian product is not the fresh shipping handoff.
    # The exact shared-origin lemma removes its independent S factor, but the
    # remaining reachable attitude/P/tuner/scheduler cover must still be built.
    fresh_entry=e['fresh_live_entry_graph']
    fresh_entry_cover=bool(
      fresh_entry['source_uniform_reachable_other_entry_coordinates_closed']
      and adaptive_source_cover_closed)

    correction=bool(rb['source_uniform_same_graph_correction_domain_closed'])
    reset_iqc=bool(rb['homogeneous_same_cell_reset_IQC_consumed'])
    reset_gain=bool(rb['reset_gain_uniform_over_zero_to_correction_ceiling'])
    endpoint=bool(g['source_uniform_exact_graph_endpoint_augmented_LDLT_closed'])
    prefixes=bool(g['source_uniform_exact_graph_every_prefix_augmented_LDLT_closed'])
    hard_prefix=bool(g['same_graph_every_prefix_hard_domain_retention_closed'])
    arithmetic=bool(fp['full_shipping_Kalman_reset_finite_precision_enclosure_closed_conditionally'] and fp['additive_ISS_channel_complete_for_conditional_P4'])
    platform_qualified=bool(fp['deployment_finite_precision_qualification_closed'])

    # A short-window increment bound is not an indefinite S/forcing bound.
    # The exact quiet-source obstruction is a mathematical scope failure, not
    # an interval failure or a deployment qualification requirement.
    infinite_source_compatible=not infinite['bounded_all18_indefinite_target_refuted_under_finite_window_definition']

    motion=bool(infinite_source_compatible and fresh_entry_cover and admissions and backbone and graph_ready and all_bias_families_closed and adaptive_source_cover_closed
                and correction and reset_iqc and reset_gain and endpoint and prefixes and hard_prefix and arithmetic)
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','declared_domain_shrunk':False,
      'filter_changed':False,'quality_gates_changed':False,'P3_delta':1e-18,
      'indefinite_source_necessary_condition':infinite,
      'indefinite_source_target_compatible':infinite_source_compatible,
      'required_bias_families':list(REQUIRED_BIAS_FAMILIES),
      'bias_family_source_admission':family_admission,
      'bias_family_joint_ISS_supply_materialized':family_supply_materialized,
      'bias_family_same_history_graph_ready':family_graph_ready,
      'bias_family_joint_supply_norm_upper_mps2':{x:bs['family_supply'][x]['joint_supply_norm_upper_per_prediction_mps2'] for x in REQUIRED_BIAS_FAMILIES},
      'BIAS2_uniform_separation_sector_closed':bias2_separation_closed,
      'BIAS2_separation_required_for_bounded_bias_objective':bias2_separation_required,
      'bias_error_compactness_upper_mps2':g['bias_error_compactness_upper_mps2'],
      'projection_and_Joseph_prerequisites_are_family_parametric':bool(g['projection_and_Joseph_prerequisites_are_family_parametric']),
      'BIAS0_or_BIAS2_inferred_from_BIAS1':False,
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
      'fresh_live_entry_graph':fresh_entry,
      'fresh_live_entry_cover_closed':fresh_entry_cover,
      'legacy_independent_S_box_can_promote_fresh_entry':False,
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
        ('B: finite-window BRMM primitives admit unbounded centered-S ambiguity; an additional S/forcing qualification would be E until proved',infinite_source_compatible),
        ('reachable fresh-Live correlated entry including prior-frequency timeout branch',fresh_entry_cover),
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
    for table in ('bias_family_source_uniform_same_history_closed','bias_family_source_admission','bias_family_joint_ISS_supply_materialized','bias_family_same_history_graph_ready'):
        if set(d.get(table,{}))!=set(REQUIRED_BIAS_FAMILIES):f.append(table+' incomplete')
    if not all(bool(v) for v in d.get('bias_family_source_admission',{}).values()):f.append('every bias family must be admitted from its own module')
    if not all(bool(v) for v in d.get('bias_family_joint_ISS_supply_materialized',{}).values()):f.append('every bias family must have a materialized joint ISS supply')
    for k in ('BRMM_required_for_any_P4_proof','signal_generated_tau_sigma_required',
              'tau_sigma_TS_to_RS_interdependence_required','adaptation_schedule_candidate_active_commit_scheduler_required',
              'adaptive_same_signal_source_contract_ready','hard_entry_set_admitted_without_covariance_membership',
              'universal_full_entry_finite_angle_differential_backbone_closed','exact_chord_projection_BIAS1_same_history_graph_ready',
              'same_cell_Joseph_reset_domain_binding_consumed','conditional_full_shipping_finite_precision_additive_ISS_closed',
              'projection_and_Joseph_prerequisites_are_family_parametric'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('declared_domain_shrunk','filter_changed','quality_gates_changed','point_capture_can_promote',
              'rowwise_coefficient_boxes_can_promote','differential_backbone_alone_can_promote',
              'independent_tau_sigma_TS_RS_boxes_can_promote','independent_frequency_sigma_coordinates_can_promote',
              'BIAS0_or_BIAS2_inferred_from_BIAS1','BIAS2_separation_required_for_bounded_bias_objective',
              'legacy_independent_S_box_can_promote_fresh_entry'):
        if d.get(k) is not False:f.append(k+' not false')
    for k in ('H18_finite_angle_worst_LDLT_pivot_lower','A21_first_active_ba_margin_lower'):
        if float(d.get(k,0))<=0:f.append(k+' not positive')

    infinite=d.get('indefinite_source_necessary_condition',{})
    f.extend('indefinite source: '+x for x in INFINITE.validate(infinite))
    compatible=not infinite.get('bounded_all18_indefinite_target_refuted_under_finite_window_definition',True)
    if d.get('indefinite_source_target_compatible') is not compatible:f.append('indefinite source compatibility changed')

    required=all(bool(d[k]) for k in (
      'indefinite_source_target_compatible',
      'fresh_live_entry_cover_closed','hard_entry_set_admitted_without_covariance_membership','universal_full_entry_finite_angle_differential_backbone_closed',
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
