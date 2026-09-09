#!/usr/bin/env python3
"""Bind exact chord/reset/projection/source/roundoff graph to the P4 master.

The bridge keeps finite-angle vector geometry, accelerometer a_w/b_a coupling,
physical BIAS1 recurrence, active radial projection and exact finite reset in one
same-history architecture. Reset correction-domain validity is consumed only
from the production same-cell Joseph binding; an independent rowwise K box is
never a theorem input.

The nonlinear/reset ledger now also consumes the exact reduced same-cell
Joseph+reset identity. This eliminates S^-1 and posterior precision from the
nonlinear/reset part and leaves K only through the same-cell graph d=Kq. The
exact covariance-frame reset identity is retained as an independent algebraic
cross-check and a smaller lifted representation for the terminal certificate.

This remains a prerequisite bridge. The production certificate must separately
close endpoint augmented LDLT, every literal-prefix augmented LDLT, and every-
prefix hard/working-domain retention. No prerequisite boolean can promote P4.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path

import ou3_p4_complete_brmm_exact_chord_joint_coordinate as CHORD
import ou3_p4_complete_brmm_exact_chord_iqc as IQC
import ou3_p4_complete_brmm_signed_information_ledger as SIGNED
import ou3_p4_complete_brmm_joint_sector_master as MASTER
import ou3_p4_joint_iss_augmented_master as ISSMASTER
import ou3_p4_bias1_joint_iss_supply as BIASISS
import ou3_p4_bias1_family as BIAS1
import ou3_p4_bias_family_joint_iss_supply as ALLBIAS
import ou3_p4_projection_sector as PROJ
import ou3_p4_rowwise_coefficient_enclosure_fast as COEFF
import ou3_p4_finite_angle_h18_universal as H18FA
import ou3_p4_finite_angle_a21_universal as A21FA
import ou3_p4_reset_graph_iqc as RESETIQC
import ou3_p4_reset_domain_binding as RESETBIND
import ou3_p4_affine_hard_tube_iqc as HARDIQC
import ou3_p4_kalman_reset_binary32_iss as FPISS
import ou3_p4_binary32_platform_arithmetic_contract as PLATFORM
import ou3_p4_joseph_reset_reduced_signed_identity as REDUCED
import ou3_p4_exact_reset_covariance_frame as RESETFRAME
import ou3_p4_universal_finite_angle_strict_blocks as STRICT

SCHEMA=7
QUALIFICATION="OU3_P4_EXACT_CHORD_SIGNED_MASTER_BRIDGE_V7"

def build()->dict:
    chord=CHORD.build();iqc=IQC.build();signed=SIGNED.build();master=MASTER.build();iss=ISSMASTER.build()
    bias=BIAS1.build();biasiss=BIASISS.build();proj=PROJ.build();coeff=COEFF.build();h18=H18FA.build();a21=A21FA.build()
    reset=RESETIQC.build();resetbind=RESETBIND.build();hard=HARDIQC.build();fp=FPISS.build();platform=PLATFORM.build()
    reduced=REDUCED.build();resetframe=RESETFRAME.build();strict=STRICT.build()
    allbias=ALLBIAS.build()
    bad={
      'chord':CHORD.validate(chord),'iqc':IQC.validate(iqc),'signed':SIGNED.validate(signed),
      'master':MASTER.validate(master),'iss_master':ISSMASTER.validate(iss),
      'bias1':BIAS1.validate(bias),'bias_iss':BIASISS.validate(biasiss),
      'all_bias_family_supply':ALLBIAS.validate(allbias),
      'projection':PROJ.validate(proj),'coeff':COEFF.validate(coeff),
      'H18_finite_angle':H18FA.validate(h18),'A21_finite_angle':A21FA.validate(a21),
      'reset_iqc':RESETIQC.validate(reset),'reset_binding':RESETBIND.validate(resetbind),
      'hard_iqc':HARDIQC.validate(hard),'finite_precision_iss':FPISS.validate(fp),'platform':PLATFORM.validate(platform),
      'reduced_joseph_reset':REDUCED.validate(reduced),'reset_covariance_frame':RESETFRAME.validate(resetframe),
      'strict_blocks':STRICT.validate(strict)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('exact-chord signed-master prerequisites failed: '+repr(bad))

    same_source=all(x=='COMPLETE_BRMM_NORMAL_LIVE_WORD' for x in (
      chord['canonical_source'],iqc['canonical_source'],signed['canonical_source'],master['canonical_source'],
      iss['canonical_source'],biasiss['canonical_source'],h18['canonical_source'],a21['canonical_source'],
      reset['canonical_source'],resetbind['canonical_source'],resetframe['canonical_source'],strict['canonical_source']))
    reduced_ready=bool(
      reduced['same_shipping_P_H_R_K_cell_required']
      and reduced['same_cell_cross_identity_consumed']
      and reduced['S_inverse_eliminated_from_reduced_nonlinear_reset_ledger']
      and reduced['posterior_precision_eliminated_from_reduced_nonlinear_reset_ledger']
      and reduced['prior_precision_is_same_moving_shipping_metric']
      and reduced['K_survives_only_through_same_cell_d_equals_Kq_graph']
      and reduced['exact_rational_total_identity_residual_zero']
      and not reduced['independent_K_box_used']
      and not reduced['independent_Sinverse_box_used']
      and not reduced['posterior_precision_box_used'])
    reset_frame_ready=bool(
      resetframe['deployed_correction_cayley_parallel_to_d']
      and resetframe['same_Joseph_information_frame_retained']
      and resetframe['homogeneous_lift_available']
      and not resetframe['reset_defect_independent_port_used'])
    strict_blocks_ready=bool(
      strict['finite_angle_H18_actual_interval_matrices_materialized']
      and strict['A21_exact_first_active_direct_sum_matrices_materialized']
      and strict['same_x_cover_as_canonical_H18_completion']
      and not strict['scalar_margin_substituted_for_matrix'])
    reset_primitive=bool(
      reset['parameterized_reset_defect_dense_IQC_available']
      and reset['reset_IQC_keeps_residual_direction_via_Etheta_K_Y']
      and reset['production_delta_must_be_certified_from_same_Dmap']
      and reset['affine_correction_domain_target_available']
      and reset['rowwise_K_correction_domain_forbidden']
      and reset['parameterized_homogeneous_reset_sector_consumed'])
    reset_binding_clean=bool(
      resetbind['same_cell_Joseph_correction_domain_consumed']
      and resetbind['same_cell_parameterized_exact_reset_lift_builder_used']
      and resetbind['same_cell_d_equals_Etheta_K_y_retained']
      and not resetbind['rowwise_K_bound_used']
      and not resetbind['independent_K_box_used']
      and not resetbind['historical_rowwise_reset_lift_path_consumed']
      and not resetbind['scalar_delta_used_as_independent_storage_port'])
    hard_ready=bool(
      hard['hard_entry_ball_IQC_available']
      and hard['same_graph_correction_domain_target_available']
      and hard['same_graph_every_prefix_ball_target_available']
      and hard['nonnegative_multiplier_Sprocedure_available']
      and hard['strict_outward_LDLT_checker_available'])
    graph_ready=bool(
      chord['exact_accelerometer_joint_identity_closed']
      and chord['mixed_c_cross_aw_retained_inside_chord_coordinate']
      and chord['accelerometer_bias_retained_linearly_in_same_residual']
      and iqc['exact_q_dot_q_minus_p_equality_encoded_as_two_IQCs']
      and iqc['q_u_Joseph_cross_terms_preserved_by_coordinate_not_scalarized']
      and iqc['mixed_c_cross_aw_has_explicit_augmented_coordinate']
      and signed['joint_complete_word_signed_information_composition_available']
      and signed['finite_reset_defect_remains_explicit']
      and reduced_ready and reset_frame_ready and strict_blocks_ready
      and reset_primitive and reset_binding_clean and hard_ready
      and master['terminal_full_augmented_interval_LDLT_available']
      and master['same_history_quadratic_graph_sector_assembler_available']
      and iss['same_augmented_coordinate_for_state_graph_source_and_roundoff'])
    # The projection/Joseph prerequisites see an accelerometer-bias family only
    # through four things: that it is admitted, that it carries ONE physical
    # bias history across the word, that its driver enters e_b and b_true
    # through the same w column, and through |b_true|.  The radial projection
    # map F_R(e,beta)=beta-Pi_R(beta-e) contains no driver term at all, and the
    # Joseph/reset coefficient family comes from the reachable P/H/R cell, not
    # from the bias model.  The three declared families share one true-bias
    # envelope, so the fourth dependence is one number for all of them and the
    # prerequisites are family-parametric rather than BIAS1-specific.
    graph_prerequisites=bool(
      proj['global_joint_sector_closed']
      and coeff['Joseph_gain_family_outwardly_bounded']
      and coeff['rowwise_K_reset_correction_domain_forbidden']
      and coeff['reset_coefficient_family_requires_same_graph_correction_domain']
      and allbias['radial_projection_sector_closed_for_every_family']
      and allbias['projection_sector_sees_no_bias_driver']
      and allbias['families_share_one_true_bias_envelope']
      and allbias['bias_compactness_is_family_uniform'])
    physical_ready_family={
      name:bool(
        graph_prerequisites
        and allbias['family_supply'][name]['event_lift']['prediction_lift_accepts_family_interval']
        and allbias['family_supply'][name]['event_lift']['supply_injection_available']
        and allbias['family_supply'][name]['event_lift']['same_w_column_shared_by_error_and_truth']
        and allbias['family_supply'][name]['event_lift']['physical_factor_carried_in_truth_block']
        and not allbias['independent_per_sample_bias_state_slots_used'])
      for name in ALLBIAS.REQUIRED_BIAS_FAMILIES}
    physical_ready=bool(
      bias['BIAS1_SOURCE_ADMISSION_PASS'] and bias['one_root_one_parameter_history_required']
      and bias['independent_per_sample_bias_slots_forbidden']
      and biasiss['same_w_enters_error_and_true_bias']
      and biasiss['one_physical_beta_state_carried_across_word']
      and physical_ready_family['BIAS1'])
    backbone=bool(
      h18['full_declared_45deg_entry_covered']
      and h18['universal_H18_finite_angle_prior_free_LDLT_closed']
      and a21['universal_A21_finite_angle_first_active_full_matrix_LDLT_closed'])
    fp_ready=bool(
      fp['full_shipping_Kalman_reset_finite_precision_enclosure_closed_conditionally']
      and fp['additive_ISS_channel_complete_for_conditional_P4']
      and not fp['rowwise_K_reset_domain_used']
      and platform['conditional_mathematical_execution_premise']
      and platform['P4_may_consume_as_explicit_execution_premise'])

    correction_domain=bool(resetbind['source_uniform_same_graph_correction_domain_closed'])
    endpoint=False;prefix=False;retention=False
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_canonical_source_across_chord_signed_master':same_source,
      'reduced_same_cell_Joseph_reset_identity_ready':reduced_ready,
      'reset_covariance_frame_identity_ready':reset_frame_ready,
      'universal_strict_interval_blocks_materialized':strict_blocks_ready,
      'reduced_nonlinear_reset_ledger':'-q^T R^-1 q+eta^T R^-1 eta+q^T R^-1 H d+2 e^T J b-2 eta^T R^-1 H b+b^T J b+(Hb)^T R^-1(Hb)',
      'full_declared_entry_chord_information_retention_lower':chord['information_retention_factor_lower_full_entry'],
      'finite_angle_differential_information_retention_lower':h18['finite_angle_vector_information_retention_lower'],
      'universal_H18_finite_angle_differential_backbone_closed':h18['universal_H18_finite_angle_prior_free_LDLT_closed'],
      'universal_H18_worst_LDLT_pivot_lower':h18['worst_full_H18_LDLT_pivot_lower'],
      'universal_A21_finite_angle_differential_backbone_closed':a21['universal_A21_finite_angle_first_active_full_matrix_LDLT_closed'],
      'universal_A21_first_active_ba_margin_lower':a21['first_active_ba_margin_lower'],
      'universal_full_entry_finite_angle_differential_backbone_closed':backbone,
      'finite_source_enumeration_required_for_differential_backbone':False,
      'exact_chord_parameterized_reset_affine_graph_ready_for_augmented_master':graph_ready,
      'parameterized_same_cell_reset_IQC_primitive_ready':reset_primitive,
      'same_cell_reset_domain_binding_clean':reset_binding_clean,
      'affine_hard_entry_correction_and_prefix_IQC_ready':hard_ready,
      'same_graph_correction_domain_certificate_required':True,
      'rowwise_K_reset_domain_forbidden':True,
      'joint_ISS_master_ready_for_BIAS1_and_roundoff':iss['smoke_ISS_LDLT_closed'],
      'conditional_full_binary32_additive_ISS_ready':fp_ready,
      'target_toolchain_qualified_here':platform['target_toolchain_qualified_here'],
      'physical_BIAS1_projection_and_Joseph_prerequisites_ready':physical_ready,
      'physical_projection_and_Joseph_prerequisites_ready_per_family':physical_ready_family,
      'projection_and_Joseph_prerequisites_are_family_parametric':graph_prerequisites,
      'all_bias_families_reach_the_same_history_graph':all(physical_ready_family.values()),
      'bias_error_compactness_upper_mps2':allbias['bias_error_compactness_upper_mps2'],
      'accelerometer_residual_coordinate':'y=q+u; p=[c]x(f_hat+R_hat*delta_a_w); u=R_hat*delta_a_w+delta_b_a',
      'joseph_favorable_q_u_cross_term_retained':True,
      'mixed_c_cross_aw_retained_as_explicit_coordinate':True,
      'finite_reset_cross_and_defect_terms_retained':True,
      'reset_correction_map':'D=E_theta*K*Y in the same Joseph cell',
      'physical_BIAS1_one_history_retained':True,'active_radial_projection_sector_retained':True,
      'actual_same_history_K_required_not_independent_row_box':True,
      'rowwise_K_enclosure_used_only_as_finite_precision_magnitude_ceiling':True,
      'packet_count_multiplier_used':False,'standalone_eta_Rinv_budget_used':False,
      'scalar_correction_radius_used_for_storage':False,
      'source_uniform_same_graph_correction_domain_closed':correction_domain,
      'same_cell_reset_domain_modes':resetbind['modes'],
      'source_uniform_exact_graph_endpoint_augmented_LDLT_closed':endpoint,
      'source_uniform_exact_graph_every_prefix_augmented_LDLT_closed':prefix,
      'same_graph_every_prefix_hard_domain_retention_closed':retention,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':(
        'form the production event/word augmented quadratic using the materialized strict H18/A21 interval blocks and '
        'the reduced same-cell Joseph/reset identity; bind chord/reset/projection/BIAS1/binary32 QCs, then close endpoint '
        'and every-prefix outward LDLT and every-prefix hard-domain targets without source replay or rowwise K domains')}

def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_canonical_source_across_chord_signed_master','reduced_same_cell_Joseph_reset_identity_ready',
              'reset_covariance_frame_identity_ready','universal_strict_interval_blocks_materialized',
              'universal_H18_finite_angle_differential_backbone_closed','universal_A21_finite_angle_differential_backbone_closed',
              'universal_full_entry_finite_angle_differential_backbone_closed','exact_chord_parameterized_reset_affine_graph_ready_for_augmented_master',
              'parameterized_same_cell_reset_IQC_primitive_ready','same_cell_reset_domain_binding_clean',
              'affine_hard_entry_correction_and_prefix_IQC_ready','same_graph_correction_domain_certificate_required',
              'rowwise_K_reset_domain_forbidden','joint_ISS_master_ready_for_BIAS1_and_roundoff',
              'conditional_full_binary32_additive_ISS_ready','physical_BIAS1_projection_and_Joseph_prerequisites_ready',
              'projection_and_Joseph_prerequisites_are_family_parametric','all_bias_families_reach_the_same_history_graph',
              'joseph_favorable_q_u_cross_term_retained','mixed_c_cross_aw_retained_as_explicit_coordinate',
              'finite_reset_cross_and_defect_terms_retained','physical_BIAS1_one_history_retained','active_radial_projection_sector_retained',
              'actual_same_history_K_required_not_independent_row_box','rowwise_K_enclosure_used_only_as_finite_precision_magnitude_ceiling'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('finite_source_enumeration_required_for_differential_backbone','packet_count_multiplier_used','standalone_eta_Rinv_budget_used',
              'scalar_correction_radius_used_for_storage','source_uniform_exact_graph_endpoint_augmented_LDLT_closed',
              'source_uniform_exact_graph_every_prefix_augmented_LDLT_closed','same_graph_every_prefix_hard_domain_retention_closed',
              'P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if d.get('source_uniform_same_graph_correction_domain_closed') not in (False,True):f.append('correction-domain flag not boolean')
    if float(d.get('full_declared_entry_chord_information_retention_lower',0))<=0.85:f.append('chord retention lost')
    if float(d.get('finite_angle_differential_information_retention_lower',0))<=0.64:f.append('differential retention lost')
    for k in ('universal_H18_worst_LDLT_pivot_lower','universal_A21_first_active_ba_margin_lower'):
        if float(d.get(k,0))<=0:f.append(k+' not positive')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'backbone':d['universal_full_entry_finite_angle_differential_backbone_closed'],
      'strict_blocks':d['universal_strict_interval_blocks_materialized'],'reduced_identity':d['reduced_same_cell_Joseph_reset_identity_ready'],
      'reset_frame':d['reset_covariance_frame_identity_ready'],'correction_domain':d['source_uniform_same_graph_correction_domain_closed'],
      'endpoint':d['source_uniform_exact_graph_endpoint_augmented_LDLT_closed'],'prefix':d['source_uniform_exact_graph_every_prefix_augmented_LDLT_closed'],
      'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
