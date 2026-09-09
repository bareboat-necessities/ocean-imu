#!/usr/bin/env python3
"""Bind exact chord/reset/projection/source/roundoff graph to the P4 master.

The bridge keeps finite-angle vector geometry, accelerometer a_w/b_a coupling,
physical BIAS1 recurrence, active radial projection and exact finite reset in one
same-history architecture.  Reset is deliberately parameterized: its correction
domain must be proved from the SAME affine augmented map D=E_theta*K*Y, never
from an independent rowwise K box.

This remains a prerequisite bridge.  The production certificate must separately
close (1) every same-cell correction domain needed by reset, (2) endpoint
augmented LDLT, (3) every literal-prefix augmented LDLT, and (4) every-prefix
hard/working-domain retention.  No prerequisite boolean can promote P4.
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
import ou3_p4_projection_sector as PROJ
import ou3_p4_rowwise_coefficient_enclosure_fast as COEFF
import ou3_p4_finite_angle_h18_universal as H18FA
import ou3_p4_finite_angle_a21_universal as A21FA
import ou3_p4_reset_graph_iqc as RESETIQC
import ou3_p4_affine_hard_tube_iqc as HARDIQC
import ou3_p4_kalman_reset_binary32_iss as FPISS
import ou3_p4_binary32_platform_arithmetic_contract as PLATFORM

SCHEMA=4
QUALIFICATION="OU3_P4_EXACT_CHORD_SIGNED_MASTER_BRIDGE_V4"

def build()->dict:
    chord=CHORD.build();iqc=IQC.build();signed=SIGNED.build();master=MASTER.build();iss=ISSMASTER.build()
    bias=BIAS1.build();biasiss=BIASISS.build();proj=PROJ.build();coeff=COEFF.build();h18=H18FA.build();a21=A21FA.build()
    reset=RESETIQC.build();hard=HARDIQC.build();fp=FPISS.build();platform=PLATFORM.build()
    bad={
      'chord':CHORD.validate(chord),'iqc':IQC.validate(iqc),'signed':SIGNED.validate(signed),
      'master':MASTER.validate(master),'iss_master':ISSMASTER.validate(iss),
      'bias1':BIAS1.validate(bias),'bias_iss':BIASISS.validate(biasiss),
      'projection':PROJ.validate(proj),'coeff':COEFF.validate(coeff),
      'H18_finite_angle':H18FA.validate(h18),'A21_finite_angle':A21FA.validate(a21),
      'reset_iqc':RESETIQC.validate(reset),'hard_iqc':HARDIQC.validate(hard),
      'finite_precision_iss':FPISS.validate(fp),'platform':PLATFORM.validate(platform)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('exact-chord signed-master prerequisites failed: '+repr(bad))

    same_source=all(x=='COMPLETE_BRMM_NORMAL_LIVE_WORD' for x in (
      chord['canonical_source'],iqc['canonical_source'],signed['canonical_source'],master['canonical_source'],
      iss['canonical_source'],biasiss['canonical_source'],h18['canonical_source'],a21['canonical_source'],reset['canonical_source']))
    reset_primitive=bool(
      reset['parameterized_reset_defect_dense_IQC_available']
      and reset['reset_IQC_keeps_residual_direction_via_Etheta_K_Y']
      and reset['production_delta_must_be_certified_from_same_Dmap']
      and reset['affine_correction_domain_target_available']
      and reset['rowwise_K_correction_domain_forbidden']
      and reset['parameterized_homogeneous_reset_sector_consumed'])
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
      and reset_primitive and hard_ready
      and master['terminal_full_augmented_interval_LDLT_available']
      and master['same_history_quadratic_graph_sector_assembler_available']
      and iss['same_augmented_coordinate_for_state_graph_source_and_roundoff'])
    physical_ready=bool(
      bias['BIAS1_SOURCE_ADMISSION_PASS'] and bias['one_root_one_parameter_history_required']
      and bias['independent_per_sample_bias_slots_forbidden']
      and biasiss['same_w_enters_error_and_true_bias']
      and biasiss['one_physical_beta_state_carried_across_word']
      and proj['global_joint_sector_closed']
      and coeff['Joseph_gain_family_outwardly_bounded']
      and coeff['rowwise_K_reset_correction_domain_forbidden']
      and coeff['reset_coefficient_family_requires_same_graph_correction_domain'])
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

    # Fail closed until a production source-correlated augmented certificate
    # supplies actual outward matrices and pivots for all four obligations.
    correction_domain=False;endpoint=False;prefix=False;retention=False
    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_canonical_source_across_chord_signed_master':same_source,
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
      'affine_hard_entry_correction_and_prefix_IQC_ready':hard_ready,
      'same_graph_correction_domain_certificate_required':True,
      'rowwise_K_reset_domain_forbidden':True,
      'joint_ISS_master_ready_for_BIAS1_and_roundoff':iss['smoke_ISS_LDLT_closed'],
      'conditional_full_binary32_additive_ISS_ready':fp_ready,
      'target_toolchain_qualified_here':platform['target_toolchain_qualified_here'],
      'physical_BIAS1_projection_and_Joseph_prerequisites_ready':physical_ready,
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
      'source_uniform_exact_graph_endpoint_augmented_LDLT_closed':endpoint,
      'source_uniform_exact_graph_every_prefix_augmented_LDLT_closed':prefix,
      'same_graph_every_prefix_hard_domain_retention_closed':retention,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':(
        'materialize source-correlated endpoint/literal-prefix augmented matrices; first certify each reset correction '
        'domain delta^2*h^2-||Etheta*K*Y*z||^2>=0 from the same hard-entry/graph IQCs, then instantiate mu_R(delta), '
        'close endpoint and every-prefix outward LDLT with chord/projection/BIAS1/binary32 supplies, and certify every '
        'prefix physical-coordinate target from those same affine matrices')}

def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_canonical_source_across_chord_signed_master','universal_H18_finite_angle_differential_backbone_closed',
              'universal_A21_finite_angle_differential_backbone_closed','universal_full_entry_finite_angle_differential_backbone_closed',
              'exact_chord_parameterized_reset_affine_graph_ready_for_augmented_master','parameterized_same_cell_reset_IQC_primitive_ready',
              'affine_hard_entry_correction_and_prefix_IQC_ready','same_graph_correction_domain_certificate_required',
              'rowwise_K_reset_domain_forbidden','joint_ISS_master_ready_for_BIAS1_and_roundoff',
              'conditional_full_binary32_additive_ISS_ready','physical_BIAS1_projection_and_Joseph_prerequisites_ready',
              'joseph_favorable_q_u_cross_term_retained','mixed_c_cross_aw_retained_as_explicit_coordinate',
              'finite_reset_cross_and_defect_terms_retained','physical_BIAS1_one_history_retained','active_radial_projection_sector_retained',
              'actual_same_history_K_required_not_independent_row_box','rowwise_K_enclosure_used_only_as_finite_precision_magnitude_ceiling'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('finite_source_enumeration_required_for_differential_backbone','packet_count_multiplier_used','standalone_eta_Rinv_budget_used',
              'scalar_correction_radius_used_for_storage','source_uniform_same_graph_correction_domain_closed',
              'source_uniform_exact_graph_endpoint_augmented_LDLT_closed','source_uniform_exact_graph_every_prefix_augmented_LDLT_closed',
              'same_graph_every_prefix_hard_domain_retention_closed','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('full_declared_entry_chord_information_retention_lower',0))<=0.85:f.append('chord retention lost')
    if float(d.get('finite_angle_differential_information_retention_lower',0))<=0.64:f.append('differential retention lost')
    for k in ('universal_H18_worst_LDLT_pivot_lower','universal_A21_first_active_ba_margin_lower'):
        if float(d.get(k,0))<=0:f.append(k+' not positive')
    return f

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'backbone':d['universal_full_entry_finite_angle_differential_backbone_closed'],
      'reset_primitive':d['parameterized_same_cell_reset_IQC_primitive_ready'],'hard_iqc':d['affine_hard_entry_correction_and_prefix_IQC_ready'],
      'binary32_iss':d['conditional_full_binary32_additive_ISS_ready'],'H18_pivot':d['universal_H18_worst_LDLT_pivot_lower'],
      'A21_ba_margin':d['universal_A21_first_active_ba_margin_lower'],'correction_domain':d['source_uniform_same_graph_correction_domain_closed'],
      'endpoint':d['source_uniform_exact_graph_endpoint_augmented_LDLT_closed'],'prefix':d['source_uniform_exact_graph_every_prefix_augmented_LDLT_closed'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
