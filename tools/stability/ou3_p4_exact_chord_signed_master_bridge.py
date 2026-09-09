#!/usr/bin/env python3
"""Bind exact finite-angle chord/reset graph to the signed Joseph P4 master.

This bridge keeps nonlinear vector geometry, accelerometer a_w/b_a coupling,
physical BIAS1 recurrence, active radial projection and finite reset transport in
one architecture. It consumes the universal full-entry finite-angle H18/A21
differential backbones, the homogeneous same-cell reset IQC D=E_theta*K*Y, and
the conditional binary32 platform/ISS arithmetic contract.

The differential backbone is NOT itself the finite-map endpoint certificate.
Endpoint/every-prefix promotion still requires the exact graph to close by the
augmented outward LDLT plus hard-domain retention at every literal prefix.
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
import ou3_p4_kalman_reset_binary32_iss as FPISS
import ou3_p4_binary32_platform_arithmetic_contract as PLATFORM

SCHEMA=3
QUALIFICATION="OU3_P4_EXACT_CHORD_SIGNED_MASTER_BRIDGE_V3"

def build()->dict:
    chord=CHORD.build();iqc=IQC.build();signed=SIGNED.build();master=MASTER.build();iss=ISSMASTER.build()
    bias=BIAS1.build();biasiss=BIASISS.build();proj=PROJ.build();coeff=COEFF.build();h18=H18FA.build();a21=A21FA.build()
    reset=RESETIQC.build();fp=FPISS.build();platform=PLATFORM.build()
    bad={
      'chord':CHORD.validate(chord),'iqc':IQC.validate(iqc),'signed':SIGNED.validate(signed),
      'master':MASTER.validate(master),'iss_master':ISSMASTER.validate(iss),
      'bias1':BIAS1.validate(bias),'bias_iss':BIASISS.validate(biasiss),
      'projection':PROJ.validate(proj),'coeff':COEFF.validate(coeff),
      'H18_finite_angle':H18FA.validate(h18),'A21_finite_angle':A21FA.validate(a21),
      'reset_iqc':RESETIQC.validate(reset),'finite_precision_iss':FPISS.validate(fp),'platform':PLATFORM.validate(platform)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('exact-chord signed-master prerequisites failed: '+repr(bad))

    same_source=all(x=='COMPLETE_BRMM_NORMAL_LIVE_WORD' for x in (
      chord['canonical_source'],iqc['canonical_source'],signed['canonical_source'],master['canonical_source'],
      iss['canonical_source'],biasiss['canonical_source'],h18['canonical_source'],a21['canonical_source'],reset['canonical_source']))
    reset_ready=bool(
      reset['reset_defect_dense_IQC_available'] and reset['reset_IQC_keeps_residual_direction_via_Etheta_K_Y']
      and reset['uniform_gain_from_homogeneous_exact_reset_sector']
      and reset['rowwise_K_used_only_for_uniform_correction_ceiling']
      and reset['rowwise_K_may_not_replace_same_cell_Dmap_in_storage']
      and all(m['same_cell_Dtheta_equals_Etheta_K_Y_required'] and m['uniform_gain_from_homogeneous_exact_reset_sector']
              and not m['endpoint_ratio_used_as_uniform_gain'] for m in reset['modes'].values()))
    graph_ready=bool(
      chord['exact_accelerometer_joint_identity_closed']
      and chord['mixed_c_cross_aw_retained_inside_chord_coordinate']
      and chord['accelerometer_bias_retained_linearly_in_same_residual']
      and iqc['exact_q_dot_q_minus_p_equality_encoded_as_two_IQCs']
      and iqc['q_u_Joseph_cross_terms_preserved_by_coordinate_not_scalarized']
      and iqc['mixed_c_cross_aw_has_explicit_augmented_coordinate']
      and signed['joint_complete_word_signed_information_composition_available']
      and signed['finite_reset_defect_remains_explicit'] and reset_ready
      and master['terminal_full_augmented_interval_LDLT_available']
      and master['same_history_quadratic_graph_sector_assembler_available']
      and iss['same_augmented_coordinate_for_state_graph_source_and_roundoff'])
    physical_ready=bool(
      bias['BIAS1_SOURCE_ADMISSION_PASS'] and bias['one_root_one_parameter_history_required']
      and bias['independent_per_sample_bias_slots_forbidden']
      and biasiss['same_w_enters_error_and_true_bias']
      and biasiss['one_physical_beta_state_carried_across_word']
      and proj['global_joint_sector_closed']
      and coeff['Kalman_and_reset_coefficient_family_outwardly_bounded'])
    backbone=bool(
      h18['full_declared_45deg_entry_covered']
      and h18['universal_H18_finite_angle_prior_free_LDLT_closed']
      and a21['universal_A21_finite_angle_first_active_full_matrix_LDLT_closed'])
    fp_ready=bool(
      fp['conditional_full_shipping_finite_precision_additive_ISS_closed']
      and platform['conditional_mathematical_execution_premise']
      and platform['P4_may_consume_as_explicit_execution_premise'])

    # Fail closed: production exact-graph endpoint/prefix certificates must set
    # these from actual outward LDLT results, never from prerequisite booleans.
    endpoint=False;prefix=False;retention=False
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
      'exact_chord_and_homogeneous_reset_graph_ready_for_augmented_master':graph_ready,
      'homogeneous_same_cell_reset_IQC_consumed':reset_ready,
      'reset_gain_uniform_over_zero_to_correction_ceiling':reset_ready,
      'joint_ISS_master_ready_for_BIAS1_and_roundoff':iss['smoke_ISS_LDLT_closed'],
      'conditional_full_binary32_additive_ISS_ready':fp_ready,
      'target_toolchain_qualified_here':platform['target_toolchain_qualified_here'],
      'physical_BIAS1_projection_and_coefficient_prerequisites_ready':physical_ready,
      'accelerometer_residual_coordinate':'y=q+u; p=[c]x(f_hat+R_hat*delta_a_w); u=R_hat*delta_a_w+delta_b_a',
      'joseph_favorable_q_u_cross_term_retained':True,
      'mixed_c_cross_aw_retained_as_explicit_coordinate':True,
      'finite_reset_cross_and_defect_terms_retained':True,
      'reset_correction_map':'D=E_theta*K*Y in the same Joseph cell',
      'physical_BIAS1_one_history_retained':True,'active_radial_projection_sector_retained':True,
      'actual_same_history_K_required_not_independent_row_box':True,
      'rowwise_K_enclosure_used_only_as_magnitude_ceiling':True,
      'packet_count_multiplier_used':False,'standalone_eta_Rinv_budget_used':False,'scalar_correction_radius_used_for_storage':False,
      'source_uniform_exact_graph_endpoint_augmented_LDLT_closed':endpoint,
      'source_uniform_exact_graph_every_prefix_augmented_LDLT_closed':prefix,
      'same_graph_every_prefix_hard_domain_retention_closed':retention,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':(
        'materialize the same-source endpoint and literal-prefix augmented matrices using the universal H18/A21 strict '
        'block, exact chord/cross-product IQCs, homogeneous reset D=E_theta*K*Y IQC, active projection sector and '
        'BIAS1/binary32 ISS supplies; run outward LDLT and derive hard-prefix retention from those same matrices')}

def validate(d):
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_canonical_source_across_chord_signed_master','universal_H18_finite_angle_differential_backbone_closed',
              'universal_A21_finite_angle_differential_backbone_closed','universal_full_entry_finite_angle_differential_backbone_closed',
              'exact_chord_and_homogeneous_reset_graph_ready_for_augmented_master','homogeneous_same_cell_reset_IQC_consumed',
              'reset_gain_uniform_over_zero_to_correction_ceiling','joint_ISS_master_ready_for_BIAS1_and_roundoff',
              'conditional_full_binary32_additive_ISS_ready','physical_BIAS1_projection_and_coefficient_prerequisites_ready',
              'joseph_favorable_q_u_cross_term_retained','mixed_c_cross_aw_retained_as_explicit_coordinate',
              'finite_reset_cross_and_defect_terms_retained','physical_BIAS1_one_history_retained','active_radial_projection_sector_retained',
              'actual_same_history_K_required_not_independent_row_box','rowwise_K_enclosure_used_only_as_magnitude_ceiling'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('finite_source_enumeration_required_for_differential_backbone','packet_count_multiplier_used','standalone_eta_Rinv_budget_used',
              'scalar_correction_radius_used_for_storage','source_uniform_exact_graph_endpoint_augmented_LDLT_closed',
              'source_uniform_exact_graph_every_prefix_augmented_LDLT_closed','same_graph_every_prefix_hard_domain_retention_closed',
              'P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
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
      'reset_iqc':d['homogeneous_same_cell_reset_IQC_consumed'],'binary32_iss':d['conditional_full_binary32_additive_ISS_ready'],
      'H18_pivot':d['universal_H18_worst_LDLT_pivot_lower'],'A21_ba_margin':d['universal_A21_first_active_ba_margin_lower'],
      'endpoint':d['source_uniform_exact_graph_endpoint_augmented_LDLT_closed'],'prefix':d['source_uniform_exact_graph_every_prefix_augmented_LDLT_closed'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
