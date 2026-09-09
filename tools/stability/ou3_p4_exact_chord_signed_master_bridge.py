#!/usr/bin/env python3
"""Bind the exact finite-angle chord graph to the signed Joseph/reset P4 master.

This is the non-scalarized bridge requested by the P4 proof architecture.  It
consumes, on one declared COMPLETE_BRMM_NORMAL_LIVE_WORD family:

* the exact Cayley chord coordinate q for vector residuals;
* the accelerometer mixed coordinate r=c x (R_hat delta_a_w);
* the linear companion u=R_hat delta_a_w+delta_b_a;
* the Joseph signed identity with the favorable -(q+u)'S^-1(q+u) term intact;
* the finite reset cross/defect terms in the same moving information metric;
* the one-root/one-driver BIAS1 recurrence family; and
* the radial accelerometer-bias projection sector.

The module deliberately does NOT convert these into an event-count eta budget,
independent K boxes, or a scalar remainder radius.  It exposes the exact graph
pieces in the form required by the existing full augmented LDLT builder and
fails closed until one source-correlated outward word/prefix cell supplies the
actual endpoint/prefix matrices.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import ou3_p4_complete_brmm_exact_chord_joint_coordinate as CHORD
import ou3_p4_complete_brmm_exact_chord_iqc as IQC
import ou3_p4_complete_brmm_signed_information_ledger as SIGNED
import ou3_p4_complete_brmm_joint_sector_master as MASTER
import ou3_p4_bias1_family as BIAS1
import ou3_p4_projection_sector as PROJ
import ou3_p4_rowwise_coefficient_enclosure_fast as COEFF

SCHEMA = 1
QUALIFICATION = "OU3_P4_EXACT_CHORD_SIGNED_MASTER_BRIDGE_V1"


def build() -> dict:
    chord=CHORD.build(); iqc=IQC.build(); signed=SIGNED.build(); master=MASTER.build()
    bias=BIAS1.build(); proj=PROJ.build(); coeff=COEFF.build()
    bad={
      'chord':CHORD.validate(chord), 'iqc':IQC.validate(iqc),
      'signed':SIGNED.validate(signed), 'master':MASTER.validate(master),
      'bias1':BIAS1.validate(bias), 'projection':PROJ.validate(proj),
      'coeff':COEFF.validate(coeff),
    }
    bad={k:v for k,v in bad.items() if v}
    if bad: raise RuntimeError('exact-chord signed-master prerequisites failed: '+repr(bad))

    same_source=all(x=='COMPLETE_BRMM_NORMAL_LIVE_WORD' for x in (
      chord['canonical_source'], iqc['canonical_source'], signed['canonical_source'],
      master['canonical_source']))
    graph_ready=bool(
      chord['exact_accelerometer_joint_identity_closed']
      and chord['mixed_c_cross_aw_retained_inside_chord_coordinate']
      and chord['accelerometer_bias_retained_linearly_in_same_residual']
      and iqc['exact_q_dot_q_minus_p_equality_encoded_as_two_IQCs']
      and iqc['q_u_Joseph_cross_terms_preserved_by_coordinate_not_scalarized']
      and iqc['mixed_c_cross_aw_has_explicit_augmented_coordinate']
      and signed['joint_complete_word_signed_information_composition_available']
      and signed['finite_reset_defect_remains_explicit']
      and master['terminal_full_augmented_interval_LDLT_available']
      and master['same_history_quadratic_graph_sector_assembler_available'])
    physical_ready=bool(
      bias['BIAS1_SOURCE_ADMISSION_PASS']
      and bias['one_root_one_parameter_history_required']
      and bias['independent_per_sample_bias_slots_forbidden']
      and proj['global_joint_sector_closed']
      and coeff['Kalman_and_reset_coefficient_family_outwardly_bounded'])

    # These remain false until a single outward cell carries the actual
    # source-dependent P,H,R,K/reset/tuner continuation through every literal
    # prefix.  A norm enclosure of K is intentionally insufficient.
    source_correlated_word_cell=False
    endpoint_ldlt=False
    every_prefix_ldlt=False

    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,
      'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'same_canonical_source_across_chord_signed_master':same_source,
      'full_declared_entry_information_retention_lower':chord['information_retention_factor_lower_full_entry'],
      'exact_chord_graph_ready_for_augmented_master':graph_ready,
      'physical_BIAS1_projection_and_coefficient_prerequisites_ready':physical_ready,
      'dependency_reduced_smallx_tube_used_for_same_bounds':True,
      'accelerometer_residual_coordinate':'y=q+u; p=[c]x(f_hat+R_hat*delta_a_w); u=R_hat*delta_a_w+delta_b_a',
      'joseph_favorable_q_u_cross_term_retained':True,
      'mixed_c_cross_aw_retained_as_explicit_coordinate':True,
      'finite_reset_cross_and_defect_terms_retained':True,
      'physical_BIAS1_one_history_retained':True,
      'active_radial_projection_sector_retained':True,
      'actual_same_history_K_required_not_independent_row_box':True,
      'rowwise_K_enclosure_used_only_as_magnitude_ceiling':True,
      'packet_count_multiplier_used':False,
      'standalone_eta_Rinv_budget_used':False,
      'scalar_correction_radius_used_for_storage':False,
      'source_correlated_outward_word_and_prefix_cell_materialized':source_correlated_word_cell,
      'endpoint_augmented_LDLT_closed':endpoint_ldlt,
      'every_literal_prefix_augmented_LDLT_closed':every_prefix_ldlt,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'next_obligation':(
        'materialize one outward source-correlated COMPLETE_BRMM_NORMAL_LIVE_WORD cell carrying '
        'P,H,R,K,tuner,R_S,event timing, exact chord/cross-product coordinates, BIAS1 beta/w and '
        'reset transport; feed each literal prefix to the existing full augmented interval LDLT'
      ),
    }


def validate(d: dict) -> list[str]:
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('same_canonical_source_across_chord_signed_master','exact_chord_graph_ready_for_augmented_master',
              'physical_BIAS1_projection_and_coefficient_prerequisites_ready','dependency_reduced_smallx_tube_used_for_same_bounds',
              'joseph_favorable_q_u_cross_term_retained','mixed_c_cross_aw_retained_as_explicit_coordinate',
              'finite_reset_cross_and_defect_terms_retained','physical_BIAS1_one_history_retained',
              'active_radial_projection_sector_retained','actual_same_history_K_required_not_independent_row_box',
              'rowwise_K_enclosure_used_only_as_magnitude_ceiling'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('packet_count_multiplier_used','standalone_eta_Rinv_budget_used','scalar_correction_radius_used_for_storage',
              'source_correlated_outward_word_and_prefix_cell_materialized','endpoint_augmented_LDLT_closed',
              'every_literal_prefix_augmented_LDLT_closed','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('full_declared_entry_information_retention_lower',0))<=0.85:f.append('full-entry chord retention lost')
    return f


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'chord_master_ready':d['exact_chord_graph_ready_for_augmented_master'],
                      'physical_prereqs':d['physical_BIAS1_projection_and_coefficient_prerequisites_ready'],
                      'endpoint':d['endpoint_augmented_LDLT_closed'],'prefixes':d['every_literal_prefix_augmented_LDLT_closed'],
                      'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
