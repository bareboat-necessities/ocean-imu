#!/usr/bin/env python3
"""Logical first-exit reduction for the corrected fresh-Live P4 entry set.

The fresh outer-wrapper handoff has centered integral-displacement error exactly
zero, so the old independent S/integral-displacement entry factor is not part of
the reachable entry graph.  The remaining product of admitted balls with S=0 is
star-shaped and is covered by the exact [0,1] radial construction.

For a hypothetical first hard-domain exit, every earlier prefix lies inside the
working domain.  Consequently a source-uniform coefficient/event cover over the
full hard domain is sufficient locally; no 600-sample nonlinear state-box
propagation is required.  The current kernel working cells validate this local
construction for captured same-cell coefficients, but they do NOT yet prove the
source-uniform P/H/R coefficient attachment.  That distinction remains
fail-closed here.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import ou3_p4_qualified_fresh_entry as FRESH
import ou3_p4_hard_entry_set as ENTRY
import ou3_p4_hard_entry_radial_ad_boxes as RADIAL
import ou3_p4_brmm_kernel_radial_working_cells as WORK
import ou3_p4_affine_hard_tube_iqc as HARD
import ou3_p4_complete_word_prefix_transport as PREFIX
import ou3_p4_complete_brmm_source_cover_transition as COVTRANS

SCHEMA=2
QUALIFICATION='OU3_P4_FRESH_ENTRY_FIRST_EXIT_REDUCTION_V2'
S_GROUP='integral_displacement_norm_m_s'

def build()->dict:
    fresh=FRESH.build();entry=ENTRY.build();radial=RADIAL.build();work=WORK.build();hard=HARD.build();prefix=PREFIX.build();cov=COVTRANS.build()
    bad={'fresh':FRESH.validate(fresh),'entry':ENTRY.validate(entry),'radial':RADIAL.validate(radial),
         'working_cells':WORK.validate(work),'hard_iqc':HARD.validate(hard),'prefix':PREFIX.validate(prefix),
         'covariance_transition':COVTRANS.validate(cov)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('first-exit prerequisites failed: '+repr(bad))
    centered=float(fresh['source_bounds']['centered_S_norm_m_s']);membership=fresh['membership_checks']
    reduced=[g for g in entry['H18_groups'] if g!=S_GROUP]
    if entry.get('A21_additional_group'):reduced.append(entry['A21_additional_group'])
    star=bool(centered==0.0 and all(membership.values()))
    captured_local=bool(radial['recursive_partition_exactly_covers_unit_interval'] and radial['rectangular_box_is_outer_AD_enclosure_only']
      and work['radial_partition_exactly_covers_unit_interval'] and work['first_exit_semantics_required_for_working_cells']
      and work['all_working_cells_and_real_event_masters_valid'] and work['all_Szero_cells_use_actual_applied_RS']
      and work['all_nonlinear_cells_retain_structured_graph_sectors'])
    transition_api=bool(cov['shipping_covariance_transition_operator_available'] and cov['theorem_transition_requires_joint_estimator_provenance'])
    target_api=bool(hard['same_graph_every_prefix_ball_target_available'] and hard['nonnegative_multiplier_Sprocedure_available'] and hard['strict_outward_LDLT_checker_available'])
    prefix_api=bool(prefix['every_prefix_transport_family_constructor_available'] and prefix['endpoint_cannot_substitute_for_prefixes'])
    logical=bool(star and radial['recursive_partition_exactly_covers_unit_interval'] and transition_api and target_api and prefix_api)
    return {'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'fresh_entry_cover_closed':bool(fresh['fresh_live_entry_cover_closed']),
      'fresh_centered_integral_displacement_norm_m_s':centered,
      'retired_independent_integral_displacement_entry_ball_reintroduced':False,
      'fresh_entry_reduced_group_order':reduced,'fresh_entry_reduced_set_star_shaped_about_zero':star,
      'continuous_zero_to_full_radial_cover_consumed':bool(radial['recursive_partition_exactly_covers_unit_interval']),
      'radial_boxes_used_only_for_outward_AD':bool(radial['rectangular_box_is_outer_AD_enclosure_only']),
      'captured_same_cell_coefficients_full_hard_domain_local_event_cover_closed':captured_local,
      'captured_kernel_working_cells_promote_source_uniform_coefficient_cover':False,
      'source_uniform_local_event_coefficient_cover_closed_here':False,
      'same_history_covariance_transition_operator_available':transition_api,
      'same_graph_every_prefix_hard_ball_target_API_available':target_api,
      'literal_every_prefix_transport_API_available':prefix_api,
      'logical_first_exit_reduction_closed':logical,
      'first_exit_argument':'before a first exit every prior prefix is inside the hard domain; a source-uniform local event cover plus strict same-graph prefix hard-ball targets contradicts the first exit',
      'explicit_600_sample_nonlinear_state_box_propagation_required':False,
      'source_uniform_endpoint_augmented_LDLT_closed_here':False,'source_uniform_every_prefix_augmented_LDLT_closed_here':False,
      'source_uniform_every_prefix_hard_domain_retention_closed_here':False,'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'remaining_obligation':'close source-uniform estimator-owned P/H/R event attachment, then endpoint/every-prefix augmented LDLT and every hard-ball target'}

def validate(d)->list[str]:
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    for k in ('fresh_entry_cover_closed','fresh_entry_reduced_set_star_shaped_about_zero','continuous_zero_to_full_radial_cover_consumed','radial_boxes_used_only_for_outward_AD','captured_same_cell_coefficients_full_hard_domain_local_event_cover_closed','same_history_covariance_transition_operator_available','same_graph_every_prefix_hard_ball_target_API_available','literal_every_prefix_transport_API_available','logical_first_exit_reduction_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('retired_independent_integral_displacement_entry_ball_reintroduced','captured_kernel_working_cells_promote_source_uniform_coefficient_cover','source_uniform_local_event_coefficient_cover_closed_here','explicit_600_sample_nonlinear_state_box_propagation_required','source_uniform_endpoint_augmented_LDLT_closed_here','source_uniform_every_prefix_augmented_LDLT_closed_here','source_uniform_every_prefix_hard_domain_retention_closed_here','P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('fresh_centered_integral_displacement_norm_m_s',1))!=0.0:f.append('fresh centered S is not zero')
    if S_GROUP in d.get('fresh_entry_reduced_group_order',[]):f.append('independent S factor remains')
    return list(dict.fromkeys(f))

def main()->int:
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'logical_first_exit':d['logical_first_exit_reduction_closed'],'captured_local':d['captured_same_cell_coefficients_full_hard_domain_local_event_cover_closed'],'source_uniform_coefficients':d['source_uniform_local_event_coefficient_cover_closed_here'],'prefix_LDLT':d['source_uniform_every_prefix_augmented_LDLT_closed_here'],'failures':f},sort_keys=True));return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
