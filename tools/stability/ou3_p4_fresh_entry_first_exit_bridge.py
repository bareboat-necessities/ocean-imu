#!/usr/bin/env python3
"""First-exit reduction from the corrected fresh-Live entry set to P4 prefix storage.

This certificate closes a logical gap, not the final numerical inequality.
For the fresh outer wrapper the qualified entry set is the product of the
admitted H18/A21 error balls with the centered integral-displacement coordinate
fixed to zero.  That reduced set is star-shaped about the origin.  Hence every
radial segment required by the generalized mean-value theorem is covered by the
existing exact [0,1] radial partition without reintroducing the retired
independent S/integral-displacement entry factor.

Assume, for contradiction, a first literal prefix that exits the declared hard
working domain.  Every earlier prefix is inside that domain.  Therefore the
kernel-backed radial working cells are valid outward local AD/event covers for
the event that could cause the first exit.  If the production augmented-prefix
certificate proves every hard-ball target strictly from the SAME source/event/
reset/projection/bias/binary32 graph, the alleged first exit is impossible.

Thus explicit interval propagation of the nonlinear physical state through 600
samples is neither required nor authorized.  The remaining mathematical task is
exactly the source-uniform endpoint/every-prefix augmented LDLT plus every-prefix
hard-face targets.  This module does not set those booleans true.
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

SCHEMA=1
QUALIFICATION='OU3_P4_FRESH_ENTRY_FIRST_EXIT_REDUCTION_V1'
S_GROUP='integral_displacement_norm_m_s'


def build()->dict:
    fresh=FRESH.build();entry=ENTRY.build();radial=RADIAL.build();work=WORK.build();hard=HARD.build();prefix=PREFIX.build();cov=COVTRANS.build()
    bad={
      'fresh':FRESH.validate(fresh),'entry':ENTRY.validate(entry),'radial':RADIAL.validate(radial),
      'working_cells':WORK.validate(work),'hard_iqc':HARD.validate(hard),'prefix':PREFIX.validate(prefix),
      'covariance_transition':COVTRANS.validate(cov)}
    bad={k:v for k,v in bad.items() if v}
    if bad:raise RuntimeError('first-exit prerequisites failed: '+repr(bad))

    centered=float(fresh['source_bounds']['centered_S_norm_m_s'])
    membership=fresh['membership_checks']
    reduced_groups=[g for g in entry['H18_groups'] if g!=S_GROUP]
    if entry.get('A21_additional_group'):
        reduced_groups.append(entry['A21_additional_group'])
    star=bool(centered==0.0 and all(membership.values()))
    local_cover=bool(
      radial['recursive_partition_exactly_covers_unit_interval']
      and radial['rectangular_box_is_outer_AD_enclosure_only']
      and work['radial_partition_exactly_covers_unit_interval']
      and work['first_exit_semantics_required_for_working_cells']
      and work['all_working_cells_and_real_event_masters_valid']
      and work['all_Szero_cells_use_actual_applied_RS']
      and work['all_nonlinear_cells_retain_structured_graph_sectors'])
    transition_api=bool(cov['shipping_covariance_transition_operator_available'] and cov['theorem_transition_requires_joint_estimator_provenance'])
    target_api=bool(hard['same_graph_every_prefix_ball_target_available'] and hard['nonnegative_multiplier_Sprocedure_available'] and hard['strict_outward_LDLT_checker_available'])
    prefix_api=bool(prefix['every_prefix_transport_family_constructor_available'] and prefix['endpoint_cannot_substitute_for_prefixes'])
    reduction=bool(star and local_cover and transition_api and target_api and prefix_api)

    return {
      'schema':SCHEMA,'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'fresh_entry_cover_closed':bool(fresh['fresh_live_entry_cover_closed']),
      'fresh_centered_integral_displacement_norm_m_s':centered,
      'retired_independent_integral_displacement_entry_ball_reintroduced':False,
      'fresh_entry_reduced_group_order':reduced_groups,
      'fresh_entry_reduced_set_star_shaped_about_zero':star,
      'continuous_zero_to_full_radial_cover_consumed':bool(radial['recursive_partition_exactly_covers_unit_interval']),
      'radial_boxes_used_only_for_outward_AD':bool(radial['rectangular_box_is_outer_AD_enclosure_only']),
      'kernel_first_exit_working_cells_cover_full_hard_domain_locally':local_cover,
      'kernel_working_cells_claimed_reachable_without_prefix_proof':False,
      'same_history_covariance_transition_operator_available':transition_api,
      'same_graph_every_prefix_hard_ball_target_API_available':target_api,
      'literal_every_prefix_transport_API_available':prefix_api,
      'first_exit_reduction_closed':reduction,
      'first_exit_argument':(
        'if a first hard-domain exit exists, all earlier prefixes lie in the hard domain; '
        'therefore the full-hard-domain radial working cell covers the candidate exit event; '
        'a strict same-graph prefix hard-ball implication contradicts that exit'),
      'explicit_600_sample_nonlinear_state_box_propagation_required':False,
      'source_uniform_endpoint_augmented_LDLT_closed_here':False,
      'source_uniform_every_prefix_augmented_LDLT_closed_here':False,
      'source_uniform_every_prefix_hard_domain_retention_closed_here':False,
      'P4_MOTION_PASS':False,'P4_PASS':False,'P5_MAY_START':False,
      'remaining_obligation':(
        'assemble production estimator-owned literal prefix graphs, including each physical bias-family and '
        'binary32 ISS channel, and close endpoint/every-prefix outward augmented LDLT plus every hard-ball target')}


def validate(d)->list[str]:
    f=[]
    if d.get('schema')!=SCHEMA or d.get('qualification')!=QUALIFICATION:f.append('schema/qualification mismatch')
    if d.get('canonical_source')!='COMPLETE_BRMM_NORMAL_LIVE_WORD':f.append('source changed')
    for k in ('fresh_entry_cover_closed','fresh_entry_reduced_set_star_shaped_about_zero','continuous_zero_to_full_radial_cover_consumed',
              'radial_boxes_used_only_for_outward_AD','kernel_first_exit_working_cells_cover_full_hard_domain_locally',
              'same_history_covariance_transition_operator_available','same_graph_every_prefix_hard_ball_target_API_available',
              'literal_every_prefix_transport_API_available','first_exit_reduction_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('retired_independent_integral_displacement_entry_ball_reintroduced','kernel_working_cells_claimed_reachable_without_prefix_proof',
              'explicit_600_sample_nonlinear_state_box_propagation_required','source_uniform_endpoint_augmented_LDLT_closed_here',
              'source_uniform_every_prefix_augmented_LDLT_closed_here','source_uniform_every_prefix_hard_domain_retention_closed_here',
              'P4_MOTION_PASS','P4_PASS','P5_MAY_START'):
        if d.get(k) is not False:f.append(k+' not false')
    if float(d.get('fresh_centered_integral_displacement_norm_m_s',1.0))!=0.0:f.append('fresh centered S is not zero')
    if S_GROUP in d.get('fresh_entry_reduced_group_order',[]):f.append('retired independent S factor remains in reduced entry')
    return list(dict.fromkeys(f))


def main()->int:
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    d=build();f=validate(d);d['validation_pass']=not f;d['validation_failures']=f
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'first_exit_reduction':d['first_exit_reduction_closed'],'centered_S':d['fresh_centered_integral_displacement_norm_m_s'],
      'local_working_cover':d['kernel_first_exit_working_cells_cover_full_hard_domain_locally'],'prefix_LDLT':d['source_uniform_every_prefix_augmented_LDLT_closed_here'],
      'hard_retention':d['source_uniform_every_prefix_hard_domain_retention_closed_here'],'failures':f},sort_keys=True))
    return int(bool(f))
if __name__=='__main__':raise SystemExit(main())
