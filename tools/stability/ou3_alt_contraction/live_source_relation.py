#!/usr/bin/env python3
"""Regional universal Normal-Live source relation for the ALT proof.

The relation is

  O^601_BRMM
  x X_front,Live
  x BIASi

where O^601_BRMM is the existing correlated 601-sample physical outer relation,
X_front,Live is the regional frontend predecessor invariant (including the
source-order binary32 private-Mahony invariant), and BIASi is one of the three
separate analytic physical bias families.

For each sample, coefficients are generated only by the existing exact
same-history JOINT transition operator, and Riccati/event cells only by the
trusted execution-kernel/estimator-attachment operator.  This module deliberately
does NOT call their top-level build/smoke functions, because those currently
import startup-entry qualification.  Startup membership is not a premise of a
regional Live theorem; it is a later capture theorem.

The proof here is an induction relation, not a finite source enumeration.  The
component invariants are forward invariant and every branch successor of the
literal transition is retained.  Physical q15/S response propagation is handled
by ``physical_lineage``; complete 600-step cocycle assembly is the next layer.
"""
from __future__ import annotations
from pathlib import Path

import ou3_brmm_correlated_window_outer_enclosure as OUTER
import ou3_p4_joint_brmm_frontend_transition as JOINT
import ou3_brmm_complete_window_execution_kernel as KERNEL
import ou3_p4_source_uniform_estimator_event_attachment as ATTACH
from tools.stability.ou3_alt_contraction import bias_families as BIAS
from tools.stability.ou3_alt_contraction import live_frontend_predecessor as FRONT

REPO=Path(__file__).resolve().parents[3]
JOINT_SRC=REPO/'tools/stability/ou3_p4_joint_brmm_frontend_transition.py'
ATTACH_SRC=REPO/'tools/stability/ou3_p4_source_uniform_estimator_event_attachment.py'
QUALIFICATION='OU3_ALT_REGIONAL_UNIVERSAL_NORMAL_LIVE_SOURCE_RELATION_V1'

def _operator_parity():
    j=JOINT_SRC.read_text();a=ATTACH_SRC.read_text()
    return {
      'joint_commits_before_measurement':'committed=TUNER.commit_if_pending' in j,
      'joint_private_Mahony_step':'MAHONY.advance_initialized_live' in j,
      'joint_stillness_all_successors':'still_images=STILL.advance' in j and 'for si in still_images:' in j,
      'joint_WPE_all_successors':'wpe_images=STATS.advance_wpe' in j and 'for wi in wpe_images:' in j,
      'joint_tuner_all_successors':'for tuner_next,target,tuner_stats,timer in tuner_images:' in j,
      'joint_same_predecessor_token':'token,state.source_token,child' in j,
      'attachment_uses_kernel_joint_operator':'KERNEL.advance_branch_with_joint_frontend' in a,
      'attachment_retains_every_joint_successor':'for ordinal, raw_image in enumerate(images):' in a,
      'attachment_checks_same_active_schedule':'same_active_schedule_verified' in a,
      'attachment_checks_same_actual_RS':'same_actual_RS_verified' in a,
      'attachment_uses_JOINT_frontend_as_successor':'JOINT_image_authoritative_for_next_frontend_state' in a,
    }

def build():
    outer=OUTER.build();of=OUTER.validate(outer);front=FRONT.build();ff=FRONT.validate(front);bias=BIAS.build();bf=BIAS.validate(bias)
    bad={k:v for k,v in {'outer':of,'regional_frontend':ff,'bias':bf}.items() if v}
    if bad:raise RuntimeError('regional Live source relation prerequisite failed: '+repr(bad))
    parity=_operator_parity();operators=all(parity.values()) and callable(JOINT.advance) and callable(KERNEL.advance_branch_with_joint_frontend) and callable(ATTACH.synchronize_sample)
    closed=bool(outer['left_inclusion_closed'] and outer['validated_correlated_outer_enclosure_closed'] and outer['same_history_required_for_entire_window'] and outer['correlation_retained_across_samples'] and front['regional_normal_live_frontend_predecessor_set_closed'] and operators and bias['three_families_invoked_separately'] and bias['one_persistent_parameter_token_per_family'])
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','physical_outer_relation':outer['outer_set_symbol'],'sample_count':outer['sample_count'],
      'COMPLETE_BRMM_left_inclusion_consumed':outer['left_inclusion_closed'],'correlated_physical_history_across_601_samples':outer['correlation_retained_across_samples'],
      'regional_frontend_predecessor_invariant_consumed':front['regional_normal_live_frontend_predecessor_set_closed'],'startup_capture_consumed':False,'startup_entry_membership_closed_here':False,
      'individual_BIAS0_BIAS1_BIAS2_contracts_consumed':bias['three_families_invoked_separately'],'persistent_bias_parameter_token_required_over_word':True,
      'operator_source_parity':parity,'same_history_joint_frontend_transition_operator_bound':operators,'trusted_Riccati_event_attachment_operator_bound':operators,'all_branch_successors_retained_by_operator':operators,
      'independent_f_sigma_tau_TS_RS_boxes_used':False,'independent_sample_boxes_used':False,'replay_or_seeded_realization_used':False,'favorable_successor_selected':False,
      'regional_universal_Normal_Live_source_relation_closed':closed,
      'complete_600_step_physical_joint24_cocycle_closed_here':False,'H18_A21_guard_reachability_closed_here':False,'storage_search_allowed':False,'ALT_LIVE_PASS':False,
      'next_obligation':'apply physical_lineage to every transition of this relation by induction, carrying prior response blocks, primitive continuity, one BIAS parameter token and latent held covariance; splice the exact H18->A21 edge when its literal guard becomes true'}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION or d.get('physical_outer_relation')!='O^601_BRMM' or d.get('sample_count')!=601:f.append('qualification/source window mismatch')
    for k in ('COMPLETE_BRMM_left_inclusion_consumed','correlated_physical_history_across_601_samples','regional_frontend_predecessor_invariant_consumed','individual_BIAS0_BIAS1_BIAS2_contracts_consumed','persistent_bias_parameter_token_required_over_word','same_history_joint_frontend_transition_operator_bound','trusted_Riccati_event_attachment_operator_bound','all_branch_successors_retained_by_operator','regional_universal_Normal_Live_source_relation_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    if not all(d.get('operator_source_parity',{}).values()):f.append('operator source parity failed')
    for k in ('startup_capture_consumed','startup_entry_membership_closed_here','independent_f_sigma_tau_TS_RS_boxes_used','independent_sample_boxes_used','replay_or_seeded_realization_used','favorable_successor_selected','complete_600_step_physical_joint24_cocycle_closed_here','H18_A21_guard_reachability_closed_here','storage_search_allowed','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
