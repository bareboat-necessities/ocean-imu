#!/usr/bin/env python3
"""Configured Normal-Live hybrid composition for ALT.

The declared operating domain says Normal-Live same-mode words contain no
hybrid events and lists ``held_to_active`` separately.  Therefore the theorem
factorization is

    H18 word* -> [optional held_to_active edge] -> A21 word*

rather than pretending the parallel A21 covariance generated during held mode
is the actual release covariance.

For the configured certificate, the public external bias-hold override is not a
source coordinate.  The shipping default is ``acc_bias_hold_=false``.  If an
external caller changes that runtime control, it is outside this configured
certificate and requires a separate hybrid-control theorem.

The held->active guard has two literal outcomes:
  * false: remain H18 with the same physical/covariance lineage;
  * true: joint24 mean state is continuous and A21 covariance is formed from the
    H18 covariance plus the latent held BA marginal, then the shipping enable
    floor is applied.  The held invariant is a fixed point of that floor.

This closes compatibility/coverage of the configured edge.  It does NOT prove
that startup or a particular 3 s word necessarily reaches the release guard;
finite-time release is a later capture/progress obligation.
"""
from __future__ import annotations
import json
from pathlib import Path

from tools.stability.ou3_alt_contraction import h18_a21_edge as EDGE
from tools.stability.ou3_alt_contraction import physical_word as WORD

REPO=Path(__file__).resolve().parents[3]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
WRAPPER=REPO/'src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h'
QUALIFICATION='OU3_ALT_CONFIGURED_H18_A21_HYBRID_WORD_V1'


def build():
    domain=json.loads(DOMAIN.read_text());normal=domain['normal_live'];edge=EDGE.build();ef=EDGE.validate(edge);word=WORD.induction_theorem();wf=WORD.validate(word)
    if ef or wf:raise RuntimeError(f'hybrid prerequisites failed edge={ef} word={wf}')
    text=WRAPPER.read_text()
    config={
      'same_mode_has_no_internal_hybrid_events':normal.get('hybrid_events_inside_same_mode_word')==[],
      'held_to_active_declared_separate':'held_to_active' in normal.get('hybrid_events_separate_from_P3_word',[]),
      'shipping_external_hold_default_false':'bool acc_bias_hold_ = false;' in text,
      'shipping_release_calls_enable':'if (!acc_bias_hold_)' in text and 'mekf_->set_acc_bias_updates_enabled(true);' in text,
    }
    closed=bool(all(config.values()) and edge['conditional_H18_A21_state_and_covariance_edge_closed'] and edge['release_guard_literal_relation_closed'] and word['H18_complete_word_construction_closed'] and word['A21_complete_word_construction_closed'])
    return {
      'qualification':QUALIFICATION,'canonical_source':'DECLARED_OU3_DEPLOYMENT_THEOREM_OPERATING_DOMAIN',
      'configured_scope':config,
      'H18_same_mode_word_closed':word['H18_complete_word_construction_closed'],'A21_same_mode_word_closed':word['A21_complete_word_construction_closed'],
      'guard_false_branch':'remain_H18_identity_mode_edge','guard_true_branch':'splice_exact_held_covariance_then_enable_A21',
      'both_configured_guard_outcomes_covered':closed,'joint24_mean_continuous_across_release':edge['physical_joint24_mean_edge_identity'],
      'one_time_Live_S_origin_preserved_across_release':edge['one_time_Live_S_origin_preserved'],'true_bias_history_preserved_across_release':edge['true_bias_history_preserved'],
      'actual_A21_release_covariance_comes_from_H18_held_covariance':True,'parallel_hypothetical_A21_covariance_used_at_release':False,
      'enable_covariance_floor_fixed_on_initial_held_invariant':edge['enable_floor_fixed_on_held_invariant'],
      'H18_A21_edge_attached':closed,'all_literal_configured_hybrid_branches_attached':closed,'hybrid_H18_A21_word_closed':closed,
      'release_guard_eventual_reachability_closed_here':False,'startup_capture_consumed':False,
      'external_runtime_bias_hold_override_admitted_by_certificate':False,
      'storage_search_allowed_here':False,'ALT_LIVE_PASS':False,
      'next_obligation':'assemble the Phase-1 closure ledger from the regional source relation, physical word induction and this configured hybrid edge; if all physical-word guards close, begin the first common joint24 storage feasibility search while keeping release progress/startup separate'}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    if not all(d.get('configured_scope',{}).values()):f.append('configured hybrid scope parity failed')
    for k in ('H18_same_mode_word_closed','A21_same_mode_word_closed','both_configured_guard_outcomes_covered','joint24_mean_continuous_across_release','one_time_Live_S_origin_preserved_across_release','true_bias_history_preserved_across_release','actual_A21_release_covariance_comes_from_H18_held_covariance','enable_covariance_floor_fixed_on_initial_held_invariant','H18_A21_edge_attached','all_literal_configured_hybrid_branches_attached','hybrid_H18_A21_word_closed'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('parallel_hypothetical_A21_covariance_used_at_release','release_guard_eventual_reachability_closed_here','startup_capture_consumed','external_runtime_bias_hold_override_admitted_by_certificate','storage_search_allowed_here','ALT_LIVE_PASS'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
