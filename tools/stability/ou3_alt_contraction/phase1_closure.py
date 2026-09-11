#!/usr/bin/env python3
"""ALT Phase-1 closure ledger: physical whole-word before storage.

This module has one job: decide whether the theorem-only plan permits a storage
search.  It consumes the regional universal Normal-Live source relation, the
physical joint24 word induction, all three analytic bias families, and the
configured H18->A21 hybrid factorization.  It does not compute a metric or rho.

A TRUE ``storage_search_allowed`` therefore means only that the *master graph is
complete enough to ask the storage question*: same-history COMPLETE-BRMM source,
physical BRMM-vs-OU forcing, finite physical S residual, BIAS0/1/2, every
configured Normal-Live event/branch, and the held->active edge are represented.
It is NOT a stability PASS, basin proof, release-progress proof, startup proof or
finite-precision certificate.
"""
from __future__ import annotations
import json
from pathlib import Path

from tools.stability.ou3_alt_contraction import bias_families as BIAS
from tools.stability.ou3_alt_contraction import hybrid_word as HYBRID
from tools.stability.ou3_alt_contraction import live_source_relation as SOURCE
from tools.stability.ou3_alt_contraction import physical_word as WORD
from tools.stability.ou3_alt_contraction import proof_plan as PLAN

REPO=Path(__file__).resolve().parents[3]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
QUALIFICATION='OU3_ALT_PHASE1_PHYSICAL_WHOLE_WORD_CLOSURE_V1'

def build():
    source=SOURCE.build();sf=SOURCE.validate(source);word=WORD.induction_theorem();wf=WORD.validate(word);hybrid=HYBRID.build();hf=HYBRID.validate(hybrid);bias=BIAS.build();bf=BIAS.validate(bias)
    bad={k:v for k,v in {'source':sf,'word':wf,'hybrid':hf,'bias':bf}.items() if v}
    if bad:raise RuntimeError('ALT Phase-1 prerequisite failed: '+repr(bad))
    domain=json.loads(DOMAIN.read_text());normal=domain['normal_live']
    branch_scope=bool(
        normal['accelerometer_update_required_each_valid_imu_sample_after_live_entry']
        and normal['accelerometer_rejection_in_normal_live_scope'] is False
        and normal['hybrid_events_inside_same_mode_word']==[]
        and source['all_branch_successors_retained_by_operator']
        and hybrid['all_literal_configured_hybrid_branches_attached'])
    status={
      'same_history_complete_BRMM_word':bool(word['same_history_COMPLETE_BRMM_single_mode_word_closed'] and hybrid['hybrid_H18_A21_word_closed'] and source['regional_universal_Normal_Live_source_relation_closed']),
      'physical_prediction_forcing_attached':word['physical_prediction_forcing_attached'],
      'physical_S_residual_attached':word['physical_S_residual_attached'],
      'all_bias_families_attached':bool(word['all_bias_families_attached_to_single_mode_word'] and bias['three_families_invoked_separately']),
      'all_literal_branches_attached':branch_scope,
      'H18_A21_edge_attached':hybrid['H18_A21_edge_attached'],
    }
    allowed=PLAN.storage_search_allowed(status)
    if allowed: PLAN.assert_storage_search_allowed(status)
    return {
      'qualification':QUALIFICATION,'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD','phase1_status':status,
      'regional_Normal_Live_source_relation_closed':source['regional_universal_Normal_Live_source_relation_closed'],
      'physical_single_mode_word_induction_closed':word['same_history_COMPLETE_BRMM_single_mode_word_closed'],
      'configured_hybrid_word_factorization_closed':hybrid['hybrid_H18_A21_word_closed'],
      'storage_search_allowed':allowed,
      'common_joint24_storage_must_be_searched_first':True,'parameter_dependent_storage_search_allowed_initially':False,
      'release_guard_eventual_reachability_closed':False,'startup_capture_closed':False,'every_prefix_chart_retention_closed':False,'deployment_finite_precision_closed':False,
      'ALT_LIVE_PASS':False,'ALT_STARTUP_PASS':False,'ALT_END_TO_END_PASS':False,
      'no_replay_or_unreachable_perturbation_used':True,'old_P2_P3_P4_P5_route_modified':False,
      'next_obligation':'run the first source-uniform COMMON joint24 storage feasibility problem on this physical master, with bounded neutral/source supply; do not switch to source-dependent/piecewise storage unless the common search fails under the plan stop rule'}
def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    status=d.get('phase1_status',{})
    for k in ('same_history_complete_BRMM_word','physical_prediction_forcing_attached','physical_S_residual_attached','all_bias_families_attached','all_literal_branches_attached','H18_A21_edge_attached'):
        if status.get(k) is not True:f.append('phase1 '+k+' not true')
    for k in ('regional_Normal_Live_source_relation_closed','physical_single_mode_word_induction_closed','configured_hybrid_word_factorization_closed','storage_search_allowed','common_joint24_storage_must_be_searched_first','no_replay_or_unreachable_perturbation_used'):
        if d.get(k) is not True:f.append(k+' not true')
    for k in ('parameter_dependent_storage_search_allowed_initially','release_guard_eventual_reachability_closed','startup_capture_closed','every_prefix_chart_retention_closed','deployment_finite_precision_closed','ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS','old_P2_P3_P4_P5_route_modified'):
        if d.get(k) is not False:f.append(k+' not false')
    return f
