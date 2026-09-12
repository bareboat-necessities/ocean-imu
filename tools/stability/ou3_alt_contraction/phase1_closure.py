#!/usr/bin/env python3
"""ALT Phase-1 source/Jacobian ledger; NOT the finite storage master.

This module predates the finite shipping-word continuation.  It still consumes
the regional universal Normal-Live source relation, physical joint24 Jacobian
cocycle, all three analytic bias families, and configured H18->A21 factorization.
Those are useful source-ancestry facts, but ``physical_word`` now explicitly
classifies its map as a *pointwise state Jacobian cocycle*, not an anchored finite
error identity.

Consequently this ledger MUST remain fail-closed for storage.  The controlling
gate is ``proof_plan.assert_finite_storage_master`` and cannot be replaced by
legacy Phase-1 labels.  Storage/rho/high-precision work remains forbidden until
the finite event composer has physical reference forcing, all coefficient
product graphs and every configured branch bound in one source-uniform finite
shipping word.
"""
from __future__ import annotations
import json
from pathlib import Path

from tools.stability.ou3_alt_contraction import bias_families as BIAS
from tools.stability.ou3_alt_contraction import hybrid_word as HYBRID
from tools.stability.ou3_alt_contraction import live_source_relation as SOURCE
from tools.stability.ou3_alt_contraction import physical_word as WORD
from tools.stability.ou3_alt_contraction import proof_plan as PLAN
from tools.stability.ou3_alt_contraction import deployment_scope as SCOPE

REPO=Path(__file__).resolve().parents[3]
DOMAIN=REPO/'tools/stability/ou3_proof_operating_domain.json'
QUALIFICATION='OU3_ALT_PHASE1_SOURCE_JACOBIAN_LEDGER_V2'


def build():
    source=SOURCE.build(); sf=SOURCE.validate(source)
    word=WORD.induction_theorem(); wf=WORD.validate(word)
    hybrid=HYBRID.build(); hf=HYBRID.validate(hybrid)
    bias=BIAS.build(); bf=BIAS.validate(bias)
    bad={k:v for k,v in {'source':sf,'word':wf,'hybrid':hf,'bias':bf}.items() if v}
    if bad: raise RuntimeError('ALT Phase-1 prerequisite failed: '+repr(bad))

    domain=json.loads(DOMAIN.read_text()); normal=domain['normal_live']
    scope=SCOPE.certified_scope(); SCOPE.assert_certified_scope(scope)
    branch_scope=bool(
        normal['accelerometer_update_required_each_valid_imu_sample_after_live_entry']
        and normal['accelerometer_rejection_in_normal_live_scope'] is False
        and normal['hybrid_events_inside_same_mode_word']==[]
        and source['all_branch_successors_retained_by_operator']
        and hybrid['all_literal_configured_hybrid_branches_attached'])

    legacy_source_status={
      'same_history_complete_BRMM_word':bool(
          word['same_history_COMPLETE_BRMM_single_mode_word_closed']
          and hybrid['hybrid_H18_A21_word_closed']
          and source['regional_universal_Normal_Live_source_relation_closed']),
      'physical_prediction_forcing_attached':word['physical_prediction_forcing_attached'],
      'physical_S_residual_attached':word['physical_S_residual_attached'],
      'all_bias_families_attached':bool(
          word['all_bias_families_attached_to_single_mode_word']
          and bias['three_families_invoked_separately']),
      'all_literal_branches_attached':branch_scope,
      'H18_A21_edge_attached':hybrid['H18_A21_edge_attached'],
      'zero_wind_heel_scope_enforced':True,
    }

    finite=WORD.finite_storage_readiness()
    finite_status={
      'map_representation':finite['map_representation'],
      'finite_error_identity_for_every_event':finite['finite_error_identity_for_every_event'],
      'physical_reference_forcing_retained':finite['physical_reference_forcing_retained'],
      'all_coefficient_product_graphs_retained':finite['all_coefficient_product_graphs_retained'],
      'all_configured_branches_bound_to_finite_graph':finite['all_configured_branches_bound_to_finite_graph'],
      'zero_wind_heel_scope_enforced':True,
    }
    finite_gate_error=None
    try:
        PLAN.assert_finite_storage_master(finite_status)
    except RuntimeError as exc:
        finite_gate_error=str(exc)
    else:
        raise AssertionError('Jacobian-only Phase-1 ledger unexpectedly passed finite storage gate')

    return {
      'qualification':QUALIFICATION,
      'canonical_source':'COMPLETE_BRMM_NORMAL_LIVE_WORD',
      'legacy_source_jacobian_status':legacy_source_status,
      'regional_Normal_Live_source_relation_closed':source['regional_universal_Normal_Live_source_relation_closed'],
      'physical_single_mode_Jacobian_induction_closed':word['same_history_COMPLETE_BRMM_single_mode_word_closed'],
      'configured_hybrid_Jacobian_factorization_closed':hybrid['hybrid_H18_A21_word_closed'],
      'legacy_source_ancestry_ledger_closed':all(legacy_source_status.values()),
      'finite_storage_status':finite_status,
      'finite_storage_gate_error':finite_gate_error,
      'finite_storage_master_closed':False,
      'storage_search_allowed':False,
      'common_joint24_storage_must_wait_for_finite_master':True,
      'parameter_dependent_storage_search_allowed_initially':False,
      'release_guard_eventual_reachability_closed':False,
      'startup_capture_closed':False,
      'every_prefix_chart_retention_closed':False,
      'deployment_finite_precision_closed':False,
      'ALT_LIVE_PASS':False,
      'ALT_STARTUP_PASS':False,
      'ALT_END_TO_END_PASS':False,
      'no_replay_or_unreachable_perturbation_used':True,
      'old_P2_P3_P4_P5_route_modified':False,
      'next_obligation':(
          'bind the finite startup-rooted IMU/magnetic/hybrid composer to the universal '
          'COMPLETE-BRMM/BIAS continuation and retain physical forcing plus every '
          'coefficient-product graph; only assert_finite_storage_master may unlock storage'),
    }


def validate(d):
    f=[]
    if d.get('qualification')!=QUALIFICATION: f.append('qualification mismatch')
    legacy=d.get('legacy_source_jacobian_status',{})
    for k in ('same_history_complete_BRMM_word','physical_prediction_forcing_attached',
              'physical_S_residual_attached','all_bias_families_attached',
              'all_literal_branches_attached','H18_A21_edge_attached',
              'zero_wind_heel_scope_enforced'):
        if legacy.get(k) is not True: f.append('legacy source/Jacobian '+k+' not true')
    for k in ('regional_Normal_Live_source_relation_closed',
              'physical_single_mode_Jacobian_induction_closed',
              'configured_hybrid_Jacobian_factorization_closed',
              'legacy_source_ancestry_ledger_closed',
              'common_joint24_storage_must_wait_for_finite_master',
              'no_replay_or_unreachable_perturbation_used'):
        if d.get(k) is not True: f.append(k+' not true')
    for k in ('finite_storage_master_closed','storage_search_allowed',
              'parameter_dependent_storage_search_allowed_initially',
              'release_guard_eventual_reachability_closed','startup_capture_closed',
              'every_prefix_chart_retention_closed','deployment_finite_precision_closed',
              'ALT_LIVE_PASS','ALT_STARTUP_PASS','ALT_END_TO_END_PASS',
              'old_P2_P3_P4_P5_route_modified'):
        if d.get(k) is not False: f.append(k+' not false')
    finite=d.get('finite_storage_status',{})
    if finite.get('map_representation')!='pointwise_state_Jacobian_cocycle':
        f.append('legacy map representation unexpectedly changed')
    for k in ('finite_error_identity_for_every_event','physical_reference_forcing_retained',
              'all_coefficient_product_graphs_retained','all_configured_branches_bound_to_finite_graph'):
        if finite.get(k) is not False: f.append('finite blocker '+k+' not false')
    if not d.get('finite_storage_gate_error'):
        f.append('finite storage guard did not fail closed')
    return f
