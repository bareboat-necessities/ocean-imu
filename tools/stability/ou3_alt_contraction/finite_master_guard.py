"""Executable readiness guard for the current ALT finite shipping master.

This module does not prove contraction and never searches storage.  It aggregates
the strongest currently composed startup-rooted, admitted-source, machine/deployment
relations and feeds only a conservative status to ``proof_plan.assert_finite_storage_master``.
Closed sub-obligations are named explicitly so stale handoffs cannot keep them as
blockers; any missing source-uniform arithmetic or literal branch keeps the guard
fail-closed.
"""
from __future__ import annotations

from tools.stability.ou3_alt_contraction import proof_plan as PLAN
from tools.stability.ou3_alt_contraction import finite_admitted_startup_machine_tunestate_word as START
from tools.stability.ou3_alt_contraction import finite_admitted_machine_clock_qualified_interleaved_prefix as LIVE
from tools.stability.ou3_alt_contraction import finite_scheduler_nextafter_binary32 as NEXT
from tools.stability.ou3_alt_contraction import finite_aw_sync_clock_binary64 as AWCLOCK
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as WPELOG
from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as WPEFREQ
from tools.stability.ou3_alt_contraction import finite_qaxis_exp_binary32 as QEXP

QUALIFICATION='OU3_ALT_FINITE_MASTER_GUARD_V1'


def build():
    startup=START.readiness()
    live=LIVE.readiness(); nxt=NEXT.readiness(); aw=AWCLOCK.readiness()
    wlog=WPELOG.readiness(); wfreq=WPEFREQ.readiness(); qexp=QEXP.readiness()

    closed={
      'startup_to_Live_machine_history_attachment': startup['startup_to_Live_machine_history_attachment_closed'],
      'admitted_source_private_Mahony_startup_to_Live_invariant': startup['admitted_source_private_Mahony_startup_to_Live_invariant_closed'],
      'guard_Mahony_frontend_Racc_accelerometer_same_event_join': live['strong_joined_guard_frontend_Racc_word_consumed'],
      'scheduler_nextafter_binary32': nxt['machine_scheduler_nextafter_binary32_correspondence_closed'],
      'bounded_default_aw_sync_binary64_clock': aw['canonical_aw_sync_binary64_predicate_closed'],
      'all_600_IMU_edges_require_aw_sync_clock_qualification': live['complete_word_requires_clock_qualification_on_all_600_IMU_edges'],
    }

    open_obligations={
      'universal_startup_source_and_branch_reachability': not startup['every_admitted_startup_history_reaches_this_boundary_product'],
      'startup_source_uniform_deployment_supplies': not startup['source_uniform_startup_deployment_supply_bounds_closed'],
      'target_libm_and_compiler_profile_correspondence': not startup['all_target_libm_and_compiler_profile_correspondence_closed'],
      'WPE_log_and_exp_libm_correspondence': not (wlog['WPE_log_std_log_target_libm_correspondence_closed'] and wlog['WPE_log_exp_target_libm_correspondence_closed']),
      'WPE_raw_period_binary32_production': not wlog['WPE_raw_period_binary32_production_closed'],
      'WPE_frequency_source_uniform_supply': not wfreq['source_uniform_WPE_frequency_supply_bound_closed'],
      'Qaxis_exp_libm_correspondence': not qexp['Qaxis_exp_libm_binary32_correspondence_closed'],
      'all_event_arithmetic_witnesses_source_uniform': not live['all_event_arithmetic_witnesses_source_uniformly_qualified'],
      'complete_source_uniform_600_step_word': not live['source_uniform_complete_600_step_word_qualified'],
      'mutable_aw_sync_cadence_setter_ancestry': not live['mutable_adapt_every_setter_ancestry_closed'],
      'arbitrary_dt_inner_clock': not live['arbitrary_dt_inner_clock_closed'],
    }

    # The current graph is a genuine finite descriptor product, but not yet a
    # complete source-uniform master.  Do not claim representation completion
    # until every literal branch and arithmetic supply is qualified.
    finite_status={
      'map_representation':'finite_physical_descriptor_partial',
      'finite_error_identity_for_every_event':False,
      'physical_reference_forcing_retained':True,
      'all_coefficient_product_graphs_retained':False,
      'all_configured_branches_bound_to_finite_graph':False,
      'zero_wind_heel_scope_enforced':True,
    }
    error=None
    try:
        PLAN.assert_finite_storage_master(finite_status)
    except RuntimeError as exc:
        error=str(exc)
    else:
        raise AssertionError('incomplete finite master unexpectedly unlocked storage')

    return {
      'qualification':QUALIFICATION,
      'closed_subobligations':closed,
      'open_obligations':open_obligations,
      'finite_storage_status':finite_status,
      'finite_storage_guard_error':error,
      'finite_master_guard_closed':False,
      'storage_search_allowed':False,
      'ALT_STARTUP_PASS':False,'ALT_LIVE_PASS':False,'ALT_END_TO_END_PASS':False,
    }


def validate(x):
    f=[]
    if x.get('qualification')!=QUALIFICATION:f.append('qualification mismatch')
    for k,v in x.get('closed_subobligations',{}).items():
        if v is not True:f.append('closed subobligation regressed: '+k)
    if not x.get('open_obligations') or not all(v is True for v in x['open_obligations'].values()):
        f.append('open-obligation polarity mismatch')
    if not x.get('finite_storage_guard_error'):f.append('finite storage guard did not fail closed')
    for k in ('finite_master_guard_closed','storage_search_allowed','ALT_STARTUP_PASS','ALT_LIVE_PASS','ALT_END_TO_END_PASS'):
        if x.get(k) is not False:f.append(k+' not false')
    return f
