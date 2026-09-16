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
from tools.stability.ou3_alt_contraction import finite_admitted_startup_wpe_machine_word as STARTWPEADMIT
from tools.stability.ou3_alt_contraction import finite_admitted_machine_clock_qualified_interleaved_prefix as LIVE
from tools.stability.ou3_alt_contraction import finite_admitted_wpe_machine_clock_interleaved_prefix as LIVEWPE
from tools.stability.ou3_alt_contraction import finite_binary32_mahony_startup as STARTSEED
from tools.stability.ou3_alt_contraction import finite_startup_wpe_machine_history as STARTWPE
from tools.stability.ou3_alt_contraction import finite_scheduler_nextafter_binary32 as NEXT
from tools.stability.ou3_alt_contraction import finite_aw_sync_clock_binary64 as AWCLOCK
from tools.stability.ou3_alt_contraction import finite_wpe_moment_binary32 as WPEMOM
from tools.stability.ou3_alt_contraction import finite_wpe_machine_binary32 as WPEMACHINE
from tools.stability.ou3_alt_contraction import finite_wpe_log_binary32 as WPELOG
from tools.stability.ou3_alt_contraction import finite_wpe_frequency_binary32 as WPEFREQ
from tools.stability.ou3_alt_contraction import finite_qaxis_exp_binary32 as QEXP
from tools.stability.ou3_alt_contraction import finite_startup_entry_obstruction as ENTRY
from tools.stability.ou3_alt_contraction import finite_startup_disturbance_obstruction as STARTDIST
from tools.stability.ou3_alt_contraction import finite_startup_handoff_control as CONTROL
from tools.stability.ou3_alt_contraction import finite_wrapper_clock_binary32 as CLOCK
from tools.stability.ou3_alt_contraction import finite_mag_startup_capture_bound as CAPTURE
from tools.stability.ou3_alt_contraction import finite_mag_call_schedule as CALLS
from tools.stability.ou3_alt_contraction import finite_attitude_atlas as ATLAS
from tools.stability.ou3_alt_contraction import finite_fresh_joint24_entry as FRESH
from tools.stability.ou3_alt_contraction import finite_live_interleave as INTERLEAVE
from tools.stability.ou3_alt_contraction import finite_mag_counter_saturation as COUNTER
from tools.stability.ou3_alt_contraction import finite_admitted_source_live_word as ADLIVE
from tools.stability.ou3_alt_contraction import finite_source_bound_mag_dual_clock as MAGWORD
from tools.stability.ou3_alt_contraction import finite_startup_sensor_contract as SENSOR_CONTRACT
from tools.stability.ou3_alt_contraction import finite_wpe_uniform_bounds as WPEBOUNDS
from tools.stability.ou3_alt_contraction import finite_mahony_prefix_totality as MAHONYTOTAL
from tools.stability.ou3_alt_contraction import finite_seed_svd_axis_reduction as SVDAXIS
from tools.stability.ou3_alt_contraction import finite_live_disturbance_overflow as LIVEOVERFLOW
from tools.stability.ou3_alt_contraction import finite_startup_direction_sampling as SAMPLING
from tools.stability.ou3_alt_contraction import finite_live_input_contract as INPUT
from tools.stability.ou3_alt_contraction import target_wpe_libm as TARGETLIBM
from tools.stability.ou3_alt_contraction import target_wpe_compiler as TARGETWPE
from tools.stability.ou3_alt_contraction import finite_seed_eigen_svd as EIGENSVD
from tools.stability.ou3_alt_contraction import finite_frontend_uniform_bounds as FRONTBOUNDS
from tools.stability.ou3_alt_contraction import finite_startup_timeout_alignment_obstruction as TIMEOUT
from tools.stability.ou3_alt_contraction import finite_candidate_uniform_bounds as CANDIDATE
from tools.stability.ou3_alt_contraction import finite_machine_startup_core as CORESTART

QUALIFICATION='OU3_ALT_FINITE_MASTER_GUARD_V1'
OPEN_QUALIFICATIONS=(
    'source_uniform_timeout_aligned_branch_reachability',
    'near_antiparallel_Eigen_JacobiSVD_solver_correspondence',
    'universal_startup_source_and_branch_reachability',
    'startup_source_uniform_deployment_supplies',
    'target_libm_and_compiler_profile_correspondence',
    'WPE_machine_target_libm_and_compiler_selection',
    'WPE_source_uniform_machine_supply_bounds',
    'WPE_log_and_exp_libm_correspondence',
    'Qaxis_exp_libm_correspondence',
    'all_event_arithmetic_witnesses_source_uniform',
    'complete_source_uniform_600_step_word',
)


def qualification_status():
    startup=START.readiness(); seed=STARTSEED.readiness(); clock=CLOCK.readiness()
    eigen=EIGENSVD.source_uniform_certificate()
    live=LIVE.readiness(); lwpe=LIVEWPE.readiness()
    wmach=WPEMACHINE.readiness(); wlog=WPELOG.readiness(); qexp=QEXP.readiness()
    target=TARGETLIBM.profile_correspondence()
    compiler=TARGETWPE.readiness()
    error_profile=WPEBOUNDS.target_error_profile_certificate()
    return {
      'source_uniform_timeout_aligned_branch_reachability': not clock['universal_startup_deadline_closed'],
      # The source-order QR/Jacobi proof is independent of the later target
      # compiler/FCR gate. Keep that deployment qualification separate (#5).
      'near_antiparallel_Eigen_JacobiSVD_solver_correspondence': not (
          eigen['source_uniform_QR_and_Jacobi_totality_closed'] and
          eigen['axis_is_computed_QR_column_two'] and
          eigen['both_pivots_rank_one_rank_two_and_tiny_tails_covered']),
      'universal_startup_source_and_branch_reachability': not startup['every_admitted_startup_history_reaches_this_boundary_product'],
      'startup_source_uniform_deployment_supplies': not startup['source_uniform_startup_deployment_supply_bounds_closed'],
      'target_libm_and_compiler_profile_correspondence': not startup['all_target_libm_and_compiler_profile_correspondence_closed'],
      'WPE_machine_target_libm_and_compiler_selection': not (
          compiler['pinned_WPE_compiler_profile_selection_closed'] and
          target['pinned_WPE_exp_log_approximation_correspondence_closed'] and
          target['pinned_WPE_sqrt_approximation_correspondence_closed'] and
          error_profile['target_exp_log_satisfy_WPE_error_profile_under_scalar_and_link_premises']),
      'WPE_source_uniform_machine_supply_bounds': not lwpe['source_uniform_WPE_machine_supply_bounds_closed'],
      'WPE_log_and_exp_libm_correspondence': not (
          target['pinned_WPE_exp_log_approximation_correspondence_closed'] and
          target['pinned_WPE_sqrt_approximation_correspondence_closed'] and
          error_profile['target_exp_log_satisfy_WPE_error_profile_under_scalar_and_link_premises']),
      'Qaxis_exp_libm_correspondence': not qexp['Qaxis_pinned_libm_approximation_correspondence_closed'],
      'all_event_arithmetic_witnesses_source_uniform': not live['all_event_arithmetic_witnesses_source_uniformly_qualified'],
      'complete_source_uniform_600_step_word': not (
          lwpe['source_uniform_complete_600_step_word_qualified'] and
          lwpe['machine_CORE_and_control_successor_algorithms_attached']),
    }


def build():
    startup=START.readiness(); seed=STARTSEED.readiness(); swpe=STARTWPE.readiness(); saw=STARTWPEADMIT.readiness()
    live=LIVE.readiness(); lwpe=LIVEWPE.readiness(); nxt=NEXT.readiness(); aw=AWCLOCK.readiness()
    wm=WPEMOM.readiness(); wmach=WPEMACHINE.readiness(); wlog=WPELOG.readiness(); wfreq=WPEFREQ.readiness(); qexp=QEXP.readiness()
    entry=ENTRY.build(); startup_disturbance=STARTDIST.build()
    control=CONTROL.readiness(); capture=CAPTURE.readiness()
    clock=CLOCK.readiness(); counter=CALLS.counter_lifetime(); atlas=ATLAS.readiness()
    counter_certificate=COUNTER.build()
    sensor_contract=SENSOR_CONTRACT.build()
    wpe_bounds=WPEBOUNDS.build()
    mahony_total=MAHONYTOTAL.build(); svd_axis=SVDAXIS.build()
    eigen_svd=EIGENSVD.source_uniform_certificate()
    frontend_bounds=FRONTBOUNDS.build()
    candidate_bounds=CANDIDATE.build()
    startup_core=CORESTART.readiness()

    closed={
      'physical_MEMS_all_initialized_scalar_finite_prefix_totality':INPUT.prefix_certificate()['all_initialized_finite_prefixes_totality_closed'],
      'near_antiparallel_source_uniform_QR_and_Jacobi_totality':eigen_svd['source_uniform_QR_and_Jacobi_totality_closed'],
      'configured_frontend_all_finite_prefix_arithmetic_supplies':frontend_bounds['source_uniform_configured_frontend_arithmetic_totality_closed'],
      'configured_candidate_tau_sigma_RS_powf_finite_supplies':candidate_bounds['candidate_arithmetic_totality_under_named_pow_range_closed'],
      'source_owned_ungauged_startup_CORE_constructor_is_sealed':bool(
          startup_core['default_retained_construction_from_selected_wrapper_and_app'] and
          startup_core['ungauged_machine_proxy_to_full_CORE_producer_available'] and
          not startup_core['caller_replacement_quaternion_or_covariance_port']),
      'conditional_ordinary_seed_scalar_Mahony_prefix_totality':mahony_total['ordinary_seed_scalar_prefix_totality_closed'],
      'returning_pinned_Eigen_SVD_seed_axis_equals_QR_third_column':svd_axis['returning_solver_axis_equals_QR_Q_column_2'],
      'local_CORE_successors_retained_and_following_predecessors_checked':lwpe['machine_CORE_local_successors_persist_and_next_predecessors_checked'],
      'complete_word_requires_CORE_and_control_successor_continuation':lwpe['complete_word_requires_machine_CORE_and_control_continuation'],
      'Qaxis_general_branch_sufficient_exp_error_budget':qexp['Qaxis_general_branch_sufficient_libm_error_budget']['source_uniform_sufficient_error_budget_proved'],
      'commissioned_startup_source_uniform_seed_norm_margin':sensor_contract['source_uniform_seed_norm_margin_closed'],
      'bounded_vertical_input_uniform_WPE_moment_raw_period_log_supplies':wpe_bounds['reset_to_every_finite_prefix_bounded_input_induction_closed'],
      'full_magnetic_frame_represented_without_small_angle_capture': bool(
          capture['unrestricted_full_frame_chord_bound_closed'] and
          atlas['all_nonzero_relative_quaternions_covered'] and
          FRESH.readiness()['independent_fresh_entry_error_box_removed'] and
          not FRESH.readiness()['small_magnetic_capture_radius_required']),
      'independent_WPE_machine_and_exact_branches_composed': bool(
          swpe['independent_machine_WPE_production_and_frequency_branches_composed'] and
          lwpe['independent_machine_WPE_production_and_frequency_branches_composed'] and
          wfreq['eager_frequency_getter_and_invalid_post_latch_fallback_composed']),
      'source_uniform_clamped_WPE_frequency_discrepancy':wfreq['source_uniform_WPE_frequency_supply_bound_closed'],
      'literal_per_compiler_WPE_usable_latch':wmach['source_produced_per_compiler_usable_latches_retained'],
      'signed_magnetic_counter_safety_on_every_finite_prefix':counter_certificate['counter_lifetime_closed'],
      'accepted_and_rejected_magnetic_sample_counter_safety':counter_certificate['accepted_and_rejected_magnetic_sample_counters_closed'],
      'certified_empty_startup_and_later_north_in_admitted_source_word': bool(
          ADLIVE.readiness()['certified_empty_magnetic_prefix_admitted_at_ungauged_timeout'] and
          MAGWORD.readiness()['ungauged_north_and_saturated_counter_use_shared_event_composer']),
      'ungauged_Live_waiting_and_later_north_runtime_branch': INTERLEAVE.readiness()['ungauged_timeout_entry_and_later_north_acquisition_composed'],
      'all_nonzero_fresh_attitudes_represented_by_joint24_atlas': bool(
          atlas['all_nonzero_relative_quaternions_covered'] and
          FRESH.readiness()['all_nonzero_fresh_relative_attitudes_represented']),
      'finite_attitude_left_right_and_chart_transport': atlas['finite_left_right_and_chart_transport_identity'],
      'inverse_free_measurement_with_state_bound_chart_offset': atlas['inverse_free_measurement_descriptor_with_chart_transport'],
      'literal_quality_and_timeout_handoff_control':control['quality_and_timeout_predicates_materialized'],
      'exact_default_timeout_clock_crossing':clock['default_timeout_first_crossing_proved'],
      'startup_to_Live_machine_history_attachment': startup['startup_to_Live_machine_history_attachment_closed'],
      'conditional_legacy_private_Mahony_invariant_under_its_original_sensor_premises': startup['admitted_source_private_Mahony_startup_to_Live_invariant_closed'],
      'near_antiparallel_seed_branch_topology': seed['near_antiparallel_JacobiSVD_branch_topology_materialized_with_solver_witness'],
      'guard_Mahony_frontend_Racc_accelerometer_same_event_join': live['strong_joined_guard_frontend_Racc_word_consumed'],
      'scheduler_nextafter_binary32': nxt['machine_scheduler_nextafter_binary32_correspondence_closed'],
      'bounded_default_aw_sync_binary64_clock': aw['canonical_aw_sync_binary64_predicate_closed'],
      'all_600_IMU_edges_require_aw_sync_clock_qualification': live['complete_word_requires_clock_qualification_on_all_600_IMU_edges'],
      'aw_sync_cadence_bound_to_persistent_runtime_config': live['adapt_every_runtime_value_ancestry_closed_for_current_word'],
      'canonical_5ms_source_dt_bound_to_clock_word': live['canonical_5ms_source_dt_ancestry_closed_for_current_word'],
      'WPE_raw_period_branch_topology_and_same_moment_arithmetic': bool(
          wm['all_raw_period_early_return_branches_materialized'] and
          wm['raw_period_has_no_independent_input_port'] and
          wm['raw_period_derived_from_same_machine_moments_and_sqrt_result']),
      'WPE_dual_compiler_moment_raw_log_persistent_product': bool(
          wmach['dual_compiler_persistent_moment_history_attached'] and
          wmach['WPE_raw_period_binary32_production_attached_to_log_history'] and
          wmach['raw_period_bound_to_same_mode_std_log_input']),
      'WPE_repeated_period_exp_bound_to_same_stored_log_state': bool(
          wmach['moment_horizon_period_bound_to_same_stored_log_state'] and
          wmach['same_exp_log_period_value_bound_across_horizon_and_log_smoothing']),
      'WPE_full_machine_history_attached_to_startup_Mahony_vertical': bool(
          swpe['startup_frontend_vertical_ancestry_attached'] and
          swpe['full_WPE_machine_log_state_equals_existing_lower_WPE_ledger_after_each_attached_step']),
      'WPE_full_machine_history_preserved_across_goLive': swpe['goLive_preserves_full_WPE_moment_log_machine_history_by_identity'],
      'WPE_startup_history_substituted_into_admitted_Live_without_snapshot': saw['startup_guard_Mahony_WPE_machine_to_admitted_Live_attachment_closed'],
      'WPE_full_machine_history_required_on_all_600_Live_IMU_edges': bool(
          lwpe['Live_600_step_WPE_machine_history_attached'] and
          lwpe['complete_word_requires_full_WPE_machine_step_on_all_600_IMU_edges']),
      'WPE_MAG_and_HOLD_preserve_full_machine_history': lwpe['MAG_and_HOLD_preserve_full_WPE_machine_history_by_identity'],
    }

    open_obligations=qualification_status()

    finite_status={
      'map_representation':'finite_physical_descriptor_partial',
      'finite_error_identity_for_every_event':False,
      'physical_reference_forcing_retained':True,
      'all_coefficient_product_graphs_retained':False,
      'all_configured_branches_bound_to_finite_graph':False,
      'all_admitted_startup_entries_represented':entry['multi_chart_or_quotient_entry_closed'],
      'source_uniform_deployment_arithmetic_closed':False,
      'finite_word_counter_safety_closed':counter.no_signed_overflow_proved,
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
      'research_outcome':'finite_master_qualification_incomplete',
      'closed_subobligations':closed,
      'open_obligations':open_obligations,
      'falsified_prerequisites':{},
      'counter_saturation_certificate':counter_certificate,
      'startup_entry_obstruction':entry,
      'startup_disturbance_obstruction':startup_disturbance,
      'startup_sensor_contract':sensor_contract,
      'bounded_input_WPE_supplies':wpe_bounds,
      'conditional_Mahony_prefix_totality':mahony_total,
      'pinned_Eigen_seed_axis_reduction':svd_axis,
      'pinned_Eigen_source_uniform_SVD':eigen_svd,
      'configured_frontend_uniform_supplies':frontend_bounds,
      'configured_candidate_uniform_supplies':candidate_bounds,
      'source_owned_startup_CORE_constructor':startup_core,
      'commissioned_Live_MEMS_input_contract':INPUT.build(),
      'bounded_raw_Live_input_totality_obstruction':LIVEOVERFLOW.build(),
      'continuous_to_sampled_startup_budget_audit':SAMPLING.build(),
      'physical_timeout_alignment_witness':TIMEOUT.build(),
      'pinned_target_library_profile':TARGETLIBM.profile_correspondence(),
      'pinned_WPE_compiler_profile':TARGETWPE.readiness(),
      'library_closure_requires_separate_firmware_compiler_qualification':True,
      'magnetic_frame_bounds':capture,
      'attitude_atlas':atlas,
      'conditional_timeout_plus_word_last_sample':CLOCK.MAX_STEPS,
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
    if x.get('open_obligations') != qualification_status():
        f.append('open-obligation status differs from supplying proof components')
    if set(x.get('open_obligations',{})) != set(OPEN_QUALIFICATIONS):
        f.append('open qualification inventory changed without proof')
    if x.get('falsified_prerequisites') != {}:
        f.append('obsolete counter falsification retained')
    f.extend(STARTDIST.validate(x.get('startup_disturbance_obstruction',{})))
    f.extend(SENSOR_CONTRACT.validate(x.get('startup_sensor_contract',{})))
    f.extend(WPEBOUNDS.validate(x.get('bounded_input_WPE_supplies',{})))
    f.extend(MAHONYTOTAL.validate(x.get('conditional_Mahony_prefix_totality',{})))
    f.extend(SVDAXIS.validate(x.get('pinned_Eigen_seed_axis_reduction',{})))
    if x.get('pinned_Eigen_source_uniform_SVD') != EIGENSVD.source_uniform_certificate():
        f.append('Eigen source-uniform certificate differs from supplying proof')
    if x.get('configured_frontend_uniform_supplies') != FRONTBOUNDS.build():
        f.append('frontend supply certificate differs from supplying proof')
    if x.get('configured_candidate_uniform_supplies') != CANDIDATE.build():
        f.append('candidate supply certificate differs from supplying proof')
    if x.get('source_owned_startup_CORE_constructor') != CORESTART.readiness():
        f.append('startup CORE construction certificate differs from supplying proof')
    f.extend(INPUT.validate(x.get('commissioned_Live_MEMS_input_contract',{})))
    f.extend(LIVEOVERFLOW.validate(x.get('bounded_raw_Live_input_totality_obstruction',{})))
    f.extend(SAMPLING.validate(x.get('continuous_to_sampled_startup_budget_audit',{})))
    f.extend(TIMEOUT.validate(x.get('physical_timeout_alignment_witness',{})))
    if x.get('pinned_target_library_profile') != TARGETLIBM.profile_correspondence():
        f.append('pinned library profile differs from target proof')
    if x.get('pinned_WPE_compiler_profile') != TARGETWPE.readiness():
        f.append('pinned WPE compiler profile differs from target proof')
    if x.get('library_closure_requires_separate_firmware_compiler_qualification') is not True:
        f.append('library approximation closure erased firmware compiler prerequisite')
    if x.get('magnetic_frame_bounds')!=CAPTURE.readiness():
        f.append('magnetic frame bound or accuracy qualification changed')
    if x.get('research_outcome') != 'finite_master_qualification_incomplete':
        f.append('research outcome differs from current qualification')
    if x.get('counter_saturation_certificate') != COUNTER.build():
        f.append('counter certificate differs from audited shipping recurrence')
    if x.get('finite_storage_status',{}).get('finite_word_counter_safety_closed') is not True:
        f.append('proved counter safety not consumed by storage guard')
    if not x.get('finite_storage_guard_error'):f.append('finite storage guard did not fail closed')
    for k in ('finite_master_guard_closed','storage_search_allowed','ALT_STARTUP_PASS','ALT_LIVE_PASS','ALT_END_TO_END_PASS'):
        if x.get(k) is not False:f.append(k+' not false')
    return f


def main():
    import argparse
    import json
    from pathlib import Path

    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    report=build()
    failures=validate(report)
    report['validation_failures']=failures
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,default=str)+'\n')
    print(json.dumps({'storage_search_allowed':report['storage_search_allowed'],
                      'research_outcome':report['research_outcome'],
                      'falsified_prerequisites':list(report['falsified_prerequisites']),
                      'open_obligations':[k for k,v in report['open_obligations'].items() if v],
                      'closed_qualifications':[k for k,v in report['open_obligations'].items() if not v],
                      'validation_failures':failures},sort_keys=True))
    return int(bool(failures))


if __name__=='__main__':
    raise SystemExit(main())
