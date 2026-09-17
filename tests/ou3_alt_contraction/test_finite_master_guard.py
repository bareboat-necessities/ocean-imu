import copy
import unittest
from tools.stability.ou3_alt_contraction import finite_master_guard as G

class Tests(unittest.TestCase):
    def test_closed_progress_is_consumed_but_storage_stays_blocked(self):
        x=G.build()
        self.assertEqual(G.validate(x),[])
        for v in x['closed_subobligations'].values(): self.assertTrue(v)
        self.assertIn('WPE_log_and_exp_libm_correspondence',x['open_obligations'])
        self.assertIn('complete_source_uniform_600_step_word',x['open_obligations'])
        self.assertNotIn('all_admitted_ungauged_entries_have_a_finite_attitude_representation',x['open_obligations'])
        self.assertTrue(x['finite_storage_status']['all_admitted_startup_entries_represented'])
        self.assertTrue(x['closed_subobligations']['all_nonzero_fresh_attitudes_represented_by_joint24_atlas'])
        self.assertNotIn('signed_magnetic_counter_safety_on_the_finite_word',x['open_obligations'])
        self.assertEqual(set(x['falsified_prerequisites']),{'declared_joint24_common_storage_contraction'})
        self.assertNotIn('WPE_frequency_source_uniform_supply',x['open_obligations'])
        self.assertNotIn('WPE_machine_vs_exact_period_branch_robustness',x['open_obligations'])
        self.assertTrue(x['closed_subobligations']['independent_WPE_machine_and_exact_branches_composed'])
        self.assertTrue(x['closed_subobligations']['source_uniform_clamped_WPE_frequency_discrepancy'])
        self.assertTrue(x['closed_subobligations']['literal_per_compiler_WPE_usable_latch'])
        self.assertTrue(x['finite_storage_status']['finite_word_counter_safety_closed'])
        self.assertEqual(x['research_outcome'],'declared_joint24_contraction_falsified_before_master_completion')
        self.assertNotIn('startup_accumulation_to_handoff_full_frame_bound',x['open_obligations'])
        self.assertTrue(x['closed_subobligations']['full_magnetic_frame_represented_without_small_angle_capture'])
        self.assertFalse(x['magnetic_frame_bounds']['small_frame_accuracy_required_by_finite_atlas_word'])
        self.assertFalse(x['magnetic_frame_bounds']['full_accumulation_to_handoff_frame_bound_source_qualified'])
        self.assertEqual(x['conditional_timeout_plus_word_last_sample'],30602)
        self.assertEqual(x['startup_entry_obstruction']['Cayley_denominator_at_south'],0)
        self.assertTrue(x['startup_disturbance_obstruction']['finite_prefix_no_initialization_induction'])
        self.assertFalse(x['startup_disturbance_obstruction']['post_Live_ISS_bound_supplies_startup_raw_norm_premise'])
        self.assertTrue(x['startup_sensor_contract']['sensor_profile_selected'])
        self.assertTrue(x['startup_sensor_contract']['source_uniform_seed_norm_margin_closed'])
        self.assertFalse(x['startup_sensor_contract']['old_Mahony_invariant_covers_new_profiles'])
        self.assertTrue(x['bounded_input_WPE_supplies']['reset_to_every_finite_prefix_bounded_input_induction_closed'])
        self.assertIn('WPE_source_uniform_machine_supply_bounds',x['open_obligations'])
        self.assertFalse(x['open_obligations']['WPE_source_uniform_machine_supply_bounds'])
        self.assertFalse(x['open_obligations']['WPE_machine_target_libm_and_compiler_selection'])
        self.assertFalse(x['open_obligations']['WPE_log_and_exp_libm_correspondence'])
        self.assertFalse(x['open_obligations']['Qaxis_exp_libm_correspondence'])
        self.assertTrue(x['open_obligations']['target_libm_and_compiler_profile_correspondence'])
        self.assertIn('source_uniform_timeout_aligned_branch_reachability',x['open_obligations'])
        self.assertFalse(x['open_obligations']['near_antiparallel_Eigen_JacobiSVD_solver_correspondence'])
        self.assertFalse(x['finite_master_guard_closed'])
        self.assertFalse(x['storage_search_allowed'])
        self.assertIn('finite-state storage blocked',x['finite_storage_guard_error'])
        x['startup_disturbance_obstruction']['arbitrary_bounded_residuals_imply_universal_startup']=True
        self.assertIn('arbitrary_bounded_residuals_imply_universal_startup differs from startup residual obstruction',G.validate(x))

    def test_closed_counter_prerequisite_must_reach_master_status(self):
        x=G.build()
        x['finite_storage_status']['finite_word_counter_safety_closed']=False
        self.assertIn('proved counter safety not consumed by storage guard',G.validate(x))

    def test_conditional_component_evidence_does_not_remove_universal_obligations(self):
        x=G.build()
        self.assertEqual(len(x['open_obligations']),11)
        self.assertTrue(x['commissioned_Live_MEMS_input_contract']['prefix']['all_initialized_finite_prefixes_totality_closed'])
        self.assertTrue(x['library_closure_requires_separate_firmware_compiler_qualification'])
        self.assertTrue(x['conditional_Mahony_prefix_totality']['ordinary_seed_scalar_prefix_totality_closed'])
        self.assertFalse(x['conditional_Mahony_prefix_totality']['every_startup_branch_totality_closed'])
        self.assertTrue(x['pinned_Eigen_seed_axis_reduction']['returning_solver_axis_equals_QR_Q_column_2'])
        self.assertFalse(x['pinned_Eigen_seed_axis_reduction']['source_uniform_Jacobi_loop_termination_closed'])
        self.assertTrue(x['pinned_Eigen_source_uniform_SVD']['source_uniform_QR_and_Jacobi_totality_closed'])
        self.assertTrue(x['configured_frontend_uniform_supplies']['all_finite_prefixes_from_literal_frontend_reset_covered'])
        self.assertTrue(x['configured_candidate_uniform_supplies']['candidate_arithmetic_totality_under_named_pow_range_closed'])
        self.assertTrue(x['configured_candidate_uniform_supplies']['pinned_powf_positive_finite_range_under_scalar_profile_closed'])
        self.assertTrue(x['source_owned_startup_CORE_constructor']['ungauged_machine_proxy_to_full_CORE_producer_available'])
        self.assertTrue(x['pinned_WPE_compiler_profile']['pinned_WPE_compiler_profile_selection_closed'])
        self.assertFalse(x['bounded_raw_Live_input_totality_obstruction']['nonfinite_execution_has_produced_RestrictedForcing'])
        del x['open_obligations']['complete_source_uniform_600_step_word']
        self.assertIn('open qualification inventory changed without proof',G.validate(x))

    def test_library_component_cannot_erase_global_compiler_gate(self):
        x=G.build()
        x['open_obligations']['target_libm_and_compiler_profile_correspondence']=False
        self.assertIn('open-obligation status differs from supplying proof components',G.validate(x))


class RhoFalsificationTests(unittest.TestCase):
    """One shared build; G.build() is expensive and these only read it."""
    @classmethod
    def setUpClass(cls): cls.report=G.build()
    def setUp(self): self.x=copy.deepcopy(self.report)

    def test_complete_word_rho_falsification_is_carried(self):
        rho=self.x['complete_word_rho_feasibility']
        self.assertEqual(rho['qualification'],'OU3_ALT_COMPLETE_WORD_RHO_FEASIBILITY_DIAGNOSTIC_V1')
        self.assertTrue(rho['declared_joint24_contraction_falsified'])
        self.assertEqual(rho['failure_classification'],'theorem_failure')
        self.assertEqual(rho['rho_floor_over_legal_words'],1.0)
        self.assertEqual(rho['rho_floor_source'],'exact_unipotent_subspace_algebra')
        self.assertFalse(rho['storage_search_allowed'])
        falsified=self.x['falsified_prerequisites']['declared_joint24_common_storage_contraction']
        self.assertEqual(falsified['classification'],'theorem_failure')
        self.assertIn('theta_z',falsified['limiting_state_direction'])
        self.assertIn('bg_z',falsified['limiting_state_direction'])

    def test_falsification_does_not_close_or_remove_a_qualification(self):
        # The eleven names are orthogonal to the falsification and must not be
        # laundered into a closed sub-obligation by it.
        self.assertEqual(len(self.x['open_obligations']),11)
        self.assertNotIn('declared_joint24_common_storage_contraction',self.x['closed_subobligations'])

    def test_a_dropped_falsification_is_rejected(self):
        self.x['falsified_prerequisites']={}
        self.assertIn('falsification inventory differs from the supplying rho diagnostic',G.validate(self.x))

    def test_a_promoting_rho_report_is_rejected(self):
        self.x['complete_word_rho_feasibility']['storage_search_allowed']=True
        self.assertIn('rho feasibility diagnostic attempted to unlock storage',G.validate(self.x))


if __name__=='__main__': unittest.main()
