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
        self.assertEqual(x['falsified_prerequisites'],{})
        self.assertNotIn('WPE_frequency_source_uniform_supply',x['open_obligations'])
        self.assertNotIn('WPE_machine_vs_exact_period_branch_robustness',x['open_obligations'])
        self.assertTrue(x['closed_subobligations']['independent_WPE_machine_and_exact_branches_composed'])
        self.assertTrue(x['closed_subobligations']['source_uniform_clamped_WPE_frequency_discrepancy'])
        self.assertTrue(x['closed_subobligations']['literal_per_compiler_WPE_usable_latch'])
        self.assertTrue(x['finite_storage_status']['finite_word_counter_safety_closed'])
        self.assertEqual(x['research_outcome'],'finite_master_qualification_incomplete')
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
        self.assertIn('source_uniform_timeout_aligned_branch_reachability',x['open_obligations'])
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
        self.assertTrue(x['conditional_Mahony_prefix_totality']['ordinary_seed_scalar_prefix_totality_closed'])
        self.assertFalse(x['conditional_Mahony_prefix_totality']['every_startup_branch_totality_closed'])
        self.assertTrue(x['pinned_Eigen_seed_axis_reduction']['returning_solver_axis_equals_QR_Q_column_2'])
        self.assertFalse(x['pinned_Eigen_seed_axis_reduction']['source_uniform_Jacobi_loop_termination_closed'])
        self.assertFalse(x['bounded_raw_Live_input_totality_obstruction']['nonfinite_execution_has_produced_RestrictedForcing'])
        del x['open_obligations']['complete_source_uniform_600_step_word']
        self.assertIn('open qualification inventory changed without proof',G.validate(x))

if __name__=='__main__': unittest.main()
