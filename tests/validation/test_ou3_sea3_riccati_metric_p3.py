from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_p3_full_preconditions as full  # noqa: E402
import ou3_sea3_p3_joint_source_contract as joint  # noqa: E402
import ou3_sea3_riccati_metric_p3 as mod  # noqa: E402
import ou3_sea3_windowed_vector_pe as pe  # noqa: E402


class Sea3FullNormalLiveP3Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()
        cls.f = full.build()
        cls.j = joint.build()
        cls.pe = pe.build()

    def test_complete_normal_live_word_is_canonical(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(full.validate(self.f), [])
        self.assertEqual(
            self.d["canonical_P3_architecture"],
            "SEA3_FULL_NORMAL_LIVE_RICCATI_WORD",
        )
        self.assertTrue(self.d["complete_precondition_contract_consumed"])
        self.assertTrue(self.d["one_common_event_word_required_for_H_and_A"])
        self.assertTrue(self.d["full_H18_A21_matrix_comparison_required"])
        self.assertGreaterEqual(self.d["common_word_horizon_s"], 3.0)
        self.assertTrue(self.d["same_joint_source_path_feeds_F_Q_TS_RS"])
        self.assertTrue(self.d["same_event_word_contains_accel_S_PE_and_aw_floor"])

    def test_every_machine_checkable_precondition_is_present(self):
        self.assertTrue(self.f["all_current_machine_checkable_preconditions_present"])
        self.assertTrue(self.f["mandatory_preconditions"])
        self.assertTrue(all(self.f["mandatory_preconditions"].values()))
        self.assertEqual(self.f["source_parity_failures"], [])
        self.assertEqual(self.f["paper_parity_failures"], [])
        c = self.f["final_numeric_contract"]
        self.assertEqual(c["H_dimension"], 18)
        self.assertEqual(c["A_dimension"], 21)
        self.assertEqual(c["useful_gate"], 1e-18)
        self.assertTrue(c["full_18x18_and_21x21_matrix_comparison_required"])
        self.assertTrue(c["same_event_word_contains_accel_S_PE_and_aw_floor"])
        self.assertTrue(c["actual_applied_per_axis_RS_required"])

    def test_windowed_pe_uses_recurrence_not_magnetometer_odr(self):
        self.assertEqual(pe.validate(self.pe), [])
        self.assertTrue(self.pe["paper_windowed_PE_semantics_consumed"])
        self.assertTrue(self.pe["asynchronous_magnetometer_semantics_consumed"])
        self.assertFalse(self.pe["hardware_magnetometer_ODR_used_as_PE_recurrence"])
        self.assertFalse(self.pe["two_consecutive_accepted_magnetic_packets_required"])
        self.assertTrue(
            self.pe["rejected_magnetic_packets_between_required_occurrences_allowed"]
        )
        self.assertTrue(self.pe["all_valid_accelerometer_packets_required"])
        self.assertFalse(self.pe["accelerometer_rejection_branch_present"])
        self.assertGreater(
            self.pe["eta6_information"]["alpha_6_information_lower"], 0.0
        )
        spread = self.pe["spread_occurrence_selection"]
        self.assertEqual(spread["first_occurrence_window_s"], [0.0, 1.0])
        self.assertEqual(spread["second_occurrence_window_s"], [2.0, 3.0])
        self.assertGreaterEqual(spread["separation_lower_s"], 1.0)
        self.assertGreaterEqual(spread["word_horizon_s"], 3.0)

    def test_a_mode_uses_declared_finite_bias_correlation(self):
        route = self.pe["A_mode_bias_route"]
        self.assertFalse(route["uses_eta9_pointwise_packet_shortcut"])
        self.assertTrue(route["uses_eta6_plus_finite_bias_correlation"])
        self.assertGreater(route["accel_bias_tau_s"], 0.0)
        self.assertGreater(route["homogeneous_bias_contraction_gap_lower"], 0.0)
        self.assertLess(route["homogeneous_bias_contraction_upper_over_word"], 1.0)
        self.assertTrue(self.d["A_mode_finite_bias_correlation_consumed"])

    def test_four_s_is_required_component_not_whole_architecture(self):
        self.assertTrue(self.d["R_S_translation_component_consumed"])
        self.assertTrue(self.d["R_S_is_primary_translation_correction_mechanism"])
        self.assertTrue(self.d["four_S_translation_word_consumed"])
        self.assertTrue(self.d["four_S_translation_observation_geometry_closed"])
        self.assertTrue(self.d["four_S_batch_noise_upper_closed"])
        self.assertTrue(self.d["R_S_component_is_not_the_whole_P3_architecture"])
        self.assertTrue(self.d["actual_applied_per_axis_RS_required_in_final_word"])
        self.assertTrue(self.d["all_valid_accelerometer_updates_consumed"])
        self.assertTrue(
            self.d["full_accelerometer_attitude_aw_ba_cross_block_information_required"]
        )
        self.assertTrue(self.d["windowed_asynchronous_vector_PE_consumed"])
        self.assertTrue(self.d["full_process_UCC_consumed"])
        self.assertTrue(self.d["aw_covariance_floor_PSD_events_consumed"])

    def test_joint_source_cannot_be_rectangularized(self):
        self.assertEqual(joint.validate(self.j), [])
        self.assertTrue(
            self.j["canonical_independent_tau_sigma_RS_TS_extrema_product_forbidden"]
        )
        self.assertTrue(self.j["rectangular_full_box_calculation_diagnostic_only"])
        self.assertTrue(
            self.j["rectangular_full_box_failure_may_not_reject_canonical_architecture"]
        )
        req = self.j["R_S_corrective_force_requirements"]
        self.assertTrue(req["use_actual_applied_R_S_on_selected_pseudo_updates"])
        self.assertTrue(req["retain_per_axis_R_S_factors"])
        self.assertTrue(req["retain_full_P_column_S_cross_covariance_action"])
        self.assertTrue(req["credit_guaranteed_recurrent_S_updates_as_measurement_dissipation"])
        self.assertTrue(req["do_not_replace_R_S_correction_by_process_strictness"])

    def test_physical_sea3_scope_is_used_but_not_overclaimed(self):
        phys = self.f["physical_SEA3_scope"]
        self.assertTrue(phys["parameter_height_period_partition_coupling_consumed"])
        self.assertFalse(phys["global_finite_window_realization_left_inclusion_closed"])
        self.assertFalse(phys["unqualified_RAO_coupling_used_as_hard_pruning"])
        self.assertTrue(phys["canonical_P3_is_conditional_on_admitted_Normal_Live_SEA3_word"])
        self.assertTrue(self.d["SEA3_height_period_partition_coupling_consumed"])
        self.assertFalse(self.d["unqualified_RAO_coupling_used_as_hard_pruning"])
        self.assertFalse(self.d["global_physical_SEA3_left_inclusion_claimed"])

    def test_reduced_dead_end_routes_cannot_promote_p3(self):
        for key in (
            "hardware_magnetometer_ODR_used_as_PE_recurrence",
            "two_consecutive_accepted_magnetic_packets_required",
            "aw_covariance_floor_marginal_Loewner_shortcut_used",
            "source_history_graph_consumed",
            "predecessor_path_enumeration_consumed",
            "old_P2_800_state_graph_consumed",
            "one_sample_strict_Riccati_margin_consumed",
            "per_sample_SPD_lower_required",
            "selected_process_mode_strictness_used",
            "determinant_trace_scalarization_used",
            "scalar_information_beta_used",
            "blockwise_minimum_ratio_used_for_final_gate",
            "independent_tau_sigma_RS_TS_extrema_product_used",
        ):
            self.assertFalse(self.d[key], key)
        self.assertEqual(self.d["useful_gate"], 1e-18)

    def test_gate_stays_fail_closed_until_one_full_matrix_word_exists(self):
        self.assertFalse(self.d["P3_FULL_WORD_ENCLOSED"])
        self.assertFalse(self.d["P3_FULL_MATRIX_COMPARISON_CLOSED"])
        self.assertFalse(self.d["P3_CANONICAL_PASS"])
        self.assertFalse(self.d["P4_MAY_CONSUME_P3"])
        self.assertTrue(self.d["P3_CANONICAL_FAIL_REASONS"])
        for mode, dim in (("H", 18), ("A", 21)):
            self.assertEqual(self.d["modes"][mode]["dimension"], dim)
            self.assertFalse(self.d["modes"][mode]["pass"])
            self.assertEqual(
                self.d["modes"][mode]["relative_Riccati_injection_margin_lower"],
                0.0,
            )


if __name__ == "__main__":
    unittest.main()
