from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_complete_source as complete  # noqa: E402
import ou3_sea3_p3_full_preconditions as full  # noqa: E402
import ou3_sea3_riccati_metric_p3 as mod  # noqa: E402
import ou3_sea3_windowed_vector_pe as pe  # noqa: E402


class Sea3CompleteSourceP3Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()
        cls.c = complete.build()
        cls.f = full.build()
        cls.pe = pe.build()

    def test_complete_sea3_is_the_only_canonical_source(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(complete.validate(self.c), [])
        self.assertEqual(full.validate(self.f), [])
        self.assertEqual(
            self.d["canonical_P3_architecture"],
            "COMPLETE_SEA3_FULL_NORMAL_LIVE_RICCATI_WORD",
        )
        self.assertEqual(self.d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(self.d["complete_SEA3_source_consumed"])
        self.assertTrue(self.d["complete_SEA3_response_couplings_consumed"])
        self.assertTrue(self.d["complete_SEA3_finite_horizon_good_event_consumed"])
        self.assertTrue(self.d["complete_SEA3_frontend_state_consumed"])
        self.assertTrue(self.d["complete_SEA3_adaptive_state_consumed"])
        self.assertTrue(self.d["same_complete_SEA3_word_used_for_H18_and_A21"])
        self.assertGreaterEqual(self.d["common_word_horizon_s"], 3.0)

    def test_complete_source_retains_all_sea3_couplings(self):
        sea = self.c["SEA3_surface_family"]
        self.assertEqual(sea["modes_max"], 3)
        self.assertEqual(sea["gamma_interval"], [1.0, 7.0])
        self.assertTrue(sea["independent_H_T_extrema_forbidden"])
        self.assertTrue(sea["independent_partition_height_maxima_forbidden"])
        response = self.c["SEA3_response_couplings"]
        self.assertTrue(response["independent_sea_x_RAO_cartesian_product_forbidden"])
        self.assertTrue(response["same_H_s_partition_energy_enters_translation_and_rotation_conditions"])
        self.assertTrue(response["only_jointly_admitted_response_tuples_may_generate_P3_words"])
        finite = self.c["SEA3_finite_horizon_good_event"]
        self.assertTrue(finite["combined_within_budget"])
        self.assertTrue(finite["same_response_word_must_satisfy_both_good_events"])

    def test_tuner_and_rs_are_derived_from_same_sea3_word(self):
        a = self.c["derived_adaptive_source"]
        self.assertFalse(a["primitive_independent_tau_sigma_RS_TS"])
        self.assertTrue(a["same_SEA3_frontend_path_generates_tau_sigma_RS_targets"])
        self.assertTrue(a["same_candidate_snapshot_commits_tau_sigma_RS"])
        self.assertTrue(a["T_S_is_function_of_same_committed_tau"])
        self.assertTrue(a["Q_uses_same_committed_tau_sigma"])
        self.assertTrue(a["rate_bounds_are_constraints_on_SEA3_derived_path_not_a_word_generator"])
        rs = self.c["R_S_regularizer"]
        self.assertEqual(rs["source_parity_failures"], [])
        self.assertEqual(rs["axis_std_factors"], [0.72, 0.72, 1.0])
        self.assertTrue(rs["actual_applied_R_S_required_at_every_due_S_update"])
        self.assertTrue(rs["all_due_S_updates_remain_in_full_word"])
        self.assertTrue(rs["full_P_column_S_cross_covariance_action_required"])
        self.assertTrue(rs["selected_four_S_subset_may_not_replace_full_scheduler_word"])
        self.assertTrue(rs["R_S_may_not_be_replaced_by_process_strictness"])

    def test_no_fallback_generator_or_gate_exists(self):
        self.assertTrue(self.d["no_fallback_route_enabled"])
        self.assertTrue(all(v is False for v in self.c["no_fallback_generators"].values()))
        self.assertTrue(all(v is False for v in self.f["no_fallback_generators"].values()))
        for key in (
            "independent_tau_sigma_RS_TS_extrema_product_used",
            "independent_sea_x_RAO_product_used",
            "point_source_word_used",
            "selected_four_S_word_used",
            "D_W_L_W_split_used",
            "blockwise_minimum_ratio_used",
            "scalar_information_beta_used",
            "determinant_trace_scalarization_used",
            "source_history_graph_used",
            "predecessor_path_enumeration_used",
            "arbitrary_P0_rectangle_used",
            "selected_process_mode_strictness_used",
        ):
            self.assertFalse(self.d[key], key)

    def test_full_literal_event_families_are_mandatory(self):
        self.assertTrue(self.d["actual_applied_per_axis_RS_consumed"])
        self.assertTrue(self.d["all_due_S_updates_required"])
        self.assertTrue(self.d["all_valid_accelerometer_updates_required"])
        self.assertTrue(self.d["asynchronous_vector_PE_required"])
        self.assertTrue(self.d["all_full_process_Q_required"])
        self.assertTrue(self.d["all_aw_covariance_floor_events_required"])
        self.assertTrue(self.d["joint_P_Psi_Omega_backend_consumed"])
        self.assertTrue(self.d["literal_full_word_assembler_consumed"])
        self.assertTrue(self.d["literal_shipping_event_order_pass"])

    def test_windowed_pe_is_asynchronous_and_complete(self):
        self.assertEqual(pe.validate(self.pe), [])
        self.assertFalse(self.pe["hardware_magnetometer_ODR_used_as_PE_recurrence"])
        self.assertFalse(self.pe["two_consecutive_accepted_magnetic_packets_required"])
        self.assertTrue(self.pe["rejected_magnetic_packets_between_required_occurrences_allowed"])
        self.assertTrue(self.pe["all_valid_accelerometer_packets_required"])
        self.assertFalse(self.pe["accelerometer_rejection_branch_present"])
        self.assertGreater(self.pe["eta6_information"]["alpha_6_information_lower"], 0.0)
        self.assertGreater(self.pe["A_mode_bias_route"]["homogeneous_bias_contraction_gap_lower"], 0.0)

    def test_precondition_contract_has_only_complete_sea3_word(self):
        self.assertTrue(self.f["all_current_machine_checkable_preconditions_present"])
        self.assertTrue(all(self.f["mandatory_preconditions"].values()))
        self.assertEqual(self.f["source_parity_failures"], [])
        self.assertEqual(self.f["paper_parity_failures"], [])
        self.assertEqual(self.f["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        c = self.f["final_numeric_contract"]
        self.assertEqual(c["H_dimension"], 18)
        self.assertEqual(c["A_dimension"], 21)
        self.assertEqual(c["useful_gate"], 1e-18)
        self.assertTrue(c["actual_applied_per_axis_RS_required"])
        self.assertTrue(c["full_18x18_and_21x21_matrix_comparison_required"])
        self.assertIn("Omega_W - delta*P_W", c["required_final_inequality"])

    def test_gate_remains_fail_closed_until_full_sea3_family_executes(self):
        self.assertFalse(self.d["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"])
        self.assertFalse(self.d["P3_FULL_WORD_ENCLOSED"])
        self.assertFalse(self.d["P3_FULL_MATRIX_COMPARISON_CLOSED"])
        self.assertFalse(self.d["P3_CANONICAL_PASS"])
        self.assertFalse(self.d["P4_MAY_CONSUME_P3"])
        self.assertTrue(self.d["P3_CANONICAL_FAIL_REASONS"])
        for mode, dim in (("H18", 18), ("A21", 21)):
            self.assertEqual(self.d["modes"][mode]["dimension"], dim)
            self.assertFalse(self.d["modes"][mode]["full_word_executed"])
            self.assertFalse(self.d["modes"][mode]["Omega_minus_delta_P_ldlt_closed"])
            self.assertEqual(self.d["modes"][mode]["certified_delta_lower"], 0.0)


if __name__ == "__main__":
    unittest.main()
