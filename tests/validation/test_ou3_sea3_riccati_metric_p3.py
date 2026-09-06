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
        self.assertEqual(
            self.d["canonical_P3_topology"],
            "H18_3S_PRIOR_FREE_THEN_PRESERVED_H_TO_A_HYBRID_A21",
        )
        self.assertEqual(self.d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(self.d["complete_SEA3_source_consumed"])
        self.assertTrue(self.d["complete_SEA3_response_couplings_consumed"])
        self.assertTrue(self.d["complete_SEA3_compact_parameter_domain_consumed"])
        self.assertTrue(self.d["complete_SEA3_compact_transition_relation_consumed"])
        self.assertTrue(self.d["complete_SEA3_phase_continuous_realization_required"])
        self.assertTrue(self.d["same_xs_lambda_drives_entire_execution"])
        self.assertTrue(self.d["stochastic_forcing_does_not_generate_source_words"])
        self.assertTrue(self.d["stochastic_forcing_does_not_prune_homogeneous_family"])
        self.assertTrue(self.d["complete_SEA3_frontend_state_consumed"])
        self.assertTrue(self.d["complete_SEA3_adaptive_state_consumed"])
        self.assertGreaterEqual(self.d["common_word_horizon_s"], 3.0)

    def test_complete_source_retains_all_sea3_couplings(self):
        sea = self.c["SEA3_surface_family"]
        self.assertEqual(sea["modes_max"], 3)
        self.assertEqual(sea["gamma_interval"], [1.0, 7.0])
        self.assertTrue(sea["parameter_domain_compact"])
        self.assertTrue(sea["compact_transition_relation_is_theorem_domain"])
        self.assertTrue(sea["independent_H_T_extrema_forbidden"])
        self.assertTrue(sea["independent_partition_height_maxima_forbidden"])

        dynamic = self.c["SEA3_dynamic_realization"]
        self.assertTrue(dynamic["phase_continuous"])
        self.assertTrue(
            dynamic["same_realization_drives_translation_rotation_frontend_tuner_geometry"]
        )
        self.assertFalse(dynamic["probabilistic_event_may_substitute_for_realization"])
        self.assertFalse(dynamic["arbitrary_bounded_input_may_substitute_for_realization"])

        response = self.c["SEA3_response_couplings"]
        self.assertTrue(response["independent_sea_x_RAO_cartesian_product_forbidden"])
        self.assertTrue(response["same_H_s_partition_energy_enters_translation_and_rotation_conditions"])
        self.assertTrue(response["only_same_phase_continuous_SEA3_realization_may_generate_P3_words"])
        self.assertTrue(response["moment_or_probability_bound_may_not_generate_P3_word"])

    def test_tuner_and_rs_are_derived_from_same_sea3_execution(self):
        a = self.c["derived_adaptive_source"]
        self.assertFalse(a["primitive_independent_tau_sigma_RS_TS"])
        self.assertTrue(a["same_SEA3_frontend_path_generates_tau_sigma_RS_targets"])
        self.assertTrue(a["same_candidate_snapshot_commits_tau_sigma_RS"])
        self.assertTrue(a["T_S_is_function_of_same_committed_tau"])
        self.assertTrue(a["Q_uses_same_committed_tau_sigma"])
        rs = self.c["R_S_regularizer"]
        self.assertEqual(rs["source_parity_failures"], [])
        self.assertEqual(rs["axis_std_factors"], [0.72, 0.72, 1.0])
        self.assertTrue(rs["actual_applied_R_S_required_at_every_due_S_update"])
        self.assertTrue(rs["all_due_S_updates_remain_in_full_word"])

    def test_no_fallback_generator_or_gate_exists(self):
        self.assertTrue(self.d["no_fallback_route_enabled"])
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
            "eta9_point_packet_shortcut_used",
        ):
            self.assertFalse(self.d[key], key)

    def test_literal_events_reset_and_measurements_are_all_retained(self):
        self.assertTrue(self.d["actual_applied_per_axis_RS_consumed"])
        self.assertTrue(self.d["all_due_S_updates_required"])
        self.assertTrue(self.d["all_valid_accelerometer_updates_required"])
        self.assertFalse(self.d["accelerometer_rejection_after_certified_Normal_Live_allowed"])
        self.assertTrue(self.d["asynchronous_vector_PE_required"])
        self.assertTrue(self.d["all_full_process_Q_required"])
        self.assertTrue(self.d["all_aw_covariance_floor_events_required"])
        self.assertTrue(self.d["joint_P_Psi_Omega_backend_consumed"])
        self.assertTrue(self.d["literal_full_word_assembler_consumed"])
        self.assertTrue(self.d["literal_shipping_event_order_pass"])
        self.assertTrue(self.d["reset_complete_literal_execution_consumed"])
        self.assertTrue(self.d["immediate_left_error_reset_congruence_consumed"])

    def test_hybrid_topology_matches_shipping_release(self):
        self.assertTrue(self.d["same_complete_SEA3_execution_continues_across_H_to_A"])
        self.assertFalse(self.d["same_three_second_same_mode_word_used_for_H18_and_A21"])
        self.assertTrue(self.d["H_to_A_is_separate_dimension_changing_hybrid_event"])
        self.assertTrue(self.d["shipping_H_mode_hold_guarantees_H18_before_A_release"])
        self.assertEqual(
            "PRIOR_FREE_18X18_INTERVAL_LDLT",
            self.d["modes"]["H18"]["closure_method"],
        )
        self.assertEqual(
            "EXACT_H18_TO_A21_HYBRID_DIRECT_SUM_FULL_MATRIX",
            self.d["modes"]["A21"]["closure_method"],
        )

    def test_windowed_pe_is_asynchronous_and_complete(self):
        self.assertEqual(pe.validate(self.pe), [])
        self.assertFalse(self.pe["hardware_magnetometer_ODR_used_as_PE_recurrence"])
        self.assertFalse(self.pe["two_consecutive_accepted_magnetic_packets_required"])
        self.assertTrue(self.pe["rejected_magnetic_packets_between_required_occurrences_allowed"])
        self.assertTrue(self.pe["all_valid_accelerometer_packets_required"])
        self.assertFalse(self.pe["accelerometer_rejection_branch_present"])
        self.assertGreater(self.pe["eta6_information"]["alpha_6_information_lower"], 0.0)

    def test_universal_certificate_chain_closes_without_finite_materialization(self):
        self.assertFalse(self.d["finite_source_family_materialization_required"])
        self.assertFalse(self.d["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"])
        self.assertTrue(
            self.d["universal_complete_SEA3_certificate_chain_used_instead_of_finite_materialization"]
        )
        self.assertTrue(self.d["UNIVERSAL_COMPLETE_SEA3_CERTIFICATE_CHAIN_CLOSED"])
        self.assertTrue(self.d["P3_FULL_WORD_ENCLOSED"])
        self.assertTrue(self.d["P3_FULL_MATRIX_COMPARISON_CLOSED"])
        self.assertTrue(self.d["P3_CANONICAL_PASS"])
        self.assertTrue(self.d["P4_MAY_CONSUME_P3"])
        self.assertEqual([], self.d["P3_CANONICAL_FAIL_REASONS"])
        for mode, dim in (("H18", 18), ("A21", 21)):
            self.assertEqual(self.d["modes"][mode]["dimension"], dim)
            self.assertTrue(self.d["modes"][mode]["Omega_minus_delta_P_full_matrix_closed"])
            self.assertGreaterEqual(self.d["modes"][mode]["certified_delta_lower"], 1.0e-18)
            self.assertGreaterEqual(
                self.d["modes"][mode]["relative_Riccati_injection_margin_lower"],
                1.0e-18,
            )

    def test_global_physical_left_inclusion_remains_separate(self):
        self.assertTrue(self.d["global_physical_deployment_left_inclusion_is_separate_obligation"])
        self.assertFalse(self.d["global_physical_deployment_left_inclusion_closed_here"])


if __name__ == "__main__":
    unittest.main()
