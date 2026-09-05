from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_sea3_p3_joint_source_contract as joint  # noqa: E402
import ou3_sea3_rs_tau_lag_envelope as lag  # noqa: E402
import ou3_sea3_riccati_metric_p3 as mod  # noqa: E402


class Sea3RsInnovationP3Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()
        cls.j = joint.build()
        cls.lag = lag.build()

    def test_architecture_uses_four_s_rs_word_as_translation_strictness(self):
        self.assertEqual(mod.validate(self.d), [])
        self.assertEqual(
            self.d["canonical_P3_architecture"],
            "SEA3_RS_INNOVATION_DISSIPATION_WORD",
        )
        self.assertTrue(self.d["R_S_is_primary_translation_correction_mechanism"])
        self.assertTrue(self.d["pseudo_update_recurrence_is_primary_word_structure"])
        self.assertTrue(self.d["four_S_translation_word_consumed"])
        self.assertTrue(self.d["four_S_translation_observation_geometry_closed"])
        self.assertTrue(self.d["four_S_batch_noise_upper_closed"])
        self.assertEqual(
            self.d["translation_correction_word"]["mechanism"],
            "FOUR_SEPARATED_S_ZERO_INNOVATIONS",
        )
        self.assertFalse(
            self.d["translation_correction_word"]["accelerometer_needed_to_close_translation"]
        )
        self.assertEqual(self.d["strictness_location"], "RECURRENT_SEA3_MEASUREMENT_WORD")

    def test_target_and_applied_rs_are_not_conflated(self):
        self.assertTrue(self.d["tau_active_pseudo_cadence_coupling_consumed"])
        self.assertTrue(self.d["SpectralMSE_target_tau_sigma_TS_coupling_consumed"])
        self.assertTrue(self.d["applied_RS_separate_EMA_acknowledged"])
        self.assertFalse(self.d["instantaneous_RS_target_substituted_for_applied_RS"])
        self.assertTrue(self.d["safe_applied_RS_invariant_used_until_lag_theorem"])

    def test_joint_sea3_source_state_is_mandatory_for_quantitative_p3(self):
        self.assertEqual(joint.validate(self.j), [])
        self.assertEqual(
            self.j["canonical_architecture"],
            "SEA3_RS_INNOVATION_DISSIPATION_WORD",
        )
        self.assertTrue(
            self.j["canonical_independent_tau_sigma_RS_TS_extrema_product_forbidden"]
        )
        self.assertTrue(self.j["rectangular_full_box_calculation_diagnostic_only"])
        self.assertTrue(
            self.j["rectangular_full_box_failure_may_not_reject_canonical_architecture"]
        )
        self.assertEqual(
            self.j["adaptive_state"],
            [
                "tau_applied",
                "sigma_aw_filter",
                "R_S_applied",
                "pseudo_update_period",
                "scheduler_progress",
            ],
        )
        proved = self.j["proved_source_couplings"]
        self.assertTrue(proved["tau_sigma_share_same_sample_EMA_alpha"])
        self.assertTrue(proved["SpectralMSE_target_uses_same_target_tau_sigma_TS"])
        self.assertTrue(proved["R_S_has_separate_EMA"])
        self.assertTrue(proved["pseudo_period_is_clamped_monotone_function_of_applied_tau"])
        self.assertTrue(proved["pseudo_scheduler_progress_preserving"])
        self.assertTrue(proved["physical_height_period_cartesian_extrema_forbidden"])
        self.assertFalse(
            self.j["not_yet_promotable_source_couplings"][
                "sea_RAO_acceleration_coupling_may_be_used_as_hard_P3_pruning"
            ]
        )

    def test_low_dimensional_rs_tau_lag_is_certified_without_history_graph(self):
        self.assertEqual(lag.validate(self.lag), [])
        self.assertTrue(self.lag["source_generated_not_trajectory_fit"])
        self.assertFalse(self.lag["old_P2_800_state_graph_consumed"])
        self.assertFalse(self.lag["source_history_graph_consumed"])
        self.assertFalse(self.lag["predecessor_path_enumeration_consumed"])
        self.assertTrue(self.lag["candidate_snapshot_commit_preserves_invariant"])
        self.assertTrue(
            self.lag["SpectralMSE_fractional_power_removed_by_14th_power_identity"]
        )
        self.assertFalse(self.lag["ordinary_libm_fractional_power_used_in_pass_decision"])
        self.assertTrue(self.lag["target_lower_curve"]["pass"])
        self.assertTrue(self.lag["applied_invariant_lower_curve"]["pass"])
        self.assertTrue(self.lag["applied_invariant_lower_curve"]["initial_state_inside"])
        self.assertGreater(
            self.lag["applied_invariant_lower_curve"]["R_at_tau_max_lower"],
            self.lag["R_S_hard_floor"],
        )

    def test_rs_corrective_force_cannot_be_scalarized_away(self):
        req = self.j["R_S_corrective_force_requirements"]
        self.assertTrue(req["use_actual_applied_R_S_on_selected_pseudo_updates"])
        self.assertTrue(req["retain_per_axis_R_S_factors"])
        self.assertTrue(req["retain_full_P_column_S_cross_covariance_action"])
        self.assertTrue(req["credit_guaranteed_recurrent_S_updates_as_measurement_dissipation"])
        self.assertTrue(req["do_not_replace_R_S_correction_by_process_strictness"])
        self.assertTrue(
            req[
                "do_not_use_global_R_S_100_at_every_firing_in_final_canonical_matrix_if_joint_source_enclosure_is_available"
            ]
        )

    def test_exact_innovation_identity_is_canonical(self):
        self.assertTrue(self.d["exact_measurement_dissipation_identity_consumed"])
        self.assertTrue(self.d["batch_innovation_information_identity_consumed"])
        self.assertTrue(self.d["process_UCC_used_as_metric_lower_not_primary_strictness"])
        self.assertIn("V_minus - V_plus", self.d["innovation_identity"]["identity"])
        self.assertIn("D_W", self.d["batch_identity"]["correction_information"])

    def test_dead_end_routes_cannot_reenter(self):
        self.assertFalse(self.d["source_history_graph_consumed"])
        self.assertFalse(self.d["predecessor_path_enumeration_consumed"])
        self.assertFalse(self.d["old_P2_800_state_graph_consumed"])
        self.assertFalse(self.d["one_sample_strict_Riccati_margin_consumed"])
        self.assertFalse(self.d["commit_aligned_source_word_consumed"])
        self.assertFalse(self.d["per_sample_SPD_lower_required"])
        self.assertFalse(self.d["selected_process_mode_strictness_used"])
        self.assertFalse(self.d["determinant_trace_scalarization_used"])
        self.assertFalse(self.d["scalar_information_beta_used"])
        self.assertFalse(self.j["source_history_graph_consumed"])
        self.assertFalse(self.j["predecessor_path_enumeration_consumed"])
        self.assertFalse(self.j["old_P2_800_state_graph_consumed"])
        self.assertFalse(self.lag["source_history_graph_consumed"])
        self.assertFalse(self.lag["predecessor_path_enumeration_consumed"])
        self.assertFalse(self.lag["old_P2_800_state_graph_consumed"])
        self.assertEqual(self.d["useful_gate"], 1e-18)

    def test_gate_fails_closed_until_lw_and_full_matrix_composition(self):
        self.assertTrue(self.d["P3_RS_TRANSLATION_OBSERVATION_GEOMETRY_CLOSED"])
        self.assertTrue(self.d["P3_RS_BATCH_NOISE_UPPER_CLOSED"])
        self.assertFalse(self.d["P3_RS_WEIGHTED_WORD_INFORMATION_CLOSED"])
        self.assertFalse(self.d["P3_UCC_METRIC_LOWER_CLOSED"])
        self.assertFalse(self.d["P3_FULL_MATRIX_COMPARISON_CLOSED"])
        self.assertFalse(self.d["P3_CANONICAL_PASS"])
        self.assertFalse(self.d["P4_MAY_CONSUME_P3"])
        self.assertTrue(self.d["P3_CANONICAL_FAIL_REASONS"])
        for mode in ("H", "A"):
            self.assertFalse(self.d["modes"][mode]["pass"])
            self.assertEqual(
                self.d["modes"][mode]["relative_Riccati_injection_margin_lower"], 0.0
            )


if __name__ == "__main__":
    unittest.main()
