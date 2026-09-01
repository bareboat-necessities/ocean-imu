from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_operation_matched_joseph_ledger as LEDGER


class OperationMatchedJosephLedgerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = LEDGER.build()

    def test_ledger_validates_but_does_not_promote_complete_p4(self):
        d = self.d
        self.assertEqual(LEDGER.validate(d), [])
        self.assertTrue(d["strong_route_operation_ledger_closed"])
        self.assertFalse(d["full_state_directional_word_credit_established_here"])
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE"])
        self.assertFalse(d["P4_USABLE_CERTIFICATE_PROMOTED"])
        self.assertFalse(d["P5_FINITE_INNER_CAPTURE_ESTABLISHED_HERE"])

    def test_stronger_route_preserves_deployment_entrance_and_domain(self):
        d = self.d
        self.assertEqual(d["declared_P5_entrance_angle_deg"], 45.0)
        self.assertTrue(d["P5_45DEG_ENTRANCE_PRESERVED"])
        self.assertFalse(d["declared_domain_changed"])
        self.assertFalse(d["filter_changed"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["candidate_angle_reduction_used_for_closure"])
        self.assertFalse(d["aw_sigma_consistency_assumption_used"])

    def test_sector_invariance_failure_is_reproduced_but_not_used_as_gate(self):
        d = self.d
        self.assertTrue(d["legacy_sector_invariance_obstruction_reproduced"])
        self.assertTrue(d["legacy_sector_budget_distance_only"])
        self.assertFalse(d["per_operation_sector_invariance_is_P4_promotion_gate"])
        self.assertIsNotNone(d["conditional_25deg_consistency_constant_diagnostic"])
        self.assertGreater(d["conditional_25deg_consistency_constant_diagnostic"], 0.0)
        self.assertLess(d["conditional_25deg_consistency_constant_diagnostic"], 1.0)

    def test_directional_packet_is_rank_five_and_not_scalarized(self):
        d = self.d
        self.assertEqual(d["directional_packet_rank_exact"], 5)
        self.assertFalse(d["instantaneous_scalar_full_state_packet_margin_valid"])
        self.assertTrue(d["directional_PSD_word_accumulation_required"])
        self.assertFalse(d["full_state_directional_word_credit_established_here"])

    def test_accelerometer_uses_effective_aw_range_but_keeps_signed_eta_term(self):
        d = self.d
        self.assertTrue(d["accelerometer_eta_absorbed_for_state_correction_range"])
        self.assertFalse(d["accelerometer_large_aw_error_charged_as_independent_measurement_eta"])
        self.assertFalse(d["accelerometer_finite_angle_eta_penalty_dropped_from_Joseph_identity"])
        self.assertFalse(d["accelerometer_finite_angle_eta_independent_norm_budget_used"])
        ledger = {row["operation"]: row for row in d["operation_ledger"]}
        acc = ledger["accelerometer_accepted"]
        self.assertFalse(acc["large_declared_aw_error_is_measurement_eta"])
        self.assertFalse(acc["finite_angle_eta_penalty_dropped"])
        self.assertFalse(acc["finite_angle_eta_independently_maximized"])
        self.assertTrue(acc["latent_aw_rotation_is_norm_preserving"])
        self.assertFalse(acc["instantaneous_attitude_only_credit_promoted"])
        self.assertFalse(acc["sector_invariance_required"])
        self.assertIn("z_eff", acc["effective_coordinate"])
        self.assertIn("y^T S^-1 y", acc["exact_information_decrease"])
        self.assertIn("eta_a^T R_a^-1 eta_a", acc["exact_information_decrease"])

    def test_magnetometer_and_reset_keep_exact_structure(self):
        d = self.d
        ledger = {row["operation"]: row for row in d["operation_ledger"]}
        mag = ledger["magnetometer_accepted"]
        self.assertEqual(mag["radial_gain_action"], "EXACTLY_ZERO")
        self.assertTrue(mag["effective_coordinate_nonexpansive"])
        self.assertFalse(mag["radial_eta_independent_penalty_used"])
        self.assertFalse(mag["sector_invariance_required"])
        reset = ledger["quaternion_injection_and_left_error_reset"]
        self.assertEqual(reset["covariance_transport"], "EXACT_CONGRUENCE")
        self.assertEqual(reset["reset_inverse_operator_norm_upper"], 1.0)
        self.assertFalse(reset["condition_number_multiplier_used"])
        self.assertEqual(reset["remaining_nonlinear_term"], "rho=z_exact-G_ext*t")

    def test_vector_pair_constant_is_only_attitude_geometry_prerequisite(self):
        x = self.d["finite_angle_vector_pair_attitude_geometry_vs_goLive_metric_lower"]
        self.assertTrue(math.isfinite(x))
        self.assertGreater(x, 0.0)
        self.assertTrue(self.d["finite_angle_vector_pair_attitude_geometry_strict"])
        self.assertFalse(self.d["full_state_directional_word_credit_established_here"])

    def test_next_obligation_is_word_level_directional_accumulation(self):
        text = self.d["next_obligation"]
        self.assertIn("PSD directional operation forms", text)
        self.assertIn("signed nonlinear eta forms", text)
        self.assertIn("18/21-state", text)
        self.assertIn("cross block", text)
        self.assertIn("before scalarization", text)
        self.assertIn("transient attitude-sector excursions", text)


if __name__ == "__main__":
    unittest.main()
