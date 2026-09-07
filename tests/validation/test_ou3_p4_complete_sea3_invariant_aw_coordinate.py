from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STABILITY = ROOT / "tools" / "stability"
sys.path.insert(0, str(STABILITY))

import ou3_p4_complete_sea3_invariant_aw_coordinate as INV


class CompleteSea3InvariantAwCoordinateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = INV.build()

    def test_certificate_validates(self):
        self.assertEqual([], INV.validate(self.d))

    def test_source_and_p3_remain_unchanged(self):
        self.assertEqual("COMPLETE_SEA3_NORMAL_LIVE_WORD", self.d["canonical_source"])
        self.assertTrue(self.d["P3_frozen_not_modified"])
        self.assertTrue(self.d["complete_SEA3_word_retained"])
        self.assertTrue(self.d["all_valid_accelerometer_updates_remain_in_complete_word"])
        self.assertTrue(self.d["all_due_S_updates_and_actual_RS_remain_in_complete_word"])
        self.assertFalse(self.d["trajectory_replay_used"])
        self.assertFalse(self.d["source_family_replaced"])
        self.assertFalse(self.d["declared_domain_changed"])

    def test_triangular_transform_is_exact_not_metric_shortcut(self):
        tri = self.d["linear_triangular_coordinate"]
        self.assertTrue(tri["T_B_unit_triangular"])
        self.assertTrue(tri["T_B_nonsingular"])
        self.assertTrue(tri["T_B_inverse_exact"])
        self.assertEqual(1.0, tri["T_B_determinant_exact"])
        self.assertGreater(tri["B_operator_norm_upper_mps2_per_cayley"], 0.0)
        metric = self.d["moving_metric_congruence"]
        self.assertTrue(metric["innovation_covariance_S_invariant"])
        self.assertTrue(metric["Joseph_covariance_congruence_exact"])
        self.assertTrue(metric["moving_Riccati_energy_invariant"])
        self.assertFalse(metric["condition_number_multiplier_used"])
        self.assertFalse(metric["group_isotropic_metric_assumption_used"])

    def test_wave_attitude_linear_term_moves_into_aw_coordinate(self):
        self.assertTrue(self.d["transformed_attitude_column_depends_only_on_gravity"])
        self.assertTrue(self.d["wave_acceleration_attitude_cross_term_is_linear_coordinate_coupling"])
        exact = self.d["exact_finite_angle_coordinate"]
        self.assertTrue(exact["nonlinear_displacement_is_full_aw_shift"])
        self.assertTrue(exact["first_order_wave_attitude_term_removed_from_remainder"])
        self.assertIn("R_hat*g", exact["exact_residual"])

    def test_actual_rs_is_retained_and_no_packet_budget_returns(self):
        self.assertTrue(self.d["actual_RS_information_matrix_retained_under_congruence"])
        self.assertFalse(self.d["standalone_eta_Rinv_packet_budget_used"])
        self.assertFalse(self.d["packet_count_multiplier_used"])

    def test_helper_does_not_promote_p4(self):
        self.assertFalse(self.d["complete_source_correlated_transport_defect_closed_here"])
        self.assertFalse(self.d["P4_promoted_here"])
        self.assertFalse(self.d["shipping_Joseph_binding_closed"])

    def test_nominal_force_is_derived_from_existing_sea3_and_error_domain(self):
        self.assertTrue(
            self.d["nominal_force_bound_from_complete_SEA3_plus_declared_error_domain_proved"]
        )
        self.assertFalse(self.d["nominal_force_bound_from_physical_SEA3_alone_claimed"])
        f = self.d["nominal_force_derivation"]
        self.assertTrue(f["uses_complete_SEA3_physical_source_bound"])
        self.assertTrue(f["uses_already_declared_P4_error_domain"])
        self.assertFalse(f["new_source_assumption_added"])
        self.assertFalse(f["replay_extrema_used"])
        self.assertAlmostEqual(13.80665, f["true_specific_force_upper_mps2"], places=6)
        self.assertAlmostEqual(2.941995, f["declared_delta_aw_error_upper_mps2"], places=6)
        self.assertGreaterEqual(f["nominal_specific_force_upper_mps2"], 16.748645)
        self.assertGreaterEqual(f["nominal_aw_norm_upper_mps2"], 26.555295)
        self.assertLess(f["nominal_aw_norm_upper_mps2"], 27.0)

    def test_candidate_eta_bounds_are_widened_before_reuse(self):
        self.assertFalse(self.d["measurement_linearizing_shift_bounds_reused_without_widening"])
        rows = self.d["measurement_linearizing_shift_bounds_widened_for_declared_delta_aw_domain"]
        self.assertEqual([30.0, 25.0, 20.0, 15.0],
                         [r["attitude_angle_deg"] for r in rows])
        nominal_force = float(self.d["source_nominal_specific_force_upper_mps2"])
        for row in rows:
            self.assertTrue(row["widened_for_declared_delta_aw_domain"])
            self.assertGreaterEqual(float(row["force_upper_mps2_used"]), nominal_force)
            self.assertGreater(float(row["nominal_force_widening_ratio"]), 1.0)
            self.assertGreater(float(row["e_eta_norm_upper_mps2"]), 0.0)
            self.assertTrue(math.isfinite(float(row["e_eta_norm_upper_mps2"])))


if __name__ == "__main__":
    unittest.main()
