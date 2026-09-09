#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_p4_correlated_entry_relation as REL


class CorrelatedEntryRelationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = REL.build()

    def test_validates_and_promotes_nothing(self):
        self.assertEqual(REL.validate(self.d), [])
        self.assertFalse(self.d["P4_promoted_here"])
        self.assertFalse(self.d["entry_relation_consumed_by_final_gate"])
        self.assertFalse(self.d["declared_domain_shrunk"])
        self.assertFalse(self.d["entry_radii_reduced_here"])
        self.assertFalse(self.d["trajectory_replay_used"])
        self.assertFalse(self.d["covariance_confidence_ellipsoid_used"])
        self.assertFalse(self.d["P5_capture_assumed"])

    def test_deployed_S_row_matches_shipping_factor(self):
        h = float(self.d["sample_period_s"])
        row = self.d["deployed_S_row_coefficients"]
        # Kalman3D_Wave_OU_III IntegratedOUChain<T,3>::transition sets the S row
        # to (0.5*h*h, h, 1, phi_Sa).
        self.assertLessEqual(0.5 * h * h, row["v"])
        self.assertEqual(row["p"], h)
        self.assertEqual(row["S"], 1.0)
        self.assertLessEqual(h * h * h / 6.0, row["a_w_upper"])
        self.assertTrue(self.d["error_obeys_same_row_exactly"])

    def test_declared_ball_is_neither_derived_nor_conservative(self):
        self.assertEqual(self.d["declared_independent_integral_radius_m_s"], 300.0)
        self.assertFalse(self.d["declared_independent_radius_has_operating_domain_derivation"])
        self.assertFalse(self.d["declared_independent_radius_is_conservative"])
        self.assertGreater(self.d["free_integral_reach_over_declared_live_entry_floor_m_s"],
                           self.d["declared_independent_integral_radius_m_s"])

    def test_accumulation_is_monotone_in_the_dwell_window(self):
        gap = self.d["dwell_accumulation"]["one_S_zero_cadence_gap"]
        word = self.d["dwell_accumulation"]["one_P4_word"]
        self.assertLess(gap["window_s"], word["window_s"])
        self.assertLess(gap["accumulation_upper_m_s"], word["accumulation_upper_m_s"])
        # Position dominates; velocity and latent acceleration are corrections.
        self.assertGreater(gap["from_position_error_m_s"], gap["from_velocity_error_m_s"])
        self.assertGreater(gap["from_velocity_error_m_s"], gap["from_latent_acceleration_error_m_s"])

    def test_scheduler_gap_bounds_the_correlated_radius(self):
        sched = self.d["S_zero_scheduler"]
        self.assertLessEqual(sched["pseudo_update_period_upper_s"], 0.15 + 1e-7)
        self.assertGreater(sched["consecutive_event_gap_upper_s"], sched["pseudo_update_period_upper_s"])
        self.assertLess(self.d["correlated_radius_one_cadence_m_s"], 300.0)

    def test_one_cadence_relation_retains_the_declared_chart(self):
        self.assertTrue(self.d["one_cadence_relation_meets_chart_threshold"])
        self.assertFalse(self.d["declared_radius_meets_chart_threshold"])
        self.assertFalse(self.d["one_word_relation_meets_chart_threshold"])
        self.assertEqual(self.d["chart_retention_binding_mode"], "H18")
        self.assertGreater(self.d["chart_retaining_integral_radius_upper_m_s"], 0.0)
        self.assertTrue(self.d["chart_threshold_is_frozen_map_diagnostic_not_outward_certificate"])

    def test_correction_ceiling_reference_is_reporting_only(self):
        c = self.d["correction_ceiling_threshold"]
        self.assertTrue(c["reference_is_reporting_only_not_a_proof_premise"])
        self.assertTrue(self.d["one_cadence_relation_meets_correction_ceiling"])
        self.assertFalse(self.d["declared_radius_meets_correction_ceiling"])

    def test_S_zero_regulation_is_non_expansive_but_not_a_uniform_anchor(self):
        reg = self.d["S_zero_regulation"]
        self.assertTrue(reg["weighted_non_expansion_always_holds"])
        self.assertFalse(reg["uniform_contraction_useful"])
        self.assertFalse(self.d["S_zero_gives_uniform_anchor_without_dwell"])
        self.assertGreater(reg["worst_cell_mu_lower"], 0.0)
        self.assertLess(reg["worst_cell_mu_lower"], 1e-6)
        self.assertLessEqual(reg["worst_cell_weighted_contraction_factor_upper"], 1.0)
        self.assertGreater(reg["events_per_e_fold_lower"], 1e6)
        self.assertTrue(self.d["unconditional_entry_anchor_requires_P5_or_reachable_P_SS_lower_bound"])

    def test_dwell_window_requirement_is_met_by_one_cadence(self):
        self.assertTrue(self.d["one_cadence_gap_meets_required_dwell"])
        self.assertGreater(self.d["dwell_window_required_for_chart_threshold_s"],
                           self.d["S_zero_scheduler"]["consecutive_event_gap_upper_s"])

    def test_accumulation_helper_is_outward_and_step_aligned(self):
        out = REL._accumulation_upper(0.02, 0.005, 10.0, 1.0, 1.0)
        self.assertGreaterEqual(out["steps"], 4)
        self.assertGreaterEqual(out["covered_span_s"], 0.02)
        self.assertGreaterEqual(out["accumulation_upper_m_s"], 0.02 * 10.0)
        with self.assertRaises(ValueError):
            REL._accumulation_upper(0.0, 0.005, 1.0, 1.0, 1.0)

    def test_contraction_helper_rejects_degenerate_inputs(self):
        with self.assertRaises(ValueError):
            REL._s_zero_contraction(0.0, 1.0, 0.005, 1.0)
        good = REL._s_zero_contraction(0.05, 12.0, 0.005, 0.108)
        self.assertTrue(math.isfinite(good["process_floor_q_SS_per_step"]))
        self.assertGreater(good["process_floor_q_SS_per_step"], 0.0)

    def test_falsification_class_names_entry_and_qualification(self):
        self.assertIn("D_entry_set_modeling", self.d["falsification_class"])
        self.assertIn("E_missing_source_qualification", self.d["falsification_class"])


if __name__ == "__main__":
    unittest.main()
