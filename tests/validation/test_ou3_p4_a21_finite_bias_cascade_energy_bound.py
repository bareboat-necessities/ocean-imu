from __future__ import annotations

import math
import unittest

import ou3_p4_a21_finite_bias_cascade_energy_bound as ENERGY


class TestA21FiniteBiasCascadeEnergyBound(unittest.TestCase):
    def test_geometric_decay_energy_is_strictly_below_packet_count(self):
        n = 600
        s = ENERGY.geometric_decay_energy_upper(n, 0.005, 5000.0)
        self.assertGreater(s, 0.0)
        self.assertLess(s, float(n))
        # Slow 5000 s GM decay over 3 s is deliberately close to N, but still strict.
        self.assertGreater(s, 599.0)

    def test_bad_geometric_inputs_fail_closed(self):
        with self.assertRaises(ValueError):
            ENERGY.geometric_decay_energy_upper(0, 0.005, 5000.0)
        with self.assertRaises(ValueError):
            ENERGY.geometric_decay_energy_upper(600, 0.0, 5000.0)
        with self.assertRaises(ValueError):
            ENERGY.geometric_decay_energy_upper(600, 0.005, 0.0)

    def test_sufficient_power10_has_explicit_decimal_slack(self):
        p = ENERGY.sufficient_power10(26.01)
        self.assertEqual(p, 28)
        self.assertGreater(float(p), 26.01)

    def test_repository_comparison_cascade_closes_but_shipping_bridge_stays_open(self):
        d = ENERGY.build()
        self.assertEqual(ENERGY.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertEqual(d["P3_frozen_gate"], 1.0e-18)
        self.assertTrue(d["P3_frozen_not_modified"])
        self.assertEqual(d["required_accelerometer_updates"], 600)
        self.assertTrue(d["all_valid_accelerometer_updates_retained"])
        self.assertTrue(d["every_due_S_update_with_actual_RS_retained_inside_H18_word"])
        self.assertEqual(d["actual_RS_axis_factors"], [0.72, 0.72, 1.0])
        self.assertTrue(d["comparison_cascade_Schur_condition_closed"])
        self.assertTrue(d["comparison_observer_full_rank_21_state_metric_closed"])
        self.assertFalse(d["actual_A21_shipping_bias_correction_row_consumed"])
        self.assertFalse(d["actual_A21_shipping_H_b_cross_terms_closed"])
        self.assertFalse(d["actual_A21_shipping_full_cross_term_metric_closed"])
        self.assertFalse(d["P4_promoted_here"])
        self.assertFalse(d["P5_may_start"])

    def test_bound_does_not_reintroduce_retired_packet_or_suffix_routes(self):
        d = ENERGY.build()
        self.assertFalse(d["suffix_norm_product_used"])
        self.assertFalse(d["N_times_worst_map_coupling_used"])
        self.assertFalse(d["packet_count_nonlinear_remainder_multiplier_used"])
        self.assertTrue(d["closed_form_same_word_bias_energy_used"])
        self.assertFalse(d["state_elimination_used"])
        self.assertFalse(d["a_w_elimination_used"])
        self.assertGreater(d["joint_C_Hb_metric_c2_upper"], 0.0)
        self.assertTrue(math.isfinite(d["mu_log10_required_upper"]))
        self.assertGreater(d["mu_sufficient_power10"], d["mu_log10_required_upper"])
        self.assertLess(d["c2_over_mu_to_g2_budget_ratio_upper"], 1.0)


if __name__ == "__main__":
    unittest.main()
