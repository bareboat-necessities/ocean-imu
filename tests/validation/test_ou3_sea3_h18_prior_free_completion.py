#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_sea3_h18_prior_free_completion as H18


class Sea3H18PriorFreeCompletionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = H18.build()
        cls.failures = H18.validate(cls.d)

    def test_complete_sea3_source_and_same_word_endpoint(self):
        d = self.d
        self.assertEqual([], self.failures)
        self.assertEqual("COMPLETE_SEA3_NORMAL_LIVE_WORD", d["canonical_source"])
        self.assertEqual(3.0, d["complete_word_horizon_s"])
        pbar = d["same_word_diffuse_prior_covariance_upper"]
        self.assertTrue(pbar["concrete_estimator_uses_only_measurements_inside_same_word"])
        self.assertTrue(pbar["diffuse_prior_endpoint_covariance_Loewner_upper"])
        self.assertTrue(pbar["TD_inverse_T_transpose_Loewner_upper"])
        self.assertTrue(pbar["marginal_bounds_not_misused_as_diagonal_Loewner_matrix"])
        self.assertEqual(3.0, pbar["timing"]["same_word_endpoint_s"])
        self.assertTrue(pbar["timing"]["selected_S_windows_fit_same_word"])

    def test_full_18x18_prior_free_ldlt_closes_at_useful_gate(self):
        d = self.d
        self.assertEqual(1.0e-18, d["useful_gate"])
        self.assertTrue(d["H18_D_is_strictly_positive_definite"])
        self.assertTrue(d["prior_free_exact_identity_consumed"])
        self.assertTrue(d["full_18x18_interval_LDLT_used"])
        self.assertTrue(d["full_H18_prior_free_matrix_condition_closed"])
        self.assertTrue(d["H18_prior_free_completion_closed"])
        self.assertGreater(d["x_cells_certified"], 0)
        self.assertEqual([], d["x_cell_failures"])
        pivot = float(d["worst_full_H18_LDLT_pivot_lower"])
        self.assertTrue(math.isfinite(pivot))
        self.assertGreater(pivot, 0.0)

    def test_actual_rs_and_stable_shipping_q_are_retained(self):
        d = self.d
        self.assertTrue(d["actual_applied_SpectralMSE_R_S_consumed_in_same_word_covariance_upper"])
        self.assertTrue(d["all_due_S_updates_remain_in_complete_word"])
        self.assertTrue(d["stable_factored_shipping_integrated_OU_Q_consumed"])
        self.assertTrue(d["event_algebra_preserves_margin_after_closure"])

    def test_no_retired_contraction_shortcuts(self):
        d = self.d
        for key in (
            "old_one_step_Euclidean_Q_min_used",
            "marginal_Pbar_bounds_used_as_Loewner_diagonal_directly",
            "blockwise_minimum_contraction_used",
            "D_W_L_W_product_used",
            "scalar_beta_contraction_used",
            "determinant_trace_final_matrix_gate_used",
            "source_family_replaced",
            "trajectory_replay_used",
            "independent_tau_sigma_RS_source_created",
            "P3_promoted",
        ):
            self.assertFalse(d[key], key)


if __name__ == "__main__":
    unittest.main()
