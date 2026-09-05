#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_sea3_a21_prior_free_completion as A21


class Sea3A21PriorFreeCompletionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = A21.build()
        cls.failures = A21.validate(cls.d)

    def test_complete_sea3_and_same_word_A21_covariance_upper(self):
        d = self.d
        self.assertEqual([], self.failures)
        self.assertEqual("COMPLETE_SEA3_NORMAL_LIVE_WORD", d["canonical_source"])
        self.assertEqual(3.0, d["complete_word_horizon_s"])
        pbar = d["same_word_diffuse_prior_covariance_upper"]
        self.assertEqual(21, pbar["dimension"])
        self.assertEqual(21, len(pbar["Pbar_diagonal_variance_upper"]))
        self.assertTrue(pbar["active_ba_nuisance_paid_in_H_subblock_upper"])
        self.assertTrue(pbar["diffuse_prior_endpoint_covariance_Loewner_upper"])
        self.assertTrue(pbar["TD_inverse_T_transpose_Loewner_upper"])
        self.assertTrue(pbar["marginal_bounds_not_misused_as_diagonal_Loewner_matrix"])
        self.assertEqual(3.0, pbar["timing"]["same_word_endpoint_s"])
        self.assertTrue(pbar["timing"]["selected_S_windows_fit_same_word"])

    def test_full_21x21_prior_free_ldlt_closes_at_useful_gate(self):
        d = self.d
        self.assertEqual(1.0e-18, d["useful_gate"])
        self.assertTrue(d["prior_free_exact_identity_consumed"])
        self.assertTrue(d["full_21x21_interval_LDLT_used"])
        self.assertTrue(d["full_21x21_Omega_minus_delta_P_LDLT_closed"])
        self.assertTrue(d["A21_prior_free_completion_closed"])
        self.assertGreater(d["x_cells_certified"], 0)
        self.assertEqual([], d["x_cell_failures"])
        pivot = float(d["worst_full_A21_LDLT_pivot_lower"])
        self.assertTrue(math.isfinite(pivot))
        self.assertGreater(pivot, 0.0)

    def test_shipping_bias_route_and_process_are_retained(self):
        d = self.d
        self.assertEqual("ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION", d["paper_active_bias_route"])
        self.assertFalse(d["eta9_point_packet_shortcut_used"])
        self.assertTrue(d["A21_finite_bias_detectability_consumed"])
        self.assertTrue(d["shipping_active_ba_GM_Q_consumed"])
        self.assertTrue(d["stable_factored_shipping_integrated_OU_Q_consumed"])
        self.assertTrue(d["event_algebra_preserves_margin_after_closure"])

    def test_actual_rs_and_no_retired_shortcuts(self):
        d = self.d
        self.assertTrue(d["actual_applied_SpectralMSE_R_S_consumed_in_same_word_covariance_upper"])
        self.assertTrue(d["all_due_S_updates_remain_in_complete_word"])
        for key in (
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
