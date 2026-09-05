#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_sea3_a21_accumulated_word_matrix as CERT


class Sea3A21AccumulatedWordMatrixTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = CERT.build()
        cls.failures = CERT.validate(cls.d)

    def test_complete_source_and_reachable_covariance_upper(self):
        d = self.d
        self.assertEqual([], self.failures)
        self.assertEqual("COMPLETE_SEA3_NORMAL_LIVE_WORD", d["canonical_source"])
        self.assertTrue(d["same_complete_SEA3_word_used"])
        self.assertEqual("P<=21*diag(u)", d["full_covariance_Loewner_upper_rule"])
        self.assertTrue(d["unknown_covariance_cross_terms_covered"])
        self.assertFalse(d["marginal_bounds_misused_as_diagonal_Loewner_upper"])

    def test_full_21x21_restartable_matrix_closes(self):
        d = self.d
        self.assertEqual(1.0e-18, d["useful_gate"])
        self.assertTrue(d["full_21x21_interval_LDLT_used"])
        self.assertTrue(d["full_21x21_generalized_Riccati_matrix_closed"])
        self.assertTrue(d["A21_restartable_reachable_tube_contraction_closed"])
        self.assertGreater(d["translation_word_X_cells_certified"], 0)
        self.assertEqual([], d["translation_word_X_failures"])
        p = float(d["worst_full_A21_inverse_space_LDLT_pivot_lower"])
        self.assertTrue(math.isfinite(p))
        self.assertGreater(p, 0.0)

    def test_measurement_cross_information_is_paid(self):
        m = self.d["measurement_information_upper"]
        self.assertGreaterEqual(m["accelerometer_partition_factor"], 3.0)
        self.assertTrue(m["accelerometer_theta_aw_ba_cross_terms_paid"])
        self.assertTrue(m["magnetometer_only_theta"])
        self.assertTrue(m["S_zero_only_S"])

    def test_no_rejected_promotion_route(self):
        d = self.d
        for key in (
            "source_family_replaced",
            "trajectory_replay_used",
            "independent_tau_sigma_RS_TS_product_used_as_source",
            "eta9_point_packet_shortcut_used",
            "blockwise_minimum_ratio_used",
            "scalar_beta_contraction_used",
            "D_W_L_W_product_used",
            "determinant_trace_final_gate_used",
            "P3_promoted",
        ):
            self.assertFalse(d[key], key)
        self.assertTrue(d["actual_applied_SpectralMSE_R_S_source_contract_retained"])
        self.assertTrue(d["all_due_S_updates_remain_in_complete_word"])
        self.assertTrue(d["event_algebra_preserves_margin_after_word_closure"])


if __name__ == "__main__":
    unittest.main()
