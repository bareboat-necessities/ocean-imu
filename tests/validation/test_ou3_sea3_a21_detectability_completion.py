#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_sea3_a21_detectability_completion as A21


class Sea3A21DetectabilityCompletionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = A21.build()
        cls.failures = A21.validate(cls.d)

    def test_complete_source_and_paper_finite_bias_route(self):
        d = self.d
        self.assertEqual([], self.failures)
        self.assertEqual("COMPLETE_SEA3_NORMAL_LIVE_WORD", d["canonical_source"])
        self.assertEqual("ETA6_PLUS_FINITE_RESIDUAL_BIAS_CORRELATION", d["paper_active_bias_route"])
        self.assertFalse(d["eta9_point_packet_shortcut_used"])
        self.assertTrue(d["H18_complete_word_contraction_consumed"])

    def test_detectability_gap_clears_useful_gate(self):
        d = self.d
        self.assertEqual(1.0e-18, d["A21_detectability_useful_gate"])
        self.assertGreater(d["bias_homogeneous_energy_gap_lower"], 1.0e-18)
        self.assertGreaterEqual(d["A21_detectability_asymptotic_word_energy_gap_lower"], 1.0e-18)
        self.assertTrue(d["A21_detectability_useful_gate_pass"])
        self.assertTrue(d["A21_finite_bias_detectability_closed"])
        self.assertTrue(d["A21_paper_UES_hypotheses_closed"])

    def test_shipping_release_and_process_are_bound(self):
        d = self.d
        self.assertTrue(d["H_to_A_release_is_bounded_one_time_mode_jump"])
        self.assertFalse(d["H_to_A_release_requires_preceding_three_second_H_interval"])
        p = d["active_bias_process"]
        self.assertGreater(p["tau_ba_s"], 0.0)
        self.assertGreater(p["continuous_driving_variance_density"], 0.0)
        self.assertGreater(p["one_sample_Q_ba_lambda_min_lower"], 0.0)
        self.assertGreater(p["uniform_release_and_active_variance_upper"], 0.0)
        self.assertTrue(p["full_A21_process_UCC_pass"])

    def test_comparison_observer_is_proof_only_and_coupling_is_finite(self):
        d = self.d
        tri = d["triangular_detectability_observer"]
        self.assertTrue(tri["H_diagonal_block_uses_complete_SEA3_H18_certificate"])
        self.assertTrue(tri["bias_diagonal_block_uses_shipping_GM_decay"])
        self.assertTrue(tri["upper_right_coupling_finite_on_compact_word"])
        self.assertTrue(tri["finite_coupling_changes_prefactor_not_asymptotic_rate"])
        self.assertTrue(tri["comparison_observer_only_not_alternate_estimator"])
        self.assertFalse(tri["shipping_filter_changed"])
        self.assertTrue(math.isfinite(float(tri["finite_coupling_log10_upper"])))

    def test_stronger_canonical_A21_matrix_gate_remains_fail_closed(self):
        d = self.d
        self.assertFalse(d["full_21x21_Omega_minus_delta_P_LDLT_closed_here"])
        self.assertFalse(d["P3_CANONICAL_PASS"])
        self.assertFalse(d["P4_MAY_CONSUME_P3"])
        self.assertFalse(d["source_family_replaced"])
        self.assertFalse(d["trajectory_replay_used"])
        self.assertFalse(d["independent_tau_sigma_RS_source_created"])


if __name__ == "__main__":
    unittest.main()
