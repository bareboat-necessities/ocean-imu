from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_nonlinear_word_certificate as P4


class Ou3P4NonlinearWordCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = P4.build()

    def test_P4_exact_H_A_words_pass_with_positive_mu(self):
        d = self.d
        self.assertEqual(P4.validate(d), [])
        self.assertEqual(d["P4_EXACT_NONLINEAR_WORD_CERTIFICATE"], "PASS")
        self.assertEqual(d["metric_route"], "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC")
        self.assertFalse(d["block_diagonal_metric_fallback"])
        self.assertFalse(d["source_replay_used"])
        for mode in ("H", "A"):
            m = d["modes"][mode]
            self.assertTrue(m["exact_nonlinear_word_pass"])
            self.assertGreater(m["certified_level_W"], 0.0)
            self.assertGreater(m["endpoint_relative_W_decrease_lower"], 0.0)
            self.assertGreater(m["mu_W_lower"], 0.0)
            self.assertGreater(m["theta_star"], 0.0)
            self.assertLess(m["theta_star"], math.pi)
            self.assertTrue(m["all_word_prefixes_safe"])
            self.assertLess(m["prefix_canonical_error_norm_upper"], m["cayley_norm_limit"])
            self.assertLess(m["accepted_correction_norm_prefix_upper"], 1.0e-2)
            self.assertTrue(m["accepted_correction_uses_source_series_branch"])

    def test_complete_source_branch_families_are_covered_without_enumeration(self):
        c = self.d["word_branch_coverage"]
        self.assertEqual(
            c["method"],
            "uniform exact-operation quadratic defect envelope plus P3 unit segment information transport",
        )
        self.assertFalse(c["explicit_exponential_branch_enumeration_required"])
        self.assertTrue(c["all_admissible_accelerometer_accept_reject_branches"])
        self.assertTrue(c["all_admissible_magnetometer_not_due_accept_reject_branches"])
        self.assertTrue(c["all_admissible_S_not_due_due_branches"])
        self.assertTrue(c["all_admissible_aw_sync_not_due_due_branches"])

    def test_endpoint_budget_is_strict_and_not_one_sample_contraction(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            self.assertFalse(m["one_sample_decrease_used"])
            self.assertEqual(m["P3_homogeneous_prefix_information_gain_upper"], 1.0)
            self.assertLessEqual(m["nonlinear_sqrt_budget_fraction_of_delta_upper"], 0.125000000000001)
            self.assertGreaterEqual(
                m["P3_word_endpoint_delta_lower"],
                2.0 * m["endpoint_relative_W_decrease_lower"],
            )

    def test_A_inner_funnel_does_not_reach_bias_projection_surface(self):
        p = self.d["modes"]["A"]["active_bias_projection"]
        self.assertFalse(p["projection_surface_reached_in_certified_funnel"])
        self.assertEqual(p["exact_projection_branch_in_certified_funnel"], "identity_interior_branch")
        self.assertLess(p["certified_error_norm_prefix_upper"], p["interior_margin_lower_mps2"])

    def test_no_replay_or_old_block_metric_is_used_by_producer(self):
        text = (ROOT / "tools" / "ou3_p4_nonlinear_word_certificate.py").read_text(encoding="utf-8").casefold()
        self.assertNotIn("numpy", text)
        self.assertNotIn("monte carlo", text)
        self.assertNotIn("sampled trajectory", text)
        self.assertIn("block_diagonal_metric_fallback", text)
        self.assertIn("cayley", text)


if __name__ == "__main__":
    unittest.main()
