import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_prefix as G


class Ou3P5Sample1PrefixTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build(source_pieces=2)

    def test_stage_is_executable_and_fail_closed(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertIn(d["P5_SAMPLE1_S_ACCEL_PREFIX_CERTIFICATE"], ("PASS", "NOT_ESTABLISHED"))
        if d["P5_SAMPLE1_S_ACCEL_PREFIX_CERTIFICATE"] == "NOT_ESTABLISHED":
            self.assertIsNotNone(d["first_failure"])
        else:
            self.assertIsNone(d["first_failure"])

    def test_estimator_mean_and_timer_semantics_are_explicit(self):
        d = self.d
        self.assertTrue(d["estimator_mean_tracked_separately_from_physical_error"])
        self.assertTrue(d["sample1_S_residual_uses_estimator_mean_not_physical_error"])
        self.assertTrue(d["sample1_pseudo_phase_conditioned_on_sample0_timer_branch"])
        self.assertTrue(d["periodic_update_due_float_tolerance_retained"])
        counts = d["sample1_phase_counts"]
        self.assertGreaterEqual(counts["due"], 0)
        self.assertGreaterEqual(counts["not_due"], 0)

    def test_accelerometer_diagnostic_keeps_exact_finite_rotation_norm(self):
        d = self.d
        self.assertTrue(d["sample1_accel_exact_rotation_difference_norm_used"])
        self.assertTrue(d["sample1_accel_generic_post_reset_orientation_box_used"])
        self.assertTrue(math.isfinite(float(d["sample1_entry_cayley_norm_upper"])))
        self.assertLess(float(d["sample1_entry_cayley_norm_upper"]), 8.0)

    def test_no_gate_relaxation_or_word_promotion(self):
        d = self.d
        self.assertEqual(d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(d["deployed_correction_limit_increased"])
        self.assertFalse(d["sample1_magnetometer_evaluated_here"])
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])


if __name__ == "__main__":
    unittest.main()
