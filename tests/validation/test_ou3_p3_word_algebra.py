import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p3_word_algebra as ALG


class Ou3P3WordAlgebraTests(unittest.TestCase):
    def test_source_bound_word_algebra_closes_every_live_covariance_class(self):
        d = ALG.build()
        self.assertEqual(ALG.validate(d), [])
        self.assertTrue(d["pass"])
        self.assertEqual(d["fixed_dimensions"], {"H": 18, "A": 21})
        self.assertFalse(d["dimension_change_inside_word"])
        self.assertEqual(
            set(d["operation_classes"]),
            {"prediction", "accepted_joseph", "rejected_or_not_due", "left_error_reset", "aw_covariance_sync"},
        )
        self.assertTrue(d["covariance_decomposition_invariant"]["Omega_s_psd"])
        self.assertTrue(d["strict_margin_preservation"]["covers_every_operation_class"])
        self.assertFalse(d["strict_margin_preservation"]["one_step_contraction_assumed"])

    def test_prefix_information_gain_is_derived_as_one(self):
        d = ALG.build()
        p = d["prefix_information_bound"]
        self.assertEqual(p["information_gain_upper"], 1.0)
        self.assertTrue(p["source_uniform"])
        self.assertFalse(p["sampled_evidence_used"])
        self.assertIn("Omega_s>=0", p["reason"])

    def test_left_error_reset_is_globally_nonsingular_for_finite_injection(self):
        d = ALG.build()
        self.assertEqual(d["reset"]["determinant_formula"], "1+||dtheta||^2/4")
        self.assertEqual(d["reset"]["determinant_lower"], 1.0)
        self.assertFalse(d["reset"]["requires_small_angle_for_nonsingularity"])
        for theta in (0.0, 1e-9, 0.2, 1.0, math.pi, 100.0):
            self.assertGreaterEqual(ALG.reset_det(theta), 1.0)

    def test_affine_psd_suffix_cannot_destroy_a_strict_margin(self):
        delta = 0.2
        self.assertEqual(ALG.margin_after_affine_psd(delta, 0.0, 0.0), 0.0)
        self.assertGreater(ALG.margin_after_affine_psd(delta, 0.01, 0.0), 0.0)
        self.assertGreater(ALG.margin_after_affine_psd(delta, 0.0, 0.01), 0.0)
        self.assertGreater(
            ALG.margin_after_affine_psd(delta, 0.01, 0.02),
            0.01,
        )

    def test_psd_aw_sync_is_nonworsening_for_p3_delta(self):
        d = ALG.build()
        aw = d["operation_classes"]["aw_covariance_sync"]
        self.assertTrue(aw["B_psd"])
        self.assertFalse(aw["can_reduce_existing_delta"])
        self.assertIn("(1-delta)", aw["margin_increment"])

    def test_rejected_and_not_due_branches_are_identity_covariance_maps(self):
        d = ALG.build()
        r = d["operation_classes"]["rejected_or_not_due"]
        self.assertEqual(r["A"], "I")
        self.assertEqual(r["B"], "0")
        self.assertIn("accelerometer rejected", r["covers"])
        self.assertIn("magnetometer rejected", r["covers"])
        self.assertIn("S pseudo not due", r["covers"])


if __name__ == "__main__":
    unittest.main()
