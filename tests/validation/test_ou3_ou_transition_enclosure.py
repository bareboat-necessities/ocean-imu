import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_ou_transition_enclosure as mod


class OUTransitionEnclosureTests(unittest.TestCase):
    def test_source_derived_scalar_transition_is_validated_but_not_promoted(self):
        d = mod.build(mod.SOURCE.DEFAULT_HEADER)
        self.assertEqual(d["qualification"], mod.QUALIFICATION)
        self.assertTrue(d["validated_arithmetic"])
        self.assertTrue(d["outward_rounded"])
        self.assertTrue(d["continuous_ou_scalar_transition_enclosed"])
        self.assertFalse(d["continuous_matrix_word_enclosed"])
        self.assertFalse(d["nonlinear_word_enclosed"])
        self.assertEqual(d["theorem_promotion"], "NOT_ESTABLISHED")
        self.assertFalse(d["runtime_timing_contract"]["arbitrary_positive_api_dt_covered"])
        self.assertLessEqual(d["x_h_over_tau"][1], 1.0)

    def test_coefficients_cover_source_tau_endpoints(self):
        d = mod.build(mod.SOURCE.DEFAULT_HEADER)
        h = 0.005
        for tau in (0.02, 12.0):
            x = h / tau
            values = {
                "alpha": math.exp(-x),
                "expm1_neg_x": math.expm1(-x),
                "phi_pa": tau * tau * (x + math.expm1(-x)),
                "phi_Sa": tau ** 3 * (0.5 * x * x - x - math.expm1(-x)),
            }
            for name, value in values.items():
                lo, hi = d["coefficients"][name]
                self.assertLessEqual(lo, value, name)
                self.assertGreaterEqual(hi, value, name)

    def test_validation_rejects_claim_upgrade(self):
        d = mod.build(mod.SOURCE.DEFAULT_HEADER)
        d["theorem_promotion"] = "PASS"
        self.assertTrue(mod.validate(d, mod.SOURCE.DEFAULT_HEADER))


if __name__ == "__main__":
    unittest.main()
