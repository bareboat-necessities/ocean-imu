#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p4_reset_directional_transport as R


class ResetDirectionalTransportTests(unittest.TestCase):
    def test_source_bound_reset_congruence_preserves_directional_rank(self):
        d = R.build()
        self.assertEqual(R.validate(d), [])
        self.assertTrue(d["P4_RESET_DIRECTIONAL_TRANSPORT_ESTABLISHED"])
        self.assertTrue(d["reset_invertible_for_every_finite_correction"])
        self.assertTrue(d["homogeneous_information_congruence_exact"])
        self.assertTrue(d["directional_form_congruence_exact"])
        self.assertEqual(d["reset_inverse_operator_norm_exact"], 1.0)
        self.assertFalse(d["condition_number_conversion_used"])
        self.assertFalse(d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"])
        self.assertFalse(d["P5_FINITE_CAPTURE_ESTABLISHED"])
        for mode, n in (("H", 18), ("A", 21)):
            m = d["modes"][mode]
            self.assertEqual(m["dimension"], n)
            self.assertEqual(m["homogeneous_information_quadratic_gain_exact"], 1.0)
            self.assertTrue(m["directional_rank_preserved_exactly"])
            self.assertTrue(m["directional_nullity_preserved_exactly"])
            self.assertFalse(m["condition_number_multiplier_required"])
            self.assertFalse(m["nonlinear_Cayley_injection_defect_closed_here"])
            self.assertFalse(m["P4_PROMOTED"])

    def test_reset_spectral_facts_have_exact_inverse_norm_one(self):
        for delta in (0.0, 1.0e-3, 0.1, 0.8, 3.0):
            s = R.reset_spectral_facts(delta)
            self.assertEqual(s["reset_min_singular_value_lower"], 1.0)
            self.assertEqual(s["reset_inverse_operator_norm_upper"], 1.0)
            self.assertGreaterEqual(s["reset_max_singular_value_upper"], 1.0)
            self.assertGreaterEqual(s["reset_determinant_lower"], 1.0)
            self.assertGreaterEqual(s["reset_determinant_upper"], 1.0)
            self.assertTrue(math.isfinite(s["reset_max_singular_value_upper"]))

    def test_validation_rejects_condition_number_or_promotion(self):
        d = R.build()
        d["condition_number_conversion_used"] = True
        d["P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED"] = True
        f = R.validate(d)
        self.assertTrue(any("condition_number_conversion_used" in x for x in f))
        self.assertTrue(any("P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED" in x for x in f))


if __name__ == "__main__":
    unittest.main()
