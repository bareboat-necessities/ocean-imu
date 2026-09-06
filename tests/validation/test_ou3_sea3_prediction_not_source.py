from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_sea3_shipping_prediction_primitives as primitive  # noqa: E402


class Sea3PredictionNotSourceTest(unittest.TestCase):
    def test_prediction_math_has_no_promotion_or_source_domain(self):
        d = primitive.build()
        self.assertEqual(primitive.validate(d), [])
        self.assertFalse(d["source_domain_created_here"])
        self.assertFalse(d["arbitrary_bounded_input_source_created_here"])
        self.assertFalse(d["independent_tau_sigma_source_created_here"])
        self.assertTrue(d["consumes_only_SEA3_derived_sample_coordinates"])
        self.assertFalse(d["P3_promoted"])


if __name__ == "__main__":
    unittest.main()
