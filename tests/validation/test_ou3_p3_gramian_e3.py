#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p3_gramian_e3 as E3
import ou3_p3_ltv_translation_ucc_probe as LTV


class GramianE3Tests(unittest.TestCase):
    def test_e3_bound_is_strict_and_stronger_than_trace3(self):
        upper = [1.0, 2.0, 3.0, 4.0]
        old = LTV.ltv_relative_process_floor(
            upper, 0.5, 1.0 / 3.0, 12.0, 0.05
        )
        sharp = E3.sharpen_probe(old, upper)
        self.assertEqual(E3.validate(sharp), [])
        self.assertGreater(
            sharp["relative_process_floor_lower"],
            old["relative_process_floor_lower"],
        )
        self.assertGreater(sharp["improvement_over_trace3"], 1.0)

    def test_common_covariance_scaling_still_has_inverse_effect(self):
        base = LTV.ltv_relative_process_floor(
            [1.0, 2.0, 3.0, 4.0], 0.5, 1.0 / 3.0, 12.0, 0.05
        )
        a = E3.sharpen_probe(base, [1.0, 2.0, 3.0, 4.0])
        b = E3.sharpen_probe(base, [10.0, 20.0, 30.0, 40.0])
        ratio = a["relative_process_floor_lower"] / b["relative_process_floor_lower"]
        self.assertGreater(ratio, 9.999999999)
        self.assertLess(ratio, 10.000000001)

    def test_no_numerical_eigendecomposition_is_used(self):
        d = E3.sharpen(
            [1.0, 1.0, 1.0, 1.0],
            0.5,
            1.0e-12,
            1.0e-3,
        )
        self.assertFalse(d["numerical_eigendecomposition_used"])
        self.assertTrue(math.isfinite(d["relative_process_floor_lower"]))
        self.assertGreater(d["relative_process_floor_lower"], 0.0)


if __name__ == "__main__":
    unittest.main()
