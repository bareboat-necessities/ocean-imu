from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_private_mahony_discrete_comparison as mod  # noqa: E402


class BrmmPrivateMahonyDiscreteComparisonTest(unittest.TestCase):
    def test_ideal_shipping_period_step_preserves_invariant(self):
        d = mod.build()
        self.assertEqual(mod.validate(d), [])
        self.assertTrue(d["ideal_5ms_discrete_PI_invariant_closed"])
        self.assertGreater(d["discrete_metric_decrease_lower"], 0.0)
        self.assertTrue(d["same_BRMM_forcing_as_continuous_invariant"])

    def test_binary32_source_order_composition_is_closed_but_toolchain_generalization_remains_fail_closed(self):
        d = mod.build()
        self.assertTrue(d["shipping_binary32_quaternion_map_error_composed"])
        self.assertTrue(d["shipping_source_order_binary32_discrete_invariant_closed"])
        self.assertFalse(d["toolchain_independent_binary32_invariant_closed"])
        self.assertFalse(d["complete_BRMM_family_materialized_here"])
        self.assertFalse(d["P3_promoted"])
        self.assertFalse(d["source_generator"])


if __name__ == "__main__":
    unittest.main()
