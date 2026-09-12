from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_private_mahony_discrete_comparison as mod  # noqa: E402
import ou3_brmm_private_mahony_discrete_invariant as DISC  # noqa: E402


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

    def test_binary32_charge_reads_the_metric_instead_of_hardcoding_it(self):
        d = DISC.build()
        self.assertEqual(DISC.validate(d), [])
        self.assertTrue(d["metric_read_from_continuous_certificate"])
        # dot(V) <= -sqrt(C)[sqrt(C) q - 2 sup] carries ONE factor of sqrt(C);
        # the retired form used two and claimed twice the proved decrease.
        self.assertTrue(d["first_order_decrease_carries_single_sqrt_C"])
        self.assertAlmostEqual(
            d["guaranteed_first_order_V_decrease_lower"],
            d["continuous_boundary_margin_after_binary32"] * d["dt_s"] * 1.33,
            places=9,
        )
        self.assertGreater(d["discrete_V_margin_lower"], 0.0)
        self.assertGreater(d["one_step_chart_round_margin_rad"], 0.0)


if __name__ == "__main__":
    unittest.main()
