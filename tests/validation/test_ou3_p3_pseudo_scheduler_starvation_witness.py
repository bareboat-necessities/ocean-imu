from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

import ou3_p3_pseudo_scheduler_starvation_witness as W


class PseudoSchedulerStarvationWitnessTests(unittest.TestCase):
    def test_binary32_cycle_repeats_without_firing(self):
        h = W._f32(0.005)
        ratio = W._f32(W._f32(0.015) / W._f32(1.1))
        high, low, tau_high, tau_low = W._cycle_periods(h, ratio, 0.005, 0.25)
        self.assertGreater(high, low)
        self.assertEqual(W._period_from_tau(tau_high, ratio, 0.005, 0.25), high)
        self.assertEqual(W._period_from_tau(tau_low, ratio, 0.005, 0.25), low)

        e0 = W._set_period(low, low)
        self.assertEqual(e0, 0.0)
        ok, e1, first = W._accumulate_no_fire(e0, low, h, W.GAP)
        self.assertTrue(ok)
        self.assertIsNone(first)
        e1 = W._set_period(e1, high)
        ok, e2, first = W._accumulate_no_fire(e1, high, h, W.GAP)
        self.assertTrue(ok)
        self.assertIsNone(first)
        self.assertEqual(e2, low)

    def test_635_sample_cycle_has_no_pseudo_firing(self):
        h = W._f32(0.005)
        ratio = W._f32(W._f32(0.015) / W._f32(1.1))
        high, low, _, _ = W._cycle_periods(h, ratio, 0.005, 0.25)
        fires, boundaries, _ = W._simulate_word(h, high, low)
        self.assertEqual(fires, [])
        self.assertEqual(sum(int(x["samples"]) for x in boundaries), W.TARGET_SAMPLES)
        self.assertEqual(boundaries[0]["period_tag"], "H")
        self.assertEqual(boundaries[1]["period_tag"], "H")
        self.assertEqual(boundaries[2]["period_tag"], "L")


if __name__ == "__main__":
    unittest.main()
