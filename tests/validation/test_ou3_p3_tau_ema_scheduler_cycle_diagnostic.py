from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))

import ou3_p4_source_path_reachability as PATH
import ou3_source_reachable_matrix_p3 as BASE
import ou3_p3_pseudo_scheduler_starvation_witness as STARVE
import ou3_p3_tau_ema_scheduler_cycle_diagnostic as D


class TauEmaSchedulerCycleDiagnosticTests(unittest.TestCase):
    def _setup(self):
        c = PATH._constants()
        sched = BASE.source_schedule()
        h = D._f32(c["dt"])
        ratio = D._f32(sched["pseudo_ratio"])
        pmin = D._f32(sched["pseudo_min_s"])
        pmax = D._f32(sched["pseudo_max_s"])
        _, low, _, tau_low = STARVE._cycle_periods(h, ratio, pmin, pmax)
        tau_high, high = D._minimal_high_tau(low, ratio, pmin, pmax)
        return c, sched, low, tau_low, high, tau_high

    def test_adjacent_high_period_clears_scheduler_tolerance(self):
        c, sched, low, tau_low, high, tau_high = self._setup()
        self.assertGreater(high, D._f32(low + D._scheduler_tolerance(low)))
        self.assertGreater(tau_high, tau_low)
        self.assertLess(tau_high - tau_low, 2.0e-4)

    def test_exact_gap_tau_images_exist_through_legal_frequencies(self):
        c, _, _, tau_low, _, tau_high = self._setup()
        up = D._find_exact_target(tau_low, tau_high, c)
        down = D._find_exact_target(tau_high, tau_low, c)
        self.assertEqual(D._tau_ema_samples(tau_low, up, c), tau_high)
        self.assertEqual(D._tau_ema_samples(tau_high, down, c), tau_low)

        f_up = D._find_frequency_for_target(up, c)
        f_down = D._find_frequency_for_target(down, c)
        lo, hi = D._effective_frequency_bounds(c)
        self.assertTrue(lo <= f_up <= hi)
        self.assertTrue(lo <= f_down <= hi)
        self.assertEqual(D._tau_target_from_frequency(f_up, c), up)
        self.assertEqual(D._tau_target_from_frequency(f_down, c), down)

    def test_tau_ema_cycle_has_no_pseudo_firing(self):
        c, sched, _, tau_low, _, tau_high = self._setup()
        up = D._find_exact_target(tau_low, tau_high, c)
        down = D._find_exact_target(tau_high, tau_low, c)
        fires, trace, _, _, _, _ = D._simulate(c, sched, tau_low, tau_high, up, down)
        self.assertEqual(fires, [])
        self.assertEqual(sum(int(x["samples"]) for x in trace), D.TARGET_SAMPLES)
        self.assertEqual(trace[0]["period_tag"], "H")
        self.assertEqual(trace[1]["period_tag"], "H")
        self.assertEqual(trace[2]["period_tag"], "L")
        self.assertEqual(trace[2]["elapsed_after_fmod_s"], 0.0)


if __name__ == "__main__":
    unittest.main()
