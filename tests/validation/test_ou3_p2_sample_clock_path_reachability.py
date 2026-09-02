from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p2_sample_clock_path_reachability as P2
import ou3_p4_source_path_reachability as LEGACY


class P2SampleClockPathReachabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = P2.build()

    def test_validates(self):
        self.assertEqual(P2.validate(self.d), [])

    def test_shipping_schedule_is_one_sample_lag(self):
        self.assertTrue(self.d["shipping_online_tune_commit_at_sample_start"])
        self.assertTrue(self.d["shipping_online_tune_stage_after_current_sample"])
        self.assertTrue(self.d["aw_covariance_maintenance_cadence_not_tuner_commit_cadence"])
        self.assertFalse(self.d["legacy_0p1s_online_tune_commit_model_used"])
        self.assertFalse(self.d["arbitrary_late_online_tune_commit_modeled"])
        dt = self.d["configured_sample_dt_s"]
        self.assertEqual(self.d["edge_elapsed_interval_s"], [dt, dt])

    def test_graph_is_strictly_sparser_than_cartesian(self):
        self.assertEqual(self.d["partition"]["states"], 800)
        self.assertGreater(self.d["transition_edges"], 0)
        self.assertLess(self.d["transition_edges"], self.d["full_cartesian_transition_edges"])
        self.assertEqual(self.d["full_cartesian_transition_edges"], 800 * 800)
        self.assertLess(self.d["transition_density"], 1.0)

    def test_one_sample_ema_does_not_teleport_to_target(self):
        dt = self.d["configured_sample_dt_s"]
        lo, hi = LEGACY._ema_image(
            (1.0, 1.0), (2.0, 2.0), (1.0, 1.0), dt, max_elapsed=dt
        )
        self.assertGreater(lo, 1.0)
        self.assertLess(hi, 1.01)
        self.assertLess(hi, 2.0)

    def test_raw_tuner_and_filter_sigma_remain_distinct(self):
        self.assertTrue(self.d["raw_tuner_sigma_subfloor_states_included"])
        self.assertLess(self.d["raw_tuner_sigma_partition_lower"], 0.05)
        self.assertEqual(self.d["filter_sigma_floor_mps2"], 0.05)
        self.assertTrue(self.d["filter_sigma_floor_separate_from_tuner_state"])

    def test_fail_closed_if_timing_semantics_are_relaxed(self):
        for key, value in (
            ("arbitrary_late_online_tune_commit_modeled", True),
            ("legacy_0p1s_online_tune_commit_model_used", True),
            ("shipping_online_tune_commit_at_sample_start", False),
            ("shipping_online_tune_stage_after_current_sample", False),
        ):
            with self.subTest(key=key):
                d = deepcopy(self.d)
                d[key] = value
                self.assertNotEqual(P2.validate(d), [])


if __name__ == "__main__":
    unittest.main()
