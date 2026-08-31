from copy import deepcopy
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_p4_source_path_reachability as R


class P2SourcePathReachabilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = R.build()

    def test_validates(self):
        self.assertEqual(R.validate(self.d), [])

    def test_source_only_and_not_promoted(self):
        self.assertTrue(self.d["source_only"])
        self.assertFalse(self.d["trajectory_replay_used"])
        self.assertFalse(self.d["usable_P4_promoted"])
        self.assertEqual(self.d["P2_SOURCE_PATH_CERTIFICATE"], "PASS")

    def test_graph_is_nonempty(self):
        self.assertGreater(self.d["partition"]["states"], 0)
        self.assertGreater(self.d["transition_edges"], 0)
        self.assertGreater(self.d["strongly_connected_components"], 0)

    def test_raw_tuner_sigma_is_not_conflated_with_filter_floor(self):
        self.assertTrue(self.d["raw_tuner_sigma_subfloor_states_included"])
        self.assertLess(self.d["raw_tuner_sigma_partition_lower"], 0.05)
        self.assertEqual(self.d["filter_sigma_floor_mps2"], 0.05)
        self.assertTrue(self.d["filter_sigma_floor_separate_from_tuner_state"])
        lo = R._filter_sigma_box((0.001, 0.01))
        self.assertLessEqual(lo[0], 0.05)
        self.assertGreaterEqual(lo[1], 0.05)
        hi = R._filter_sigma_box((0.10, 0.20))
        self.assertLessEqual(hi[0], 0.10)
        self.assertGreaterEqual(hi[1], 0.20)

    def test_path_arithmetic_fails_closed_on_unqualified_source_shortcuts(self):
        self.assertTrue(self.d["source_float_literals_rounded_as_binary32"])
        self.assertTrue(self.d["validated_exponential_used_for_ema"])
        self.assertTrue(self.d["arbitrary_late_commit_overapproximated"])
        self.assertIsNone(self.d["inter_commit_elapsed_upper_assumed_s"])
        self.assertTrue(self.d["RS_discrepancy_slew_horizon_covered"])
        self.assertTrue(self.d["RS_target_full_deployed_clamp_overapprox"])
        self.assertFalse(self.d["RS_target_powf_tightening_used"])

    def test_arbitrarily_late_commit_image_contains_target(self):
        image = R._ema_image((1.0, 1.2), (2.0, 2.2), (0.5, 1.0), 0.1)
        self.assertLessEqual(image[0], 2.0)
        self.assertGreaterEqual(image[1], 2.2)

    def test_old_worst_corner_is_explicit(self):
        self.assertGreater(self.d["old_worst_corner_state_count"], 0)
        self.assertGreaterEqual(self.d["old_worst_corner_states_in_any_recurrent_SCC"], 0)
        self.assertLessEqual(
            self.d["old_worst_corner_states_in_any_recurrent_SCC"],
            self.d["old_worst_corner_state_count"],
        )
        self.assertIn(self.d["old_worst_corner_has_internal_recurrent_cycle"], (True, False))

    def test_mutations_of_source_completeness_fail_validation(self):
        mutations = (
            ("raw_tuner_sigma_subfloor_states_included", False),
            ("filter_sigma_floor_separate_from_tuner_state", False),
            ("validated_exponential_used_for_ema", False),
            ("arbitrary_late_commit_overapproximated", False),
            ("RS_discrepancy_slew_horizon_covered", False),
            ("RS_target_full_deployed_clamp_overapprox", False),
            ("RS_target_powf_tightening_used", True),
        )
        for key, value in mutations:
            with self.subTest(key=key):
                d = deepcopy(self.d)
                d[key] = value
                self.assertNotEqual(R.validate(d), [])


if __name__ == "__main__":
    unittest.main()
