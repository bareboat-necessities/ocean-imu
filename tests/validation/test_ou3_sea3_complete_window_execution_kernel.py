#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_sea3_complete_window_execution_kernel as KERNEL


class Sea3CompleteWindowExecutionKernelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = KERNEL.build()

    def test_kernel_status_is_non_promoting(self):
        self.assertEqual(KERNEL.validate(self.payload), [])
        self.assertTrue(self.payload["typed_execution_kernel_ready"])
        self.assertFalse(self.payload["source_generator"])
        self.assertFalse(self.payload["trajectory_replay_used"])
        self.assertFalse(self.payload["P3_promoted"])

    def test_raw_gyro_and_corrected_rate_remain_distinct_coordinates(self):
        self.assertTrue(
            self.payload["raw_gyro_and_bias_corrected_rate_are_distinct_same_witness_coordinates"]
        )

    def test_floor_is_computed_from_each_mode_covariance(self):
        self.assertTrue(self.payload["covariance_floor_request_not_increment_is_source_event"])
        self.assertTrue(self.payload["covariance_floor_increment_computed_from_current_mode_P"])
        self.assertTrue(self.payload["H18_A21_floor_increments_not_forced_equal"])

    def test_one_sample_smoke_executes_every_shipping_event_type(self):
        s = self.payload["smoke"]
        self.assertEqual(s["samples_executed"], 1)
        self.assertGreaterEqual(s["endpoint_branches"], 1)
        self.assertTrue(s["all_endpoint_H_events_present"])
        self.assertTrue(s["all_endpoint_A_events_present"])
        self.assertTrue(s["decomposition_identity_H_enclosed"])
        self.assertTrue(s["decomposition_identity_A_enclosed"])
        self.assertFalse(s["favorable_frontend_successor_selected"])


if __name__ == "__main__":
    unittest.main()
