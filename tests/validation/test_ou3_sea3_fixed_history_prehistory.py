#!/usr/bin/env python3
from __future__ import annotations

import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools" / "stability"
sys.path.insert(0, str(TOOLS))

import ou3_sea3_fixed_history_prehistory as PRE


class CompleteSea3SameDriverPrehistoryTest(unittest.TestCase):
    def test_prehistory_is_same_continuum_driver_and_joins_word(self) -> None:
        d = PRE.build()
        self.assertEqual(PRE.validate(d), [])
        self.assertEqual(d["canonical_source"], "COMPLETE_SEA3_NORMAL_LIVE_WORD")
        self.assertTrue(d["same_driver_field_as_word"])
        self.assertTrue(d["stationary_fixed_history_correlation_extension"])
        self.assertTrue(d["word_join"]["same_t0_value"])
        self.assertTrue(d["prehistory_inside_declared_acceleration_cap"])
        self.assertGreater(d["prehistory_duration_s"], d["wpe_usable_time_floor_s"])
        self.assertGreater(d["prehistory_duration_s"], d["outer_online_tune_warmup_s"])

    def test_no_startup_surrogate_is_used(self) -> None:
        d = PRE.build()
        for key in (
            "trajectory_replay_used",
            "finite_harmonic_source_used",
            "independent_startup_signal_used",
            "independent_sample_boxes_used",
            "P3_changed",
            "P4_promoted",
        ):
            self.assertFalse(d[key], key)
        self.assertFalse(d["actual_cpp_frontend_TunerReady_verified_here"])


if __name__ == "__main__":
    unittest.main()
