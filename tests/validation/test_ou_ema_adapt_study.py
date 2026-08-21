#!/usr/bin/env python3

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou_ema_adapt_study as study
import ou_validation as validation


class RoundTripTransitionTests(unittest.TestCase):
    def test_default_schedule_preserves_120_second_transition_duration(self):
        self.assertEqual(study.TRANSITION_START_SEC, 400.0)
        self.assertEqual(study.TRANSITION_END_SEC, 520.0)
        self.assertEqual(study.TRANSITION_RETURN_START_SEC, 800.0)
        self.assertEqual(study.TRANSITION_RETURN_END_SEC, 920.0)
        self.assertEqual(
            study.TRANSITION_END_SEC - study.TRANSITION_START_SEC,
            120.0,
        )
        self.assertEqual(
            study.TRANSITION_RETURN_END_SEC - study.TRANSITION_RETURN_START_SEC,
            120.0,
        )

    def test_roundtrip_high_sea_weight_rises_holds_and_returns_to_zero(self):
        times = np.asarray(
            [300.0, 400.0, 460.0, 520.0, 700.0, 800.0, 860.0, 920.0, 1040.0]
        )
        weight, rate, acceleration = study.roundtrip_profile(validation, times)
        np.testing.assert_allclose(
            weight,
            [0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 0.5, 0.0, 0.0],
            atol=1e-12,
        )
        boundary = np.asarray([400.0, 520.0, 800.0, 920.0])
        _, boundary_rate, boundary_acceleration = study.roundtrip_profile(
            validation, boundary
        )
        np.testing.assert_allclose(boundary_rate, 0.0, atol=1e-12)
        np.testing.assert_allclose(boundary_acceleration, 0.0, atol=1e-12)
        self.assertTrue(np.all(np.isfinite(rate)))
        self.assertTrue(np.all(np.isfinite(acceleration)))

    def test_roundtrip_segments_expose_both_transients_and_both_recoveries(self):
        segments, report = study.transition_segment_spec(["roundtrip"])
        self.assertEqual(
            segments,
            "low_start:300:400,rise:400:520,high_recover:520:640,"
            "high:640:800,fall:800:920,low_recover:920:1040,"
            "low_return:1040:1200",
        )
        self.assertEqual(
            report,
            [
                "window",
                "low_start",
                "rise",
                "high_recover",
                "high",
                "fall",
                "low_recover",
                "low_return",
            ],
        )

    def test_roundtrip_cannot_be_mixed_with_legacy_one_way_records(self):
        with self.assertRaises(SystemExit):
            study.transition_segment_spec(["roundtrip", "up"])


if __name__ == "__main__":
    unittest.main()
