#!/usr/bin/env python3

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou_roundtrip_transition as transition
import ou_validation as validation


class FakeValidation:
    smoothstep_profile = staticmethod(validation.smoothstep_profile)

    @staticmethod
    def phase_randomize_wave(columns, data, seed):
        del columns, seed
        return data.copy()

    @staticmethod
    def scale_wave_motion(columns, data, scale):
        result = data.copy()
        for name in (
            "disp_x", "disp_y", "disp_z",
            "vel_x", "vel_y", "vel_z",
            "acc_x", "acc_y", "acc_z",
        ):
            result[:, columns.index(name)] *= scale
        return result

    @staticmethod
    def rebuild_body_imu(columns, data):
        del columns
        return data.copy()


class RoundTripTransitionTests(unittest.TestCase):
    def test_default_schedule_keeps_both_crossfades_at_120_seconds(self):
        rise_start, rise_end, fall_start, fall_end = transition.transition_bounds()
        self.assertEqual((rise_start, rise_end), (400.0, 520.0))
        self.assertEqual((fall_start, fall_end), (800.0, 920.0))
        self.assertEqual(rise_end - rise_start, 120.0)
        self.assertEqual(fall_end - fall_start, 120.0)

    def test_roundtrip_weight_is_c2_and_returns_to_low(self):
        times = np.asarray(
            [300.0, 400.0, 460.0, 520.0, 700.0, 800.0, 860.0, 920.0, 1000.0]
        )
        weight, rate, acceleration = transition.roundtrip_profile(validation, times)
        np.testing.assert_allclose(
            weight,
            [0.0, 0.0, 0.5, 1.0, 1.0, 1.0, 0.5, 0.0, 0.0],
            atol=1e-12,
        )
        boundary = np.asarray([400.0, 520.0, 800.0, 920.0])
        _, boundary_rate, boundary_acceleration = transition.roundtrip_profile(
            validation, boundary
        )
        np.testing.assert_allclose(boundary_rate, 0.0, atol=1e-12)
        np.testing.assert_allclose(boundary_acceleration, 0.0, atol=1e-12)
        self.assertTrue(np.all(np.isfinite(rate)))
        self.assertTrue(np.all(np.isfinite(acceleration)))

    def test_generated_record_uses_same_low_realization_after_return(self):
        columns = [
            "time",
            "disp_x", "disp_y", "disp_z",
            "vel_x", "vel_y", "vel_z",
            "acc_x", "acc_y", "acc_z",
            "roll_deg", "pitch_deg", "yaw_deg",
        ]
        times = np.asarray([300.0, 520.0, 700.0, 920.0, 1000.0])
        low = np.zeros((len(times), len(columns)), dtype=float)
        high = np.zeros_like(low)
        low[:, 0] = times
        high[:, 0] = times
        low[:, 1:] = 1.0
        high[:, 1:] = 3.0

        result = transition.make_roundtrip_wave(
            FakeValidation, columns, low, high, seed=11, high_scale=1.0
        )
        np.testing.assert_allclose(result[2, 1:], high[2, 1:], atol=1e-12)
        np.testing.assert_allclose(result[3, 1:], low[3, 1:], atol=1e-12)
        np.testing.assert_allclose(result[4, 1:], low[4, 1:], atol=1e-12)

    def test_scoring_reports_every_segment_and_both_crossfades(self):
        segments = transition.roundtrip_segments()
        self.assertEqual(
            segments,
            (
                ("low_start", 300.0, 400.0),
                ("rise", 400.0, 520.0),
                ("high_recover", 520.0, 640.0),
                ("high", 640.0, 800.0),
                ("fall", 800.0, 920.0),
                ("low_recover", 920.0, 1040.0),
                ("low_return", 1040.0, 1200.0),
            ),
        )
        metrics = {}
        for index, (name, _, _) in enumerate(segments, start=1):
            for metric in transition.SEGMENT_METRICS:
                metrics[f"seg_{name}_{metric}"] = float(index)
        rows = transition.roundtrip_score_rows(metrics, segments)
        self.assertEqual([row["segment"] for row in rows], [s[0] for s in segments])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "scores.tex"
            transition.write_roundtrip_score_table(path, rows)
            text = path.read_text(encoding="utf-8")
        self.assertIn("Rise crossfade", text)
        self.assertIn("Fall crossfade", text)
        self.assertIn(r"\label{tab:ou-roundtrip-scores}", text)
        self.assertNotIn("one-way", text)

    def test_return_must_leave_room_for_unchanged_transition_duration(self):
        with self.assertRaises(ValueError):
            transition.transition_bounds(400.0, 1120.0)


if __name__ == "__main__":
    unittest.main()
