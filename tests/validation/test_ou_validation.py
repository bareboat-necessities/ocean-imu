import math
import sys
import unittest
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))

import ou_validation as validation  # noqa: E402


class WaveSurrogateTests(unittest.TestCase):
    def setUp(self):
        self.columns = [
            "time",
            "disp_x", "disp_y", "disp_z",
            "vel_x", "vel_y", "vel_z",
            "acc_x", "acc_y", "acc_z",
            "acc_bx", "acc_by", "acc_bz",
            "gyro_x", "gyro_y", "gyro_z",
            "roll_deg", "pitch_deg", "yaw_deg",
            "q_wb_zu_w", "q_wb_zu_x", "q_wb_zu_y", "q_wb_zu_z",
        ]
        count = 4096
        time = np.arange(count) * validation.DT_SECONDS
        self.data = np.zeros((count, len(self.columns)), dtype=np.float64)
        self.data[:, 0] = time
        frequencies = (0.25, 0.45, 0.75)
        for offset, name in enumerate(("disp_x", "disp_y", "disp_z")):
            self.data[:, self.columns.index(name)] = np.sin(
                2.0 * np.pi * frequencies[offset] * time + 0.2 * offset
            )
        for group, derivative in (("vel", 1.0), ("acc", 2.0)):
            for offset, axis in enumerate("xyz"):
                self.data[:, self.columns.index(f"{group}_{axis}")] = (
                    derivative
                    * np.cos(2.0 * np.pi * frequencies[offset] * time + 0.2 * offset)
                )
        self.data[:, self.columns.index("roll_deg")] = 3.0 * np.sin(0.7 * time)
        self.data[:, self.columns.index("pitch_deg")] = 2.0 * np.cos(0.9 * time)
        self.data[:, self.columns.index("yaw_deg")] = 5.0
        self.data = validation.rebuild_body_imu(self.columns, self.data)

    def test_common_phase_preserves_auto_and_cross_spectra(self):
        randomized = validation.phase_randomize_wave(self.columns, self.data, seed=73)
        names = (
            "disp_x", "disp_y", "disp_z",
            "vel_x", "vel_y", "vel_z",
            "acc_x", "acc_y", "acc_z",
            "roll_deg", "pitch_deg", "yaw_deg",
        )
        indices = validation._column_indices(self.columns, names)
        original_fft = np.fft.rfft(self.data[:, indices], axis=0)
        randomized_fft = np.fft.rfft(randomized[:, indices], axis=0)
        np.testing.assert_allclose(
            np.abs(randomized_fft), np.abs(original_fft), rtol=2e-11, atol=2e-9
        )
        np.testing.assert_allclose(
            randomized_fft[:, 0] * np.conj(randomized_fft[:, 1]),
            original_fft[:, 0] * np.conj(original_fft[:, 1]),
            rtol=2e-11,
            atol=2e-8,
        )
        quaternion_indices = validation._column_indices(
            self.columns,
            ("q_wb_zu_w", "q_wb_zu_x", "q_wb_zu_y", "q_wb_zu_z"),
        )
        norms = np.linalg.norm(randomized[:, quaternion_indices], axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=2e-12)

    def test_phase_seed_changes_realization(self):
        first = validation.phase_randomize_wave(self.columns, self.data, seed=1)
        second = validation.phase_randomize_wave(self.columns, self.data, seed=2)
        index = self.columns.index("disp_z")
        self.assertGreater(np.max(np.abs(first[:, index] - second[:, index])), 0.1)

    def test_smoothstep_has_exact_endpoints(self):
        values = validation.smoothstep_weight(np.array([0.0, 2.0, 3.0, 4.0, 8.0]), 2.0, 4.0)
        np.testing.assert_allclose(values, (0.0, 0.0, 0.5, 1.0, 1.0))


class StatisticsTests(unittest.TestCase):
    @staticmethod
    def row(family, mode, repetition, value):
        row = {
            "scenario": "sea",
            "family": family,
            "mode": mode,
            "repetition": repetition,
            "wave_phase_seed": repetition,
            "imu_noise_seed": 10 + repetition,
            "initialization_seed": 20 + repetition,
        }
        row.update({metric: value for metric in validation.METRIC_NAMES})
        return row

    def test_seed_broadcasting_is_paired(self):
        seeds = validation.broadcast_seed_triplets([1, 2], [10], [20, 21])
        self.assertEqual(
            seeds,
            [
                validation.SeedTriplet(1, 10, 20),
                validation.SeedTriplet(2, 10, 21),
            ],
        )

    def test_summary_and_paired_effect_use_all_repetitions(self):
        rows = []
        for repetition, ou2, ou3 in ((1, 3.0, 2.0), (2, 5.0, 3.0), (3, 7.0, 4.0)):
            rows.append(self.row("OU_II", "Adaptive", repetition, ou2))
            rows.append(self.row("OU_III", "Adaptive", repetition, ou3))
        summary = validation.summarize_rows(rows, bootstrap_resamples=500, stats_seed=7)
        selected = next(
            row for row in summary
            if row["family"] == "OU_II"
            and row["mode"] == "Adaptive"
            and row["metric"] == "disp_3d_rms_m"
        )
        self.assertEqual(selected["n"], 3)
        self.assertAlmostEqual(selected["mean"], 5.0)
        self.assertAlmostEqual(selected["std"], 2.0)

        effects = validation.paired_effect_rows(rows, 500, 7)
        selected_effect = next(
            row for row in effects
            if row["comparison"] == "OU_III_minus_OU_II"
            and row["metric"] == "disp_3d_rms_m"
        )
        self.assertEqual(selected_effect["n_pairs"], 3)
        self.assertAlmostEqual(selected_effect["mean_paired_difference"], -2.0)
        self.assertTrue(math.isfinite(selected_effect["cohen_dz"]))
        self.assertTrue(math.isfinite(selected_effect["hedges_gz"]))

    def test_machine_readable_metric_parser(self):
        parsed = validation.parse_validation_metrics(
            "VALIDATION_METRICS family=OU_II tuning_mode=adaptive "
            "input=wave.csv window_s=900 samples=180000 disp_z_rms_m=0.25"
        )
        self.assertEqual(parsed["family"], "OU_II")
        self.assertEqual(parsed["samples"], 180000)
        self.assertEqual(parsed["disp_z_rms_m"], 0.25)


if __name__ == "__main__":
    unittest.main()
