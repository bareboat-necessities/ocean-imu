"""The article's statistical units and plotted geometry must stay explicit."""
from __future__ import annotations

import json
import hashlib
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import carma_shared_modes_validation as study
import carma_shared_modes_validation_plots as plots


class PublicationEvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.df = pd.read_csv(study.OUT / "scores.csv")
        cls.protocol = json.loads((study.OUT / "protocol.json").read_text())

    def test_full_rate_complete_paired_matrix(self):
        d = self.df
        self.assertEqual(len(d), 1038)
        self.assertFalse(d.key.duplicated().any())
        self.assertTrue((d.score_samples == 180000).all())
        expected_seeds = {s[0] for s in self.protocol["seed_triplets"]}
        expected_heights = {.27, 1.5, 4., 8.5}
        stationary = d[(d.study == "stationary") & d.method.isin(study.KALMAN)]
        self.assertEqual(len(stationary), 2*4*len(expected_seeds)*5)
        for _, g in stationary.groupby(["method", "spectrum"]):
            self.assertEqual(set(g.wave_seed), expected_seeds)
            for _, draw in g.groupby("wave_seed"):
                self.assertEqual(set(draw.hs_m), expected_heights)
                self.assertEqual(len(draw), 4)
        for _, g in d[d.study == "roundtrip"].groupby("method"):
            self.assertTrue((g.hs_m == 4.0).all())
            self.assertTrue((g.source_record_hs_m == 8.5).all())
            self.assertEqual(set(g.wave_seed), expected_seeds)
            for segment, _, _ in study.rt.roundtrip_segments():
                self.assertTrue(np.isfinite(g[f"seg_{segment}_disp_z_rms_m"]).all())

    def test_retained_bundle_hashes(self):
        manifest=json.loads((study.OUT/"publication_summary.json").read_text())
        for name,digest in manifest["files_sha256"].items():
            with self.subTest(file=name):
                self.assertEqual(hashlib.sha256((study.OUT/name).read_bytes()).hexdigest(),digest)

    def test_unavailable_observer_outputs_are_unscored(self):
        d = self.df[self.df.method.isin(["PII", "TVG-NLO"])]
        self.assertEqual(len(d), 16)
        for field in ("disp_x_rms_m", "disp_y_rms_m", "disp_3d_rms_m",
                      "accel_bias_3d_rms_mps2", "gyro_bias_3d_rms_radps"):
            self.assertTrue(d[field].isna().all(), field)
        self.assertTrue(d[d.method == "TVG-NLO"].yaw_rms_deg.isna().all())

    def test_matched_power_control_and_pooling(self):
        d = self.df[self.df.study == "engine"]
        reference = d[(d.method == "CSM") & (d.condition == "rpm2400") & (d.guard == "off")]
        target = reference.recorded_rms_mps2.iloc[0]
        matched = d[d.condition.str.startswith("matched")]
        np.testing.assert_allclose(matched.recorded_rms_mps2, target, rtol=2e-5)
        summaries = plots.engine_summary(self.df)
        self.assertTrue((summaries.records == 8).all())
        self.assertTrue((summaries.disp_offset_m <= summaries.disp_3d_rms_m + 1e-7).all())

    def test_guard_is_transparent_on_quiet_records(self):
        d = self.df[(self.df.study == "engine") & (self.df.condition == "off")]
        for method in ("CSM", "OU-III"):
            g = d[d.method == method]
            for field in ("disp_z_rms_m", "disp_3d_rms_m", "pitch_rms_deg", "yaw_rms_deg"):
                table = g.pivot(index="input", columns="guard", values=field)
                np.testing.assert_allclose(table["off"], table["guard"], rtol=1e-6, atol=1e-9)
                np.testing.assert_allclose(table["off"], table["deployed"], rtol=1e-6, atol=1e-9)

    def test_drawing_traces_are_bounded_and_keep_endpoints(self):
        paths = list((study.OUT / "traces").glob("*.csv.gz"))
        self.assertTrue(paths)
        for path in paths:
            with self.subTest(path=path.name):
                d = pd.read_csv(path)
                self.assertLessEqual(len(d), 2000)
                # The shared C++ harness serializes six significant digits:
                # adjacent 200 Hz timestamps can coincide beyond 1000 s.
                # Drawing must retain their order and the last sample within
                # that one-sample clock quantization, rather than drop peaks.
                self.assertTrue((np.diff(d.time) >= 0).all())
                self.assertLessEqual(abs(d.time.iloc[-1]-1199.995), .005001)
                if "nominal" in path.name and ("bias" in path.name or "applied" in path.name or "tuner" in path.name):
                    self.assertAlmostEqual(d.time.iloc[0], 0., places=3)


class StatisticalUnitTests(unittest.TestCase):
    def test_interval_resamples_paired_differences_not_separate_arms(self):
        # Large shared nuisance variation cancels in a paired contrast.
        baseline = np.array([1., 100., 2., 200., 3., 300., 4., 400., 5., 500.])
        candidate = baseline + .125
        mean, low, high = plots.mean_interval(candidate - baseline)
        np.testing.assert_allclose([mean, low, high], .125, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
