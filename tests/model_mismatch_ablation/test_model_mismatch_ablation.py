import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import model_mismatch_ablation as mma  # noqa: E402


class ModelMismatchAblationTests(unittest.TestCase):
    def test_pooled_rms_is_concatenated_equal_duration_rms(self):
        rows = [{"x": 3.0}, {"x": 4.0}]
        self.assertAlmostEqual(mma.pooled_rms(rows, "x"), math.sqrt(12.5))

    def test_summary_covers_all_three_families_and_reference_normalization(self):
        rows = []
        for family in mma.FAMILIES:
            for index, record in enumerate(mma.RECORDS):
                base = float(index + 1)
                row = {
                    "family": family,
                    "disp_z_ref_rms_m": 2.0 * base,
                }
                for field in mma.RMS_FIELDS:
                    row[field] = base
                rows.append(row)

        summary = mma.summarize(rows)
        self.assertEqual([row["family"] for row in summary], list(mma.FAMILIES))
        expected = math.sqrt(sum(float(i * i) for i in range(1, 9)) / 8.0)
        for row in summary:
            self.assertEqual(row["records"], 8)
            self.assertAlmostEqual(row["disp_z_rms_m"], expected)
            self.assertAlmostEqual(row["disp_z_ref_rms_m"], 2.0 * expected)
            self.assertAlmostEqual(row["disp_z_pct_refrms"], 50.0)

    def test_report_states_interpretation_boundary(self):
        rows = []
        for family in mma.FAMILIES:
            for record in mma.RECORDS:
                row = {
                    "family": family,
                    "spectrum": record.spectrum,
                    "hs_m": record.hs_m,
                    "input": record.filename,
                    "disp_z_ref_rms_m": 1.0,
                    "disp_z_pct_refrms": 10.0,
                    "disp_z_pct_hs": 1.0,
                    "wave_period_s": 5.0,
                    "tau_applied_s": 2.5,
                    "sigma_applied_mps2": 0.5,
                }
                for field in mma.RMS_FIELDS:
                    row[field] = 0.1
                rows.append(row)
        text = mma.markdown_report(rows, mma.summarize(rows), 900.0, "abc123")
        self.assertIn("--no-noise", text)
        self.assertIn("not a claim of pure plant-model mismatch", text)
        self.assertIn("OU-II", text)
        self.assertIn("OU-III", text)
        self.assertIn("TFG", text)


if __name__ == "__main__":
    unittest.main()
