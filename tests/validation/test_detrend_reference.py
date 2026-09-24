"""The committed 1D fixture follows a separate binary64 mathematical oracle."""

import csv
import hashlib
import importlib.util
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("detrend_reference", ROOT / "tools/detrend_reference.py")
MODEL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODEL)


class DetrendReferenceTests(unittest.TestCase):
    def test_committed_fixture_matches_independent_equations(self):
        subprocess.run([sys.executable, str(ROOT / "tools/detrend_reference.py")], check=True)

    def test_reference_replay_preserves_input_text(self):
        with MODEL.FIXTURE.open(newline="", encoding="utf-8") as source:
            rows = list(csv.DictReader(source))
        inputs = "".join(row["x_axis"] + "," + row["original_cm"] + "\n" for row in rows)
        self.assertEqual(hashlib.sha256(inputs.encode()).hexdigest(),
                         "1883a669a3739d402a8d8fd3e2583b7632aa721553bf62bb7a7af150b3811fab")
        self.assertEqual(len(rows), 1830)
        regenerated = MODEL.reference_rows(rows)
        self.assertEqual(len(rows), len(regenerated))
        for before, after in zip(rows, regenerated):
            self.assertEqual((before["x_axis"], before["original_cm"]),
                             (after["x_axis"], after["original_cm"]))
        # 0.35 is the current scalar default; the old fixture used 0.25.
        self.assertAlmostEqual(regenerated[0]["baseline_cutoff_hz"], 0.35 * 0.12)

    def test_wave_accuracy_benchmark_has_an_explicit_operating_point(self):
        source = (ROOT / "tests/detrend/detrend-wave-test.cpp").read_text()
        self.assertIn("kRmsGateFractionOfHeight = 0.16;", source)
        self.assertIn("cfg.baseline_cutoff_fraction = 0.25f;", source)
        self.assertIn("cfg.enable_wave_cleanup = true;", source)
        self.assertIn("cfg.cleanup_stages = 1;", source)
        self.assertIn("AdaptiveWaveDetrender detrender(cfg);", source)
        self.assertIn("improvement_ratio,baseline_cutoff_fraction,cleanup_stages", source)

    def test_invalid_input_is_not_silently_rebaselined(self):
        for rows in ([], [{"x_axis": "0", "original_cm": "nan"}],
                     [{"x_axis": "0", "original_cm": "1"}] * 2):
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                MODEL.reference_rows(rows)


if __name__ == "__main__":
    unittest.main()
