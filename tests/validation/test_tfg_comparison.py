"""Pairing, publication and provenance regressions for the TFG companion."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import tfg_comparison as study  # noqa: E402


def fixture():
    p = study.protocol("full")
    rows = []
    for scenario in p["scenarios"]:
        for seed in p["seed_triplets"]:
            for family in study.BINARY:
                # Large PM error must not leak into the JONSWAP aggregate.
                value = 100 if "pmstokes" in scenario else 2
                difference = -0.5 if family == "TFG" else 0
                rows.append({"scenario": scenario, **seed, "family": family,
                    "input_sha256": "a" * 64, "samples": 180000, "window_s": 900.0,
                    **{m: value + difference for m in study.METRICS}})
    return p, rows


class TFGComparisonTests(unittest.TestCase):
    def test_complete_protocol_and_seed_level_effect(self):
        p, rows = fixture()
        self.assertEqual(len(rows), 180)
        summary = study.summarize(rows, p)
        primary = next(r for r in summary if r["group"] == "jonswap" and r["metric"] == "disp_z_pct_hs")
        self.assertEqual(primary["n_pairs"], 10)
        self.assertEqual(primary["ou3_mean"], 2)
        self.assertEqual(primary["difference"], -0.5)
        self.assertEqual(primary["ci95_high"], -0.5)
        self.assertEqual(primary["randomization_patterns"], 1024)
        self.assertIn("TFG has lower", study.publication(summary))
        for row in rows:
            if row["family"] == "TFG":
                for m in study.METRICS:
                    row[m] += 1
        self.assertIn("OU--III has lower", study.publication(study.summarize(rows, p)))

    def test_missing_duplicate_mismatched_and_invalid_rows_rejected(self):
        p, rows = fixture()
        for bad in (rows[:-1], rows + [rows[0]]):
            with self.assertRaises(ValueError):
                study.validate_rows(bad, p)
        for field, value in (("input_sha256", "b" * 64), ("samples", 12000),
                             ("window_s", 60.0), ("disp_z_rms_m", float("nan")),
                             ("imu_noise_seed", 999999)):
            bad = copy.deepcopy(rows)
            bad[0][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                study.validate_rows(bad, p)

    def test_ambient_tuning_does_not_override_shipping_defaults(self):
        p = study.protocol("full")
        with patch.dict("os.environ", {"TFG_SIGMA_COEFF": "99", "OU_III_R_S_COEFF": "99",
                                      "W3D_FIXED_TAU_S": "99", "SF_SIGMA_A_SCALE": "99"}):
            env = study.environment("TFG", p["seed_triplets"][0], p)
        self.assertNotIn("TFG_SIGMA_COEFF", env)
        self.assertNotIn("SF_SIGMA_A_SCALE", env)
        self.assertNotIn("W3D_FIXED_TAU_S", env)
        self.assertEqual(env["TFG_TUNING"], "adaptive")
        self.assertEqual(env["TFG_AW_COV_SYNC"], "1")

    def test_smoke_cannot_be_published(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / "manifest.json").write_text(json.dumps({"protocol": study.protocol("smoke")}))
            with self.assertRaisesRegex(ValueError, "complete declared full"):
                study.check_bundle(out)

    def test_source_dependency_closure_contains_both_estimators(self):
        hashes = study.source_hashes()
        for path in ("src/kalman_tfg/SeaStateFusionFilter_TFG.h",
                     "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
                     "src/util/W3dSimCommon.cpp", "tools/ou_validation.py"):
            self.assertIn(path, hashes)

    def test_all_five_deterministic_methods_are_required_and_rendered(self):
        spec = importlib.util.spec_from_file_location("baseline_comparison",
            ROOT / "plots/kalman_ou_iii/baseline-comparison.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertEqual(set(module.PATTERNS), {"ou3", "ou2", "pii", "nlo", "tfg"})
        regex, _ = module.PATTERNS["tfg"]
        self.assertIsNotNone(regex.match("w3d_jonswap_H1.500_L50.710_fusion_tfg.csv"))
        self.assertIsNone(regex.match("w3d_jonswap_H1.500_L50.710_fusion_tfg_nomag.csv"))
        rows = [{"wave": wave, "hs": hs, **{m: {"z_rms": 0.1, "z_pct": 3,
                    "roll_rms": 0.2, "pitch_rms": 0.3} for m in module.PATTERNS}}
                for wave in ("jonswap", "pmstokes") for hs in module.HEIGHTS]
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "baseline.tex"
            module.write_tables(rows, out)
            self.assertIn(r"\renewcommand{\OUBaselineTFGMeanVerticalPercent}{3.0}", out.read_text())
            self.assertIn("TFG &", out.read_text())
        draw = (ROOT / "plots/kalman_ou_iii/draw_plots.sh").read_text()
        self.assertIn('run_comparison "../../tests/kalman_tfg" "*_fusion_tfg.csv"', draw)
        self.assertIn("! -name '*_fusion_tfg.csv'", draw)


if __name__ == "__main__":
    unittest.main()
