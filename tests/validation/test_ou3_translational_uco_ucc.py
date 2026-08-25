import importlib.util
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "tools" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mod = load("ou3_translational_uco_ucc")
spread = load("ou3_translational_spread_uco")


class TranslationalUcoUccTests(unittest.TestCase):
    def test_source_uniform_translation_constants_are_strict(self):
        d = mod.build(mod.DEFAULT_HEADER)
        self.assertEqual(mod.validate(d), [])
        self.assertTrue(d["translation_source_complete"])
        self.assertTrue(d["process_ucc"]["pass"])
        self.assertTrue(d["S_observation_uco"]["pass"])
        self.assertGreater(d["process_ucc"]["Q_axis_lambda_min_lower"], 0.0)
        self.assertGreater(d["S_observation_uco"]["information_gramian_lambda_min_lower"], 0.0)
        self.assertFalse(d["continuous_word_enclosed"])
        self.assertEqual(d["theorem_promotion"], "NOT_ESTABLISHED")

    def test_process_bound_covers_complete_tau_sigma_source_box(self):
        d = mod.build(mod.DEFAULT_HEADER)
        p = d["process_ucc"]
        source = mod.SOURCE.build(mod.DEFAULT_HEADER)
        box = source["validated_parameter_box"]["continuous_parameters"]
        self.assertEqual(p["tau_s"], box["tau_aw_s"])
        self.assertEqual(p["sigma_aw_mps2"], box["sigma_aw_mps2"])
        expected = 2.0 * box["sigma_aw_mps2"][0] ** 2 / box["tau_aw_s"][1]
        self.assertLessEqual(p["ou_driving_intensity_lower"], expected)
        self.assertGreater(p["ou_driving_intensity_lower"], 0.0)

    def test_scheduler_bound_uses_source_cadence_not_replay_spacing(self):
        d = mod.build(mod.DEFAULT_HEADER)
        s = d["S_observation_uco"]
        source = mod.SOURCE.build(mod.DEFAULT_HEADER)
        h = source["configured_runtime_assumption"]["imu_dt_outward_interval_s"]
        ts = source["validated_parameter_box"]["continuous_parameters"]["pseudo_update_period_s"]
        self.assertLessEqual(s["pseudo_gap_min_s"], h[0])
        self.assertGreaterEqual(s["pseudo_gap_max_s"], ts[1] + h[1])
        self.assertGreaterEqual(s["aligned_window_s"], 3.0 * s["pseudo_gap_max_s"])
        self.assertEqual(s["aligned_firing_count"], 4)

    def test_spread_search_selects_best_validated_information_bound(self):
        # Deployment theorem currently declares a 1 s recurring-PE word.
        d = spread.build(1.0, mod.DEFAULT_HEADER)
        self.assertEqual(spread.validate(d), [])
        self.assertGreaterEqual(d["admissible_q_max"], 1)
        self.assertGreaterEqual(d["best"]["q"], 1)
        self.assertLessEqual(d["best"]["q"], d["admissible_q_max"])
        self.assertGreater(d["best"]["information_gramian_lambda_min_lower"], 0.0)
        self.assertGreaterEqual(d["information_widening_factor_vs_adjacent_lower"], 1.0)
        self.assertFalse(d["three_S_detectability_used_for_this_UCO"])

    def test_spread_search_never_worse_than_adjacent_for_longer_word(self):
        d = spread.build(4.0, mod.DEFAULT_HEADER)
        self.assertEqual(spread.validate(d), [])
        self.assertGreaterEqual(
            d["best"]["information_gramian_lambda_min_lower"],
            d["adjacent_q1"]["information_gramian_lambda_min_lower"],
        )

    def test_wide_exp_reduction_contains_libm_reference_for_diagnostics(self):
        I = mod.Interval.outward_bounds(38.25, 38.25)
        E = mod._exp_negative_wide(I)
        ref = math.exp(-38.25)
        self.assertLessEqual(E.lo, ref)
        self.assertGreaterEqual(E.hi, ref)
        self.assertGreater(E.lo, 0.0)

    def test_certificate_does_not_use_replay_modules(self):
        for filename in ("ou3_translational_uco_ucc.py", "ou3_translational_spread_uco.py"):
            text = (ROOT / "tools" / filename).read_text()
            for forbidden in ("ou3_exact_replay", "ou3_numerical_certificate", "ou_sweep_common",
                              "path_metrics.npz", "neighborhood_radius_search"):
                self.assertNotIn(forbidden, text)


if __name__ == "__main__":
    unittest.main()
