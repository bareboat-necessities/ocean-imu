import importlib.util
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location(
    "ou3_scalar_ou_enclosure", ROOT / "tools" / "ou3_scalar_ou_enclosure.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def contains(bounds, value):
    return float(bounds[0]) <= float(value) <= float(bounds[1])


def direct_source_formula(h, tau):
    x = h / tau
    alpha = math.exp(-x)
    em1 = math.expm1(-x)
    if abs(x) < 1.0e-2:
        x2 = x * x
        x3 = x2 * x
        x4 = x3 * x
        x5 = x4 * x
        pa = tau * tau * (0.5 * x2 - (1.0 / 6.0) * x3 + (1.0 / 24.0) * x4)
        Sa = tau * tau * tau * (
            (1.0 / 6.0) * x3 - (1.0 / 24.0) * x4 + (1.0 / 120.0) * x5
        )
    else:
        pa = tau * tau * (x + em1)
        Sa = tau * tau * tau * (0.5 * x * x - x - em1)
    return x, alpha, em1, pa, Sa


class ScalarOUEnclosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = mod.build(mod.DEFAULT_HEADER.resolve(), cells_per_branch=24)

    def test_generated_enclosure_validates_and_stays_unpromoted(self):
        failures = mod.validate(self.payload)
        self.assertEqual(failures, [])
        self.assertTrue(self.payload["one_step_scalar_ou_enclosed"])
        self.assertFalse(self.payload["deployment_timing_complete"])
        self.assertFalse(self.payload["continuous_word_enclosed"])
        self.assertFalse(self.payload["nonlinear_word_enclosed"])
        self.assertEqual(self.payload["theorem_promotion"], "NOT_ESTABLISHED")

    def test_nominal_source_schedule_is_inside_audited_transcendental_range(self):
        x = self.payload["x_h_over_tau_source_box"]
        self.assertGreater(x[0], 0.0)
        self.assertLessEqual(x[1], mod.VT.MAX_ABS_ARGUMENT)
        self.assertLessEqual(x[1], 0.26)

    def test_source_tau_endpoints_and_branch_boundary_are_enclosed(self):
        h = float(self.payload["nominal_imu_dt_source_value_s"])
        tau_box = self.payload["tau_aw_source_box_s"]
        samples = [
            max(float(tau_box[0]), 0.02),
            0.05,
            0.2,
            0.5,
            1.0,
            min(float(tau_box[1]), 12.0),
        ]
        for tau in samples:
            with self.subTest(tau=tau):
                x, alpha, em1, pa, Sa = direct_source_formula(h, tau)
                candidates = [
                    c for c in self.payload["cells"]
                    if contains(c["x_h_over_tau"], x)
                ]
                self.assertTrue(candidates, f"no cell covers x={x} for tau={tau}")
                self.assertTrue(
                    any(
                        contains(c["alpha"], alpha)
                        and contains(c["em1"], em1)
                        and contains(c["phi_pa_s2"], pa)
                        and contains(c["phi_Sa_s3"], Sa)
                        for c in candidates
                    ),
                    f"source formula not enclosed for tau={tau}, x={x}",
                )

    def test_global_bounds_contain_every_cell_and_keep_expected_signs(self):
        g = self.payload["global_bounds"]
        self.assertGreater(g["alpha"][0], 0.0)
        self.assertLessEqual(g["alpha"][1], 1.0)
        self.assertGreater(g["em1"][0], -1.0)
        self.assertLess(g["em1"][1], 0.0)
        self.assertGreaterEqual(g["phi_pa_s2"][0], 0.0)
        self.assertGreaterEqual(g["phi_Sa_s3"][0], 0.0)
        for c in self.payload["cells"]:
            for key in ("alpha", "em1", "phi_pa_s2", "phi_Sa_s3"):
                self.assertLessEqual(g[key][0], c[key][0])
                self.assertGreaterEqual(g[key][1], c[key][1])

    def test_result_names_external_dt_as_an_open_deployment_obligation(self):
        msg = self.payload["deployment_timing_open_obligation"]
        self.assertIn("updateTime", msg)
        self.assertIn("dt", msg)
        self.assertIn("admissible", msg)


if __name__ == "__main__":
    unittest.main()
