import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_scalar_ou_enclosure as mod


def contains(bounds, value):
    return float(bounds[0]) <= float(value) <= float(bounds[1])


def direct_source_formula(h, tau):
    x = h / tau
    alpha = math.exp(-x)
    em1 = math.expm1(-x)
    va = -tau * em1
    if abs(x) < 1.0e-2:
        x2 = x * x
        x3 = x2 * x
        x4 = x3 * x
        x5 = x4 * x
        pa = tau * tau * (0.5 * x2 - (1.0 / 6.0) * x3 + (1.0 / 24.0) * x4)
        Sa = tau * tau * tau * ((1.0 / 6.0) * x3 - (1.0 / 24.0) * x4 + (1.0 / 120.0) * x5)
    else:
        pa = tau * tau * (x + em1)
        Sa = tau * tau * tau * (0.5 * x * x - x - em1)
    return x, alpha, em1, va, pa, Sa


def direct_axis_transition(h, tau):
    _, alpha, _, va, pa, Sa = direct_source_formula(h, tau)
    return [
        [1.0, 0.0, 0.0, va],
        [h, 1.0, 0.0, pa],
        [0.5 * h * h, h, 1.0, Sa],
        [0.0, 0.0, 0.0, alpha],
    ]


class ScalarOUEnclosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = mod.build(mod.DEFAULT_HEADER.resolve(), cells_per_branch=24)

    def test_generated_enclosure_validates_and_stays_unpromoted(self):
        self.assertEqual(mod.validate(self.payload), [])
        self.assertEqual(self.payload["schema"], 2)
        self.assertTrue(self.payload["one_step_scalar_ou_enclosed"])
        self.assertTrue(self.payload["one_step_axis_transition_enclosed"])
        self.assertTrue(self.payload["deployment_timing_complete"])
        self.assertFalse(self.payload["continuous_word_enclosed"])
        self.assertFalse(self.payload["nonlinear_word_enclosed"])
        self.assertEqual(self.payload["theorem_promotion"], "NOT_ESTABLISHED")
        self.assertFalse(self.payload["configured_runtime_assumption"]["api_enforces_this_bound"])

    def test_nominal_source_schedule_is_inside_audited_transcendental_range(self):
        x = self.payload["x_h_over_tau_source_box"]
        self.assertGreater(x[0], 0.0)
        self.assertLessEqual(x[1], mod.VT.MAX_ABS_ARGUMENT)
        self.assertLessEqual(x[1], 0.26)

    def test_source_tau_endpoints_and_branch_boundary_are_enclosed(self):
        h = float(self.payload["nominal_imu_dt_source_value_s"])
        tau_box = self.payload["tau_aw_source_box_s"]
        samples = [max(float(tau_box[0]), 0.02), 0.05, 0.2, 0.5, 1.0, min(float(tau_box[1]), 12.0)]
        samples.extend([h / math.nextafter(1.0e-2, 0.0), h / 1.0e-2])
        for tau in samples:
            with self.subTest(tau=tau):
                x, alpha, em1, va, pa, Sa = direct_source_formula(h, tau)
                candidates = [c for c in self.payload["cells"] if contains(c["x_h_over_tau"], x)]
                self.assertTrue(candidates, f"no cell covers x={x} for tau={tau}")
                self.assertTrue(any(
                    contains(c["alpha"], alpha)
                    and contains(c["em1"], em1)
                    and contains(c["phi_va_s"], va)
                    and contains(c["phi_pa_s2"], pa)
                    and contains(c["phi_Sa_s3"], Sa)
                    for c in candidates
                ), f"source formula not enclosed for tau={tau}, x={x}")

    def test_axis_transition_matrix_contains_source_formula(self):
        h = float(self.payload["nominal_imu_dt_source_value_s"])
        for tau in (0.02, 0.5, 12.0):
            x = h / tau
            exact = direct_axis_transition(h, tau)
            candidates = [c for c in self.payload["cells"] if contains(c["x_h_over_tau"], x)]
            self.assertTrue(candidates)
            self.assertTrue(any(
                all(
                    contains(c["transition_axis_v_p_S_aw"][i][j], exact[i][j])
                    for i in range(4) for j in range(4)
                )
                for c in candidates
            ))

    def test_global_bounds_contain_every_cell_and_keep_expected_signs(self):
        g = self.payload["global_bounds"]
        self.assertGreater(g["alpha"][0], 0.0)
        self.assertLessEqual(g["alpha"][1], 1.0)
        self.assertGreater(g["em1"][0], -1.0)
        self.assertLess(g["em1"][1], 0.0)
        self.assertGreaterEqual(g["phi_va_s"][0], 0.0)
        self.assertGreaterEqual(g["phi_pa_s2"][0], 0.0)
        self.assertGreaterEqual(g["phi_Sa_s3"][0], 0.0)
        self.assertGreaterEqual(g["transition_axis_spectral_norm_upper"], 1.0)
        for c in self.payload["cells"]:
            for key in ("alpha", "em1", "phi_va_s", "phi_pa_s2", "phi_Sa_s3"):
                self.assertLessEqual(g[key][0], c[key][0])
                self.assertGreaterEqual(g[key][1], c[key][1])
            self.assertGreaterEqual(
                g["transition_axis_spectral_norm_upper"],
                c["transition_axis_spectral_norm_upper"],
            )

    def test_timing_scope_is_explicit_not_an_api_guard_claim(self):
        scope = self.payload["configured_runtime_assumption"]
        self.assertEqual(scope["sample_period_contract"], "FIXED_SOURCE_NOMINAL")
        self.assertFalse(scope["api_enforces_this_bound"])
        self.assertIn("arbitrary positive finite caller dt", scope["theorem_scope_note"])


if __name__ == "__main__":
    unittest.main()
