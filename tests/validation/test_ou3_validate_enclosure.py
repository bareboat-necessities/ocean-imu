import importlib.util
from fractions import Fraction
import pathlib
import random
import sys
import types
import unittest


ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "ou3_validate_enclosure.py"
sys.path.insert(0, str(ROOT / "tools"))

import ou3_source_domain_contract as SOURCE  # noqa: E402
import ou3_interval as IA  # noqa: E402
from ou3_interval import Interval  # noqa: E402
import ou3_source_interval_box as SOURCE_BOX  # noqa: E402


def load_tool():
    base = types.ModuleType("ou3_numerical_certificate")
    base.DEFAULT_OUT = ROOT / "reports" / "results" / "ou3_numerical_certificate"
    sys.modules.setdefault("ou3_numerical_certificate", base)
    spec = importlib.util.spec_from_file_location("ou3_validate_enclosure", TOOL)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def contract_mode():
    return {
        "recommended_word_horizon_s": 16.0,
        "executed_reference_only": {
            "relative_Riccati_injection_margin_worst": 0.05,
            "Sigma_endpoint_lambda_min": 0.01,
            "Sigma_endpoint_lambda_max": 100.0,
        },
    }


def mode_payload(level=10.0, ratio=0.8):
    return {
        "source_complete": True,
        "outward_rounded": True,
        "word_horizon_s": 16.0,
        "relative_Riccati_injection_margin_lower": 0.04,
        "Sigma_lambda_min_lower": 0.005,
        "Sigma_lambda_max_upper": 120.0,
        "prefix_information_gain_upper": 3.0,
        "theta_star": 1.0,
        "endpoint_W_ratio_upper": ratio,
        "certified_level_W": level,
        "all_word_prefixes_safe": True,
    }


def exact(x: float) -> Fraction:
    return Fraction.from_float(float(x))


def assert_contains_exact(test: unittest.TestCase, interval: Interval, value: Fraction):
    test.assertLessEqual(exact(interval.lo), value)
    test.assertLessEqual(value, exact(interval.hi))


class Ou3IntervalArithmeticTests(unittest.TestCase):
    def test_basic_operations_enclose_exact_binary64_endpoint_arithmetic(self):
        rng = random.Random(403)
        for _ in range(300):
            a0, a1 = sorted((rng.uniform(-20.0, 20.0), rng.uniform(-20.0, 20.0)))
            b0, b1 = sorted((rng.uniform(-20.0, 20.0), rng.uniform(-20.0, 20.0)))
            A = Interval(a0, a1)
            B = Interval(b0, b1)

            C = A + B
            assert_contains_exact(self, C, exact(a0) + exact(b0))
            assert_contains_exact(self, C, exact(a1) + exact(b1))

            C = A - B
            assert_contains_exact(self, C, exact(a0) - exact(b1))
            assert_contains_exact(self, C, exact(a1) - exact(b0))

            products = [exact(x) * exact(y) for x in (a0, a1) for y in (b0, b1)]
            C = A * B
            assert_contains_exact(self, C, min(products))
            assert_contains_exact(self, C, max(products))

    def test_division_encloses_exact_values_when_denominator_avoids_zero(self):
        rng = random.Random(404)
        for _ in range(200):
            a0, a1 = sorted((rng.uniform(-10.0, 10.0), rng.uniform(-10.0, 10.0)))
            b0, b1 = sorted((rng.uniform(0.2, 10.0), rng.uniform(0.2, 10.0)))
            A = Interval(a0, a1)
            B = Interval(b0, b1)
            quotients = [exact(x) / exact(y) for x in (a0, a1) for y in (b0, b1)]
            C = A / B
            assert_contains_exact(self, C, min(quotients))
            assert_contains_exact(self, C, max(quotients))

    def test_square_handles_zero_crossing(self):
        C = Interval(-3.0, 2.0).square()
        self.assertEqual(C.lo, 0.0)
        assert_contains_exact(self, C, Fraction(9, 1))

    def test_interval_matrix_product_contains_exact_product(self):
        A = [
            [Interval.outward_bounds(0.9, 1.1), Interval.outward_bounds(-0.2, -0.1)],
            [Interval.outward_bounds(0.3, 0.4), Interval.outward_bounds(1.9, 2.1)],
        ]
        B = IA.matrix_point([[2.0, -1.0], [0.5, 3.0]])
        C = IA.matrix_mul(A, B)
        for a00, a01, a10, a11 in (
            (0.9, -0.2, 0.3, 1.9),
            (1.1, -0.1, 0.4, 2.1),
            (1.0, -0.15, 0.35, 2.0),
        ):
            exact_c = (
                (2.0 * a00 + 0.5 * a01, -a00 + 3.0 * a01),
                (2.0 * a10 + 0.5 * a11, -a10 + 3.0 * a11),
            )
            for i in range(2):
                for j in range(2):
                    self.assertTrue(C[i][j].contains(exact_c[i][j]))

    def test_gershgorin_certifies_uniform_spd_without_eigensolver(self):
        A = [
            [Interval.outward_bounds(3.9, 4.1), Interval.outward_bounds(-0.2, 0.2)],
            [Interval.outward_bounds(-0.2, 0.2), Interval.outward_bounds(2.9, 3.1)],
        ]
        ok, lower = IA.symmetric_positive_definite_gershgorin(A)
        self.assertTrue(ok)
        self.assertGreater(lower, 2.6)
        self.assertGreaterEqual(IA.symmetric_gershgorin_upper(A), 4.1)

    def test_gershgorin_refuses_to_claim_spd_when_box_is_too_wide(self):
        A = [
            [Interval.outward_bounds(0.9, 1.1), Interval.outward_bounds(-2.0, 2.0)],
            [Interval.outward_bounds(-2.0, 2.0), Interval.outward_bounds(0.9, 1.1)],
        ]
        ok, lower = IA.symmetric_positive_definite_gershgorin(A)
        self.assertFalse(ok)
        self.assertLess(lower, 0.0)

    def test_validated_spectral_norm_bound_uses_absolute_sums(self):
        A = IA.matrix_point([[1.0, -2.0], [3.0, 4.0]])
        bound = IA.matrix_spectral_norm_upper(A)
        self.assertGreaterEqual(bound, 7.0)

    def test_source_box_encloses_every_source_domain_endpoint_without_promoting_theorem(self):
        source = SOURCE.build(SOURCE.DEFAULT_HEADER.resolve())
        box = SOURCE_BOX.build(SOURCE.DEFAULT_HEADER.resolve())
        failures = SOURCE_BOX.validate(box, SOURCE.DEFAULT_HEADER.resolve())
        self.assertEqual(failures, [])
        self.assertTrue(box["validated_arithmetic"])
        self.assertTrue(box["outward_rounded"])
        self.assertEqual(box["theorem_promotion"], "NOT_ESTABLISHED")
        self.assertFalse(box["continuous_word_enclosed"])
        self.assertFalse(box["nonlinear_word_enclosed"])
        for name, endpoints in source["continuous_parameters"].items():
            I = Interval(*box["continuous_parameters"][name])
            self.assertTrue(I.contains(endpoints[0]), name)
            self.assertTrue(I.contains(endpoints[1]), name)


class Ou3ValidatedEnclosureTests(unittest.TestCase):
    def test_nonlinear_margin_is_derived_from_endpoint_ratio(self):
        tool = load_tool()
        out = tool.validate_mode("H", mode_payload(ratio=0.8), contract_mode())
        self.assertTrue(out["pass"], out["failures"])
        self.assertAlmostEqual(out["mu_W_lower"], 0.2)

    def test_sampled_counterexample_caps_validated_level(self):
        tool = load_tool()
        out = tool.validate_mode(
            "H", mode_payload(level=12.0), contract_mode(), {"first_fail_W": 11.0}
        )
        self.assertFalse(out["pass"])
        self.assertTrue(any("sampled failure" in x for x in out["failures"]))

    def test_hybrid_margin_is_recomputed_and_dimension_change_explicit(self):
        tool = load_tool()
        modes = {
            "H": {"certified_level_W": 10.0},
            "A": {"certified_level_W": 10.0},
        }
        common = dict(
            source_complete=True,
            outward_rounded=True,
            source_level_W_upper=2.0,
            jump_gain_upper=0.5,
            additive_W_upper=0.2,
            destination_level_W=5.0,
            destination_mode="H",
        )
        rows = [
            {"kind": "startup_handoff", **common},
            {
                "kind": "held_to_active",
                **{**common, "destination_mode": "A"},
                "source_dimension": 18,
                "destination_dimension": 21,
                "dimension_change_handled_by_embedding": True,
                "new_coordinate_W_upper": 0.3,
                "inward_margin_lower": 999.0,
            },
            {"kind": "magnetic_regauge", **common},
            {"kind": "tilt_reset", **common},
            {"kind": "cooldown", **common},
        ]
        out = tool.validate_hybrid(rows, modes)
        self.assertTrue(out["pass"], out["failures"])
        held = next(x for x in out["bounds"] if x["kind"] == "held_to_active")
        self.assertAlmostEqual(held["post_jump_W_upper"], 1.5)
        self.assertAlmostEqual(held["inward_margin_lower"], 3.5)

    def test_stochastic_probability_is_recomputed_from_source_moments(self):
        tool = load_tool()
        modes = {
            "H": {"endpoint_W_ratio_upper": 0.5, "certified_level_W": 100.0, "word_horizon_s": 0.01},
            "A": {"endpoint_W_ratio_upper": 0.6, "certified_level_W": 100.0, "word_horizon_s": 0.01},
        }
        source_noise = {
            "schema": 1,
            "source_generated_not_trajectory_fit": True,
            "standardized_increment": {"dimension": 18, "covariance_upper_identity": True},
            "physical_scales": {"imu_dt_s": 0.005},
        }
        payload = {
            "source_complete": True,
            "outward_rounded": True,
            "localization_prefix_safe": True,
            "localization_radius_standardized": 20.0,
            "word_samples_upper": 2,
            "finite_horizon_words": 1,
            "funnel_level_a": 100.0,
            "W0_upper": 1.0,
            "L_X_upper": 0.1,
            "G_bar_upper": 1e-4,
            "c_zw_upper": 0.0,
            "r_star_upper": 1.0,
            "c_ww_upper": 0.0,
            "g_W_upper": 1e-4,
            "h_W_upper": 1e-4,
            "b_W_upper": 0.0,
            "v_W_upper": 0.0,
            "finite_horizon_failure_probability_upper": 0.0,
        }
        out = tool.validate_stochastic(payload, modes, source_noise)
        self.assertTrue(out["pass"], out["failures"])
        self.assertEqual(out["source_noise_moments"]["s2_upper"], 18.0)
        self.assertEqual(out["source_noise_moments"]["s4_upper"], 360.0)
        self.assertEqual(out["b_W_upper"], 100.0)
        self.assertEqual(out["v_W_upper"], 2500.0)
        self.assertLess(out["finite_horizon_failure_probability_upper"], 1.0)

    def test_stochastic_localization_radius_must_exceed_gaussian_rms(self):
        tool = load_tool()
        modes = {
            "H": {"endpoint_W_ratio_upper": 0.5, "certified_level_W": 10.0, "word_horizon_s": 0.01},
            "A": {"endpoint_W_ratio_upper": 0.5, "certified_level_W": 10.0, "word_horizon_s": 0.01},
        }
        source_noise = {
            "schema": 1,
            "source_generated_not_trajectory_fit": True,
            "standardized_increment": {"dimension": 18, "covariance_upper_identity": True},
            "physical_scales": {"imu_dt_s": 0.005},
        }
        payload = {
            "source_complete": True, "outward_rounded": True, "localization_prefix_safe": True,
            "localization_radius_standardized": 4.0,
            "word_samples_upper": 2, "finite_horizon_words": 1,
            "funnel_level_a": 10.0, "W0_upper": 1.0,
            "L_X_upper": 0.1, "G_bar_upper": 0.0, "c_zw_upper": 0.0,
            "r_star_upper": 1.0, "c_ww_upper": 0.0, "g_W_upper": 0.0, "h_W_upper": 0.0,
        }
        out = tool.validate_stochastic(payload, modes, source_noise)
        self.assertFalse(out["pass"])
        self.assertTrue(any("Gaussian RMS" in x for x in out["failures"]))


if __name__ == "__main__":
    unittest.main()
