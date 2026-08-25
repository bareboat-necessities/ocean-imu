import copy
import importlib.util
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location("ou3_deployment_gate", ROOT / "tools" / "ou3_deployment_gate.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def hybrid_row(kind):
    row = {
        "kind": kind,
        "destination_mode": "H",
        "source_level_W_upper": 2.0,
        "jump_gain_upper": 0.5,
        "additive_W_upper": 0.2,
        "new_coordinate_W_upper": 0.0,
        "destination_level_W": 5.0,
        "inward_margin_lower": 3.8,
        "pass": True,
    }
    if kind == "tilt_reset":
        row.update({
            "discarded_pre_reset_tilt_excluded_from_multiplicative_gain": True,
            "reset_to_funnel_exact_map": True,
        })
    if kind == "cooldown_reentry":
        row.update({
            "reachable_word_product_used": True,
            "global_worst_word_power_used": False,
        })
    return row


def complete_check():
    return {
        "validated_enclosure_provenance_pass": True,
        "modes": {
            "H": {"linear_pass": True, "nonlinear_pass": True},
            "A": {"linear_pass": True, "nonlinear_pass": True},
        },
        "hybrid": {
            "pass": False,
            "bounds": [
                hybrid_row("startup_handoff"),
                hybrid_row("held_to_active"),
                hybrid_row("magnetic_lock"),
                hybrid_row("magnetic_regauge_refinement"),
                hybrid_row("tilt_reset"),
                hybrid_row("tilt_relock"),
                hybrid_row("cooldown_reentry"),
            ],
        },
    }


def primitive_payload():
    return {
        "stochastic_primitives": {
            "lambda_W_upper": 0.5,
            "word_length_samples_upper": 2,
            "L_X_upper": 0.5,
            "G_bar_upper": 0.001,
            "c_zw_upper": 0.0,
            "r_star_upper": 1.0,
            "c_ww_upper": 0.0,
            "s2_upper": 1e-6,
            "s4_upper": 3e-12,
            "g_W_upper": 0.01,
            "h_W_upper": 0.01,
            "Sigma_trace_upper": 1e-6,
            "Sigma_trace_square_upper": 1e-12,
            "Sigma_norm_upper": 1e-6,
            "localization_radius_lower": 0.1,
            "word_horizon_count": 10,
            "b_W_upper": 1e-3,
            "v_W_upper": 1e-8,
            "funnel_level_a_lower": 0.1,
            "initial_W_upper": 0.01,
            "failure_probability_budget": 0.1,
        },
        "capture_primitives": {
            mode: {
                "lambda_upper": 0.5,
                "gamma_upper": 0.25,
                "initial_level_upper": 8.0,
                "strict_superlevel_factor": 1.001,
                "strict_superlevel_absolute": 1e-12,
                "word_horizon_s_upper": 16.0,
            }
            for mode in ("H", "A")
        },
    }


class DeploymentGateTests(unittest.TestCase):
    def test_gaussian_threshold_is_recomputed(self):
        t = mod.t_star_for_radius(100.0, 1.0, 1.0, 1.0)
        self.assertGreater(t, 0.0)
        lhs = 1.0 + 2.0 * math.sqrt(t) + 2.0 * t
        self.assertLessEqual(lhs, 100.0 + 1e-9)

    def test_capture_uses_strict_superlevel(self):
        c = mod.derive_capture({
            "lambda_upper": 0.5, "gamma_upper": 0.25,
            "initial_level_upper": 8.0, "strict_superlevel_factor": 1.001,
            "strict_superlevel_absolute": 1e-12, "word_horizon_s_upper": 16.0,
        })
        self.assertTrue(c["pass"])
        self.assertGreater(c["strict_capture_level_b_eta"], c["asymptotic_level_b_star_upper"])
        self.assertGreater(c["capture_words_upper"], 0)
        self.assertEqual(c["capture_time_s_upper"], 16.0 * c["capture_words_upper"])

    def test_stochastic_probability_is_derived_from_primitives(self):
        s = mod.derive_stochastic(primitive_payload()["stochastic_primitives"])
        self.assertTrue(s["pass"])
        self.assertLessEqual(s["finite_horizon_failure_probability_upper"], 0.1)
        self.assertLess(s["finite_horizon_failure_probability_upper"], 1.0)
        self.assertGreater(s["gaussian_t_star_lower"], 0.0)

    def test_source_domain_is_regenerated_and_compared_to_current_source(self):
        expected = mod.SOURCE_DOMAIN.build(mod.SOURCE_DOMAIN.DEFAULT_HEADER.resolve())
        out = mod.validate_source_domain(expected)
        self.assertTrue(out["pass"], out["failures"])
        self.assertIn("configured_runtime_assumption", out)
        self.assertFalse(out["configured_runtime_api_enforced"])
        stale = copy.deepcopy(expected)
        stale["continuous_parameters"].pop("tau_aw_s")
        out = mod.validate_source_domain(stale)
        self.assertFalse(out["pass"])

    def test_source_domain_rejects_stale_runtime_scope(self):
        expected = mod.SOURCE_DOMAIN.build(mod.SOURCE_DOMAIN.DEFAULT_HEADER.resolve())
        stale = copy.deepcopy(expected)
        stale["configured_runtime_assumption"]["imu_dt_s"] *= 2.0
        self.assertFalse(mod.validate_source_domain(stale)["pass"])

    def test_required_hybrid_set_is_complete(self):
        self.assertEqual(mod.REQUIRED_HYBRID, {
            "startup_handoff", "held_to_active", "magnetic_lock",
            "magnetic_regauge_refinement", "tilt_reset", "tilt_relock",
            "cooldown_reentry", "periodic_aw_covariance_sync",
        })

    def test_final_gate_recomputes_current_hybrid_closure(self):
        source = mod.SOURCE_DOMAIN.build(mod.SOURCE_DOMAIN.DEFAULT_HEADER.resolve())
        out = mod.compose(complete_check(), source, primitive_payload())
        self.assertTrue(out["hybrid_pass"], out["hybrid"]["failures"])
        self.assertEqual(out["hybrid"]["missing"], [])
        self.assertFalse(out["hybrid"]["legacy_name_aliases_used"])
        self.assertTrue(out["hybrid"]["periodic_aw_covariance_sync"]["pass"])
        self.assertEqual(out["deployment_theorem_certificate"], "PASS")

    def test_final_gate_fails_if_current_hybrid_obligation_is_missing(self):
        check = complete_check()
        check["hybrid"]["bounds"] = [x for x in check["hybrid"]["bounds"] if x["kind"] != "magnetic_lock"]
        source = mod.SOURCE_DOMAIN.build(mod.SOURCE_DOMAIN.DEFAULT_HEADER.resolve())
        out = mod.compose(check, source, primitive_payload())
        self.assertFalse(out["hybrid_pass"])
        self.assertEqual(out["deployment_theorem_certificate"], "FAIL")
        self.assertIn("magnetic_lock", out["hybrid"]["missing"])


if __name__ == "__main__":
    unittest.main()
