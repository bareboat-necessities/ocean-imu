import importlib.util
import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
sys.path.insert(0, str(TOOLS))


def load(name: str):
    spec = importlib.util.spec_from_file_location(name, TOOLS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


COMP = load("ou3_certificate_completion")
CONTRACT = load("ou3_information_enclosure_contract")
ENC = load("ou3_validate_enclosure")


def cayley_metric():
    return {
        "kind": "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC",
        "chart_coordinate": "c(R)=2*tan(theta/2)*u=4*e_R/(1+tr(R))",
        "chart_domain": "theta<pi",
        "exact_group_metric": "W_g=[c(R);xi]^T Sigma_KF(g)^-1 [c(R);xi]",
        "source_covariance_inverse": True,
        "node_dependent": True,
        "full_attitude_linear_cross_terms_retained": True,
        "block_diagonal_metric_used": False,
        "common_Euclidean_metric_used": False,
        "local_coordinate_matches_P3_delta_theta": True,
        "local_quadratic_equals_P3_information_metric": True,
        "endpoint_metric_must_match_endpoint_source_covariance": True,
        "joint_source_reachability_required": True,
        "metric_lambda_min_lower": 1.0 / 120.0,
        "metric_lambda_max_upper": 2.0e6,
    }


class Ou3CertificateCompletionTests(unittest.TestCase):
    def test_recurrence_has_finite_capture_when_lambda_lt_one(self):
        r = COMP.recurrence_capture(c0=10.0, lam=0.8, gamma=0.1, b=0.5)
        self.assertEqual(r["status"], "PASS")
        self.assertGreater(r["N_H_words"], 0)

    def test_recurrence_refuses_noncontracting_word_family(self):
        r = COMP.recurrence_capture(c0=10.0, lam=1.01, gamma=0.1, b=1.0)
        self.assertEqual(r["status"], "NOT_AVAILABLE")

    @staticmethod
    def contract_mode():
        return {
            "mode": "A",
            "required_path_metric": "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC",
            "endpoint_metric_source_correlation_required": True,
            "full_attitude_linear_cross_terms_retained": True,
            "executed_reference_only": {
                "relative_Riccati_injection_margin_worst": 0.01,
                "Sigma_endpoint_lambda_min": 1.0e-6,
                "Sigma_endpoint_lambda_max": 100.0,
            },
        }

    @staticmethod
    def valid_mode_payload():
        mmin = 1.0 / 120.0
        decrease = 1.0e-4
        return {
            "source_complete": True,
            "outward_rounded": True,
            "joint_source_reachability": True,
            "one_sample_decrease_used": False,
            "source_replay_used": False,
            "word_horizon_s": 1.0,
            "word_endpoint_relative_Riccati_injection_margin_lower": 0.005,
            "Sigma_lambda_min_lower": 5.0e-7,
            "Sigma_lambda_max_upper": 120.0,
            "prefix_information_gain_upper": 1.0,
            "path_metric": cayley_metric(),
            "theta_star": 1.0,
            "endpoint_relative_W_decrease_lower": decrease,
            "mu_W_lower": math.nextafter(decrease * mmin, -math.inf),
            "certified_level_W": 0.05,
            "all_word_prefixes_safe": True,
            "accepted_correction_uses_source_series_branch": True,
            "prefix_canonical_error_norm_upper": 1.0e-6,
            "cayley_norm_limit": 1.0,
            "accepted_correction_norm_prefix_upper": 1.0e-7,
        }

    def test_information_mode_accepts_strict_cayley_bounds(self):
        ans = ENC.validate_mode("A", self.valid_mode_payload(), self.contract_mode())
        self.assertTrue(ans["linear_pass"], ans["failures"])
        self.assertTrue(ans["nonlinear_pass"], ans["failures"])
        self.assertAlmostEqual(ans["P3_inverse_covariance_conditioning_lambda_min_lower"], 1.0 / 120.0)
        self.assertAlmostEqual(ans["P3_inverse_covariance_conditioning_lambda_max_upper"], 2.0e6)
        self.assertGreater(ans["mu_W_lower"], 0.0)
        self.assertEqual(ans["path_metric"]["kind"], "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC")
        self.assertTrue(ans["path_metric"]["full_attitude_linear_cross_terms_retained"])

    def test_information_mode_rejects_zero_endpoint_injection(self):
        payload = self.valid_mode_payload()
        payload["word_endpoint_relative_Riccati_injection_margin_lower"] = 0.0
        ans = ENC.validate_mode("A", payload, self.contract_mode())
        self.assertFalse(ans["linear_pass"])

    def test_information_mode_rejects_old_field_alias(self):
        payload = self.valid_mode_payload()
        payload.pop("word_endpoint_relative_Riccati_injection_margin_lower")
        payload["relative_Riccati_injection_margin_lower"] = 0.005
        ans = ENC.validate_mode("A", payload, self.contract_mode())
        self.assertFalse(ans["linear_pass"])

    def test_information_mode_rejects_retired_block_metric(self):
        payload = self.valid_mode_payload()
        payload["path_metric"]["kind"] = "GROUP_COMPATIBLE_NODE_METRIC"
        payload["path_metric"]["block_diagonal_metric_used"] = True
        ans = ENC.validate_mode("A", payload, self.contract_mode())
        self.assertFalse(ans["nonlinear_pass"])

    def test_contract_keeps_executed_horizons_diagnostic_only(self):
        info = {
            "status": "PASS",
            "held": {"selected": {"horizon_s": 2.0, "lambda_worst_information": 0.9999},
                     "strongest_executed_margin": {"horizon_s": 16.0, "lambda_worst_information": 0.95,
                        "relative_Riccati_injection_margin_worst": 0.05,
                        "Sigma_endpoint_lambda_min": 1e-7, "Sigma_endpoint_lambda_max": 10.0}},
            "active": {"selected": {"horizon_s": 2.0, "lambda_worst_information": 0.999},
                       "strongest_executed_margin": {"horizon_s": 4.0, "lambda_worst_information": 0.99,
                        "relative_Riccati_injection_margin_worst": 0.01,
                        "Sigma_endpoint_lambda_min": 1e-7, "Sigma_endpoint_lambda_max": 10.0}},
        }
        completion = {"status": "PASS_EXECUTED_REPLAY",
                      "held": {"asymptotic_floor_b_star_replay": 10.0},
                      "active": {"asymptotic_floor_b_star_replay": 20.0}}
        c = CONTRACT.build_contract(info, completion)
        self.assertEqual(c["schema"], 4)
        self.assertEqual(c["modes"]["H"]["executed_reference_only"]["strongest_executed_horizon_s"], 16.0)
        self.assertEqual(c["modes"]["A"]["executed_reference_only"]["strongest_executed_horizon_s"], 4.0)
        self.assertIn("DIAGNOSTIC_ONLY", c["modes"]["H"]["executed_reference_only"]["qualification"])
        policy = c["metric_policy"]
        self.assertTrue(policy["P3_conditioning_coordinate_invariant"])
        self.assertTrue(policy["P4_local_quadratic_equals_P3_information_metric"])
        self.assertTrue(policy["P4_exact_Cayley_group_lift_required"])
        self.assertTrue(policy["P4_endpoint_metric_matches_source_covariance"])
        self.assertTrue(policy["P4_full_attitude_linear_cross_terms_retained"])
        self.assertFalse(policy["block_diagonal_group_metric_fallback_allowed"])
        self.assertFalse(policy["common_Euclidean_or_common_quadratic_fallback_allowed"])

    def test_completion_does_not_promote_sampled_replay(self):
        text = (TOOLS / "ou3_information_completion.py").read_text()
        self.assertIn('"numerical_neighborhood_certificate": "NOT_ESTABLISHED"', text)
        self.assertIn('"deployment_theorem_certificate": "NOT_ESTABLISHED"', text)

    def test_enclosure_gate_uses_endpoint_information_and_cayley_metric(self):
        text = (TOOLS / "ou3_validate_enclosure.py").read_text()
        self.assertIn("word_endpoint_relative_Riccati_injection_margin_lower", text)
        self.assertIn("CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC", text)
        self.assertIn("full_attitude_linear_cross_terms_retained", text)
        self.assertIn("endpoint_relative_W_decrease_lower", text)
        self.assertNotIn("GROUP_COMPATIBLE_NODE_METRIC", text)
        self.assertNotIn("np.load(", text)
        self.assertNotIn("load_metrics(", text)
        self.assertNotIn("robust_box_lmi_upper(", text)


if __name__ == "__main__":
    unittest.main()
