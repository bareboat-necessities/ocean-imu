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
    scale = 120.0
    sigma_lo = 5.0e-7
    return {
        "kind": "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC",
        "chart_coordinate": "c(R)=2*tan(theta/2)*u=4*e_R/(1+tr(R))",
        "chart_domain": "theta<pi",
        "exact_group_metric": "W_g=s_mode*[c(R);xi]^T Sigma_KF(g)^-1 [c(R);xi]",
        "source_covariance_inverse": True,
        "mode_global_positive_scale": scale,
        "same_scale_on_every_source_node_in_mode": True,
        "node_dependent": True,
        "full_attitude_linear_cross_terms_retained": True,
        "block_diagonal_metric_used": False,
        "common_Euclidean_metric_used": False,
        "local_coordinate_matches_P3_delta_theta": True,
        "local_quadratic_is_positive_scalar_multiple_of_P3_information_metric": True,
        "endpoint_metric_must_match_endpoint_source_covariance": True,
        "joint_source_reachability_required": True,
        "metric_lambda_min_lower": 1.0,
        "metric_lambda_max_upper": scale / sigma_lo,
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
        return {"mode":"A","required_path_metric":"CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC",
                "endpoint_metric_source_correlation_required":True,
                "full_attitude_linear_cross_terms_retained":True}

    @staticmethod
    def valid_mode_payload():
        decrease=1e-4
        return {"source_complete":True,"outward_rounded":True,"joint_source_reachability":True,
                "one_sample_decrease_used":False,"source_replay_used":False,"word_horizon_s":1.0,
                "word_endpoint_relative_Riccati_injection_margin_lower":0.005,
                "Sigma_lambda_min_lower":5e-7,"Sigma_lambda_max_upper":120.0,
                "prefix_information_gain_upper":1.0,"path_metric":cayley_metric(),"theta_star":1.0,
                "endpoint_relative_W_decrease_lower":decrease,"mu_W_lower":math.nextafter(decrease,-math.inf),
                "certified_level_W":0.05,"all_word_prefixes_safe":True,
                "accepted_correction_uses_source_series_branch":True,
                "prefix_canonical_error_norm_upper":1e-6,"cayley_norm_limit":1.0,
                "accepted_correction_norm_prefix_upper":1e-7,
                "active_bias_projection":{"projection_surface_reached_in_certified_funnel":False,
                                          "exact_projection_branch_in_certified_funnel":"identity_interior_branch"}}

    def test_information_mode_accepts_normalized_cayley_bounds(self):
        ans=ENC.validate_mode("A",self.valid_mode_payload(),self.contract_mode())
        self.assertTrue(ans["linear_pass"],ans["failures"]); self.assertTrue(ans["nonlinear_pass"],ans["failures"])
        self.assertAlmostEqual(ans["P3_inverse_covariance_conditioning_lambda_min_lower"],1/120)
        self.assertAlmostEqual(ans["P3_inverse_covariance_conditioning_lambda_max_upper"],2e6)
        self.assertGreater(ans["mu_W_lower"],0); self.assertEqual(ans["path_metric"]["kind"],"CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC")
        self.assertTrue(ans["path_metric"]["same_scale_on_every_source_node_in_mode"])

    def test_zero_endpoint_injection_old_alias_and_block_metric_are_rejected(self):
        p=self.valid_mode_payload(); p["word_endpoint_relative_Riccati_injection_margin_lower"]=0.0
        self.assertFalse(ENC.validate_mode("A",p,self.contract_mode())["linear_pass"])
        p=self.valid_mode_payload(); p.pop("word_endpoint_relative_Riccati_injection_margin_lower"); p["relative_Riccati_injection_margin_lower"]=0.005
        self.assertFalse(ENC.validate_mode("A",p,self.contract_mode())["linear_pass"])
        p=self.valid_mode_payload(); p["path_metric"]["kind"]="GROUP_COMPATIBLE_NODE_METRIC"; p["path_metric"]["block_diagonal_metric_used"]=True
        self.assertFalse(ENC.validate_mode("A",p,self.contract_mode())["nonlinear_pass"])

    def test_contract_keeps_executed_horizons_diagnostic_and_requires_cayley(self):
        info={"status":"PASS","held":{"selected":{"horizon_s":2.0,"lambda_worst_information":.9999},"strongest_executed_margin":{"horizon_s":16.0,"lambda_worst_information":.95,"relative_Riccati_injection_margin_worst":.05,"Sigma_endpoint_lambda_min":1e-7,"Sigma_endpoint_lambda_max":10}},"active":{"selected":{"horizon_s":2.0,"lambda_worst_information":.999},"strongest_executed_margin":{"horizon_s":4.0,"lambda_worst_information":.99,"relative_Riccati_injection_margin_worst":.01,"Sigma_endpoint_lambda_min":1e-7,"Sigma_endpoint_lambda_max":10}}}
        completion={"status":"PASS_EXECUTED_REPLAY","held":{"asymptotic_floor_b_star_replay":10},"active":{"asymptotic_floor_b_star_replay":20}}
        c=CONTRACT.build_contract(info,completion); self.assertEqual(c["schema"],4)
        self.assertIn("DIAGNOSTIC_ONLY",c["modes"]["H"]["executed_reference_only"]["qualification"])
        policy=c["metric_policy"]
        self.assertTrue(policy["P4_exact_Cayley_group_lift_required"]); self.assertTrue(policy["P4_endpoint_metric_matches_source_covariance"]); self.assertTrue(policy["P4_full_attitude_linear_cross_terms_retained"])
        self.assertFalse(policy["block_diagonal_group_metric_fallback_allowed"]); self.assertFalse(policy["common_Euclidean_or_common_quadratic_fallback_allowed"])

    def test_completion_does_not_promote_sampled_replay(self):
        text=(TOOLS/"ou3_information_completion.py").read_text(); self.assertIn('"numerical_neighborhood_certificate": "NOT_ESTABLISHED"',text); self.assertIn('"deployment_theorem_certificate": "NOT_ESTABLISHED"',text)

    def test_enclosure_gate_uses_cayley_cross_terms_and_direct_gap(self):
        text=(TOOLS/"ou3_validate_enclosure.py").read_text(); self.assertIn("CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC",text); self.assertIn("same_scale_on_every_source_node_in_mode",text); self.assertIn("full_attitude_linear_cross_terms_retained",text); self.assertIn("endpoint_relative_W_decrease_lower",text); self.assertNotIn("np.load(",text)


if __name__=="__main__":
    unittest.main()
