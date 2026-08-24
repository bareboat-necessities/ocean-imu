import importlib.util
from fractions import Fraction
import math
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
        "required_path_metric": "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC",
        "endpoint_metric_source_correlation_required": True,
        "full_attitude_linear_cross_terms_retained": True,
    }


def path_metric():
    sigma_hi = 120.0
    sigma_lo = 0.005
    scale = sigma_hi
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


def mode_payload(level=10.0, decrease=0.2):
    return {
        "source_complete": True,
        "outward_rounded": True,
        "joint_source_reachability": True,
        "one_sample_decrease_used": False,
        "source_replay_used": False,
        "word_horizon_s": 1.0,
        "word_endpoint_relative_Riccati_injection_margin_lower": 0.04,
        "Sigma_lambda_min_lower": 0.005,
        "Sigma_lambda_max_upper": 120.0,
        "prefix_information_gain_upper": 1.0,
        "path_metric": path_metric(),
        "theta_star": 1.0,
        "endpoint_relative_W_decrease_lower": decrease,
        "mu_W_lower": decrease,
        "certified_level_W": level,
        "all_word_prefixes_safe": True,
        "accepted_correction_uses_source_series_branch": True,
        "prefix_canonical_error_norm_upper": 1.0e-4,
        "cayley_norm_limit": 1.0,
        "accepted_correction_norm_prefix_upper": 1.0e-5,
    }


def exact(x: float) -> Fraction:
    return Fraction.from_float(float(x))


def assert_contains_exact(test: unittest.TestCase, interval: Interval, value: Fraction):
    test.assertLessEqual(exact(interval.lo), value)
    test.assertLessEqual(value, exact(interval.hi))


class Ou3IntervalArithmeticTests(unittest.TestCase):
    def test_basic_operations_enclose_exact_binary64_endpoint_arithmetic(self):
        rng = random.Random(403)
        for _ in range(100):
            a0, a1 = sorted((rng.uniform(-20,20), rng.uniform(-20,20)))
            b0, b1 = sorted((rng.uniform(-20,20), rng.uniform(-20,20)))
            A, B = Interval(a0,a1), Interval(b0,b1)
            assert_contains_exact(self, A+B, exact(a0)+exact(b0))
            assert_contains_exact(self, A+B, exact(a1)+exact(b1))
            products=[exact(x)*exact(y) for x in (a0,a1) for y in (b0,b1)]
            C=A*B; assert_contains_exact(self,C,min(products)); assert_contains_exact(self,C,max(products))

    def test_division_encloses_exact_values_when_denominator_avoids_zero(self):
        A=Interval(-2.0,3.0); B=Interval(0.5,2.0); C=A/B
        for x in (-2.0,3.0):
            for y in (0.5,2.0): assert_contains_exact(self,C,exact(x)/exact(y))

    def test_square_handles_zero_crossing(self):
        C=Interval(-3.0,2.0).square(); self.assertEqual(C.lo,0.0); assert_contains_exact(self,C,Fraction(9,1))

    def test_interval_matrix_product_contains_exact_product(self):
        A=[[Interval.outward_bounds(0.9,1.1),Interval.outward_bounds(-0.2,-0.1)],
           [Interval.outward_bounds(0.3,0.4),Interval.outward_bounds(1.9,2.1)]]
        B=IA.matrix_point([[2.0,-1.0],[0.5,3.0]]); C=IA.matrix_mul(A,B)
        self.assertTrue(C[0][0].contains(2*1.0+0.5*(-0.15)))

    def test_gershgorin_certifies_and_refuses_as_expected(self):
        A=[[Interval.outward_bounds(3.9,4.1),Interval.outward_bounds(-0.2,0.2)],
           [Interval.outward_bounds(-0.2,0.2),Interval.outward_bounds(2.9,3.1)]]
        self.assertTrue(IA.symmetric_positive_definite_gershgorin(A)[0])
        B=[[Interval.outward_bounds(0.9,1.1),Interval.outward_bounds(-2,2)],
           [Interval.outward_bounds(-2,2),Interval.outward_bounds(0.9,1.1)]]
        self.assertFalse(IA.symmetric_positive_definite_gershgorin(B)[0])

    def test_source_box_is_outward_and_nonpromoting(self):
        source=SOURCE.build(SOURCE.DEFAULT_HEADER.resolve()); box=SOURCE_BOX.build(SOURCE.DEFAULT_HEADER.resolve())
        self.assertEqual(SOURCE_BOX.validate(box,SOURCE.DEFAULT_HEADER.resolve()),[])
        self.assertEqual(box["theorem_promotion"],"NOT_ESTABLISHED")
        for name,endpoints in source["continuous_parameters"].items():
            I=Interval(*box["continuous_parameters"][name]); self.assertTrue(I.contains(endpoints[0])); self.assertTrue(I.contains(endpoints[1]))


class Ou3ValidatedEnclosureTests(unittest.TestCase):
    def test_cayley_mode_accepts_direct_positive_mu(self):
        out=load_tool().validate_mode("H",mode_payload(),contract_mode())
        self.assertTrue(out["pass"],out["failures"])
        self.assertEqual(out["path_metric"]["kind"],"CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC")
        self.assertTrue(out["path_metric"]["full_attitude_linear_cross_terms_retained"])
        self.assertAlmostEqual(out["endpoint_W_ratio_upper"],0.8)

    def test_tiny_positive_gap_is_not_erased_by_one_minus(self):
        p=mode_payload(decrease=1e-30); p["mu_W_lower"]=math.nextafter(1e-30,-math.inf)
        out=load_tool().validate_mode("H",p,contract_mode())
        self.assertTrue(out["pass"],out["failures"]); self.assertIsNone(out["endpoint_W_ratio_upper"]); self.assertGreater(out["mu_W_lower"],0)

    def test_retired_block_metric_and_cross_term_loss_are_rejected(self):
        tool=load_tool(); p=mode_payload(); p["path_metric"]["kind"]="GROUP_COMPATIBLE_NODE_METRIC"; p["path_metric"]["block_diagonal_metric_used"]=True
        self.assertFalse(tool.validate_mode("H",p,contract_mode())["nonlinear_pass"])
        p=mode_payload(); p["path_metric"]["full_attitude_linear_cross_terms_retained"]=False
        self.assertFalse(tool.validate_mode("H",p,contract_mode())["nonlinear_pass"])

    def test_old_injection_alias_is_not_schema4_fallback(self):
        p=mode_payload(); p.pop("word_endpoint_relative_Riccati_injection_margin_lower"); p["relative_Riccati_injection_margin_lower"]=0.04
        self.assertFalse(load_tool().validate_mode("H",p,contract_mode())["linear_pass"])

    def test_sampled_counterexample_can_only_falsify_level(self):
        out=load_tool().validate_mode("H",mode_payload(level=12),contract_mode(),{"first_fail_W":11})
        self.assertFalse(out["pass"]); self.assertTrue(any("sampled failure" in x for x in out["failures"]))

    def test_hybrid_widening_semantics_are_required(self):
        tool=load_tool(); modes={"H":{"certified_level_W":10},"A":{"certified_level_W":10}}
        common=dict(source_complete=True,outward_rounded=True,source_level_W_upper=2.0,jump_gain_upper=0.5,additive_W_upper=0.2,destination_level_W=5.0,destination_mode="H")
        rows=[{"kind":"startup_handoff",**common},
              {"kind":"held_to_active",**{**common,"destination_mode":"A"},"source_dimension":18,"destination_dimension":21,"dimension_change_handled_by_embedding":True,"new_coordinate_W_upper":0.3},
              {"kind":"magnetic_lock",**common},{"kind":"magnetic_regauge_refinement",**common},
              {"kind":"tilt_reset",**common,"discarded_pre_reset_tilt_excluded_from_multiplicative_gain":True,"reset_to_funnel_exact_map":True},
              {"kind":"tilt_relock",**common},{"kind":"cooldown_reentry",**common,"reachable_word_product_used":True,"global_worst_word_power_used":False}]
        out=tool.validate_hybrid(rows,modes); self.assertTrue(out["pass"],out["failures"]); self.assertEqual(out["analytic_obligation_separate"],"periodic_aw_covariance_sync")

    @staticmethod
    def stochastic_payload(radius=20.0,level=100.0):
        return {"source_complete":True,"outward_rounded":True,"localization_prefix_safe":True,"gaussian_localization_used":True,"freedman_excursion_used":True,"markov_union_fallback_used":False,"localization_radius_standardized":radius,"word_samples_upper":2,"finite_horizon_words":1,"funnel_level_a":level,"W0_upper":1.0,"L_X_upper":0.1,"G_bar_upper":1e-4,"c_zw_upper":0.0,"r_star_upper":1.0,"c_ww_upper":0.0,"g_W_upper":1e-4,"h_W_upper":1e-4}

    def test_stochastic_probability_uses_gaussian_and_freedman(self):
        tool=load_tool(); modes={"H":{"endpoint_W_ratio_upper":0.5,"certified_level_W":100,"word_horizon_s":0.01},"A":{"endpoint_W_ratio_upper":0.6,"certified_level_W":100,"word_horizon_s":0.01}}
        noise={"schema":1,"source_generated_not_trajectory_fit":True,"standardized_increment":{"dimension":18,"covariance_upper_identity":True},"physical_scales":{"imu_dt_s":0.005}}
        out=tool.validate_stochastic(self.stochastic_payload(),modes,noise); self.assertTrue(out["pass"],out["failures"]); self.assertLess(out["finite_horizon_failure_probability_upper"],1)


if __name__=="__main__":
    unittest.main()
