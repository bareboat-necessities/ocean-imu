import json
import sys
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_explicit_information_word_certificate as CERT


class Ou3ExplicitInformationWordCertificateTests(unittest.TestCase):
    def test_source_uniform_H_A_information_margins_are_word_endpoint_certified(self):
        d = CERT.build()
        diagnostics = {}
        for mode in ("H", "A"):
            row = d["modes"][mode]
            c = row["matrix_comparison"]
            g = c.get("generalized_matrix_inequality", {})
            diagnostics[mode] = {
                "delta": row["word_endpoint_relative_Riccati_injection_margin_lower"],
                "gate": row["useful_margin_gate"],
                "tau_s": c["tau_s"],
                "sigma_aw_mps2": c["sigma_aw_mps2"],
                "R_S_filter_std": c["R_S_filter_std"],
                "direct_translation_delta": c.get("direct_translation_generalized_margin_lower"),
                "direct_nontranslation_delta": c.get("direct_nontranslation_margin_lower"),
                "limiting_block": g.get("limiting_block"),
            }
        print("P3_MATRIX_MARGINS=" + json.dumps(diagnostics, sort_keys=True), flush=True)

        self.assertEqual(CERT.validate(d), [])
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertTrue(d["validated_arithmetic"])
        self.assertTrue(d["outward_rounded"])
        self.assertEqual(d["p3_backend"], "SOURCE_CELL_GENERALIZED_MATRIX_COMPARISON")
        self.assertEqual(d["generalized_matrix_backend"], "DIRECT_RL_INVERSE_CONGRUENCE_INTERVAL_LDLT")
        self.assertEqual(d["p3_window_backend"], "SOURCE_COMPLETE_WORD_ENDPOINT_GENERALIZED_INFORMATION")
        self.assertFalse(d["old_scalar_min_Q_route_used"])
        self.assertFalse(d["old_scalar_generalized_ratio_used"])
        self.assertEqual(d["continuous_linear_information_certificate"], "PASS")
        self.assertFalse(d["nonlinear_word_enclosed"])
        self.assertEqual(d["theorem_promotion"], "LINEAR_ONLY")

        requirement = d["theorem_margin_requirement"]
        self.assertEqual(requirement["predicate"], "> 0")
        self.assertEqual(requirement["numeric_boundary"], 0.0)
        self.assertFalse(requirement["old_fixed_1e_minus_18_gate_is_theorem_requirement"])
        self.assertEqual(requirement["numerical_search_seed_only"], 1.0e-18)

        binding = d["source_word_binding"]
        self.assertTrue(binding["source_complete_relative_to_declared_theorem_hypotheses"])
        self.assertTrue(binding["arbitrary_source_branches_between_required_PE_events_remain_admissible"])
        self.assertTrue(binding["joint_source_reachability_required"])
        self.assertFalse(binding["one_sample_decrease_required"])
        self.assertEqual(binding["translation_primary_route"], "FOUR_S_SPREAD_COMPLETE_V_P_S_AW_UCO")
        self.assertEqual(binding["translation_aligned_firing_count"], 4)
        self.assertGreaterEqual(binding["translation_spread_index_q_W"], 1)
        self.assertGreater(binding["translation_spread_selected_spacing_lower_s"], 0.0)
        self.assertGreaterEqual(binding["translation_det_spacing_widening_factor_vs_adjacent"], 1)

        for mode in ("H", "A"):
            row = d["modes"][mode]
            self.assertGreater(row["Sigma_lambda_min_lower"], 0.0)
            self.assertGreater(row["Sigma_lambda_max_upper"], row["Sigma_lambda_min_lower"])
            self.assertGreater(row["word_noise_Omega_lambda_min_lower"], 0.0)
            self.assertGreater(row["word_endpoint_relative_Riccati_injection_margin_lower"], 0.0)
            self.assertEqual(
                row["relative_Riccati_injection_margin_lower"],
                row["word_endpoint_relative_Riccati_injection_margin_lower"],
            )
            self.assertEqual(row["useful_margin_gate"], 0.0)
            self.assertEqual(row["useful_margin_gate_predicate"], "> 0")
            self.assertTrue(row["useful_margin_pass"])
            self.assertEqual(row["prefix_information_gain_upper"], 1.0)

            matrix = row["matrix_comparison"]
            self.assertGreater(matrix["direct_translation_generalized_margin_lower"], 0.0)
            g = matrix["generalized_matrix_inequality"]
            self.assertTrue(g["validated_interval_ldlt"])
            self.assertTrue(g["reported_delta_recertified"])
            self.assertFalse(g["old_scalar_rho_over_max_scaled_upper_used"])
            self.assertEqual(g["full_mode_dimension"], row["dimension"])

            arg = matrix["word_endpoint_information_argument"]
            self.assertFalse(arg["repeated_one_step_contraction_used"])
            self.assertFalse(arg["source_replay_used"])
            self.assertTrue(arg["four_S_spread_translation_route"])
            self.assertEqual(arg["coupled_translation_block"], "[v,p,S,a_w]")

            metric = matrix["state_information_metric"]
            self.assertEqual(metric["kind"], "COMPUTATIONAL_CONGRUENCE_FOR_GENERALIZED_MATRIX_INEQUALITY")
            self.assertEqual(len(metric["D_diagonal_squared"]), row["dimension"])
            self.assertTrue(metric["same_congruence_applied_to_noise_and_covariance"])
            self.assertFalse(metric["raw_Euclidean_eigenvalue_gate_used"])
            self.assertFalse(metric["is_nonlinear_Lyapunov_metric"])
            self.assertEqual(metric["translation_coupling_retained"], "full 4x4 [v,p,S,a_w] block")

    def test_no_repeated_step_contraction_shortcut_remains(self):
        text = (ROOT / "tools" / "ou3_explicit_information_word_certificate.py").read_text()
        self.assertNotIn("_window_margin_lower", text)
        self.assertNotIn("1-(1-delta_1)^N", text)
        self.assertNotIn("complete_steps_lower_used_for_information_composition", text)
        self.assertIn("repeated_one_step_contraction_used", text)

    def test_exact_RL_inverse_process_congruence_is_bound_to_certificate(self):
        d = CERT.build()
        p = d["process_congruence_preconditioner"]
        self.assertEqual(p["form"], "C=R*L_inverse")
        self.assertTrue(p["exact_rational"])
        self.assertEqual(p["R_diagonal"], ["1", "10", "100", "2"])
        self.assertEqual(p["limiting_transformed_diagonal"], ["2/3", "5/8", "200/567", "1/2"])
        self.assertEqual(p["norm_2_squared_upper_exact"], 37310)
        self.assertEqual(p["certified_object"], "C*(Q_scaled/x)*C^T")

    def test_pseudo_cadence_is_coupled_to_tau_inside_worst_cells(self):
        d = CERT.build()
        ratio = d["source_schedule"]["pseudo_ratio"]
        lo_guard = d["source_schedule"]["pseudo_min_s"]
        hi_guard = d["source_schedule"]["pseudo_max_s"]
        for mode in ("H", "A"):
            c = d["modes"][mode]["matrix_comparison"]
            tau_lo, tau_hi = c["tau_s"]
            expected_lo = min(max(ratio * tau_lo, lo_guard), hi_guard)
            expected_hi = min(max(ratio * tau_hi, lo_guard), hi_guard)
            cadence_lo, cadence_hi = c["cadence_s"]
            self.assertLessEqual(cadence_lo, expected_lo)
            self.assertGreaterEqual(cadence_hi, expected_hi)

    def test_certificate_depends_on_declared_physical_upper_bounds(self):
        base = json.loads(CERT.DEFAULT_DOMAIN.read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "domain.json"
            broken = json.loads(json.dumps(base))
            broken["normal_live"].pop("specific_force_norm_upper_mps2")
            p.write_text(json.dumps(broken), encoding="utf-8")
            with self.assertRaises(KeyError):
                CERT.build(p)

    def test_compatibility_entry_point_contains_no_retired_scalar_floor(self):
        text = (ROOT / "tools" / "ou3_explicit_information_word_certificate.py").read_text()
        self.assertNotIn("prediction_Q_lambda_min_lower", text)
        self.assertNotIn("posterior_floor(", text)
        self.assertIn("ou3_source_reachable_matrix_p3_direct", text)
        self.assertIn("THEOREM_REQUIRED_MARGIN_PREDICATE = \"> 0\"", text)


if __name__ == "__main__":
    unittest.main()
