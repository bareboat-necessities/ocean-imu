import json
import sys
from pathlib import Path
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_explicit_information_word_certificate as CERT


class Ou3ExplicitInformationWordCertificateTests(unittest.TestCase):
    def test_source_uniform_H_A_information_margins_are_matrix_certified(self):
        d = CERT.build()
        diagnostics = {}
        for mode in ("H", "A"):
            row = d["modes"][mode]
            c = row["matrix_comparison"]
            diagnostics[mode] = {
                "delta": row["relative_Riccati_injection_margin_lower"],
                "gate": row["useful_margin_gate"],
                "tau_s": c["tau_s"],
                "sigma_aw_mps2": c["sigma_aw_mps2"],
                "R_S_filter_std": c["R_S_filter_std"],
                "cadence_s": c["cadence_s"],
                "x_h_over_tau": c["x_h_over_tau"],
                "process_scaled_lambda_min_lower": c["process_scaled_lambda_min_lower"],
                "post_measurement_scaled_Omega_lambda_min_lower": c["post_measurement_scaled_Omega_lambda_min_lower"],
                "Sigma_scaled_lambda_max_upper": c["Sigma_scaled_lambda_max_upper"],
                "Sigma_lambda_max_upper": row["Sigma_lambda_max_upper"],
            }
        print("P3_MATRIX_MARGINS=" + json.dumps(diagnostics, sort_keys=True), flush=True)
        self.assertEqual(CERT.validate(d), [])
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertTrue(d["validated_arithmetic"])
        self.assertTrue(d["outward_rounded"])
        self.assertEqual(d["p3_backend"], "SOURCE_CELL_GENERALIZED_MATRIX_COMPARISON")
        self.assertFalse(d["old_scalar_min_Q_route_used"])
        self.assertGreater(d["cell_partition"]["joint_cells"], 0)
        self.assertGreater(d["source_schedule"]["tau_applied_invariant_s"][0], 0.3)
        self.assertEqual(d["continuous_linear_information_certificate"], "PASS")
        self.assertFalse(d["nonlinear_word_enclosed"])
        self.assertEqual(d["theorem_promotion"], "LINEAR_ONLY")
        for mode in ("H", "A"):
            row = d["modes"][mode]
            self.assertGreater(row["Sigma_lambda_min_lower"], 0.0)
            self.assertGreater(row["Sigma_lambda_max_upper"], row["Sigma_lambda_min_lower"])
            self.assertGreater(row["word_noise_Omega_lambda_min_lower"], 0.0)
            self.assertGreaterEqual(row["relative_Riccati_injection_margin_lower"], row["useful_margin_gate"])
            self.assertLess(row["relative_Riccati_injection_margin_lower"], 1.0)
            self.assertTrue(row["useful_margin_pass"])
            self.assertEqual(row["prefix_information_gain_upper"], 1.0)
            matrix = row["matrix_comparison"]
            self.assertEqual(len(matrix["comparison_scale_diagonal_squared"]), row["dimension"])
            self.assertEqual(len(matrix["Sigma_diagonal_upper"]), row["dimension"])
            self.assertGreater(matrix["post_measurement_scaled_Omega_lambda_min_lower"], 0.0)
            self.assertGreater(matrix["Sigma_scaled_lambda_max_upper"], 0.0)

    def test_exact_RL_inverse_process_congruence_is_bound_to_certificate(self):
        d = CERT.build()
        p = d["process_congruence_preconditioner"]
        self.assertEqual(p["form"], "C=R*L_inverse")
        self.assertTrue(p["exact_rational"])
        self.assertEqual(p["R_diagonal"], ["1", "10", "100", "2"])
        self.assertEqual(
            p["L_inverse"],
            [
                ["1", "0", "0", "0"],
                ["-3/8", "1", "0", "0"],
                ["1/15", "-4/9", "1", "0"],
                ["-15/2", "30", "-105/2", "1"],
            ],
        )
        self.assertEqual(p["limiting_transformed_diagonal"], ["2/3", "5/8", "200/567", "1/2"])
        self.assertEqual(p["norm_inf_exact"], 182)
        self.assertEqual(p["norm_1_exact"], 205)
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

    def test_compatibility_entry_point_does_not_contain_retired_scalar_floor(self):
        text = (ROOT / "tools" / "ou3_explicit_information_word_certificate.py").read_text()
        self.assertNotIn("prediction_Q_lambda_min_lower", text)
        self.assertNotIn("posterior_floor(", text)
        self.assertIn("ou3_source_reachable_matrix_p3", text)


if __name__ == "__main__":
    unittest.main()
