from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_radius_continuation_certificate as P4R
import ou3_validate_enclosure as ENC


class Ou3P4RadiusContinuationCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = P4R.build()

    def test_radius_continuation_passes_and_is_monotone_vs_legacy(self):
        d = self.d
        self.assertEqual(P4R.validate(d), [])
        self.assertEqual(d["P4_RADIUS_CONTINUATION_WORD_CERTIFICATE"], "PASS")
        self.assertTrue(d["radius_continuation_source_only"])
        self.assertFalse(d["radius_continuation_trajectory_sampling_used"])
        self.assertFalse(d["radius_continuation_changes_filter"])
        self.assertFalse(d["radius_continuation_changes_P3_margin"])
        for mode in ("H", "A"):
            m = d["modes"][mode]
            self.assertTrue(m["radius_continuation_used"])
            self.assertTrue(m["radius_exact_prefix_bootstrap_used"])
            self.assertTrue(m["radius_exact_endpoint_budget_used"])
            self.assertGreaterEqual(m["certified_level_W"], m["certified_level_W_legacy"])
            self.assertGreaterEqual(m["radius_W_widening_factor_lower"], 1.0)
            self.assertGreaterEqual(m["radius_sqrtW_widening_factor_lower"], 1.0)
            self.assertGreater(m["radius_candidate_count_certified"], 0)
            self.assertGreater(m["radius_selected_design_norm"], 0.0)
            self.assertLess(m["prefix_canonical_error_norm_upper"], m["cayley_norm_limit"])
            self.assertLess(m["accepted_correction_norm_prefix_upper"], 1.0e-2)

    def test_selected_candidate_is_present_and_maximal_on_certified_grid(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            rows = m["radius_continuation_candidates"]
            self.assertEqual(len(rows), m["radius_candidate_count_certified"])
            best = max(r["W"] for r in rows)
            self.assertEqual(m["certified_level_W"], best)
            selected = [r for r in rows if r["q_design"] == m["radius_selected_design_norm"]]
            self.assertTrue(selected)
            self.assertEqual(selected[0]["W"], best)

    def test_exact_prefix_factor_is_used_not_fixed_factor_two(self):
        for mode in ("H", "A"):
            m = self.d["modes"][mode]
            self.assertGreaterEqual(m["radius_prefix_factor_upper"], 1.0)
            self.assertLessEqual(m["radius_prefix_factor_upper"], 2.0)
            self.assertAlmostEqual(
                m["prefix_W_factor_upper"],
                m["radius_prefix_factor_upper"] ** 2,
                delta=8.0 * math.ulp(max(1.0, m["prefix_W_factor_upper"])),
            )

    def test_schema4_validator_accepts_radius_continuation_rows(self):
        for mode in ("H", "A"):
            out = ENC.validate_mode(
                mode,
                self.d["modes"][mode],
                {"required_path_metric": "CAYLEY_LIFTED_SOURCE_INFORMATION_METRIC"},
            )
            self.assertTrue(out["linear_pass"], out["failures"])
            self.assertTrue(out["nonlinear_pass"], out["failures"])
            self.assertGreater(out["mu_W_lower"], 0.0)

    def test_A_projection_remains_exact_interior_identity(self):
        p = self.d["modes"]["A"]["active_bias_projection"]
        self.assertFalse(p["projection_surface_reached_in_certified_funnel"])
        self.assertEqual(p["exact_projection_branch_in_certified_funnel"], "identity_interior_branch")
        self.assertLess(p["certified_error_norm_prefix_upper"], p["interior_margin_lower_mps2"])

    def test_no_replay_numpy_or_sampled_state_search(self):
        text = (ROOT / "tools" / "ou3_p4_radius_continuation_certificate.py").read_text(
            encoding="utf-8"
        ).casefold()
        self.assertNotIn("numpy", text)
        self.assertNotIn("monte carlo", text)
        self.assertNotIn("sampled trajectory", text)
        self.assertIn("radius-continuation", text)


if __name__ == "__main__":
    unittest.main()
