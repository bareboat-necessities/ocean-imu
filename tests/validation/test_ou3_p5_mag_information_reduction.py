import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_mag_information_reduction as M


class Ou3P5MagInformationReductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = M.build()

    def test_reduction_passes_without_promoting_p5(self):
        d = self.d
        self.assertEqual(M.validate(d), [])
        self.assertEqual(d["P5_MAGNETOMETER_INFORMATION_REDUCTION_CERTIFICATE"], "PASS")
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertFalse(d["source_replay_used"])
        self.assertFalse(d["filter_changed"])
        self.assertFalse(d["complete_word_numerical_certificate_closed_here"])

    def test_radial_joseph_information_cancels_exactly(self):
        d = self.d
        self.assertTrue(d["configured_R_isotropic"])
        self.assertTrue(d["radial_subspace_invariant_under_S_and_R"])
        self.assertTrue(d["radial_residual_equals_radial_eta"])
        self.assertTrue(d["radial_Joseph_information_cancels_exactly"])
        self.assertFalse(d["standalone_radial_eta_information_budget_used"])
        self.assertFalse(d["standalone_full_vector_eta_information_budget_used"])
        self.assertIn("H d_eff", d["exact_reduced_Joseph_identity"])
        self.assertIn("d_eff-c_perp", d["exact_reduced_Joseph_identity"])

    def test_finite_cayley_cells_have_strict_tangent_penalty_below_one(self):
        d = self.d
        cells = d["annular_information_cells"]
        self.assertEqual(len(cells), d["subdivision_cell_count"])
        self.assertGreater(len(cells), 1)
        last = -1.0
        for row in cells:
            self.assertTrue(row["radial_Joseph_information_contribution_exact_zero"])
            self.assertTrue(row["strict_tangent_penalty_ratio_below_one"])
            p = row["effective_vs_linear_tangent_penalty_information_ratio_upper"]
            self.assertGreaterEqual(p, 0.0)
            self.assertLess(p, 1.0)
            self.assertGreaterEqual(p, last)
            last = p
            self.assertGreater(row["effective_tangent_gain_lower"], 0.0)
            self.assertLessEqual(row["effective_tangent_gain_lower"], row["effective_tangent_gain_upper"])
            self.assertLessEqual(row["effective_tangent_gain_upper"], 1.0)
        self.assertLess(d["widened_tangent_penalty_ratio_upper"], 1.0)
        self.assertGreater(d["widened_effective_tangent_gain_lower"], 0.0)

    def test_useful_information_still_requires_source_correlated_interval_cell(self):
        d = self.d
        self.assertTrue(d["useful_S_inverse_term_still_requires_source_correlated_interval_cell"])
        self.assertIn("P,H,R,K,S", d["next_obligation"])


if __name__ == "__main__":
    unittest.main()
