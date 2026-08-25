import math
import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p5_sample1_rotation_gauge_refinement as G


class Ou3P5Sample1RotationGaugeRefinementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = G.build(source_pieces=2, source_cell_index=0)

    def test_refinement_validates_and_is_fail_closed(self):
        d = self.d
        self.assertEqual(G.validate(d), [])
        self.assertGreater(d["evaluated_branch_count"], 0)
        self.assertIn(
            d["P5_SAMPLE1_ROTATION_GAUGE_WITNESS_REFINEMENT"],
            ("PASS", "NOT_ESTABLISHED"),
        )
        if d["P5_SAMPLE1_ROTATION_GAUGE_WITNESS_REFINEMENT"] == "PASS":
            self.assertIsNone(d["first_unclosed_branch"])
        else:
            self.assertIsNotNone(d["first_unclosed_branch"])

    def test_exact_rotation_structure_is_retained(self):
        d = self.d
        self.assertTrue(d["sample0_gravity_SO2_equivariance_used"])
        self.assertTrue(d["sample0_accel_attitude_correction_gravity_axis_component_exact_zero"])
        self.assertTrue(d["sample0_accepted_correction_canonicalized_to_positive_tangent_axis"])
        self.assertTrue(d["world_linear_groups_congruence_transformed_into_corrected_body_gauge"])
        self.assertTrue(d["next_prediction_body_rotation_applied_to_world_linear_gauge"])
        self.assertTrue(d["sample1_J_aw_exact_identity_in_transported_body_gauge"])
        self.assertTrue(d["safe_ldlt_solver_identity_branches_kept_separate"])

    def test_numerical_outputs_are_finite_or_witnessed(self):
        d = self.d
        self.assertTrue(math.isfinite(d["max_sample1_residual_norm_upper_mps2"]))
        self.assertTrue(math.isfinite(d["max_sample1_correction_norm_upper_rad"]))
        self.assertGreaterEqual(d["fixed_pivot_inverse_count"], 0)
        self.assertGreaterEqual(d["spectral_fallback_inverse_count"], 0)
        for row in d["branches"]:
            if "exception" not in row:
                self.assertIn(
                    row["sample1_inverse_backend"],
                    (
                        "FIXED_PIVOT_INTERVAL_GAUSS_JORDAN",
                        "SPD_S_GE_R_SPECTRAL_ENTRY_ENCLOSURE",
                    ),
                )
                self.assertGreaterEqual(row["sample1_correction_norm_upper_rad"], 0.0)

    def test_no_gate_relaxation_or_premature_promotion(self):
        d = self.d
        self.assertEqual(d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(d["deployed_correction_limit_increased"])
        self.assertFalse(d["complete_source_family_refined_here"])
        self.assertFalse(d["whole_word_promoted_here"])
        self.assertFalse(d["N_H_words_set_here"])


if __name__ == "__main__":
    unittest.main()
