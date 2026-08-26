import math
import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2] / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p5_sample1_coupled_tangent_refinement as C


class Ou3P5Sample1CoupledTangentRefinementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = C.build(source_pieces=2, source_cell_index=0, delta_pieces=4)

    def test_coupling_semantics(self):
        self.assertTrue(self.d["same_first_residual_attitude_aw_coupling_used"])
        self.assertTrue(self.d["beta_source_derived_from_certified_first_gain"])
        self.assertTrue(self.d["coupling_remainder_retained"])
        self.assertTrue(self.d["forward_E_formed_before_posterior_multiplication"])
        self.assertTrue(self.d["sample1_innovation_reconstructed_from_same_forward_map"])
        self.assertTrue(self.d["sample1_S_identity_subbranch_only"])
        self.assertGreater(self.d["beta_mps2_per_rad"], 0.0)

    def test_finite_diagnostics(self):
        for k in (
            "first_tangent_combined_gain_norm_upper",
            "first_tangent_relation_remainder_norm_upper_mps2",
            "max_E_theta_norm_upper",
            "max_Pj_Et_first6_norm_upper",
            "max_actual_Ctheta_norm_upper",
            "max_sample1_acc_correction_norm_upper_rad",
        ):
            self.assertTrue(math.isfinite(float(self.d[k])), k)
            self.assertGreaterEqual(float(self.d[k]), 0.0)
        self.assertEqual(
            self.d["fixed_pivot_inverse_count"] + self.d["spectral_fallback_inverse_count"],
            self.d["evaluated_delta_cells"],
        )

    def test_no_relaxation_or_promotion(self):
        self.assertFalse(self.d["filter_changed"])
        self.assertEqual(self.d["deployed_correction_limit_rad"], 6.0)
        self.assertFalse(self.d["deployed_correction_limit_increased"])
        self.assertFalse(self.d["complete_sample1_branch_refined_here"])
        self.assertFalse(self.d["whole_word_promoted_here"])
        self.assertFalse(self.d["N_H_words_set_here"])

    def test_validates_fail_closed(self):
        self.assertEqual(C.validate(self.d), [])
        status = self.d["P5_SAMPLE1_COUPLED_TANGENT_WITNESS_REFINEMENT"]
        self.assertIn(status, ("PASS", "NOT_ESTABLISHED"))
        if status == "PASS":
            self.assertIsNone(self.d["first_unclosed_delta_cell"])
        else:
            self.assertIsNotNone(self.d["first_unclosed_delta_cell"])


if __name__ == "__main__":
    unittest.main()
