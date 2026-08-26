import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_structured_full_gain_v4 as V4


class StructuredFullGainV4Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d=V4.build(p_pieces=4,d_pieces=4,axial_pieces=4)

    def test_full_force_and_latent_semantics(self):
        for k in (
            "canonical_first_full_theta_aw_Joseph_marginal_used",
            "shipping_left_error_reset_used_in_theta_marginal",
            "corrected_body_Rx_gauge_used_in_aw_marginal",
            "source_correlated_tangent_aw_mean_relation_used",
            "sample1_nonaxial_force_included",
            "complete_3x3_accelerometer_innovation_used",
            "complete_3x3_attitude_gain_used",
            "full_propagated_aw_error_norm_used",
            "latent_finite_rotation_covered_by_orthogonal_norm_invariance",
            "sample1_body_rotation_removed_by_simultaneous_orthogonal_gauge",
        ):
            self.assertTrue(self.d[k])
        self.assertFalse(self.d["temporal_force_slew_assumed"])

    def test_fail_closed_numeric_result(self):
        self.assertGreater(self.d["evaluated_joint_cells"],0)
        self.assertTrue(math.isfinite(self.d["max_sample1_residual_norm_upper_mps2"]))
        self.assertTrue(math.isfinite(self.d["max_Ktheta_operator_norm_upper"]))
        self.assertTrue(math.isfinite(self.d["max_correction_norm_upper_rad"]))
        status=self.d["P5_SAMPLE1_STRUCTURED_FULL_GAIN_V4"]
        self.assertIn(status,("PASS","NOT_ESTABLISHED"))
        self.assertEqual(self.d["first_unclosed_joint_cell"] is None,status=="PASS")

    def test_no_promotion(self):
        self.assertFalse(self.d["first_attitude_PSD_cross_axis_remainder_included"])
        self.assertFalse(self.d["sample1_S_covariance_update_included"])
        self.assertFalse(self.d["sample1_S_attitude_injection_included"])
        self.assertFalse(self.d["complete_sample1_branch_closed_here"])
        self.assertFalse(self.d["q8_word_promoted_here"])
        self.assertFalse(self.d["whole_word_promoted_here"])
        self.assertFalse(self.d["N_H_words_set_here"])

    def test_validate(self):
        self.assertEqual(V4.validate(self.d),[])


if __name__=="__main__": unittest.main()
