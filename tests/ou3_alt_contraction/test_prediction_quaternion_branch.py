import unittest
from tools.stability.ou3_alt_contraction import prediction_quaternion_branch as B


class PredictionQuaternionBranchTests(unittest.TestCase):
    def test_regional_domain_excludes_trig_branch(self):
        d=B.build();self.assertEqual(B.validate(d),[])
        self.assertLess(d['nominal_step_angle_norm_upper_rad'],B.THRESHOLD_RAD)
        self.assertLess(d['shadow_step_angle_norm_upper_rad'],B.THRESHOLD_RAD)
        self.assertGreater(d['shadow_margin_to_threshold_rad'],0.006)
        self.assertFalse(d['trigonometric_prediction_branch_reachable_in_regional_domain'])
        self.assertFalse(d['binary32_polynomial_and_normalization_roundoff_closed'])
        self.assertFalse(d['ALT_LIVE_PASS'])

    def test_polynomial_parameters_match_branch_origin(self):
        w,k=B.small_branch_parameters(0.0)
        self.assertEqual(w,1.0);self.assertEqual(k,0.5)

    def test_threshold_is_strict(self):
        with self.assertRaises(ValueError):B.small_branch_parameters(B.THRESHOLD_RAD**2)
        with self.assertRaises(ValueError):B.small_branch_parameters(-1.0)


if __name__=='__main__': unittest.main()
