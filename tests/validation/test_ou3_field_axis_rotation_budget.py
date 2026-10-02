import unittest
from fractions import Fraction as F
from tools.stability.ou3_theorem.field_axis_rotation_budget import certificate

class T(unittest.TestCase):
    def test_raw_angle_is_not_promoted(self):
        r=certificate()
        self.assertEqual(F(r["fast_integral_rad"]),F("0.002"))
        self.assertTrue(r["slow_must_use_rate_paired_Abel"])
        self.assertFalse(r["per_reset_norm_sum_used"])
        self.assertFalse(r["source_uniform_twist_sum_closed"])

if __name__=="__main__":
    unittest.main()
