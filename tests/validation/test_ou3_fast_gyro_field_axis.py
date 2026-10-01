import unittest
from fractions import Fraction as F
from tools.stability.ou3_theorem.fast_gyro_field_axis import certificate

class T(unittest.TestCase):
    def test_candidate_cap_kills_old_fast_gyro_obstruction(self):
        r=certificate()
        self.assertTrue(r["old_sinusoid_excluded_by_candidate_fast_cap"])
        self.assertLess(F(r["fast_signed_charge_mps2"]),F("0.000648"))
        self.assertLess(F(r["slow_plus_fast_charge_mps2"]),F("0.013644"))
        self.assertFalse(r["coefficient_rotation_variation_included"])
        self.assertFalse(r["theorem_closed"])

if __name__=="__main__":
    unittest.main()
