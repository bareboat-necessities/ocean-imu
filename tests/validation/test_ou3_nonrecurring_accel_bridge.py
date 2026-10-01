import unittest
from fractions import Fraction as F
from tools.stability.ou3_theorem.nonrecurring_accel_bridge import certificate

class T(unittest.TestCase):
    def test_candidate_has_physical_sensor_mean_margin(self):
        r=certificate()
        self.assertEqual(F(r["physical_signed_mean_supply_mps2"]),F(67,80))
        self.assertEqual(F(r["candidate_fast_mean_cap_mps2"]),F(1,320))
        self.assertTrue(r["sensor_record_margin_positive"])
        self.assertGreater(F(r["sensor_record_margin_mps2"]),F(89,100))
        self.assertFalse(r["nominal_AW_mean_follows_sensor_mean"])
        self.assertFalse(r["G0_nominal_window_premise_closed"])

if __name__=="__main__":
    unittest.main()
