import math
import unittest
from tools.stability.ou3_theorem.marine_30s_qualification import certificate

class T(unittest.TestCase):
    def test_exact_profile_consequences(self):
        r=certificate()
        self.assertEqual(r["T_E_s"],30.0)
        self.assertEqual(r["T_P_s"],30.0)
        self.assertAlmostEqual(r["theta_E_rad"],math.pi/90)
        self.assertAlmostEqual(r["attitude_span_average_scale_rad_s"],math.pi/2700)
        self.assertEqual(r["displacement_velocity_chord_lower_mps"],.001)
        self.assertEqual(r["signed_acceleration_chord_lower_m"],0.0)
        self.assertTrue(r["zero_translation_MOVING_excluded"])
        self.assertFalse(r["pointwise_acceleration_floor_inferred"])
        self.assertFalse(r["full_joint_gauge_breaking_closed"])

if __name__=="__main__":
    unittest.main()
