import unittest
from fractions import Fraction as F
from tools.stability.ou3_theorem.candidate_joint_qualification import certificate

class T(unittest.TestCase):
    def test_candidate_excludes_known_oscillatory_witness(self):
        r=certificate()
        self.assertEqual(F(r["fast_30s_accumulation_cap_accel_mps"]),F(1,20))
        self.assertEqual(F(r["fast_30s_accumulation_cap_gyro_rad"]),F(1,500))
        self.assertTrue(r["old_oscillatory_witness_excluded_by_accel"])
        self.assertTrue(r["old_oscillatory_witness_excluded_by_gyro"])
        self.assertTrue(r["old_oscillatory_witness_jointly_excluded"])
        self.assertFalse(r["full_same_history_gauge_breaking_closed"])
        self.assertTrue(r["hardware_qualification_pending"])

if __name__=="__main__":
    unittest.main()
