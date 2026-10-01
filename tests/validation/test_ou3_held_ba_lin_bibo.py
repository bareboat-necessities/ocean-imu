import unittest
from tools.stability.ou3_theorem.held_ba_lin_bibo import *
class T(unittest.TestCase):
 def test_lift(self):
  r=h18_lin_covariance_upper(); self.assertTrue(r["P_LL_uniform_upper_after_17s_held_H18"]); self.assertFalse(r["held_BA_activity_used_in_first_four_blocks"])
 def test_compact(self): self.assertTrue(coefficient_compactness()["held_H18_LIN_coefficient_family_compact_after_17s"])
 def test_status(self):
  r=certificate(); self.assertTrue(r["uniform_homogeneous_rho_LIN_lt_1_exists"]); self.assertFalse(r["all_time_affine_BIBO_closed"])
if __name__=="__main__": unittest.main()
