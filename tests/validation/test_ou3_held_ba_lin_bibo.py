import unittest
from tools.stability.ou3_theorem.held_ba_lin_bibo import *
class T(unittest.TestCase):
 def test_covariance_lift(self): self.assertTrue(h18_lin_covariance_upper()["P_LL_uniform_upper_after_17s_held_H18"])
 def test_affine_compactness(self):
  r=affine_input_compactness(); self.assertTrue(r["source_uniform_fixed_word_affine_bound_exists"]); self.assertFalse(r["fast_temporal_H_C_needed_for_BIBO"])
 def test_bibo_closes_without_promoting_outer_supply(self):
  r=certificate(); self.assertTrue(r["all_time_affine_BIBO_closed"]); self.assertTrue(r["qTv_boundary_action_closed"]); self.assertTrue(r["release_LIN_mean_compactness_closed"]); self.assertFalse(r["sharp_outer_entry_supply_closed"])
if __name__=="__main__": unittest.main()
