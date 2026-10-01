import unittest
from tools.stability.ou3_theorem.h18_release import *
class T(unittest.TestCase):
 def test_gates(self):
  r=default_captured_refinement(); self.assertTrue(r["all_tuner_gates_uniformly_pass"]); self.assertGreater(r["margins"]["norm_ratio_margin"],0); self.assertGreater(r["margins"]["horizontal_fraction_margin"],0)
 def test_release(self):
  r=release_compactness(); self.assertTrue(r["A21_release_set_compact_conditional_on_captured_domain"]); self.assertFalse(r["general_capture_into_that_domain_proved"]); self.assertFalse(r["release_into_inner_0p15_ball_proved"])
if __name__=="__main__": unittest.main()
