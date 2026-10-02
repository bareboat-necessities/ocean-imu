import unittest,math
from tools.stability.ou3_theorem.structured_aw_mean_recurrence import *
class AWMeanTests(unittest.TestCase):
 def test_zero_source_decays(self):
  z=step(2,.9,1,.2,0,0);self.assertAlmostEqual(z["A_next"],1.8)
 def test_low_aw_axis_small(self):
  z=weighted_inj_defect(1,0,.1,.9);self.assertLess(z["D"],.01)
 def test_source_terms_linked(self):
  z=aw_mean_step(1,.8,.2,.3,.1,.4);self.assertAlmostEqual(z["A_next"],.9)
if __name__=="__main__":unittest.main()
