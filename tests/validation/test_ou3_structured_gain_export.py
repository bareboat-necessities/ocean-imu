import unittest,math
from tools.stability.ou3_theorem.inj_block_algebra import INJ
from tools.stability.ou3_theorem.structured_gain_export import *
class StructuredGainTests(unittest.TestCase):
 def test_norm_exact_diagonal_projector(self):
  x=INJ(2,-.5,.3)
  self.assertAlmostEqual(inj_spectral_norm(x),max(1.5,math.sqrt(4.09)))
 def test_s_gain(self):
  self.assertAlmostEqual(s_aw_gain_scalar(.2,.5,.5),.2)
 def test_acc_gain_zero_cross(self):
  z=acc_aw_gain(INJ(0,0,0),INJ(1,0,0),None,INJ(0,0,0),INJ(1,0,0),INJ(2,0,0),False)
  self.assertAlmostEqual(z.i,.5)
if __name__=="__main__":unittest.main()
