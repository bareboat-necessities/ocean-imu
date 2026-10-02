import unittest,math
from tools.stability.ou3_theorem.inj_block_algebra import INJ
from tools.stability.ou3_theorem.structured_gain_export import *
class StructuredGainTests(unittest.TestCase):
 def test_norm_exact_diagonal_projector(self):
  x=INJ(2,-.5,.3)
  self.assertAlmostEqual(inj_spectral_norm(x),max(1.5,math.sqrt(4.09)))
 def test_s_gain_scalar_special_case(self):
  self.assertAlmostEqual(s_aw_gain_scalar(.2,.5,.5),.2)
 def test_s_gain_retains_inj_structure(self):
  z=s_aw_gain(INJ(.2,.1,.03),INJ(.5,.2,.04),INJ(.5,0,0))
  self.assertTrue(math.isfinite(inj_spectral_norm(z)))
  self.assertNotEqual(z.n,0)
 def test_acc_gain_zero_cross(self):
  z=acc_aw_gain(INJ(0,0,0),INJ(1,0,0),None,INJ(0,0,0),INJ(1,0,0),INJ(2,0,0),False)
  self.assertAlmostEqual(z.i,.5)
 def test_lever_bg_column_is_retained(self):
  z=acc_aw_gain(INJ(0,0,0),INJ(0,0,0),None,INJ(0,0,0),INJ(1,0,0),
                INJ(2,0,0),False,INJ(1,0,0),INJ(.4,0,0),True)
  self.assertAlmostEqual(z.i,.2)
 def test_mag_aw_gain(self):
  z=mag_aw_gain(INJ(.4,0,0),INJ(.5,0,0),INJ(2,0,0))
  self.assertAlmostEqual(z.i,.1)
if __name__=="__main__":unittest.main()
