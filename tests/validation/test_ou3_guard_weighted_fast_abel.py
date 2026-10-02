# ruff: noqa: F401, F811
import unittest,math
from tools.stability.ou3_theorem.guard_weighted_fast_abel import *
class WeightedGuardFastAbelTests(unittest.TestCase):
 def test_zero_weight_zero_charge(self):
  z=integrated_primitive_bound(math.exp(-2*math.pi*3*.005),.005,.05,.3,60.,[0.]*100)
  self.assertEqual(z["bound"],0.)
 def test_bound_scales_with_C_not_C_over_h(self):
  a=math.exp(-2*math.pi*3*.005)
  z=integrated_primitive_bound(a,.005,.05,.3,60.,[.2]*100)
  self.assertAlmostEqual(z["bound"],.05*z["charge"])
 def test_slow_slew_finite_charge(self):
  a=math.exp(-2*math.pi*3*.005);s=1-math.exp(-.005/5.)
  w=[];x=0.
  for _ in range(1000):
   x=x+s*(1-x);w.append(x)
  z=integrated_primitive_bound(a,.005,.05,.3,60.,w)
  self.assertTrue(math.isfinite(z["bound"]))
  self.assertLess(z["weight_variation"],1.)
if __name__=="__main__":unittest.main()
