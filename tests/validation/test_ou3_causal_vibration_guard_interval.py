import unittest
from tools.stability.ou3_theorem.causal_tuner_interval import I
from tools.stability.ou3_theorem.causal_vibration_guard_interval import initial_guard
class GuardIntervalTests(unittest.TestCase):
 def test_constant_input_stays_transparent(self):
  g=initial_guard();a=(I(0,0),I(0,0),I(9.8,9.8));z=g.step(a,.005)
  for _ in range(20):z=g.step(a,.005)
  self.assertEqual(z['excess'].lo,0.);self.assertEqual(z['weight'].lo,0.)
 def test_no_midpoint_branch(self):
  g=initial_guard();g.step((I(0,0),I(0,0),I(0,0)),.005)
  z=g.step((I(-1,1),I(0,0),I(0,0)),.005)
  self.assertFalse(z['branch_midpoint_used'])
 def test_lp_raw_difference_starts_exactly_zero(self):
  g=initial_guard();a=(I(0,0),I(0,0),I(9.8,9.8))
  z=g.step(a,.005,.4,.7)
  self.assertEqual(g.lp_raw_diff_norm,[0.0]*4)
  self.assertAlmostEqual(z["conditioned_norm_lower"],.4)
 def test_conditioned_floor_uses_weighted_lp_difference(self):
  g=initial_guard();a=(I(0,0),I(0,0),I(9.8,9.8));g.step(a,.005,.4,.7)
  z=g.step((I(-1,1),I(-1,1),I(8.8,10.8)),.005,.4,.7)
  self.assertIn("lp_raw_diff_norm_upper",z)
  self.assertGreaterEqual(z["conditioned_norm_lower"],0.)
  self.assertLessEqual(z["conditioned_norm_lower"],.4+1e-15)
 def test_fast_charge_is_not_recursively_accumulated(self):
  g=initial_guard();a=(I(0,0),I(0,0),I(9.8,9.8));g.step(a,.005,.4,0.,0.)
  g.step(a,.005,.4,0.,.12)
  self.assertAlmostEqual(g.smooth_lp_raw_diff_norm[1],0.)
  self.assertAlmostEqual(g.lp_raw_diff_norm[1],.12)
  g.step(a,.005,.4,0.,.13)
  self.assertAlmostEqual(g.lp_raw_diff_norm[1],.13)
if __name__=='__main__':unittest.main()
