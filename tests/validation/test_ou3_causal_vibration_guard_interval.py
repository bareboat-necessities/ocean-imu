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
if __name__=='__main__':unittest.main()
