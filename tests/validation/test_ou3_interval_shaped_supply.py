import unittest
from tools.stability.ou3_theorem.interval_riccati_21 import diagonal_interval
from tools.stability.ou3_theorem.interval_shaped_supply import verified_precision,IJointQuadratic
from tools.stability.ou3_theorem.rank_loss_interval_factor import exact
class IntervalShapedSupplyTests(unittest.TestCase):
 def test_full_precision_contains_exact_diagonal_inverse(self):
  P=diagonal_interval((1.,2.,3.),(1.,2.,3.));Q,c=verified_precision(P);self.assertTrue(c["verified"])
  self.assertAlmostEqual(Q.mid[0][0],1.);self.assertAlmostEqual(Q.mid[1][1],.5)
 def test_exact_operation_interval_joint(self):
  Q=diagonal_interval((2.,3.),(2.,3.));A=exact([[.5,0],[0,.8]]);D=exact([[.1],[.2]])
  q=IJointQuadratic.zeros(2,1,"h");q.add_operation(Q,Q,A,D)
  self.assertEqual(q.dependency_token,"h");self.assertGreater(q.A.mid[0][0]-q.A.rad[0][0],0)
if __name__=="__main__":unittest.main()
