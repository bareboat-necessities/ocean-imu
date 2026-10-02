import unittest,numpy as np
from tools.stability.ou3_theorem.verified_linked_qcqp import Box,verified_interval_supply_bnb
from tools.stability.ou3_theorem.interval_riccati_21 import IMat
class IntervalQCQPTests(unittest.TestCase):
 def im(self,m,r):return IMat(tuple(tuple(x for x in row) for row in m),tuple(tuple(x for x in row) for row in r))
 def test_interval_coefficients_increase_safe_upper(self):
  B=self.im([[.2]],[[.01]]);C=self.im([[-.1]],[[.02]])
  e=Box(np.array([-1.]),np.array([1.]));u=Box(np.array([-.5]),np.array([.5]))
  A=np.array([[1.],[-1.]]);b=np.array([.5,.5])
  z=verified_interval_supply_bnb(B,C,e,u,A,b,lambda box:True,tol=.05,max_leaves=10000)
  self.assertTrue(z.verified);self.assertGreater(z.upper,.225)
if __name__=="__main__":unittest.main()
