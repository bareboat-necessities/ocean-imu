# ruff: noqa: F401, F811
import unittest,numpy as np
from tools.stability.ou3_theorem.interval_riccati_21 import IMat
from tools.stability.ou3_theorem.interval_kr import *
class IntervalKrTests(unittest.TestCase):
 def test_point_product(self):
  K=IMat(((1.,0.,0.),(0.,2.,0.)),((0.,0.,0.),(0.,0.,0.)))
  m,r=interval_matvec(K,np.array([3.,4.,0.]),np.zeros(3))
  self.assertTrue(np.allclose(m,[3.,8.]));self.assertTrue(np.all(r>=0))
 def test_joint_radius_contains_corner(self):
  K=IMat(((1.,0.,0.),),((.1,0.,0.),))
  m,r=interval_matvec(K,np.array([2.,0,0]),np.array([.2,0,0]))
  self.assertLessEqual(m[0]-r[0],.9*1.8);self.assertGreaterEqual(m[0]+r[0],1.1*2.2)
if __name__=="__main__":unittest.main()
