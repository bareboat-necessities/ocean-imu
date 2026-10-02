# ruff: noqa: F401, F811
import unittest
import numpy as np
from tools.stability.ou3_theorem.inj_block_algebra import INJ
from tools.stability.ou3_theorem.structured_reset_split import *
class ResetSplitTests(unittest.TestCase):
 def test_parallel_reset_stays_inj(self):
  b=INJ(2.,.3,.2);n=np.array([0.,0.,1.])
  z=right_reset_inj(b,np.array([0.,0.,.1]),n)
  self.assertLess(z["defect_norm"],1e-15)
 def test_perpendicular_reset_defect_reconstructs(self):
  b=INJ(2.,.3,.2);n=np.array([0.,0.,1.]);d=np.array([.1,.2,.03])
  z=right_reset_inj(b,d,n)
  x,y,q=d;X=np.array([[0.,-q,y],[q,0.,-x],[-y,x,0.]])
  exact=b.matrix(n)@(np.eye(3)-.5*X)
  recon=z["structured"].matrix(n)+z["defect"]
  self.assertLess(np.linalg.norm(exact-recon,2),1e-12)
  self.assertLessEqual(z["defect_norm"],reset_defect_bound(b,d,n)+1e-12)
if __name__=="__main__":unittest.main()
