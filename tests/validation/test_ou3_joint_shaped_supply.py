import unittest,numpy as np
from tools.stability.ou3_theorem.joint_shaped_supply import JointQuadratic
class JointShapedSupplyTests(unittest.TestCase):
 def test_exact_one_operation_identity(self):
  Q0=np.diag([2.,3.]);Q1=np.diag([4.,5.]);A=np.array([[.8,.1],[0,.7]]);D=np.array([[1.],[-.2]])
  q=JointQuadratic.zeros(2,1,"h");q.add_operation(Q0,Q1,A,D)
  e=np.array([.3,-.4]);u=np.array([.2]);en=A@e+D@u
  lhs=e@Q0@e-en@Q1@en
  z=np.r_[e,u];self.assertAlmostEqual(lhs,z@q.matrix()@z,12)
 def test_source_cross_term_retained(self):
  q=JointQuadratic.zeros(1,1,"h");q.add_operation(np.eye(1),np.eye(1),np.array([[.5]]),np.array([[2.]]))
  self.assertNotEqual(float(q.B[0,0]),0.)
if __name__=="__main__":unittest.main()
