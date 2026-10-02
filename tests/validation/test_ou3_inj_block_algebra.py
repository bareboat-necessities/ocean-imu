import unittest,numpy as np
from tools.stability.ou3_theorem.inj_block_algebra import INJ
from tools.stability.ou3_theorem.structured_riccati_closure import inverse_coefficients
class INJTests(unittest.TestCase):
 def test_multiplication_matches_matrix(self):
  n=np.array([.2,.3,.9327379053]);n/=np.linalg.norm(n)
  a=INJ(2,.4,-.3);b=INJ(.8,-.2,.5)
  self.assertTrue(np.allclose(a.mul(b).matrix(n),a.matrix(n)@b.matrix(n),atol=1e-10))
 def test_inverse(self):
  n=np.array([0.,0.,1.]);a=INJ(2,.5,.3);b=inverse_coefficients(a)
  self.assertTrue(np.allclose(a.matrix(n)@b.matrix(n),np.eye(3),atol=1e-12))
if __name__=="__main__":unittest.main()
