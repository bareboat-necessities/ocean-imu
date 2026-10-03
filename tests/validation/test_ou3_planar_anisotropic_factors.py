import unittest, numpy as np
from tools.stability.ou3_theorem.planar_anisotropic_factors import parity_factors,generalized_eta
class PlanarAnisotropicFactorsTests(unittest.TestCase):
 def test_spd_parity_factors(self):
  f=parity_factors()
  for x in f.values():
   B=np.array(x["lower_matrix"]);L=np.array(x["cholesky"])
   self.assertTrue(np.allclose(L@L.T,B));self.assertGreater(x["lambda_min"],0)
 def test_generalized_eta_keeps_anisotropy(self):
  L=np.diag([2.,.1]);Q=np.diag([1.,.001])
  self.assertAlmostEqual(generalized_eta(Q,L),.25)
if __name__=="__main__":unittest.main()
