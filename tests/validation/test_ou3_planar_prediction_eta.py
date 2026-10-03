import unittest
from tools.stability.ou3_theorem.planar_prediction_eta import certificate,retention_from_etas
class PlanarPredictionEtaTests(unittest.TestCase):
 def test_matrix_eta_path(self):
  c=certificate();self.assertGreater(c["synthetic_eta"],0);self.assertGreater(c["synthetic_retention"],0);self.assertFalse(c["one_second_product_verified"])
 def test_product(self):self.assertAlmostEqual(retention_from_etas([.1,.2]),1/1.1/1.2)
if __name__=="__main__":unittest.main()
