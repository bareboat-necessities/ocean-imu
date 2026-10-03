import unittest
from tools.stability.ou3_theorem.paired_prediction_loss import certificate
class PairedPredictionLossTests(unittest.TestCase):
 def test_identity_and_psd_loss(self):
  c=certificate();self.assertLess(c["identity_defect"],1e-10);self.assertGreater(c["loss_min_eigenvalue"],-1e-10);self.assertGreater(c["retention_factor"],0);self.assertFalse(c["shipping_one_second_bound_verified"])
if __name__=="__main__":unittest.main()
