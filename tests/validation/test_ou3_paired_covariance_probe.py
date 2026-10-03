import unittest
from tools.stability.ou3_theorem.paired_covariance_probe import certificate
class PairedCovarianceProbeTests(unittest.TestCase):
 def test_exact_correction_storage_identity(self):
  c=certificate();self.assertLess(c["identity_defect"],1e-10);self.assertFalse(c["shipping_service_proved"])
if __name__=="__main__":unittest.main()
