import unittest
from tools.stability.ou3_theorem.lower_factor_transport import certificate
class LowerFactorTransportTests(unittest.TestCase):
 def test_matrix_correction_transport(self):
  c=certificate();self.assertGreater(c["synthetic_loewner_margin"],-1e-10);self.assertFalse(c["literal_information_ceilings_verified"])
if __name__=="__main__":unittest.main()
