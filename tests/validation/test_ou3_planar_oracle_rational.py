import unittest
from tools.stability.ou3_theorem.planar_oracle_rational import certificate
class PlanarOracleRationalTests(unittest.TestCase):
 def test_exact_oracle_closes_but_shipping_stays_open(self):
  c=certificate()
  self.assertTrue(c["one_second_covariance_self_inclusion_exact"])
  self.assertTrue(c["magnetic_information_minus_identity_spd_exact"])
  self.assertGreater(c["information_minus_identity_det_decimal"],.9)
  self.assertFalse(c["literal_shipping_service_proved"])
  self.assertFalse(c["theorem_closed"])
if __name__=="__main__":unittest.main()
