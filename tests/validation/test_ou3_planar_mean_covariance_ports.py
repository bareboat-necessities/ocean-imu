import unittest
from tools.stability.ou3_theorem.planar_mean_covariance_ports import certificate
class MeanCovariancePortTests(unittest.TestCase):
 def test_literal_planar_row_bounds(self):
  c=certificate()
  self.assertGreater(c["acc_H_lipschitz_attitude"],9.80665)
  self.assertLess(c["acc_H_lipschitz_attitude"],9.809)
  self.assertEqual(c["mag_H_lipschitz_attitude"],75.0)
  self.assertGreater(c["innovation_inverse_norm_upper_from_carried_floor"],24.0)
  self.assertLess(c["innovation_inverse_norm_upper_from_carried_floor"],24.5)
  self.assertFalse(c["joint_cell_forward_invariant"])
if __name__=="__main__":unittest.main()
