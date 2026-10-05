import unittest
from tools.stability.ou3_theorem.planar_generated_port_structure import certificate
class GeneratedPortStructureTest(unittest.TestCase):
 def test_default_tuner_has_no_mekf_feedback(self):
  c=certificate()
  self.assertEqual(c["mekf_mean_or_covariance_to_tuner_gain"],0)
  self.assertEqual(c["mekf_mean_or_covariance_to_S_period_gain"],0)
  self.assertEqual(c["mekf_mean_or_covariance_to_AW_target_gain"],0)
  self.assertTrue(c["tuner_to_covariance_is_one_way"])
  self.assertFalse(c["all_time_magnetic_service_verified"])
if __name__=="__main__":unittest.main()
