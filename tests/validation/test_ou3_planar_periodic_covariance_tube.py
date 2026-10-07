import unittest,numpy as np
from tools.stability.ou3_theorem.planar_periodic_covariance_tube import correction_monotone_check,certificate
class PlanarPeriodicTubeTests(unittest.TestCase):
 def test_correction_monotone(self):
  lo=np.diag([1.,2.]);hi=lo+.2*np.eye(2);H=np.array([[1.,.3]]);R=np.array([[.4]])
  self.assertGreaterEqual(correction_monotone_check(lo,hi,H,R),-1e-12)
 def test_fails_closed(self):
  c=certificate();self.assertFalse(c["periodic_lower_tube_verified"])
  self.assertIsNone(c["required_retention"])
  self.assertFalse(c["prediction_product_times_replay_service_is_a_service_bound"])
  self.assertFalse(c["literal_reset_AW_and_joint_coefficient_stream_attached"])
if __name__=="__main__":unittest.main()
