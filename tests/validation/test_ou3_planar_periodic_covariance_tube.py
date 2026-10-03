import unittest,numpy as np
from tools.stability.ou3_theorem.planar_periodic_covariance_tube import correction_monotone_check,certificate
class PlanarPeriodicTubeTests(unittest.TestCase):
 def test_correction_monotone(self):
  lo=np.diag([1.,2.]);hi=lo+.2*np.eye(2);H=np.array([[1.,.3]]);R=np.array([[.4]])
  self.assertGreaterEqual(correction_monotone_check(lo,hi,H,R),-1e-12)
 def test_fails_closed(self):
  c=certificate();self.assertFalse(c["periodic_lower_tube_verified"]);self.assertGreater(c["required_retention"],.14)
if __name__=="__main__":unittest.main()
