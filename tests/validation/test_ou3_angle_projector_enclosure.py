import unittest,numpy as np
from tools.stability.ou3_theorem.rank_loss_interval_factor import exact
from tools.stability.ou3_theorem.angle_projector_enclosure import projector_cone_restricted_lower
class AngleProjectorEnclosureTests(unittest.TestCase):
 def test_zero_angle_recovers_midpoint_floor(self):
  z=projector_cone_restricted_lower(exact(np.diag([0.,2.,4.]).tolist()),[1,0,0],0.)
  self.assertTrue(z["verified"]);self.assertAlmostEqual(z["restricted_eigen_lower"],2.)
 def test_small_cone_positive(self):
  z=projector_cone_restricted_lower(exact(np.diag([0.,2.,4.]).tolist()),[1,0,0],1e-3)
  self.assertTrue(z["verified"]);self.assertLess(z["restricted_eigen_lower"],2.)
 def test_large_cone_fails(self):
  z=projector_cone_restricted_lower(exact(np.diag([0.,.01,4.]).tolist()),[1,0,0],.2)
  self.assertFalse(z["verified"])
if __name__=="__main__":unittest.main()
