import unittest
from tools.stability.ou3_theorem.planar_common_norm_radius import solve
class CommonNormRadiusTest(unittest.TestCase):
 def test_fails_closed_without_rigorous_ports(self):
  o=solve(.8738212970667967,.5827057639270218,None,None,None,None,6.024764605642485)
  self.assertFalse(o["joint_cell_forward_invariant"]);self.assertFalse(o["all_time_magnetic_service_verified"])
 def test_linked_solver_algebra(self):
  o=solve(.8,.5,.1,.2,.01,.02,6.,1.)
  self.assertTrue(o["joint_cell_forward_invariant"]);self.assertGreater(o["candidate_P_radius"],0)
 def test_relaxed_failure_is_D_not_counterexample(self):
  o=solve(.8,.5,1.,1.,.01,.02,6.,1.)
  self.assertEqual(o["failure_class"],"D_SUFFICIENT_BOUND_FAILURE")
  self.assertFalse(o["all_time_magnetic_service_verified"])
if __name__=="__main__":unittest.main()
