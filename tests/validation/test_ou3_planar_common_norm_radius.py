import unittest
from tools.stability.ou3_theorem.planar_common_norm_radius import solve
class CommonNormRadiusTest(unittest.TestCase):
 def test_fails_closed_without_rigorous_ports(self):
  o=solve(.8738212970667967,None,None,None,None,None,6.024764605642485)
  self.assertFalse(o["joint_cell_forward_invariant"]);self.assertFalse(o["all_time_magnetic_service_verified"])
 def test_linked_solver_algebra(self):
  o=solve(.8,.5,.1,.2,.01,.02,6.,1.,gauge_injection=.3,gauge_amplitude=.1,covariance_gauge_injection=0.)
  self.assertTrue(o["comparison_fixed_point_exists"])
  self.assertAlmostEqual(o["candidate_P_radius"],.125)
  self.assertAlmostEqual(o["candidate_quotient_radius"],.15)
  self.assertTrue(o["candidate_information_margin_positive"])
  for key in ("joint_cell_forward_invariant","information_cell_perturbation_certified","all_time_magnetic_service_verified","theorem_closed"):
   self.assertFalse(o[key])
 def test_relaxed_failure_is_D_not_counterexample(self):
  o=solve(.8,.5,1.,1.,.01,.02,6.,1.,gauge_injection=0.,gauge_amplitude=0.,covariance_gauge_injection=0.)
  self.assertEqual(o["failure_class"],"D_SUFFICIENT_BOUND_FAILURE")
  self.assertFalse(o["all_time_magnetic_service_verified"])
 def test_missing_gauge_cannot_be_discarded(self):
  self.assertFalse(solve(.8,.5,.1,.2,.01,.02,6.,1.)["comparison_fixed_point_exists"])
 def test_invalid_comparison_and_positive_determinant_trap(self):
  for v in (-1.,float('nan'),float('inf')):
   with self.assertRaises(ValueError):solve(v,.5,.1,.2,.01,.02,6.)
  o=solve(2.,2.,0.,0.,.01,.02,6.,gauge_injection=0.,gauge_amplitude=0.,covariance_gauge_injection=0.)
  self.assertEqual(o["failure_class"],"D_SUFFICIENT_BOUND_FAILURE")
if __name__=="__main__":unittest.main()
