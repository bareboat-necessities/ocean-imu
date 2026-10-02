import unittest,numpy as np
from tools.stability.ou3_theorem.rank_loss_interval_factor import exact
from tools.stability.ou3_theorem.kernel_restricted_action import restrict_interval_action,certify_rotating_line
class KernelRestrictedActionTests(unittest.TestCase):
 def test_exact_one_dim_kernel_restriction(self):
  A=np.diag([0.,2.,3.]);R,c=restrict_interval_action(exact(A.tolist()),[1,0,0])
  self.assertTrue(c["verified"]);self.assertAlmostEqual(c["restricted_eigen_lower"],2.)
 def test_interval_line_cone_certifies_with_margin(self):
  R,c=restrict_interval_action(exact(np.diag([0.,2.,3.]).tolist()),[1,0,0],[.001,0,0])
  self.assertTrue(c["verified"]);self.assertTrue(c["whole_line_cone_covered"])
 def test_wide_line_cone_fails_closed(self):
  with self.assertRaises(ArithmeticError):restrict_interval_action(exact(np.diag([0.,.01,.01]).tolist()),[1,0,0],[.5,0,0])
 def test_later_action_breaks_first_line(self):
  A0=exact(np.diag([0.,1.,1.]).tolist());A1=exact(np.diag([1.,0.,1.]).tolist());T=exact(np.eye(3).tolist())
  R,c=certify_rotating_line(A0,[1,0,0],T,A1,[0,1,0]);self.assertTrue(c["verified"])
if __name__=="__main__":unittest.main()
