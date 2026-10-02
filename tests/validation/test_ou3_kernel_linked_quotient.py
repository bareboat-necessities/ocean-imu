import unittest
from tools.stability.ou3_theorem.interval_shaped_supply import IJointQuadratic
from tools.stability.ou3_theorem.rank_loss_interval_factor import exact
from tools.stability.ou3_theorem.kernel_linked_quotient import optimize_lambda_grid
class KernelLinkedQuotientTests(unittest.TestCase):
 def test_linked_schur_can_certify(self):
  q=IJointQuadratic.zeros(1,1,"h")
  q.A=exact([[2.]]);q.B=exact([[.1]]);q.C=exact([[-.2]])
  R=exact([[1.]])
  z=optimize_lambda_grid(q,R,[.1,.2,.5,1.,2.])
  self.assertTrue(z["verified"]);self.assertGreater(z["dissipation_lower"],0)
if __name__=="__main__":unittest.main()
