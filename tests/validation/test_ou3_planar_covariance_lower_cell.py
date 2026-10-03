import unittest
from tools.stability.ou3_theorem.planar_covariance_lower_cell import certificate,event_counts,variational_lower_from_action
class PlanarCovarianceLowerCellTests(unittest.TestCase):
 def test_scheduler_phase_counts_and_fail_closed(self):
  c=certificate();self.assertGreaterEqual(c["scheduler_phase_uniform_event_counts"]["S_max"],6);self.assertFalse(c["lower_cell_verified"])
 def test_variational_inverse(self):
  x=variational_lower_from_action([[2.,0.],[0.,4.]])
  self.assertAlmostEqual(x["covariance_scalar_floor"],.25)
if __name__=="__main__":unittest.main()
