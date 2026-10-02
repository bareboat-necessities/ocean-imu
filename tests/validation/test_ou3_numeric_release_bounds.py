import unittest
from tools.stability.ou3_theorem.numeric_release_bounds import bounds
class NumericReleaseBoundsTests(unittest.TestCase):
 def test_explicit_AG_BG_BA_numbers(self):
  b=bounds();self.assertAlmostEqual(b["release_horizon_s"],469.)
  self.assertGreater(b["AG_covariance_spectral_upper"],0)
  self.assertAlmostEqual(b["BG_mean_error_norm_upper_rad_s"],.52)
  self.assertAlmostEqual(b["BA_graph_operator_norm_upper"],18.60665,5)
  self.assertFalse(b["LIN_mean_numeric_radius_available"])
if __name__=="__main__":unittest.main()
