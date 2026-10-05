import unittest
from tools.stability.ou3_theorem.planar_linked_riccati_mean import certificate
class LinkedRiccatiMeanTest(unittest.TestCase):
 def test_identity_is_linked(self):
  c=certificate();self.assertFalse(c["independent_extrema_used"]);self.assertFalse(c["finite_word_uniform_bound_closed"])
  self.assertIn("dK*r+K*dr",c["mean_differential"])
if __name__=="__main__":unittest.main()
