import unittest
from tools.stability.ou3_theorem.reachable_graph import ReachableGraphCertificate,reject_cartesian_generated_boxes

class ReachableGraphTests(unittest.TestCase):
 def test_existence_not_constructive_entry(self):
  c=ReachableGraphCertificate(100,True,True,True,True,True,False,False,False)
  self.assertTrue(c.validates_source_uniform_existence())
  self.assertFalse(c.validates_constructive_entry())
 def test_independent_tuner_box_rejected(self):
  with self.assertRaises(ValueError):
   reject_cartesian_generated_boxes({"tau":(.02,12),"sigma_aw":(.05,4)})
 def test_history_parameterization_allowed(self):
  self.assertTrue(reject_cartesian_generated_boxes({"history_cell":"shared causal cell"}))
if __name__=="__main__": unittest.main()
