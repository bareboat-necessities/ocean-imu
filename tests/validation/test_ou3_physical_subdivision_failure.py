import unittest
from tools.stability.ou3_theorem.enclosure_failure import EnclosureFailure
from tools.stability.ou3_theorem.source_uniform_root import root_cells
from tools.stability.ou3_theorem.physical_history_subdivision import split_guard_failure,is_physical_coordinate
class PhysicalSubdivisionFailureTests(unittest.TestCase):
 def test_failure_carries_time(self):
  e=EnclosureFailure("variance","nonfinite",.13,{"physical_v_0_0.0":1.})
  self.assertAlmostEqual(e.time,.13);self.assertIn("t=0.130000",str(e))
 def test_split_is_physical(self):
  z=split_guard_failure(root_cells(60.)[0],.13)
  self.assertTrue(is_physical_coordinate(z["coordinate"]))
  self.assertFalse(z["generated_split"])
if __name__=="__main__":unittest.main()
