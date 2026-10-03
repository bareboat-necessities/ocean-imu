import unittest
from tools.stability.ou3_theorem.scheduler_branch_enumeration import certificate,representatives
class SchedulerBranchEnumerationTests(unittest.TestCase):
 def test_finite(self):
  c=certificate(0.1363605111837387,.005,4000);self.assertTrue(c["finite_branch_reduction"]);self.assertGreater(c["open_cells"],1)
if __name__=="__main__":unittest.main()
