import unittest
from tools.stability.ou3_theorem.scheduler_s_shift_bound import certificate
class SShiftBoundTests(unittest.TestCase):
 def test_identity(self):
  c=certificate();self.assertLess(c["identity_defect"],1e-10);self.assertLess(c["defect_symmetric"],1e-10);self.assertFalse(c["phase_uniform_bound_verified"])
if __name__=="__main__":unittest.main()
