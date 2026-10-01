import unittest
from tools.stability.ou3_theorem.physical_constant_qualification import audit
class T(unittest.TestCase):
 def test_fail_closed(self):
  r=audit()
  for k in ("H_a","C_a","H_g","C_g","T_E","theta_E","T_P","P_E"): self.assertEqual(r[k]["status"],"OPEN")
  self.assertIsNone(r["chi_over_delta_numeric"]); self.assertFalse(r["may_invent_values"])
if __name__=="__main__": unittest.main()
