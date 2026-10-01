import unittest
from tools.stability.ou3_theorem.six_column_candidate import evaluate
class T(unittest.TestCase):
 def test_candidate_does_not_promote_carried_G0(self):
  r=evaluate(); self.assertTrue(r["candidate_profile_physical_gate_passes"]); self.assertGreater(r["G0_injection_free_six_column_floor_squared"],0)
  self.assertFalse(r["source_uniform_nominal_AW_window_statistics"]); self.assertFalse(r["literal_injection_transport_uniform_bound"]); self.assertFalse(r["literal_six_column_J_positive_certified"])
if __name__=="__main__": unittest.main()
