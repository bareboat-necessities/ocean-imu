import unittest
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from tools.stability.ou3_theorem.six_column_candidate import evaluate
class T(unittest.TestCase):
 def test_candidate_does_not_promote_carried_G0(self):
  r=evaluate(); self.assertTrue(r["candidate_profile_physical_gate_passes"]); self.assertGreater(r["G0_injection_free_six_column_floor_squared"],0)
  self.assertFalse(r["proposed_zero_action_implies_physical_compatibility"]); self.assertFalse(r["outer_physical_to_nominal_accelerometer_bridge"]); self.assertFalse(r["literal_six_column_J_positive_certified"]); self.assertFalse(r["theorem_closed"])
if __name__=="__main__": unittest.main()
