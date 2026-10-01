import unittest
from tools.stability.ou3_theorem.candidate_profile_feasibility import evaluate
class T(unittest.TestCase):
 def test_candidate_passes_existing_physical_gate_but_not_final_theorem(self):
  r=evaluate(); self.assertTrue(r["necessary_physical_separation_gate_passes"]); self.assertGreater(r["remaining_physical_angle_margin_rad"],0)
  self.assertIsNone(r["delta_ann_numeric_source_uniform"]); self.assertIsNone(r["decisive_ratio_numeric"]); self.assertFalse(r["full_mathematical_theorem_closed_for_candidate_profile"])
if __name__=="__main__": unittest.main()
