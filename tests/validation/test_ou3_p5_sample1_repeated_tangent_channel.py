import math
import sys
import unittest
from pathlib import Path
TOOLS=Path(__file__).resolve().parents[2]/"tools"
if str(TOOLS) not in sys.path: sys.path.insert(0,str(TOOLS))
import ou3_p5_sample1_repeated_tangent_channel as T

class RepeatedTangentTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=T.build(source_pieces=2,source_cell_index=0,p_pieces=4,axial_pieces=4)
 def test_semantics(self):
  self.assertTrue(self.d["exact_first_structured_posterior_used"])
  self.assertTrue(self.d["dependency_preserving_positive_repeated_innovation_formula_used"])
  self.assertTrue(self.d["signed_aligned_sample1_force_used"])
  self.assertFalse(self.d["reset_process_and_tangent_force_perturbations_included"])
 def test_positive_finite(self):
  self.assertGreater(self.d["minimum_scalar_innovation_variance_lower"],0.0)
  for k in ("max_scalar_Ktheta_abs_upper","max_scalar_correction_norm_upper_rad"):
   self.assertTrue(math.isfinite(float(self.d[k])))
  self.assertEqual(self.d["evaluated_joint_cells"],16)
 def test_no_promotion(self):
  self.assertFalse(self.d["filter_changed"]); self.assertEqual(self.d["deployed_correction_limit_rad"],6.0); self.assertFalse(self.d["deployed_correction_limit_increased"]); self.assertFalse(self.d["complete_sample1_branch_refined_here"]); self.assertFalse(self.d["whole_word_promoted_here"]); self.assertFalse(self.d["N_H_words_set_here"])
 def test_validate(self): self.assertEqual(T.validate(self.d),[])
if __name__=="__main__": unittest.main()
