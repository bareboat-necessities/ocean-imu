import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_structured_S_prefix as S

class T(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=S.build(source_pieces=2,source_cell_index=0,d_pieces=4)
 def test_semantics(self):
  self.assertTrue(self.d["sample0_due_implies_sample1_due_for_refined_source_cell"])
  self.assertTrue(self.d["exact_first_Joseph_reset_gauge_child_used"])
  self.assertTrue(self.d["one_step_process_prediction_included"])
  self.assertTrue(self.d["actual_estimator_S_residual_used"])
  self.assertTrue(self.d["sample1_S_shipping_gain_and_Joseph_reset_used"])
  self.assertTrue(self.d["rigorous_S_ge_R_spectral_inverse_is_admissible"])
 def test_finite_enclosure(self):
  n=self.d["evaluated_d_cells"]
  self.assertGreater(n,0)
  self.assertEqual(self.d["fixed_pivot_inverse_count"]+self.d["spectral_fallback_inverse_count"],2*n)
  self.assertTrue(math.isfinite(self.d["max_sample1_S_attitude_correction_norm_upper_rad"]))
  self.assertTrue(math.isfinite(self.d["max_sample1_q_after_S_upper"]))
 def test_no_promotion(self):
  self.assertFalse(self.d["filter_changed"])
  self.assertFalse(self.d["second_accelerometer_evaluated_here"])
  self.assertFalse(self.d["complete_sample1_branch_closed_here"])
  self.assertFalse(self.d["whole_word_promoted_here"])
  self.assertFalse(self.d["N_H_words_set_here"])
 def test_validate(self): self.assertEqual(S.validate(self.d),[])

if __name__=="__main__": unittest.main()
