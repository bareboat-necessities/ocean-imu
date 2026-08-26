import math
import sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_crosscov_gain_refinement as G

class T(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=G.build(delta_pieces=8,force_pieces=4)
 def test_validates_fail_closed(self):
  d=self.d; self.assertEqual(G.validate(d),[]); self.assertGreater(d["evaluated_joint_cells"],0)
  self.assertIn(d["P5_SAMPLE1_CROSSCOV_GAIN_WITNESS_REFINEMENT"],("PASS","NOT_ESTABLISHED"))
  if d["P5_SAMPLE1_CROSSCOV_GAIN_WITNESS_REFINEMENT"]=="PASS": self.assertIsNone(d["first_unclosed_joint_cell"])
  else: self.assertIsNotNone(d["first_unclosed_joint_cell"])
 def test_crosscov_semantics(self):
  d=self.d; self.assertTrue(d["uses_actual_P_Ht_attitude_crosscovariance"]); self.assertTrue(d["uses_S_ge_R_inverse_norm_only"]); self.assertTrue(d["sample1_S_posterior_loewner_bound_used"]); self.assertTrue(d["sample1_J_aw_exact_identity"])
 def test_finite_outputs(self):
  d=self.d
  for k in ("max_Ctheta_norm_upper","max_Ktheta_norm_upper","max_sample1_acc_residual_norm_upper_mps2","max_sample1_acc_correction_norm_upper_rad"): self.assertTrue(math.isfinite(d[k])); self.assertGreaterEqual(d[k],0)
 def test_no_promotion(self):
  d=self.d; self.assertEqual(d["deployed_correction_limit_rad"],6.0); self.assertFalse(d["deployed_correction_limit_increased"]); self.assertFalse(d["complete_source_cell_refined_here"]); self.assertFalse(d["whole_word_promoted_here"]); self.assertFalse(d["N_H_words_set_here"])
if __name__=="__main__":unittest.main()
