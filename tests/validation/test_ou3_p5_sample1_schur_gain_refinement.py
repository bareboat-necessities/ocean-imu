import math
import sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_schur_gain_refinement as G

class T(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=G.build(delta_pieces=8)
 def test_validates_fail_closed(self):
  d=self.d; self.assertEqual(G.validate(d),[]); self.assertGreater(d["evaluated_delta_cells"],0)
  self.assertIn(d["P5_SAMPLE1_SCHUR_GAIN_WITNESS_REFINEMENT"],("PASS","NOT_ESTABLISHED"))
  if d["P5_SAMPLE1_SCHUR_GAIN_WITNESS_REFINEMENT"]=="PASS": self.assertIsNone(d["first_unclosed_delta_cell"])
  else: self.assertIsNotNone(d["first_unclosed_delta_cell"])
 def test_schur_semantics(self):
  d=self.d; self.assertTrue(d["accepted_gain_uses_KSK_le_P"]); self.assertTrue(d["solver_failure_identity_branch_included"]); self.assertTrue(d["no_sample1_interval_inverse_used"]); self.assertTrue(d["S_posterior_loewner_upper_by_prior_used"]); self.assertTrue(d["reset_operator_exact_norm_formula_used"])
 def test_finite_outputs(self):
  d=self.d
  for k in ("max_sample1_S_attitude_correction_norm_upper_rad","max_sample1_S_aw_correction_norm_upper_mps2","max_sample1_acc_residual_norm_upper_mps2","max_sample1_acc_correction_norm_upper_rad"): self.assertTrue(math.isfinite(d[k])); self.assertGreaterEqual(d[k],0.0)
 def test_no_promotion(self):
  d=self.d; self.assertEqual(d["deployed_correction_limit_rad"],6.0); self.assertFalse(d["deployed_correction_limit_increased"]); self.assertFalse(d["complete_source_cell_refined_here"]); self.assertFalse(d["whole_word_promoted_here"]); self.assertFalse(d["N_H_words_set_here"])
if __name__=="__main__":unittest.main()
