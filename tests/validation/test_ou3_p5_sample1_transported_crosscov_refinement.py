import math
import sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_transported_crosscov_refinement as G

class T(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=G.build(delta_pieces=8,force_pieces=4)
 def test_validates_fail_closed(self):
  d=self.d; self.assertEqual(G.validate(d),[]); self.assertGreater(d["evaluated_joint_cells"],0)
  self.assertIn(d["P5_SAMPLE1_TRANSPORTED_CROSSCOV_WITNESS_REFINEMENT"],("PASS","NOT_ESTABLISHED"))
  if d["P5_SAMPLE1_TRANSPORTED_CROSSCOV_WITNESS_REFINEMENT"]=="PASS": self.assertIsNone(d["first_unclosed_joint_cell"])
  else: self.assertIsNotNone(d["first_unclosed_joint_cell"])
 def test_identity_transport_semantics(self):
  d=self.d; self.assertTrue(d["first_posterior_identity_P_Ht_equals_KR_used"]); self.assertTrue(d["dual_reference_jacobian_satisfies_Href_A_equals_H0"]); self.assertTrue(d["process_crosscovariance_term_retained"]); self.assertTrue(d["actual_H_minus_Href_mismatch_retained"]); self.assertTrue(d["sample1_S_identity_subbranch_only"])
 def test_finite_outputs(self):
  d=self.d
  for k in ("max_transported_reference_Ctheta_norm_upper","max_reference_mismatch_Ctheta_norm_upper","max_actual_Ctheta_norm_upper","max_sample1_acc_correction_norm_upper_rad"): self.assertTrue(math.isfinite(d[k])); self.assertGreaterEqual(d[k],0)
 def test_no_promotion(self):
  d=self.d; self.assertEqual(d["deployed_correction_limit_rad"],6.0); self.assertFalse(d["deployed_correction_limit_increased"]); self.assertFalse(d["complete_sample1_branch_refined_here"]); self.assertFalse(d["whole_word_promoted_here"]); self.assertFalse(d["N_H_words_set_here"])
if __name__=="__main__":unittest.main()
