import math
import sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_force_cone_crosscov_refinement as G

class T(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=G.build(delta_pieces=8,axial_pieces=8)
 def test_validates_fail_closed(self):
  d=self.d; self.assertEqual(G.validate(d),[]); self.assertGreater(d["evaluated_joint_cells"],0)
  self.assertIn(d["P5_SAMPLE1_FORCE_CONE_CROSSCOV_WITNESS_REFINEMENT"],("PASS","NOT_ESTABLISHED"))
  if d["P5_SAMPLE1_FORCE_CONE_CROSSCOV_WITNESS_REFINEMENT"]=="PASS": self.assertIsNone(d["first_unclosed_joint_cell"])
  else: self.assertIsNotNone(d["first_unclosed_joint_cell"])
 def test_force_cone_semantics(self):
  d=self.d; self.assertTrue(d["first_aw_tangent_gain_uses_gravity_information_denominator"]); self.assertTrue(d["first_aw_axial_gain_kept_separate"]); self.assertTrue(d["reset_and_prediction_preserve_force_cone_decomposition"]); self.assertTrue(d["transported_first_posterior_crosscovariance_used"]); self.assertTrue(d["sample1_S_identity_subbranch_only"])
  self.assertLessEqual(d["first_aw_tangent_gain_norm_upper"],d["first_aw_axial_gain_abs_upper"])
 def test_finite_diagnostics(self):
  d=self.d
  for k in ("sample1_aw_tangent_mean_norm_upper_mps2","sample1_aw_axial_mean_abs_upper_mps2","max_transported_reference_Ctheta_norm_upper","max_theta_mismatch_contribution_norm_upper","max_bg_mismatch_contribution_norm_upper","max_aw_mismatch_contribution_norm_upper","max_total_mismatch_Ctheta_norm_upper","max_actual_Ctheta_norm_upper","max_sample1_acc_correction_norm_upper_rad"): self.assertTrue(math.isfinite(d[k])); self.assertGreaterEqual(d[k],0)
 def test_no_promotion(self):
  d=self.d; self.assertEqual(d["deployed_correction_limit_rad"],6.0); self.assertFalse(d["deployed_correction_limit_increased"]); self.assertFalse(d["complete_sample1_branch_refined_here"]); self.assertFalse(d["whole_word_promoted_here"]); self.assertFalse(d["N_H_words_set_here"])
if __name__=="__main__":unittest.main()
