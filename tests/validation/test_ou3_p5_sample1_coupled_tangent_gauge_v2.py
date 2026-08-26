import math
import sys
import unittest
from pathlib import Path

TOOLS=Path(__file__).resolve().parents[2]/"tools"
if str(TOOLS) not in sys.path: sys.path.insert(0,str(TOOLS))
import ou3_p5_sample1_coupled_tangent_gauge_v2 as G

class T(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=G.build(source_pieces=2,source_cell_index=0,delta_pieces=4,axial_pieces=4)
 def test_gauge_contract(self):
  self.assertTrue(self.d["sample0_canonical_specific_force_is_plus_g_e3"])
  self.assertTrue(self.d["sample1_zero_correction_reference_specific_force_is_plus_g_e3"])
  self.assertTrue(self.d["gravity_sign_consistent_with_canonical_H0"])
  self.assertTrue(self.d["same_first_residual_attitude_aw_coupling_used"])
  self.assertTrue(self.d["structured_tangent_gain_cancellation_used_before_interval_inversion"])
 def test_finite(self):
  for k in ("first_tangent_relation_remainder_norm_upper_mps2","max_E_theta_norm_upper","max_Pj_Et_first6_norm_upper","max_actual_Ctheta_norm_upper","max_sample1_acc_correction_norm_upper_rad"):
   self.assertTrue(math.isfinite(float(self.d[k])),k)
  self.assertEqual(self.d["fixed_pivot_inverse_count"]+self.d["spectral_fallback_inverse_count"],self.d["evaluated_joint_cells"])
 def test_no_promotion(self):
  self.assertFalse(self.d["filter_changed"]); self.assertEqual(self.d["deployed_correction_limit_rad"],6.0); self.assertFalse(self.d["deployed_correction_limit_increased"]); self.assertFalse(self.d["complete_sample1_branch_refined_here"]); self.assertFalse(self.d["whole_word_promoted_here"]); self.assertFalse(self.d["N_H_words_set_here"])
 def test_validate(self): self.assertEqual(G.validate(self.d),[])
if __name__=="__main__": unittest.main()
