import math
import sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_reset_perp_scalar_channel as C
class T(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=C.build(source_pieces=2,source_cell_index=0,p_pieces=4,d_pieces=4,axial_pieces=4)
 def test_semantics(self):
  self.assertTrue(self.d["exact_reset_plus_body_gauge_yaw_mixing_used"])
  self.assertTrue(self.d["one_step_attitude_PSD_remainder_included"])
  self.assertTrue(self.d["one_step_OU_process_variance_included"])
  self.assertTrue(self.d["aligned_force_core_only"])
  self.assertFalse(self.d["sample1_tangent_force_perturbation_included"])
 def test_finite(self):
  self.assertGreater(self.d["evaluated_joint_cells"],0)
  self.assertGreater(self.d["minimum_scalar_innovation_variance_lower"],0.0)
  for k in ("max_post_reset_predicted_tangent_variance","max_Ktheta_abs_upper","max_correction_norm_upper_rad"):
   self.assertTrue(math.isfinite(float(self.d[k])))
 def test_no_promotion(self):
  self.assertFalse(self.d["filter_changed"])
  self.assertFalse(self.d["complete_sample1_branch_closed_here"])
  self.assertFalse(self.d["whole_word_promoted_here"])
  self.assertFalse(self.d["N_H_words_set_here"])
 def test_validate(self): self.assertEqual(C.validate(self.d),[])
if __name__=="__main__": unittest.main()
