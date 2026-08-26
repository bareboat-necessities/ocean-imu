import math
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"tools"))
import ou3_p5_sample1_reset_tangent_structure as D

class T(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=D.build(source_pieces=2,source_cell_index=0,p_pieces=4,d_pieces=4)
 def test_semantics(self):
  self.assertTrue(self.d["first_tangent_scalar_posterior_used"])
  self.assertTrue(self.d["first_correction_axis_fixed_by_rotational_symmetry"])
  self.assertTrue(self.d["shipping_left_error_reset_exact_yz_block_used"])
  self.assertTrue(self.d["proof_gauge_exact_yz_rotation_used"])
  self.assertFalse(self.d["process_noise_included_here"])
  self.assertFalse(self.d["sample1_S_due_not_due_included_here"])
 def test_finite(self):
  self.assertGreater(self.d["evaluated_cells"],0)
  for k in ("max_post_reset_gauge_tangent_variance","max_post_reset_gauge_yaw_variance","max_post_reset_gauge_tangent_yaw_covariance_abs","max_tangent_variance_multiplier_upper"):
   self.assertTrue(math.isfinite(float(self.d[k])))
   self.assertGreaterEqual(float(self.d[k]),0.0)
 def test_no_promotion(self):
  self.assertFalse(self.d["filter_changed"])
  self.assertFalse(self.d["complete_sample1_branch_closed_here"])
  self.assertFalse(self.d["whole_word_promoted_here"])
  self.assertFalse(self.d["N_H_words_set_here"])
 def test_validate(self): self.assertEqual(D.validate(self.d),[])

if __name__=="__main__": unittest.main()
