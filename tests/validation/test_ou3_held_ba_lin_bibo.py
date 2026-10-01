from fractions import Fraction as F
import unittest
from tools.stability.ou3_theorem.held_ba_lin_bibo import (
 four_s_zero_action_kernel,covariance_metric_nonexpansion,bibo_reduction,certificate)

class HeldBALinBiboTests(unittest.TestCase):
 def test_three_distinct_S_rows_kill_neutral_chain(self):
  r=four_s_zero_action_kernel((0,1,2,3))
  self.assertEqual(r["neutral_three_row_det"],1)
  self.assertTrue(r["full_12D_kernel_zero"])
 def test_bad_knots_fail(self):
  with self.assertRaises(ValueError): four_s_zero_action_kernel((0,1,1,2))
 def test_acc_correction_is_metric_nonexpansive_not_assumed_euclidean_damping(self):
  r=covariance_metric_nonexpansion()
  self.assertTrue(r["literal_accelerometer_correction_nonexpansive"])
  self.assertFalse(r["euclidean_gain_required"])
 def test_fail_closed_at_missing_compactness(self):
  r=bibo_reduction()
  self.assertTrue(r["four_S_structural_detectability"])
  self.assertFalse(r["these_premises_currently_proved"])
  self.assertFalse(r["all_time_BIBO_closed"])
  self.assertFalse(r["qTv_boundary_action_closed"])
if __name__=="__main__": unittest.main()
