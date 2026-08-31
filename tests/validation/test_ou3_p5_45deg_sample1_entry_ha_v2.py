import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import ou3_p5_45deg_sample1_entry_ha_v2 as V2


class Ou3P545DegSample1EntryHAV2Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=V2.build(source_pieces=2)

 def test_source_reachability_is_pinned(self):
  d=self.d
  self.assertTrue(all(d['source_reachability'].values()),d['source_reachability'])
  self.assertEqual(d['deployed_default_first_live_mode'],'H')
  self.assertFalse(d['deployed_default_A_first_live_reachable'])
  self.assertEqual(d['A_sample0_evaluation_role'],'CONFIGURATION_DIAGNOSTIC_AND_LATER_H_TO_A_BOUND')

 def test_outer_projection_is_exact_branch_not_smooth_p4_promotion(self):
  d=self.d
  self.assertTrue(d['A_shipping_closed_ball_projection_used'])
  self.assertTrue(d['A_projection_nonexpansive_used'])
  self.assertFalse(d['A_projection_surface_forbidden_in_outer_P5'])
  self.assertFalse(d['A_P4_inner_045_ball_promoted_here'])

 def test_H_and_projection_aware_A_diagnostic_reach_sample1(self):
  d=self.d
  self.assertEqual(V2.validate(d),[])
  self.assertEqual(d['P5_45DEG_SAMPLE1_ENTRY_HA_V2_CERTIFICATE'],'PASS')
  self.assertTrue(d['sample1_entry_inside_q8'])
  self.assertLess(d['sample1_pre_measurement_q_upper'],8.0)
  self.assertTrue(d['modes']['H']['complete'])
  self.assertTrue(d['modes']['A']['complete'])
  self.assertIsNone(d['modes']['H']['first_failure'])
  self.assertIsNone(d['modes']['A']['first_failure'])

 def test_projection_does_not_change_filter_or_promote_capture(self):
  d=self.d
  self.assertTrue(d['source_generated_not_trajectory_fit'])
  self.assertFalse(d['source_replay_used'])
  self.assertFalse(d['filter_changed'])
  self.assertFalse(d['sample1_measurement_prefix_evaluated_here'])
  self.assertFalse(d['returned_to_30deg_P4_sector_here'])
  self.assertFalse(d['P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE'])
  self.assertFalse(d['P5_FINITE_CAPTURE_TO_P4_ESTABLISHED_HERE'])


if __name__=='__main__': unittest.main()
