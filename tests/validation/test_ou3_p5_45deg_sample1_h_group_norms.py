import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import ou3_p5_45deg_sample1_h_group_norms as G
import ou3_p5_45deg_first_accel_signed_source_bound_v2 as S


class Ou3P545DegSample1HGroupNormTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.d=G.build(source_pieces=2)
  cls.s=S.build(source_pieces=2,tangent_cells=96)

 def test_certificate_and_semantics(self):
  d=self.d
  self.assertEqual(G.validate(d),[])
  self.assertEqual(d['P5_45DEG_H_SAMPLE1_GROUP_NORM_CERTIFICATE'],'PASS')
  self.assertEqual(d['deployed_first_live_mode'],'H')
  self.assertTrue(d['norm_balls_not_reinterpreted_as_cartesian_cubes'])
  self.assertTrue(d['position_componentwise_half_Hs_converted_to_sqrt3_norm'])
  self.assertTrue(d['integrated_OU_axis_isotropy_used'])
  self.assertFalse(d['new_hard_bounds_invented'])

 def test_declared_norm_bounds_are_not_sqrt3_inflated_at_entry(self):
  i=self.d['initial_group_norm_caps']
  self.assertAlmostEqual(i['gyro_bias'],0.01)
  self.assertAlmostEqual(i['velocity'],5.0)
  self.assertAlmostEqual(i['S'],300.0)
  self.assertAlmostEqual(i['aw'],10.0)
  self.assertGreater(i['position'],7.3)
  self.assertLess(i['position'],7.4)

 def test_sample1_physical_norms_improve_cube_readback(self):
  d=self.d
  n=d['sample1_max_group_norm_caps']; b=d['old_interval_box_max_reported_group_norms']
  for k in ('velocity','S','aw'):
   self.assertLess(n[k],b[k],k)
  self.assertTrue(d['rows'])

 def test_signed_first_accel_bound_strictly_improves_scalar_q8_bridge(self):
  d=self.s
  self.assertEqual(S.validate(d),[])
  self.assertEqual(d['P5_45DEG_FIRST_ACCEL_SIGNED_SOURCE_BOUND'],'PASS')
  self.assertTrue(d['strictly_improves_sign_agnostic_bridge'])
  self.assertTrue(d['inside_q8'])
  self.assertFalse(d['favorable_correction_direction_assumed'])
  self.assertTrue(d['ideal_first_accel_yaw_injection_exact_zero'])
  self.assertLess(d['signed_source_correlated_post_update_q_upper'],d['sign_agnostic_scalar_post_update_q_upper'])
  print('SIGNED_45DEG_FIRST_ACCEL',
        'q_pre',d['pre_update_q_upper'],
        'q_scalar',d['sign_agnostic_scalar_post_update_q_upper'],
        'q_signed',d['signed_source_correlated_post_update_q_upper'],
        'factor',d['q_upper_improvement_factor'],
        'max_d',d['max_signed_decomposition_correction_norm_upper_rad'],
        'min_den',d['minimum_signed_composition_denominator_lower'],
        'returned30',d['returned_to_30deg_P4_sector_here'])

 def test_no_filter_or_replay_change(self):
  d=self.d
  self.assertTrue(d['source_generated_not_trajectory_fit'])
  self.assertFalse(d['source_replay_used'])
  self.assertFalse(d['filter_changed'])
  s=self.s
  self.assertTrue(s['source_generated_not_trajectory_fit'])
  self.assertFalse(s['source_replay_used'])
  self.assertFalse(s['filter_changed'])
  self.assertFalse(s['deployed_correction_limit_increased'])


if __name__=='__main__': unittest.main()
