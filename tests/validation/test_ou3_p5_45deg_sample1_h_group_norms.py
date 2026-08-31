import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import ou3_p5_45deg_sample1_h_group_norms as G


class Ou3P545DegSample1HGroupNormTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=G.build(source_pieces=2)

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

 def test_no_filter_or_replay_change(self):
  d=self.d
  self.assertTrue(d['source_generated_not_trajectory_fit'])
  self.assertFalse(d['source_replay_used'])
  self.assertFalse(d['filter_changed'])


if __name__=='__main__': unittest.main()
