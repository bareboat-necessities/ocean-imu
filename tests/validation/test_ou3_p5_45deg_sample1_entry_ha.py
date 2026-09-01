import sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'))
import ou3_p5_45deg_sample1_entry_ha as E

class Ou3P545DegSample1EntryHATests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.d=E.build(source_pieces=2)
 def test_semantics_and_dimensions(self):
  d=self.d
  self.assertEqual(d['H_dimension'],18); self.assertEqual(d['A_dimension'],21)
  self.assertTrue(d['starts_from_45deg_sign_complete_q8_bridge'])
  self.assertTrue(d['position_entrance_uses_half_Hs'])
  self.assertTrue(d['shipping_Joseph_update_used'])
  self.assertTrue(d['shipping_left_error_reset_congruence_used'])
  self.assertTrue(d['accepted_and_identity_branches_hulled'])
  self.assertFalse(d['source_replay_used']); self.assertFalse(d['filter_changed'])
 def test_no_promotion_shortcut(self):
  d=self.d
  self.assertFalse(d['sample1_measurement_prefix_evaluated_here'])
  self.assertFalse(d['returned_to_30deg_P4_sector_here'])
  self.assertFalse(d['P4_COMPLETE_WORD_DISSIPATION_ESTABLISHED_HERE'])
  self.assertFalse(d['P5_FINITE_CAPTURE_TO_P4_ESTABLISHED_HERE'])
 def test_state_bounds_preserve_operator_norm_structure(self):
  d=self.d
  self.assertTrue(d['state_corrections_use_group_operator_norm_caps'])
  self.assertTrue(d['raw_residual_component_cube_not_used_as_group_norm'])
 def test_fail_closed_or_pass_is_explicit(self):
  d=self.d; st=d['P5_45DEG_SAMPLE1_ENTRY_HA_CERTIFICATE']
  self.assertIn(st,('PASS','NOT_ESTABLISHED'))
  if st=='PASS':
   self.assertEqual(E.validate(d),[])
   self.assertTrue(d['sample1_entry_inside_q8'])
   for m in ('H','A'):
    self.assertTrue(d['modes'][m]['complete'])
  else:
   self.assertTrue(d['failures'])
   self.assertTrue(any(d['modes'][m]['first_failure'] is not None for m in ('H','A')))

if __name__=='__main__': unittest.main()
