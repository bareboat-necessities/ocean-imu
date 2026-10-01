import unittest
from tools.stability.ou3_theorem.displacement_excitation import (
 chord_moment_lower, velocity_chord_lower, zero_translation_excluded, certificate)

class DisplacementExcitationTests(unittest.TestCase):
 def test_zero_translation_is_excluded_once_PE_is_qualified_positive(self):
  self.assertTrue(zero_translation_excluded(.01))
 def test_mean_value_velocity_consequence(self):
  self.assertEqual(velocity_chord_lower(4,2),.5)
 def test_no_fake_acceleration_floor(self):
  self.assertEqual(chord_moment_lower(10,1,5.5),0)
  self.assertEqual(chord_moment_lower(1,8,5.5),2.5)
  self.assertFalse(certificate()["instantaneous_acceleration_floor_inferred"])
 def test_invalid(self):
  for x in (0,-1):
   with self.assertRaises(ValueError): zero_translation_excluded(x)
if __name__=="__main__": unittest.main()
