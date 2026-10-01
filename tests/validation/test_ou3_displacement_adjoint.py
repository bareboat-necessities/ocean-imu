from fractions import Fraction as F
import unittest
from tools.stability.ou3_theorem.displacement_adjoint import (
 four_s_boundary_jets, displacement_chord_identity, chord_multiplier_jet,
 chord_integration_by_parts, augmented_multiplier_identity, certificate)

class DisplacementAdjointTests(unittest.TestCase):
 def test_four_s_exterior_jets_are_zero(self):
  j=four_s_boundary_jets((0,1,2,3))
  self.assertEqual(j["left"],(F(0),F(0),F(0)))
  self.assertEqual(j["right"],(F(0),F(0),F(0)))
 def test_chord_multiplier_exact_jets(self):
  self.assertEqual(chord_multiplier_jet(4,0),(F(4),F(-1),F(0)))
  self.assertEqual(chord_multiplier_jet(4,4),(F(0),F(-1),F(0)))
 def test_integration_by_parts_boundary(self):
  r=chord_integration_by_parts(4)
  self.assertEqual(r["boundary_velocity_coefficient"],-4)
  self.assertEqual(r["boundary_position_coefficients"],(F(-1),F(1)))
  self.assertEqual(r["interior_atoms"],())
 def test_four_s_superposition_does_not_change_chord_boundary(self):
  for a in (-3,0,7):
   r=augmented_multiplier_identity(5,a)
   self.assertEqual(r["boundary_velocity_coefficient"],-5)
   self.assertEqual(r["boundary_position_coefficients"],(F(-1),F(1)))
   self.assertTrue(r["boundary_coefficients_independent_of_four_S_scale"])
 def test_status(self):
  r=certificate()
  self.assertTrue(r["boundary_velocity_completion_identity_proved"])
  self.assertFalse(r["AG_reader_can_supply_velocity_boundary"])
  self.assertFalse(r["boundary_velocity_numerical_action_bound"])
if __name__=="__main__": unittest.main()
