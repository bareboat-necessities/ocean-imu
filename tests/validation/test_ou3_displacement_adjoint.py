from fractions import Fraction as F
import unittest
from tools.stability.ou3_theorem.displacement_adjoint import (
 four_s_boundary_jets, displacement_chord_identity, certificate)

class DisplacementAdjointTests(unittest.TestCase):
 def test_four_s_exterior_jets_are_zero(self):
  j=four_s_boundary_jets((0,1,2,3))
  self.assertEqual(j["left"],(F(0),F(0),F(0)))
  self.assertEqual(j["right"],(F(0),F(0),F(0)))
 def test_displacement_chord_has_nonzero_velocity_boundary(self):
  r=displacement_chord_identity(4,2)
  self.assertEqual(r["left_velocity_coefficient"],4)
  self.assertEqual(r["chord"],2)
 def test_fail_closed_conclusion(self):
  r=certificate()
  self.assertFalse(r["displacement_chord_in_standard_four_S_row_space"])
  self.assertFalse(r["source_uniform_outer_compatibility_closed"])
if __name__=="__main__": unittest.main()
