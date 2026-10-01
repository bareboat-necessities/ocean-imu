import unittest
from fractions import Fraction as F
from tools.stability.ou3_theorem.aw_reader_uniform_enclosure import certificate, required_factor_radius

class T(unittest.TestCase):
    def test_exact_remaining_uniform_radius(self):
        r=certificate()
        self.assertTrue(r["source_uniform_finiteness_closed"])
        self.assertEqual(required_factor_radius(),F("0.52439539501604595"))
        self.assertFalse(r["composed_Lipschitz_bound_certified"])
        self.assertFalse(r["reachable_factor_cover_radius_certified"])
        self.assertFalse(r["theorem_closed"])

if __name__=="__main__":
    unittest.main()
