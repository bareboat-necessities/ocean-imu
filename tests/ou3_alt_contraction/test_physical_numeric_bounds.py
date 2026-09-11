import unittest
from tools.stability.ou3_alt_contraction import physical_numeric_bounds as P


class PhysicalNumericBoundsTest(unittest.TestCase):
    def test_repaired_physics_supplies_alt_centered_S_bound(self):
        d=P.build(); self.assertEqual(P.validate(d),[])
        self.assertTrue(d['numeric_complete_family_physical_envelope_closed'])
        self.assertEqual(d['centered_S_norm_upper_m_s'],1100.0)
        self.assertFalse(d['legacy_300m_s_used_as_physical_bound'])
        self.assertFalse(d['D_S_comes_from_P4_working_radius'])
        self.assertTrue(d['usable_by_ALT_finite_event_outer_enclosure'])


if __name__=='__main__': unittest.main()
