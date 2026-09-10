import pathlib
import sys
import unittest

ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools'/'stability'))

import ou3_p4_acceleration_moment_iqc_sector as S


class AccelerationMomentSectorTests(unittest.TestCase):
    def test_sector_status_and_boundary(self):
        d=S.build()
        self.assertEqual(S.validate(d),[])
        self.assertTrue(d['moment_sector_matrix_available'])
        self.assertTrue(d['aligned_constant_acceleration_boundary_contains_zero'])
        self.assertTrue(d['detached_box_corner_strictly_excluded'])
        self.assertTrue(d['zero_radial_zero_moment_boundary_contains_zero'])
        self.assertTrue(d['same_radial_coordinate_must_parameterize_entire_source_history'])
        self.assertFalse(d['independent_source_supply_port_used'])
        self.assertFalse(d['production_source_map_bound_here'])
        self.assertFalse(d['P4_PASS'])

    def test_bad_dimensions_are_rejected(self):
        z=S.Interval.point(0.0)
        with self.assertRaises(ValueError):
            S.moment_sector([[z]],[[z]],8.0)

    def test_false_promotion_is_rejected(self):
        d=S.build(); d['P4_PASS']=True
        self.assertIn('P4_PASS not false',S.validate(d))


if __name__=='__main__': unittest.main()
