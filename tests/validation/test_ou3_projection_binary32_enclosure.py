from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

TOOLS=Path(__file__).resolve().parents[2]/'tools'/'stability'
sys.path.insert(0,str(TOOLS))

import ou3_p4_projection_binary32_enclosure as P


class ProjectionBinary32Tests(unittest.TestCase):
    def test_enclosure_is_small_and_fail_scoped(self):
        d=P.build()
        self.assertEqual(P.validate(d),[])
        self.assertTrue(d['projection_finite_precision_enclosure_closed'])
        self.assertTrue(d['branch_ambiguity_covered'])
        self.assertTrue(d['FMA_contraction_covered_by_overcounted_rounding_model'])
        self.assertLess(d['single_evaluation_projection_value_error_norm_upper_mps2'],1e-5)
        self.assertFalse(d['full_Kalman_reset_finite_precision_enclosure_closed_here'])
        self.assertFalse(d['P4_MOTION_PASS'])
        self.assertFalse(d['P4_PASS'])

    def test_binary32_radius_literal_is_exactly_runtime_literal(self):
        d=P.build()
        self.assertEqual(d['shipping_radius_literal_binary32'],P.f32(.4))
        self.assertGreater(d['radius_literal_abs_error'],0.0)

    def test_rounding_gamma_is_monotone(self):
        self.assertEqual(P.gamma(0),math.nextafter(0.0,math.inf))
        self.assertLess(P.gamma(1),P.gamma(2))
        self.assertLess(P.gamma(2),P.gamma(6))


if __name__=='__main__':
    unittest.main()
