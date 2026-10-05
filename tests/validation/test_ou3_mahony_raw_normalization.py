"""Regression for the analytical bit-seed certificate, not a trajectory test."""
from fractions import Fraction as F
from pathlib import Path
import unittest
from tools.stability.ou3_theorem.mahony_raw_normalization import certificate, cubic_range


class RawNormalizationTests(unittest.TestCase):
    def test_exact_certificate_and_source_binding(self):
        c = certificate()
        lo, hi = map(F, c['raw_output_squared_norm_bounds'])
        self.assertGreater(lo, F(995, 1000))
        self.assertLess(hi, F(1001, 1000))
        self.assertFalse(c['all_time_pre_normalization_domain_verified'])
        source = (Path(__file__).resolve().parents[2]/'src/ahrs/Mahony_AHRS.h').read_text()
        self.assertIn('0x5f375a86', source)
        self.assertIn('y = y * (1.5f - (number * 0.5f * y * y))', source)

    def test_cubic_has_interior_maximum(self):
        self.assertEqual(cubic_range(F(1, 2), F(2), F(3), F(1)), (F(2), F(4)))


if __name__ == '__main__':
    unittest.main()
