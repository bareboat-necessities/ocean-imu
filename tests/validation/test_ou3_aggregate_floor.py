"""Tube lemma and the injection-free aggregate six-column floor (Theorem G0)."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tools.stability.ou3_theorem.aggregate_floor import (  # noqa: E402
    OMEGA_INVARIANT, OMEGA_PHYSICAL, certificate, cos_lower, optimize_eps, sqrt_lower,
    sqrt_upper, theorem_g0, tube_constant)

KW = dict(window=16, separation=64, sigma_w=F(1, 5), m_perp=F(2, 5), force_l1=F(6, 5),
          rows_per_window=2666, start_offset=2)


class AggregateFloorTest(unittest.TestCase):
    def test_committed_certificate_matches_exact_reproduction(self):
        path = ROOT/'reports/results/ou3_stability/aggregate-floor-certificate.json'
        record = json.loads(path.read_text())
        self.assertEqual(record, certificate())
        for key in ('injection_transport_included', 'nominal_window_premises_proved',
                    'source_uniform_six_column_floor', 'rho0_certified', 'theorem_closed'):
            self.assertFalse(record[key])
        self.assertFalse(record['downstream']['practical_linear_margin'])
        for row in record['synthetic_falsification_audit'].values():
            self.assertLessEqual(row['G0_floor'], row['actual_sigma_min'])

    def test_rational_cosine_and_root_bounds(self):
        import mpmath as mp
        for x in (F(0), F(1, 3), F(1), OMEGA_INVARIANT, F(3, 2), F(3)):
            self.assertLessEqual(float(cos_lower(x)), float(mp.cos(float(x)))+1e-15)
        for v in (F(2), F(1, 7), F(10**6+3)):
            self.assertLessEqual(sqrt_lower(v)**2, v)
            self.assertGreaterEqual(sqrt_upper(v)**2, v)

    def test_tube_constant_needs_short_service_gaps(self):
        self.assertGreater(tube_constant(OMEGA_INVARIANT, 1, 0), F(4, 10))
        self.assertLess(tube_constant(OMEGA_INVARIANT, 2, 0), 0)
        self.assertGreater(tube_constant(OMEGA_PHYSICAL, 2, 0), 0)

    def test_theorem_g0_positive_and_monotone_in_premises(self):
        best = optimize_eps(OMEGA_INVARIANT, 1, **KW)
        self.assertIsNotNone(best)
        s2 = best[1]['six_column_floor_squared']
        self.assertGreater(s2, F(14, 10000))
        for m_perp in (F(9, 20), F(1)):
            worse = theorem_g0(OMEGA_INVARIANT, 1, best[0], **{**KW, 'm_perp': m_perp})
            self.assertLessEqual(worse.get('six_column_floor_squared', 0), s2)
        self.assertFalse(theorem_g0(OMEGA_INVARIANT, 2, F(1, 10), **KW)['positive'])
        self.assertFalse(theorem_g0(OMEGA_INVARIANT, 1, F(1, 10), **{**KW, 'm_perp': F(2)})['positive'])
        short = theorem_g0(OMEGA_INVARIANT, 1, F(41, 200), **{**KW, 'separation': 2})
        self.assertFalse(short['positive'])


if __name__ == '__main__':
    unittest.main()
