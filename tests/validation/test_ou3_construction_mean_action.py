from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from tools.stability.ou3_theorem.construction_mean_action import (
    RationalInterval as I, correction, inverse_action, ldlt, verify_summary,
)
from tools.stability.ou3_theorem.matrix_certificates import is_psd


class ConstructionMeanActionTests(unittest.TestCase):
    def test_coupled_varying_gain_mean_gram_identity(self):
        d = [[F(0)]*6 for _ in range(6)]
        eta, u, energy = [F(0)]*6, [F(0)]*6, F(0)
        for step in range(3):
            phi = F(2, 3+step)
            f = [F(1)]*3+[phi]*3
            u = [f[i]*u[i] for i in range(6)]
            eta = [f[i]*eta[i] for i in range(6)]
            d = [[f[i]*d[i][j]*f[j] for j in range(6)] for i in range(6)]
            k = [[F((i+step+j*j)%5-2, 7) for j in range(3)] for i in range(6)]
            s = [[F(3), F(1), F(1, 2)], [F(1), F(2), F(1, 3)], [F(1, 2), F(1, 3), F(4)]]
            r = [F(1+step, 2), F(-2), F(3, 7)]
            after = [u[i]+sum(k[i][j]*r[j] for j in range(3))+F(i-2, 10**8) for i in range(6)]
            energy = correction(d, eta, energy, k, s, r, u, after)
            u = after
        v = [u[i]-eta[i] for i in range(6)]
        joint = [d[i]+[v[i]] for i in range(6)]+[v+[energy]]
        self.assertTrue(is_psd(joint))
        self.assertNotEqual(d[0][3], 0)
        self.assertNotEqual(eta[0], 0)
        # An overlarge correction cannot be hidden by a small roundoff label.
        bad = deepcopy(joint)
        bad[0][6] = bad[6][0] = F(100)
        self.assertFalse(is_psd(bad))

    def test_rational_interval_factor_encloses_full_inverse_action(self):
        a = [[I(F(2), F(21, 10)), I(F(1, 2), F(3, 5))],
             [I(F(1, 2), F(3, 5)), I(F(3), F(31, 10))]]
        cost = inverse_action(a, [I.cast(1), I.cast(2)], I.cast(0), I.cast(1))
        for x in (F(2), F(21, 10)):
            for y in (F(1, 2), F(3, 5)):
                for z in (F(3), F(31, 10)):
                    exact = inverse_action([[x, y], [y, z]], [F(1), F(2)])
                    self.assertLessEqual(cost.lo, exact)
                    self.assertGreaterEqual(cost.hi, exact)

    def test_indefinite_or_uncertain_factor_fails_closed(self):
        for a in ([[F(1), F(2)], [F(2), F(1)]], [[I(F(-1), F(1))]]):
            with self.assertRaisesRegex(ArithmeticError, 'positive factor pivot'):
                ldlt(a)
        with self.assertRaisesRegex(ArithmeticError, 'contains zero'):
            I.cast(1)/I(F(-1), F(1))

    def test_committed_gyro_barrier_and_failed_collinearity_margin(self):
        report = json.loads((ROOT/'reports/results/ou3_stability/construction-mean-action.json').read_text())
        actual = verify_summary(report['enclosure'], report['committed_field'])
        self.assertEqual(actual, report['exact_summary'])
        self.assertTrue(actual['finite_recorded_gyro_mean_norm_below_one_verified'])
        self.assertLess(F(actual['finite_recorded_angular_increment_upper']), F('.01'))
        self.assertLess(F(actual['collinearity_exclusion_margin_upper']), -817884)
        self.assertFalse(actual['source_uniform_verified'])
        self.assertFalse(actual['theorem_closed'])

    def test_gyro_barrier_rejects_increased_energy(self):
        report = json.loads((ROOT/'reports/results/ou3_stability/construction-mean-action.json').read_text())
        report['enclosure']['innovation_energy'] = ['10000000', '10000000']
        with self.assertRaisesRegex(ArithmeticError, 'positive factor pivot'):
            verify_summary(report['enclosure'], report['committed_field'])


if __name__ == '__main__':
    unittest.main()
