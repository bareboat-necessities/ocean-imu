"""Independent finite-energy checks for the native trace analyzer."""
from fractions import Fraction as F
import tempfile
import unittest
from pathlib import Path
import numpy as np
from tools.stability.ou3_alt_contraction.finite_error_storage_diagnostic import analyze, nonlinear_budget

class FiniteEnergyTests(unittest.TestCase):
    def test_conditional_budget(self):
        rho, supply = nonlinear_budget(F(989,1000), F(1,1000), F(1,500),
                                       F(1,500), F(1,200), 1)
        self.assertLess(rho, 1)
        self.assertEqual(supply, F(100701,500))
        # A large remainder must not be reported contractive.
        self.assertGreater(nonlinear_budget(F(989,1000), F(1,10), 0, 0,
                                           F(1,200), 1)[0], 1)
        with self.assertRaises(ValueError):
            nonlinear_budget(1, 0, 1, 0, 1, 1)

    def trace(self, path, asymmetric=False):
        rows = np.zeros((1801, 464))
        rows[:, 0] = np.arange(1801)
        # Include held-bias/output energy; a motion-only projection would fail.
        rows[:, 2] = 2**(-np.arange(1801)/600)
        rows[:, 22] = rows[:, 2]
        P = np.eye(21)
        P[0, 20] = P[20, 0] = .5
        if asymmetric:
            P[20, 0] = 0
        rows[:, 23:] = P.reshape(-1)
        np.savetxt(path, rows, fmt='%.17g')

    def test_full_covariance_finite_energy(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'trace.txt'
            self.trace(path)
            result = analyze(path)
        self.assertAlmostEqual(float(result['endpoints'][0]['V_60digits']), 4/3)
        for word in result['words']:
            self.assertAlmostEqual(float(word['ratio_60digits']), .25)
            self.assertAlmostEqual(word['prefix_max_over_start'], 1)
        self.assertTrue(all(e['exact_covariance_SPD'] for e in result['endpoints']))

    def test_no_symmetry_repair(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)/'trace.txt'
            self.trace(path, asymmetric=True)
            with self.assertRaisesRegex(ValueError, 'symmetric'):
                analyze(path)

if __name__ == '__main__':
    unittest.main()
