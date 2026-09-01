from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
import ou3_interval_ad as AD


class Ou3IntervalADTests(unittest.TestCase):
    def test_rational_arithmetic_derivative_contains_exact_point(self):
        x = AD.independent(Interval.point(2.0), 0, 1)
        y = (x * x + 3.0 * x) / (x + 1.0)
        # f=(x^2+3x)/(x+1); f'(2)=(7*3-10)/9=11/9.
        self.assertLessEqual(y.der[0].lo, 11.0 / 9.0)
        self.assertGreaterEqual(y.der[0].hi, 11.0 / 9.0)

    def test_cayley_rotation_is_identity_and_has_correct_local_generator(self):
        c = [AD.independent(Interval.point(0.0), i, 3) for i in range(3)]
        R = AD.rotation_from_cayley(c)
        for i in range(3):
            for j in range(3):
                self.assertEqual(R[i][j].val.lo, 1.0 if i == j else 0.0)
                self.assertEqual(R[i][j].val.hi, 1.0 if i == j else 0.0)
        # dR/dc_x at zero is [e_x]_x.
        self.assertLessEqual(R[1][2].der[0].lo, -1.0)
        self.assertGreaterEqual(R[1][2].der[0].hi, -1.0)
        self.assertLessEqual(R[2][1].der[0].lo, 1.0)
        self.assertGreaterEqual(R[2][1].der[0].hi, 1.0)

    def test_zero_deployed_correction_has_identity_cayley_derivative(self):
        n = 6
        c = [AD.independent(Interval.point(0.0), i, n) for i in range(3)]
        d = [AD.independent(Interval.point(0.0), 3 + i, n) for i in range(3)]
        cp = AD.deployed_correct_cayley(c, d)
        J = AD.jacobian(cp)
        for i in range(3):
            for j in range(3):
                self.assertLessEqual(J[i][j].lo, 1.0 if i == j else 0.0)
                self.assertGreaterEqual(J[i][j].hi, 1.0 if i == j else 0.0)
                self.assertLessEqual(J[i][3 + j].lo, 1.0 if i == j else 0.0)
                self.assertGreaterEqual(J[i][3 + j].hi, 1.0 if i == j else 0.0)

    def test_axis_branch_generalized_derivatives_are_finite_through_large_correction(self):
        n = 3
        d = [
            AD.independent(Interval.outward_bounds(3.18, 3.22), 0, n),
            AD.independent(Interval.outward_bounds(-0.002, 0.002), 1, n),
            AD.independent(Interval.outward_bounds(-0.002, 0.002), 2, n),
        ]
        w, v, dn = AD.deployed_quaternion_ad(d)
        self.assertGreater(dn, math.pi)
        for x in [w, *v]:
            self.assertTrue(math.isfinite(x.val.lo))
            self.assertTrue(math.isfinite(x.val.hi))
            for q in x.der:
                self.assertTrue(math.isfinite(q.lo))
                self.assertTrue(math.isfinite(q.hi))

    def test_interval_matrix_norm_upper_dominates_identity(self):
        z = Interval.point(0.0)
        o = Interval.point(1.0)
        A = [[o if i == j else z for j in range(4)] for i in range(4)]
        self.assertGreaterEqual(AD.interval_matrix_op2_upper(A), 1.0)
        self.assertLess(AD.interval_matrix_op2_upper(A), 1.0 + 1e-12)


if __name__ == "__main__":
    unittest.main()
