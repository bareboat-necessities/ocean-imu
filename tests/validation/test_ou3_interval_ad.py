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

    def test_square_preserves_nonnegative_range_across_zero(self):
        x = AD.independent(Interval.outward_bounds(-3.0, 2.0), 0, 1)
        y = x.square()
        self.assertGreaterEqual(y.val.lo, 0.0)
        self.assertLessEqual(y.val.lo, 0.0)
        self.assertGreaterEqual(y.val.hi, 9.0)
        self.assertLessEqual(y.der[0].lo, -6.0)
        self.assertGreaterEqual(y.der[0].hi, 4.0)

    def test_squared_norm_and_wide_cayley_denominator_do_not_spuriously_cross_zero(self):
        c = [
            AD.independent(Interval.outward_bounds(-5.0, 5.0), i, 3)
            for i in range(3)
        ]
        c2 = AD.squared_norm(c)
        self.assertGreaterEqual(c2.val.lo, 0.0)
        self.assertGreaterEqual(c2.val.hi, 75.0)
        # The exact inverse-Cayley denominator is 4+||c||^2 >= 4.  A broad
        # coordinate cell may widen the rotation entries, but it must not fail
        # merely because interval x*x invented a negative lower square bound.
        R = AD.rotation_from_cayley(c)
        for row in R:
            for x in row:
                self.assertTrue(math.isfinite(x.val.lo))
                self.assertTrue(math.isfinite(x.val.hi))

    def test_cayley_rotation_is_identity_and_has_correct_local_generator(self):
        c = [AD.independent(Interval.point(0.0), i, 3) for i in range(3)]
        R = AD.rotation_from_cayley(c)
        for i in range(3):
            for j in range(3):
                expected = 1.0 if i == j else 0.0
                # Proof arithmetic is deliberately outward rounded, so an exact
                # algebraic identity is tested by containment, not endpoint
                # equality of its binary64 enclosure.
                self.assertLessEqual(R[i][j].val.lo, expected)
                self.assertGreaterEqual(R[i][j].val.hi, expected)
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
