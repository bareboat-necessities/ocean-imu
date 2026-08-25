import importlib.util
import re
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval, matrix_point, symmetric_positive_definite_ldlt
from ou3_interval_linear_algebra import (
    IntervalPivotError,
    matrix_inverse_gauss_jordan,
)

spec = importlib.util.spec_from_file_location(
    "ou3_validated_kalman_interval",
    ROOT / "tools" / "ou3_validated_kalman_interval.py",
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def contains_matrix(A, values):
    return all(A[i][j].contains(values[i][j]) for i in range(len(values)) for j in range(len(values[0])))


class IntervalLinearAlgebraTests(unittest.TestCase):
    def test_diagonal_inverse_is_enclosed(self):
        inv = matrix_inverse_gauss_jordan(matrix_point([[2.0, 0.0], [0.0, 4.0]]))
        self.assertTrue(contains_matrix(inv, [[0.5, 0.0], [0.0, 0.25]]))

    def test_spd_two_by_two_inverse_is_enclosed(self):
        A = matrix_point([[2.0, 0.5], [0.5, 1.5]])
        inv = matrix_inverse_gauss_jordan(A)
        det = 2.0 * 1.5 - 0.25
        expected = [[1.5 / det, -0.5 / det], [-0.5 / det, 2.0 / det]]
        self.assertTrue(contains_matrix(inv, expected))

    def test_zero_crossing_pivot_is_rejected(self):
        A = [
            [Interval(-0.1, 0.1), Interval.point(0.0)],
            [Interval.point(0.0), Interval.point(1.0)],
        ]
        with self.assertRaises(IntervalPivotError):
            matrix_inverse_gauss_jordan(A)


class ValidatedKalmanIntervalTests(unittest.TestCase):
    def test_scalar_prediction(self):
        P = matrix_point([[2.0]])
        F = matrix_point([[1.0]])
        Q = matrix_point([[0.5]])
        pred = mod.covariance_predict(F, P, Q)
        self.assertTrue(pred[0][0].contains(2.5))

    def test_scalar_joseph_update_matches_closed_form(self):
        P = matrix_point([[2.0]])
        H = matrix_point([[1.0]])
        R = matrix_point([[3.0]])
        out = mod.joseph_measurement_update(P, H, R)
        self.assertTrue(out["S"][0][0].contains(5.0))
        self.assertTrue(out["K"][0][0].contains(0.4))
        self.assertTrue(out["A_correction"][0][0].contains(0.6))
        self.assertTrue(out["P_plus"][0][0].contains(1.2))

    def test_two_state_joseph_update_remains_certifiably_spd(self):
        P = matrix_point([[2.0, 0.2], [0.2, 1.0]])
        H = matrix_point([[1.0, 0.0]])
        R = matrix_point([[0.5]])
        out = mod.joseph_measurement_update(P, H, R)
        ok, pivots = symmetric_positive_definite_ldlt(out["P_plus"])
        self.assertTrue(ok, pivots)
        self.assertEqual(len(out["K"]), 2)
        self.assertEqual(len(out["K"][0]), 1)

    def test_narrow_interval_family_propagates_without_unvalidated_inverse(self):
        P = [
            [Interval.outward_bounds(1.9, 2.1), Interval.outward_bounds(0.09, 0.11)],
            [Interval.outward_bounds(0.09, 0.11), Interval.outward_bounds(0.95, 1.05)],
        ]
        H = matrix_point([[1.0, 0.0]])
        R = [[Interval.outward_bounds(0.49, 0.51)]]
        out = mod.joseph_measurement_update(P, H, R)
        self.assertTrue(out["S"][0][0].lo > 0.0)
        self.assertTrue(out["P_plus"][0][0].lo > 0.0)
        self.assertTrue(out["P_plus"][1][1].lo > 0.0)

    def test_no_unvalidated_linear_algebra_calls_in_proof_modules(self):
        for path in (
            ROOT / "tools" / "ou3_interval_linear_algebra.py",
            ROOT / "tools" / "ou3_validated_kalman_interval.py",
        ):
            text = path.read_text()
            # Match executable imports/calls rather than prose such as
            # "no NumPy eigensolver" in the module's audit documentation.
            forbidden_patterns = (
                r"^\s*import\s+numpy\b",
                r"^\s*from\s+numpy\b",
                r"\bnp\.linalg\b",
                r"\bnumpy\.linalg\b",
                r"\blinalg\.inv\s*\(",
                r"\bpinv\s*\(",
                r"\beig(?:vals?|envalues?)?\s*\(",
            )
            for pattern in forbidden_patterns:
                self.assertIsNone(re.search(pattern, text, flags=re.MULTILINE), pattern)


if __name__ == "__main__":
    unittest.main()
