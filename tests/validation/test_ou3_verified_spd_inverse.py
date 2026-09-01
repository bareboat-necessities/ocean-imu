from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan
import ou3_verified_spd_inverse as V


CERT = {"symmetric_certified": True, "spd_certified": True}


class VerifiedSPDInverseTests(unittest.TestCase):
    def test_diagonal_spd_family_is_enclosed(self):
        S = [
            [Interval(1.9, 2.1), Interval.point(0), Interval.point(0)],
            [Interval.point(0), Interval(2.9, 3.1), Interval.point(0)],
            [Interval.point(0), Interval.point(0), Interval(3.9, 4.1)],
        ]
        X, meta = V.inverse_enclosure(S, **CERT)
        self.assertTrue(meta["criterion_strict"])
        self.assertTrue(meta["source_symmetry_certified"])
        self.assertTrue(meta["source_SPD_certified"])
        self.assertLess(meta["neumann_q_inf_upper"], 1.0)
        for i, (lo, hi) in enumerate(((1.9, 2.1), (2.9, 3.1), (3.9, 4.1))):
            self.assertTrue(X[i][i].contains(1.0 / lo))
            self.assertTrue(X[i][i].contains(1.0 / hi))

    def test_point_spd_inverse_is_contained(self):
        S = [
            [Interval.point(2.0), Interval.point(0.2), Interval.point(0.1)],
            [Interval.point(0.2), Interval.point(1.5), Interval.point(0.05)],
            [Interval.point(0.1), Interval.point(0.05), Interval.point(1.0)],
        ]
        X, _ = V.inverse_enclosure(S, **CERT)
        Y = matrix_inverse_gauss_jordan(S)
        for i in range(3):
            for j in range(3):
                y = Y[i][j].lo + 0.5 * (Y[i][j].hi - Y[i][j].lo)
                self.assertTrue(X[i][j].contains(y))
                self.assertEqual(X[i][j].lo, X[j][i].lo)
                self.assertEqual(X[i][j].hi, X[j][i].hi)

    def test_missing_source_certificates_fail_closed(self):
        S = [[Interval.point(1.0)]]
        with self.assertRaisesRegex(V.VerifiedInverseFailure, "SPD certificate"):
            V.inverse_enclosure(S)
        with self.assertRaisesRegex(V.VerifiedInverseFailure, "symmetry certificate"):
            V.inverse_enclosure(S, spd_certified=True)

    def test_certified_symmetry_rejects_disjoint_paired_boxes(self):
        S = [
            [Interval.point(2.0), Interval(0.1, 0.2)],
            [Interval(0.3, 0.4), Interval.point(2.0)],
        ]
        with self.assertRaisesRegex(V.VerifiedInverseFailure, "paired enclosures became disjoint"):
            V.inverse_enclosure(S, **CERT)

    def test_noncontractive_preconditioned_family_fails_closed(self):
        S = [
            [Interval(-1.0, 3.0), Interval.point(0), Interval.point(0)],
            [Interval.point(0), Interval.point(1.0), Interval.point(0)],
            [Interval.point(0), Interval.point(0), Interval.point(1.0)],
        ]
        with self.assertRaises(V.VerifiedInverseFailure):
            V.inverse_enclosure(S, **CERT)


if __name__ == "__main__":
    unittest.main()
