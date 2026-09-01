from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval, matrix_add, matrix_mul, matrix_transpose
from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan
import ou3_correlated_kalman_gain as CG


def iv(x):
    return Interval.point(float(x))


class CorrelatedKalmanGainTests(unittest.TestCase):
    def test_point_gain_contains_direct_point_solution(self):
        P = [[iv(2.0), iv(0.3)], [iv(0.3), iv(1.0)]]
        H = [[iv(1.0), iv(0.2)]]
        R = [[iv(0.5)]]
        d = CG.gain_enclosure(P, H, R, P_psd_certified=True)
        PHt = matrix_mul(P, matrix_transpose(H))
        S = matrix_add(matrix_mul(H, PHt), R)
        Kdirect = matrix_mul(PHt, matrix_inverse_gauss_jordan(S))
        for i in range(2):
            target = 0.5 * (Kdirect[i][0].lo + Kdirect[i][0].hi)
            self.assertTrue(d["K"][i][0].contains(target))
        self.assertTrue(d["source_P_PSD_certified"])
        self.assertTrue(d["S_ge_R_premise_from_P_PSD"])
        self.assertTrue(d["KSKt_le_P_premise_from_P_PSD"])
        self.assertFalse(d["interval_S_inverse_formed"])
        self.assertTrue(d["correlated_gain_equation_used"])
        self.assertGreater(d["R_eigenvalue_lower"], 0.0)

    def test_interval_family_produces_finite_gain_enclosure(self):
        P = [
            [Interval(1.9, 2.1), Interval(0.25, 0.35)],
            [Interval(0.25, 0.35), Interval(0.9, 1.1)],
        ]
        H = [[Interval(0.98, 1.02), Interval(0.18, 0.22)]]
        R = [[Interval(0.49, 0.51)]]
        d = CG.gain_enclosure(P, H, R, P_psd_certified=True)
        self.assertEqual(len(d["K"]), 2)
        self.assertEqual(len(d["K"][0]), 1)
        self.assertTrue(all(x >= 0.0 for x in d["row_gain_radius_upper"]))
        self.assertTrue(all(x < float("inf") for x in d["row_gain_radius_upper"]))
        for i in range(2):
            z = d["residual_intersection"][i][0]
            self.assertTrue(d["residual_direct"][i][0].contains_interval(z))
            self.assertTrue(d["residual_factored"][i][0].contains_interval(z))

    def test_joseph_inverse_identity_uses_same_gain_enclosure(self):
        P = [[iv(1.4), iv(0.1)], [iv(0.1), iv(0.8)]]
        H = [[iv(1.0), iv(0.25)]]
        R = [[iv(0.4)]]
        d = CG.gain_enclosure(P, H, R, P_psd_certified=True)
        Sinv_from_gain = CG.joseph_s_inverse_from_gain(H, R, d["K"])
        PHt = matrix_mul(P, matrix_transpose(H))
        S = matrix_add(matrix_mul(H, PHt), R)
        direct = matrix_inverse_gauss_jordan(S)
        target = 0.5 * (direct[0][0].lo + direct[0][0].hi)
        self.assertTrue(Sinv_from_gain[0][0].contains(target))

    def test_missing_P_PSD_source_certificate_fails_closed(self):
        P = [[iv(1.0)]]
        H = [[iv(1.0)]]
        R = [[iv(0.5)]]
        with self.assertRaisesRegex(CG.CorrelatedGainFailure, "PSD certificate"):
            CG.gain_enclosure(P, H, R)

    def test_nonpositive_R_floor_fails_closed(self):
        P = [[iv(1.0)]]
        H = [[iv(1.0)]]
        R = [[Interval(-0.1, 0.2)]]
        with self.assertRaises(CG.CorrelatedGainFailure):
            CG.gain_enclosure(P, H, R, P_psd_certified=True)


if __name__ == "__main__":
    unittest.main()
