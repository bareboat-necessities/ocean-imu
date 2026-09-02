from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

from ou3_interval import Interval, symmetric_positive_definite_ldlt
import ou3_p3_word_process_floor as W


class P3WordProcessFloorTests(unittest.TestCase):
    def test_translation_information_is_positive_on_word_scale(self):
        for lo, hi in ((0.05, 0.051), (0.5, 0.505), (2.0, 2.01)):
            with self.subTest(x=(lo, hi)):
                info = W.translation_information_upper(Interval.outward_bounds(lo, hi))
                self.assertIsNotNone(info)
                self.assertTrue(symmetric_positive_definite_ldlt(info)[0])

    def test_translation_margin_is_strict(self):
        info = W.translation_information_upper(Interval.outward_bounds(0.5, 0.501))
        self.assertIsNotNone(info)
        delta = W.translation_margin_from_information(
            info,
            [2.0, 3.0, 4.0, 1.5],
            [0.0, 0.0, 0.1, 0.2],
        )
        self.assertTrue(math.isfinite(delta))
        self.assertGreater(delta, 0.0)

    def test_attitude_bias_doubling_floor_is_spd(self):
        Omega = W.attitude_bias_word_noise(0.5, 0.01, 7, 2.0)
        self.assertTrue(symmetric_positive_definite_ldlt(Omega)[0])
        Sigma = [[W.I(2.0), W.I(0.0)], [W.I(0.0), W.I(3.0)]]
        self.assertGreater(W.generalized_delta(Omega, Sigma), 0.0)

    def test_word_step_count_is_lower_rounded(self):
        k = W.word_step_doublings(1.0, 0.005)
        self.assertGreaterEqual(k, 0)
        self.assertLessEqual(2 ** k * 0.005, 1.0)
        self.assertGreater(2 ** (k + 1) * 0.005, 1.0)

    def test_word_series_cap_fails_closed(self):
        with self.assertRaises(ValueError):
            W.word_normalized_matrix(Interval.outward_bounds(2.5, 2.6))


if __name__ == "__main__":
    unittest.main()
