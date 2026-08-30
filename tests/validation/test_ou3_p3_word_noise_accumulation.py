from pathlib import Path
import math
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p3_word_noise_accumulation as WORDNOISE
import ou3_source_reachable_matrix_p3_direct as DIRECT
from ou3_interval import Interval, symmetric_positive_definite_ldlt


def _point(value):
    return Interval.outward_bounds(float(value), float(value))


def _identity(n):
    return [[_point(1.0 if i == j else 0.0) for j in range(n)] for i in range(n)]


class Ou3P3WordNoiseAccumulationTests(unittest.TestCase):
    def test_doubling_reproduces_the_explicit_prediction_sum(self):
        # Omega_N = sum_k Phi^k Q (Phi^k)'.  Check the doubling identity against
        # the sum written out term by term.
        Phi = [[_point(1.0), _point(-0.25)], [_point(0.0), _point(1.0)]]
        Q = [[_point(2.0), _point(0.5)], [_point(0.5), _point(3.0)]]
        for doublings in range(4):
            got = WORDNOISE.accumulate_word_noise(Q, Phi, doublings)
            power = _identity(2)
            total = [[_point(0.0), _point(0.0)], [_point(0.0), _point(0.0)]]
            for _ in range(2 ** doublings):
                term = [
                    [
                        sum(
                            (power[i][a] * Q[a][b] * power[j][b] for a in range(2) for b in range(2)),
                            _point(0.0),
                        )
                        for j in range(2)
                    ]
                    for i in range(2)
                ]
                total = [[total[i][j] + term[i][j] for j in range(2)] for i in range(2)]
                power = [
                    [sum((power[i][a] * Phi[a][j] for a in range(2)), _point(0.0)) for j in range(2)]
                    for i in range(2)
                ]
            for i in range(2):
                for j in range(2):
                    self.assertAlmostEqual(
                        0.5 * (got[i][j].lo + got[i][j].hi),
                        0.5 * (total[i][j].lo + total[i][j].hi),
                        delta=1.0e-9,
                    )

    def test_accumulation_is_monotone_in_the_horizon(self):
        Phi = [[_point(1.0), _point(-0.1)], [_point(0.0), _point(1.0)]]
        Q = [[_point(1.0), _point(0.0)], [_point(0.0), _point(1.0)]]
        previous = None
        for doublings in range(5):
            got = WORDNOISE.accumulate_word_noise(Q, Phi, doublings)
            trace = got[0][0].lo + got[1][1].lo
            if previous is not None:
                self.assertGreater(trace, previous)
            previous = trace

    def test_rebalancing_keeps_the_recursion_definite(self):
        # Without rebalancing the S row grows like N^7 and the a_w row like N,
        # so the enclosure loses definiteness long before the word horizon.
        Phi = [
            [_point(1.0), _point(0.0), _point(0.0), _point(1.0)],
            [_point(1.0), _point(1.0), _point(0.0), _point(0.5)],
            [_point(0.5), _point(1.0), _point(1.0), _point(1.0 / 6.0)],
            [_point(0.0), _point(0.0), _point(0.0), _point(0.999)],
        ]
        Q = [[Interval.outward_bounds(0.0, 0.0) for _ in range(4)] for _ in range(4)]
        for i in range(4):
            Q[i][i] = Interval.outward_bounds(1.0, 1.0000001)
        raw = WORDNOISE.accumulate_word_noise(Q, Phi, 9)
        balanced = WORDNOISE.accumulate_word_noise(Q, Phi, 9, [0.5, 0.25, 0.125, 1.0])
        raw_spread = max(abs(raw[i][i].hi) for i in range(4)) / min(
            abs(raw[i][i].lo) for i in range(4)
        )
        balanced_spread = max(abs(balanced[i][i].hi) for i in range(4)) / min(
            abs(balanced[i][i].lo) for i in range(4)
        )
        self.assertGreater(raw_spread, 1.0e6)
        self.assertLess(balanced_spread, 1.0e4)

    def test_measurement_posterior_never_exceeds_the_prior_floor(self):
        Omega = [[_point(2.0), _point(0.3)], [_point(0.3), _point(1.0)]]
        post = WORDNOISE.measurement_posterior(Omega, [0.0, 0.0])
        for i in range(2):
            for j in range(2):
                self.assertAlmostEqual(post[i][j].lo, Omega[i][j].lo, delta=1.0e-12)
        shrunk = WORDNOISE.measurement_posterior(Omega, [5.0, 5.0])
        for i in range(2):
            self.assertLess(shrunk[i][i].hi, Omega[i][i].lo)
            self.assertGreater(shrunk[i][i].lo, 0.0)

    def test_word_step_count_uses_a_lower_bound_on_the_horizon(self):
        self.assertEqual(WORDNOISE.word_step_doublings(1.0, 0.005), 7)
        self.assertEqual(WORDNOISE.word_step_doublings(3.14, 0.005), 9)
        with self.assertRaises(RuntimeError):
            WORDNOISE.word_step_doublings(0.001, 0.005)

    def test_attitude_bias_floor_is_positive_and_bias_keeps_its_accumulation(self):
        # No shipping update measures the gyro bias directly, so the bias
        # coordinate must keep the whole accumulated floor while the attitude
        # coordinate is shrunk by the vector information.
        block = WORDNOISE.attitude_bias_word_noise(0.9999, 1.0e-4, 9, 0.066)
        self.assertTrue(symmetric_positive_definite_ldlt(block)[0])
        self.assertGreater(block[1][1].lo, 400.0)
        self.assertLess(block[0][0].hi, block[1][1].lo)
        self.assertGreater(block[0][0].lo, 0.0)

    def test_normalised_sigma_margin_agrees_with_the_explicit_posterior(self):
        # (Omega^-1 + D)^-1 >= delta Sigma  and  lambda_max(G (Omega^-1 + D) G) <= 1/delta
        # with G = Sigma^{1/2} are the same statement.  The second is the one
        # that survives a strongly measured word, because it never forms the
        # posterior and never tests definiteness of a near-singular Sigma^-1.
        from ou3_interval_linear_algebra import matrix_inverse_gauss_jordan
        Omega = [
            [_point(0.80), _point(0.10), _point(0.02), _point(0.01)],
            [_point(0.10), _point(0.50), _point(0.03), _point(0.00)],
            [_point(0.02), _point(0.03), _point(0.40), _point(0.05)],
            [_point(0.01), _point(0.00), _point(0.05), _point(0.60)],
        ]
        sigma_diag = [40.0, 90.0, 25.0, 10.0]
        info = [0.0, 0.0, 2.0, 7.0]
        Sigma = [[_point(sigma_diag[i] if i == j else 0.0) for j in range(4)] for i in range(4)]
        post = WORDNOISE.measurement_posterior(Omega, info)
        explicit = DIRECT._certified_generalized_delta(post, Sigma, 1.0e-18)

        inverse = matrix_inverse_gauss_jordan([list(r) for r in Omega])
        sigma_root = [math.sqrt(v) for v in sigma_diag]
        normalised = DIRECT._word_translation_block_margin(
            inverse, sigma_root, [info[i] * sigma_diag[i] for i in range(4)]
        )
        self.assertGreater(explicit, 0.0)
        self.assertGreater(normalised, 0.0)
        # The normalised route bounds lambda_max by Gershgorin, so it can only
        # under-report, and on a 4x4 by at most the row count.
        self.assertLessEqual(normalised, explicit * (1.0 + 1.0e-9))
        self.assertGreater(normalised, explicit / 4.0)


if __name__ == "__main__":
    unittest.main()
