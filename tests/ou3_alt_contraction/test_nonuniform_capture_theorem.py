from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import nonuniform_capture_theorem as T


class NonuniformCaptureTests(unittest.TestCase):
    def test_product_sum_identity_including_amplifying_words(self):
        rhos, supplies = [F(2), F(1, 4), F(3, 4)], [F(1), F(2), F(3)]
        r = T.finite_comparison(5, rhos, supplies)
        # Explicit expansion, independent of the recurrence implementation.
        self.assertEqual(r['storage_upper_bound'], F(5)*2/F(4)*F(3, 4)
                         + F(1)/4*F(3, 4) + F(2)*F(3, 4) + 3)
        self.assertFalse(r['infinite_tail_certified'])

    def test_every_factor_strict_but_product_stays_above_half(self):
        for n in (0, 1, 2, 17, 64):
            rhos = [1-F(1, (k+2)**2) for k in range(n)]
            actual = T.finite_comparison(1, rhos, [0]*n)
            self.assertEqual(actual['storage_upper_bound'], F(n+2, 2*(n+1)))
            self.assertGreater(actual['storage_upper_bound'], F(1, 2))

    def test_vanishing_product_does_not_control_constant_noise(self):
        for n in (1, 2, 17, 64):
            rhos = [F(k+1, k+2) for k in range(n)]
            actual = T.finite_comparison(0, rhos, [1]*n)
            self.assertEqual(actual['initial_multiplier'], F(1, n+1))
            self.assertEqual(actual['forced_response'], F(n*(n+3), 2*(n+1)))

    def test_supply_scaled_to_decrement_gives_telescoping_bound(self):
        n = 31
        rhos = [F(k+1, k+2) for k in range(n)]
        r = T.normalized_supply_bound(10, rhos, 3)
        self.assertEqual(r['storage_upper_bound'], 3+F(7, n+1))
        self.assertFalse(r['product_tends_to_zero_proved'])

    def test_arbitrarily_slow_history_rate_has_finite_own_bound(self):
        r = T.fixed_history_bound(2, F(999, 1000), 1, 13)
        self.assertEqual(r['conditional_ultimate_storage'], 1000)
        self.assertFalse(r['history_uniform_rate_claimed'])
        self.assertEqual(r['storage_upper_bound'],
                         T.finite_comparison(2, [F(999, 1000)]*13, [1]*13)['storage_upper_bound'])

    def test_decaying_metric_does_not_prove_state_decay(self):
        # e_n=1, M_n=2^-n: V halves forever, while the physical error is fixed.
        for n in (1, 9, 31):
            metric, error = F(1, 2**n), F(1)
            storage = T.fixed_history_bound(1, F(1, 2), 0, n)['storage_upper_bound']
            self.assertEqual(metric*error**2, storage)
            self.assertEqual(error, 1)

    def test_unknown_capture_does_not_close_shipping_gates(self):
        status = T.selection_status()
        self.assertFalse(status['common_capture_deadline_required'])
        self.assertTrue(status['capture_is_an_explicit_assumption'])
        for key in ('finite_capture_proved_for_shipping_filter',
                    'shipping_nonlinear_word_inequality_proved',
                    'history_wise_tail_coercivity_proved',
                    'history_wise_prefix_retention_proved',
                    'uniform_ISS_claimed', 'lyapunov_stability_from_startup_claimed',
                    'storage_search_allowed', 'ALT_STARTUP_PASS', 'ALT_LIVE_PASS',
                    'ALT_END_TO_END_PASS'):
            self.assertFalse(status[key], key)

    def test_invalid_inputs_and_float_promotion_rejected(self):
        with self.assertRaises(TypeError):
            T.fixed_history_bound(1, .99, 1, 5)
        with self.assertRaises(ValueError):
            T.fixed_history_bound(1, 1, 1, 5)
        with self.assertRaises(ValueError):
            T.finite_comparison(1, [F(1, 2)], [])
        with self.assertRaises(ValueError):
            T.normalized_supply_bound(1, [2], 1)


if __name__ == '__main__':
    unittest.main()
