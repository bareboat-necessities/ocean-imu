from fractions import Fraction as F
from pathlib import Path
import tempfile
import unittest

from tools.stability.ou3_alt_contraction import finite_startup_trap as T


class StartupTrapTests(unittest.TestCase):
    def test_conditional_margin_does_not_promote_unaudited_majorants(self):
        result = T.certificate()
        self.assertGreater(result['invariant_radius_margin'], 0)
        self.assertLess(result['integral_increment_upper'],
                        result['smallest_integral_rounding_half_cell'])
        for key in ('conditional_observer_tail_invariant_closed',
                    'rounding_majorants_audit_closed',
                    'indefinite_wrapper_guard_composition_closed',
                    'eventual_finite_startup_refuted',
                    'source_uniform_rho_certified', 'storage_search_allowed',
                    'ALT_STARTUP_PASS', 'ALT_LIVE_PASS', 'ALT_END_TO_END_PASS'):
            self.assertFalse(result[key], key)

    def test_direction_is_exact_and_scale_invariant(self):
        self.assertEqual(T.dot(T.REFERENCE, T.REFERENCE), 1)
        self.assertEqual(T.down(tuple(3*x for x in T.ROOT_Q)), T.REFERENCE)
        with self.assertRaises(ValueError):
            T.down((F(0),)*4)

    def test_literal_tail_preserves_attained_integral_on_finite_prefix(self):
        state = T.V.State(q=T.ROOT_Q, integral=T.ROOT_I, initialized=True)
        for _ in range(12):
            gyro, acc = T.tail_packet(state)
            state = T.M.step_initialized(state, T.CFG, dt=T.DT,
                                         gyro=gyro, acc=acc).vertical.state
            self.assertEqual(state.integral, T.ROOT_I)
            delta = tuple(x-y for x, y in zip(T.down(state.q), T.REFERENCE))
            self.assertLess(T.dot(delta, delta), T.RADIUS**2)

    def test_incomplete_prefix_cannot_supply_reachability(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'prefix.txt'
            path.write_text('1 0x0p+0\n')
            with self.assertRaisesRegex(ValueError, 'complete reset'):
                T.verify_prefix(path)


if __name__ == '__main__':
    unittest.main()
