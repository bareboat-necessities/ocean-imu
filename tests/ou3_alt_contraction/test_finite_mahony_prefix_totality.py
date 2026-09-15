"""Check the source-uniform inequalities and the boundaries they protect."""
from dataclasses import replace
from fractions import Fraction as F
import unittest

from tools.stability.ou3_alt_contraction import finite_mahony_prefix_totality as X
from tools.stability.ou3_alt_contraction import finite_vertical_complementary_runtime as V


class Tests(unittest.TestCase):
    def test_both_commissioned_profiles_close_the_scalar_prefix(self):
        r = X.build()
        self.assertEqual(X.validate(r), [])
        for p in r['profiles'].values():
            self.assertEqual(p['max_steps'], X.CLOCK.MAX_STEPS)
            self.assertLess(p['ordinary_seed']['ordinary_seed_norm2_upper'], F(53,50))
            self.assertLess(p['feedback_component_abs_upper'], X.ERROR)
            self.assertLess(p['integral_component_abs_upper'], F(41,10))
            self.assertLess(p['Euler_norm_sum_upper'], 6)
            self.assertLess(p['vertical_abs_upper'], 32)
            self.assertTrue(p['ordinary_seed_and_scalar_state_prefix_totality_closed'])

    def test_first_seed_near_branch_boundary_retains_correlation(self):
        # These are graph checks of cancellation-sensitive branch endpoints;
        # the universal result is the rational inequality in the certificate.
        p = X.SENSOR.profile('MPU6886_COMMISSIONED_V1')
        cert = X.ordinary_seed_certificate(p)
        g = X.SEED.rn(X.SENSOR.G)
        for horizontal in (F(1,20), F(1,10), F(1), F(10)):
            acc = tuple(map(X.SEED.rn, (horizontal, F(0), g)))
            s = X.SEED.seed(acc)
            self.assertEqual(s.branch, 'ordinary-FromTwoVectors')
            self.assertLessEqual(sum(q*q for q in s.quaternion), cert['ordinary_seed_norm2_upper'])

    def test_largest_certified_integral_still_returns_to_shell(self):
        r = X.prefix_certificate(X.SENSOR.profile('MPU6886_COMMISSIONED_V1'))
        cfg = V.Config(X.SEED.rn(F(1,5)), X.SEED.rn(F(1,50)), X.SEED.rn(X.SENSOR.G), 20)
        # Four equal q components exercise the largest simultaneous products;
        # this is a robustness check, not a source-admission assertion.
        s = V.State(q=(F(1,2),)*4,
                    integral=(X.SEED.rn(F(4)),)*3,
                    initialized=True, elapsed=X.SEED.rn(F(150)))
        out = X.M.step_initialized(s,cfg,dt=X.CLOCK.DT_FLOAT,
            gyro=(X.SEED.rn(F(3,5)),)*3,
            acc=(X.SEED.rn(F(1,5)),X.SEED.rn(F(1,3)),-cfg.gravity))
        for op in out.operations:
            op.validate()
        self.assertLess(sum(q*q for q in out.vertical.state.q), X.Q2)
        self.assertLess(max(abs(i) for i in out.vertical.state.integral), F(41,10))
        self.assertLess(abs(out.vertical.state.up), r['vertical_abs_upper'])

    def test_no_promotion_of_solver_target_or_timeout(self):
        r = X.build()
        for key in ('every_startup_branch_totality_closed',
                    'startup_deployment_supply_bounds_closed', 'storage_search_allowed'):
            self.assertFalse(r[key])
        for p in r['profiles'].values():
            for key in ('startup_timeout_reachability_inferred_from_horizon',
                        'near_antiparallel_seed_totality_closed',
                        'unused_Euler_angle_library_calls_qualified',
                        'target_compiler_and_sqrt_qualified',
                        'indefinite_integral_invariant_closed'):
                self.assertFalse(p[key])
        r['every_startup_branch_totality_closed'] = True
        self.assertTrue(X.validate(r))


if __name__ == '__main__':
    unittest.main()
