"""World-frame historical rows, injection budget and the same-cell witness."""
from fractions import Fraction as F
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tools.stability.ou3_theorem.matrix_certificates import add, identity, is_psd, matmul, transpose
from tools.stability.ou3_theorem.world_frame import (
    certificate, collinear_cadence_floor, collinear_same_cell_witness, injection_budget,
    nominal_attitude_column_floor,
    quaternion_rotation, reset_angle_upper, reset_discrepancy_gram, reset_distance_upper,
    same_cell_world_floor, skew, world_factorization)


def column(v):
    return [[F(x)] for x in v]


def body_group(rotation, injection_rotation, d, force, field):
    """Literal same-cell rows at attitude R: acc, reset, then mag at R_2."""
    r2 = matmul(injection_rotation, rotation)
    acc = [[-x for x in row] for row in skew([r[0] for r in matmul(rotation, column(force))])]
    mag = [[-x for x in row] for row in skew([r[0] for r in matmul(r2, column(field))])]
    return acc + matmul(mag, add(identity(3), skew(d), F(1, 2)))


class WorldFrameRowsTest(unittest.TestCase):
    def test_rational_word_factors_exactly_and_predictions_are_world_neutral(self):
        r0 = quaternion_rotation((2, -1, 1, 3))
        step = quaternion_rotation((30, 1, -1, 0))
        r1 = matmul(step, r0)
        bs = [[F(1, 200) if i == j else F(0) for j in range(3)] for i in range(3)]
        inj = quaternion_rotation((50, F(1, 3), 1, F(-1, 2)))
        r2 = matmul(inj, r1)
        word = [('row', [F(1), F(0), F(-9)]), ('predict', r1, step, bs),
                ('reset', r2, [F(1, 30), F(1, 10), F(-1, 50)]), ('row', [F(40), F(1), F(-7)])]
        out = world_factorization(r0, word)
        self.assertEqual(out['row_factorization_residual'], 0)
        self.assertEqual(out['factors'][0][1], identity(3))
        n, x = out['factors'][1][1], [row[0] for row in out['factors'][1][2]]
        self.assertEqual(matmul(transpose(n), n), reset_discrepancy_gram(x))
        self.assertTrue(is_psd(add(reset_discrepancy_gram(x), identity(3), F(-1))))

    def test_prediction_mismatch_is_reported_not_hidden(self):
        r0 = quaternion_rotation((1, 0, 0, 0))
        mean = quaternion_rotation((10, 0, 0, 1))
        covariance = quaternion_rotation((10, 0, 1, 0))
        out = world_factorization(r0, [('predict', mean, covariance, identity(3)),
                                       ('row', [F(0), F(0), F(-9)])])
        self.assertEqual(out['row_factorization_residual'], 0)
        self.assertNotEqual(out['factors'][0][1], identity(3))

    def test_same_cell_floor_is_attitude_invariant_and_valid(self):
        # Fixed WORLD injection x and world injection rotation; the body
        # injection d=R x and mean injection R W R' follow the attitude.
        force, field, x = [F(1, 2), F(-1, 3), F(-9)], [F(30), F(5), F(60)], [F(1, 100), 0, F(-1, 80)]
        world_inj = quaternion_rotation((200, F(1, 2), 0, F(-5, 8)))
        floor = same_cell_world_floor(force, field, F(1, 50))
        self.assertGreater(floor, 0)
        grams = []
        for q in ((1, 0, 0, 0), (3, -2, 1, 5), (1, 7, -4, 2)):
            rot = quaternion_rotation(q)
            d = [r[0] for r in matmul(rot, column(x))]
            inj = matmul(rot, matmul(world_inj, transpose(rot)))
            c = body_group(rot, inj, d, force, field)
            gram = matmul(transpose(c), c)
            grams.append(matmul(transpose(rot), matmul(gram, rot)))
            self.assertTrue(is_psd(add(gram, identity(3), -floor)))
        self.assertEqual(grams[0], grams[1])
        self.assertEqual(grams[1], grams[2])

    def test_collinear_same_cell_has_zero_floor_and_bad_inputs_fail(self):
        self.assertEqual(same_cell_world_floor([F(-7), 0, F(-24)], [F(21), 0, F(72)], 0), 0)
        with self.assertRaises(ValueError):
            same_cell_world_floor([0, 0, 0], [1, 0, 0], 0)
        with self.assertRaises(ValueError):
            reset_angle_upper(4)
        with self.assertRaises(ValueError):
            reset_distance_upper(F(-1))

    def test_reset_bounds_dominate_exact_values(self):
        import mpmath as mp
        with mp.workdps(40):
            for t in (F(1, 10**6), F(1, 100), F(1, 3), F(3, 2), F(39, 10)):
                th = mp.mpf(t.numerator)/t.denominator
                rho = 1/mp.sqrt(1+th**2/4)
                exact = th-mp.atan(th/2)+mp.atan((1-rho)/(2*mp.sqrt(rho)))
                bound = reset_angle_upper(t)
                self.assertGreaterEqual(mp.mpf(bound.numerator)/bound.denominator, exact)
                if t <= 2:
                    dist = mp.sqrt(2+th**2/4-2*mp.cos(th)-th*mp.sin(th))
                    bound = reset_distance_upper(t)
                    self.assertGreaterEqual(mp.mpf(bound.numerator)/bound.denominator, dist)

    def test_injection_budget_needs_the_joseph_structure_for_the_prior_bound(self):
        p = [[F(2), F(1, 2), F(0)], [F(1, 2), F(1), F(0)], [F(0), F(0), F(3)]]
        h = [[F(1), F(0), F(1)], [F(0), F(2), F(0)], [F(1), F(1), F(0)]]
        r = [[F(1), F(0), F(0)], [F(0), F(1, 2), F(0)], [F(0), F(0), F(2)]]
        s = add(matmul(h, matmul(p, transpose(h))), r)
        from tools.stability.ou3_theorem.lin_path_certificate import inverse
        k = matmul(p, matmul(transpose(h), inverse(s)))
        ok = injection_budget(k, s, [F(1), F(-2), F(3)], p)
        self.assertTrue(ok['gain_action_dominates'] and ok['prior_dominates_gain_action'])
        self.assertTrue(ok['prior_dominates_injection'])
        # An arbitrary gain keeps the Cauchy--Schwarz bound, not the prior one.
        bad = injection_budget([[F(10), 0, 0], [0, F(10), 0], [0, 0, F(10)]], s, [F(1), 0, 0], p)
        self.assertTrue(bad['gain_action_dominates'])
        self.assertFalse(bad['prior_dominates_gain_action'])

    def test_attitude_column_transfer_threshold(self):
        base = nominal_attitude_column_floor(16, 0, F(1, 5))
        self.assertEqual(base['aw_error_threshold'], F(112383, 100000))
        at = nominal_attitude_column_floor(16, base['aw_error_threshold'], F(1, 5))
        self.assertEqual(at['normalized_gram_floor_without_injection'], 0)
        self.assertFalse(at['positive'])
        below = nominal_attitude_column_floor(16, F(1), F(1, 5))
        self.assertTrue(below['positive'])
        with self.assertRaises(ValueError):
            nominal_attitude_column_floor(16, F(-1), F(1, 5))

    def test_collinear_witness_is_exact_but_not_a_service_claim(self):
        witness = collinear_same_cell_witness()
        for key in ('acceleration_upper', 'jerk_upper', 'velocity_upper', 'displacement_upper',
                    'primitive_upper', 'horizontal_fraction_admitted',
                    'integer_time_force_parallel_field', 'ideal_same_cell_kernel_along_field',
                    'half_integer_transverse_force_squared', 'marine_motion_and_imu_bias_admitted'):
            self.assertIs(witness[key], True)
        self.assertIs(witness['magnetic_service_admission_claimed'], False)

    def test_jerk_limited_collinear_cadence_is_tight(self):
        g, v, j, length = F('9.80665'), F('5.5'), F(100), F(16)
        floor = collinear_cadence_floor(F(1, 5), length)
        self.assertEqual(floor, F(127383, 2500000))
        # Tent dips c.a=|c|-J min(t-t_k,t_(k+1)-t) at uniform gap D have mean
        # |c|-J D/4; at the floor the velocity change is exactly 2V.
        c = g/5
        self.assertEqual((c-j*floor/4)*length, 2*v)
        self.assertGreater((c-j*F(1, 25)/4)*length, 2*v)
        self.assertEqual(collinear_cadence_floor(F(1, 5), F(1, 10)), 0)
        with self.assertRaises(ValueError):
            collinear_cadence_floor(0, 16)

    def test_committed_certificate_matches_exact_reproduction(self):
        path = ROOT/'reports/results/ou3_stability/world-frame-certificate.json'
        record = json.loads(path.read_text())
        self.assertEqual(record, certificate())
        self.assertFalse(record['same_cell_floor_from_motion_and_bias_bounds_alone'])
        self.assertFalse(record['same_cell_route_refuted_under_three_assumptions'])
        self.assertTrue(record['jerk_limited_collinear_cadence']['regular_25Hz_all_collinear_excluded'])
        self.assertFalse(record['nominal_attitude_column_transfer']['gyro_columns_certified'])
        self.assertFalse(record['theorem_closed'])

    def test_native_diagnostic_provenance_and_scope_fail_closed(self):
        from tools.stability.ou3_theorem.world_frame_source_diagnostic import verify_diagnostic
        record = json.loads((ROOT/'reports/results/ou3_stability/world-frame-source-feasibility.json').read_text())
        self.assertTrue(verify_diagnostic(record))
        collinear = record['cases'][2]
        self.assertLess(float(collinear['same_cell_actual_singular_max']), 1)
        self.assertGreater(float(collinear['aggregate_six_column_singular']), 10)
        # One applied correction per 1-s window: far below mu_M=1, so the
        # 1-Hz cadence is excluded by MAGNETIC SERVICE.
        self.assertLess(float(collinear['one_correction_window_service_lambda_min_max']), 1e-3)
        for key, value in [('driver_sha256', '0'*64), ('source_uniform_verified', True),
                           ('theorem_closed', True)]:
            changed = copy.deepcopy(record)
            changed[key] = value
            with self.assertRaises(ValueError):
                verify_diagnostic(changed)
        changed = copy.deepcopy(record)
        changed['cases'][1]['same_cell_floor_to_actual_ratio_max'] = '1.5'
        with self.assertRaises(ValueError):
            verify_diagnostic(changed)


if __name__ == '__main__':
    unittest.main()
