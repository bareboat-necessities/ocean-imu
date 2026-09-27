"""Exact chronological grouping, including nonzero and noncommuting resets."""
from fractions import Fraction as F
import copy
import json
from pathlib import Path
import unittest
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.stability.ou3_theorem.lin_path_certificate import inverse
from tools.stability.ou3_theorem.matrix_certificates import identity, matmul
from tools.stability.ou3_theorem.moving_pivots import same_prediction_cell_groups
from tools.stability.ou3_theorem.moving_pivots import same_cell_geometry_squared_floor


def prediction(h=F(1, 200)):
    f = identity(6)
    for i in range(3):
        f[i][i+3] = h
    return {'kind': 'prediction', 'F_AG': f}


def observation(sensor):
    c = ([[0, -9, 0], [9, 0, 0], [0, 0, 0]] if sensor == 'acc'
         else [[0, 0, 0], [0, 0, 10], [0, -10, 0]])
    return {'kind': 'correction', 'sensor': sensor, 'H': [r+[0]*18 for r in c]}


def reset(x, y, z):
    return {'kind': 'reset', 'd': [[x], [y], [z]]}


def skew(v):
    x, y, z = map(F, v)
    return [[F(0), -z, y], [z, F(0), -x], [-y, x, F(0)]]


def scaled(a, c):
    return [[c*x for x in row] for row in a]


class SameCellGroupsTest(unittest.TestCase):
    def test_reset_pulled_field_geometry_with_noncommuting_resets(self):
        from tools.stability.ou3_theorem.matrix_certificates import add, is_psd, transpose
        g = matmul(add(identity(3), scaled(skew((0, F(1, 5), 0)), F(1, 2))),
                   add(identity(3), scaled(skew((F(1, 10), 0, 0)), F(1, 2))))
        f, b = [[F(0)], [F(0)], [F(9)]], [[F(10)], [F(0)], [F(0)]]
        k = matmul(inverse(g), b)
        dot = lambda a, b: sum(x[0]*y[0] for x, y in zip(a, b))
        self.assertLess(dot(f, k)**2, dot(f, f)*dot(k, k)/100)
        c = skew([x[0] for x in f]) + matmul(skew([x[0] for x in b]), g)
        floor = same_cell_geometry_squared_floor(9, 10, F(1, 10))
        self.assertEqual(floor, F(729, 10))
        self.assertTrue(is_psd(add(matmul(transpose(c), c), scaled(identity(3), floor), -1)))

    def test_raw_force_field_nonparallel_is_insufficient_after_reset(self):
        from tools.stability.ou3_theorem.matrix_certificates import add
        # Relaxed row group only: no shipping-reachability assertion.
        g = add(identity(3), scaled(skew((0, 2, 0)), F(1, 2)))
        f, b = [1, 0, 1], [1, 0, 0]
        c = skew(f) + matmul(skew(b), g)
        self.assertEqual(matmul(c, [[F(x)] for x in f]), [[F(0)]]*6)
        self.assertEqual(same_cell_geometry_squared_floor(1, 1, 1), 0)
        for args in ((0, 1, 0), (1, 1, -1), (1, 1, 2)):
            with self.assertRaises(ValueError):
                same_cell_geometry_squared_floor(*args)

    def test_native_diagnostic_provenance_and_scope_fail_closed(self):
        from tools.stability.ou3_theorem.moving_transport_source_diagnostic import verify_diagnostic
        root = Path(__file__).resolve().parents[2]
        record = json.loads((root/'reports/results/ou3_stability/moving-transport-source-feasibility.json').read_text())
        self.assertTrue(verify_diagnostic(record))
        for key, value in [('shipping_header_sha256', '0'*64),
                           ('source_uniform_verified', True), ('theorem_closed', True)]:
            changed = copy.deepcopy(record)
            changed[key] = value
            with self.assertRaises(ValueError):
                verify_diagnostic(changed)

    def test_literal_noncommuting_resets_and_original_root(self):
        events = [prediction(), observation('acc'), reset(F(1, 3), 0, 0),
                  {'kind': 'correction', 'sensor': 'S'}, reset(0, F(1, 5), 0),
                  observation('mag'), reset(0, 0, F(1, 7)), prediction(),
                  {'kind': 'sync'}, observation('acc'), reset(0, 0, F(1, 11)),
                  observation('mag')]
        first, last = same_prediction_cell_groups(events)
        between = matmul(last['anchor_transport'], inverse(first['anchor_transport']))
        raw = first['raw_rows'] + last['raw_rows']
        expected = [r+[F(0)]*3 for r in first['C']] + matmul(last['C'], between[:3])
        self.assertEqual(matmul(raw, inverse(first['anchor_transport'])), expected)
        self.assertNotEqual(first['C'][3:], [r[:3] for r in observation('mag')['H']])
        self.assertEqual(first['acc_event'], 1)
        self.assertEqual(last['mag_event'], 11)

    def test_prediction_prevents_false_group(self):
        self.assertEqual(same_prediction_cell_groups([
            observation('acc'), prediction(), observation('mag')]), [])

    def test_actual_applied_rows_required(self):
        self.assertEqual(same_prediction_cell_groups([
            prediction(), observation('mag'), reset(0, 0, 1)]), [])
        groups = same_prediction_cell_groups([
            observation('acc'), observation('mag'), reset(0, 1, 0), observation('mag')])
        self.assertEqual(len(groups), 2)
        self.assertNotEqual(groups[0]['C'], groups[1]['C'])

    def test_direct_gyro_and_hard_event_fail_closed(self):
        acc = observation('acc')
        acc['H'][0][3] = 1
        for events in ([acc], [{'kind': 'reference_change'}],
                       [{'kind': 'correction', 'sensor': 'unknown'}]):
            with self.assertRaises(ValueError):
                same_prediction_cell_groups(events)
        f = prediction()
        f['F_AG'][5][5] = 2
        with self.assertRaises(ValueError):
            same_prediction_cell_groups([f])


if __name__ == '__main__':
    unittest.main()
