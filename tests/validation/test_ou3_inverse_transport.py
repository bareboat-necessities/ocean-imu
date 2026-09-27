"""Exact algebra audits; supplied reset sequences are not reachable witnesses."""
from fractions import Fraction as F
import unittest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tools.stability.ou3_theorem.lin_path_certificate import inverse
from tools.stability.ou3_theorem.matrix_certificates import add, identity, is_psd, matmul, transpose
from tools.stability.ou3_theorem.moving_pivots import (
    gyro_transport_prefixes, inverse_frame_gyro_prefixes,
)


def skew(v):
    x, y, z = map(F, v)
    return [[F(0), -z, y], [z, F(0), -x], [-y, x, F(0)]]


def scaled(a, c):
    return [[c * x for x in row] for row in a]


class InverseTransportTest(unittest.TestCase):
    def test_literal_reset_inverse_is_nonexpansive_with_large_injection(self):
        for d in ((0, 0, 0), (1, 2, -3), (1000, -1, 9)):
            g = add(identity(3), scaled(skew(d), F(1, 2)))
            gi = inverse(g)
            self.assertTrue(is_psd(add(matmul(transpose(g), g), identity(3), -1)))
            self.assertTrue(is_psd(add(identity(3), matmul(transpose(gi), gi), -1)))

    def test_both_small_rate_source_inverse_bounds_exactly(self):
        for h in (F('.004'), F('.006')):
            theta = h / 10**8
            wh = skew((0, 0, theta))
            wh2 = matmul(wh, wh)
            r = add(add(identity(3), wh, -1), wh2, F(1, 2))
            d = scaled(add(add(identity(3), wh, F(-1, 2)), wh2, F(1, 6)), h)
            ri = inverse(r)
            self.assertTrue(is_psd(add(identity(3), matmul(transpose(ri), ri), -1)))
            u = theta + theta**2 / 2
            v = h * (theta / 2 + theta**2 / 3)
            for error, bound in ((add(ri, identity(3), -1), u),
                                 (add(matmul(ri, d), scaled(identity(3), h), -1), v)):
                self.assertTrue(is_psd(add(scaled(identity(3), bound**2),
                                           matmul(transpose(error), error), -1)))

    def test_terminal_reset_cannot_consume_previous_gyro_floor(self):
        h = F(1, 200)
        new = inverse_frame_gyro_prefixes([('predict', h, 0), ('reset', 4)])
        old = gyro_transport_prefixes([(h, 1, 0, 0), (0, 3, 2, 0)])
        self.assertEqual(new[-1]['gyro_singular_lower'], h)
        self.assertEqual(new[0]['normalized_gyro_defect_upper'],
                         new[1]['normalized_gyro_defect_upper'])
        self.assertEqual(old[-1]['gyro_singular_lower'], 0)

    def test_every_prefix_keeps_noncommuting_reset_action(self):
        ops = [('predict', F('.004'), 0), ('reset', F('.1')),
               ('predict', F('.006'), 0), ('reset', F('.2')),
               ('predict', F('.005'), 0)]
        a, b = identity(3), scaled(identity(3), F(0))
        axis = 0
        for op, bound in zip(ops, inverse_frame_gyro_prefixes(ops)):
            if op[0] == 'predict':
                b = add(b, scaled(identity(3), op[1]))
            else:
                d = [F(0)] * 3
                d[axis] = op[1]
                axis += 1
                g = add(identity(3), scaled(skew(d), F(1, 2)))
                a, b = matmul(g, a), matmul(g, b)
            c = matmul(inverse(a), b)
            error = add(c, scaled(identity(3), bound['elapsed_time']), -1)
            e, floor = bound['normalized_gyro_defect_upper'], bound['gyro_singular_lower']
            self.assertTrue(is_psd(add(scaled(identity(3), e**2),
                                       matmul(transpose(error), error), -1)))
            self.assertTrue(is_psd(add(matmul(transpose(b), b),
                                       scaled(identity(3), floor**2), -1)))

    def test_unbounded_reset_action_can_cancel_interanchor_gyro_transport(self):
        # Qualified predictions, zero corrected rate, three literal resets.
        # This tests a relaxed transport family, NOT shipping reachability,
        # applied service, or the rank of all intermediate sensor rows.
        b = scaled(identity(3), F(3, 625))
        for d in (F(4), F(4), F(8, 3)):
            g = add(identity(3), scaled(skew((0, 0, d)), F(1, 2)))
            b = matmul(g, b)
        for _ in range(8):
            b = add(b, scaled(identity(3), F(1, 200)))
        self.assertEqual(b, [[F(0), F(0), F(0)], [F(0), F(0), F(0)],
                             [F(0), F(0), F(28, 625)]])
        ops = [('predict', F(3, 625), 0)] + [('reset', d) for d in (4, 4, F(8, 3))]
        ops += [('predict', F(1, 200), 0)] * 8
        self.assertEqual(inverse_frame_gyro_prefixes(ops)[-1]['gyro_singular_lower'], 0)

    def test_invalid_operation_cannot_disappear(self):
        for op in (('predict', 0, 0), ('predict', 1, -1), ('reset', -1), ('correction', 1)):
            with self.assertRaises(ValueError):
                inverse_frame_gyro_prefixes([op])


if __name__ == '__main__':
    unittest.main()
