from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tools.stability.ou3_theorem.moving_pivots import (
    certificate, gyro_transport_prefixes, six_pivot_floor_valid,
    two_group_six_column_floor,
)


class MovingPivotTests(unittest.TestCase):
    def test_chronological_reset_defect_is_not_dropped(self):
        before = gyro_transport_prefixes([(1, 1, 0, 0), (0, 2, 1, 0)])[-1]
        after = gyro_transport_prefixes([(0, 2, 1, 0), (1, 1, 0, 0)])[-1]
        self.assertEqual(before["gyro_singular_lower"], 0)
        self.assertEqual(after["gyro_singular_lower"], 1)

    def test_positive_one_step_floors_do_not_imply_chronological_floor(self):
        # A pi rotation about x cancels the transverse sum of two identity
        # injections. This is a scalar-bound counterexample, not a qualified
        # shipping prediction (whose angle is <.007).
        row = gyro_transport_prefixes([(1, 1, 0, 0), (1, 1, 2, 0)])[-1]
        self.assertEqual(row["gyro_singular_lower"], 0)

    def test_source_shaped_small_prediction_and_reset_bounds(self):
        row = gyro_transport_prefixes([(F(1, 200), 1, F(7, 1000), F(7, 400000)),
                                       (0, F(1001, 1000), F(1, 1000), 0),
                                       (F(1, 200), 1, F(7, 1000), F(7, 400000))])[-1]
        self.assertGreater(row["gyro_singular_lower"], F(9, 1000))
        with self.assertRaises(ValueError):
            gyro_transport_prefixes([(0, 1, 0, F(1, 10))])

    def test_full_six_column_defect_budget_is_fail_closed(self):
        passing = two_group_six_column_floor(1, F(1, 2), 1, F(1, 100))
        self.assertEqual(passing["six_column_singular_lower"], F(19, 100))
        failing = two_group_six_column_floor(1, F(1, 2), 1, F(1, 5))
        self.assertFalse(failing["six_columns_certified_conditionally"])
        with self.assertRaises(ValueError):
            two_group_six_column_floor(1, 0, 1, 0)

    def test_all_six_pivots_and_row_count_are_charged(self):
        self.assertTrue(six_pivot_floor_valid(F(1, 10), 100, F(1, 100)))
        self.assertFalse(six_pivot_floor_valid(F(1, 10), 101, F(1, 100)))
        with self.assertRaises(ValueError):
            six_pivot_floor_valid(1, 5, F(1, 100))

    def test_quiet_nominal_floor_does_not_resolve_physical_ambiguity(self):
        row = certificate()
        self.assertEqual(F(row["quiet_nominal_two_group_singular_floor"]), F(18, 127))
        self.assertEqual(F(row["quiet_nominal_all_six_pivot_floor"]), F(9, 254))
        self.assertFalse(row["quiet_nominal_floor_identifies_physical_attitude_and_bias"])
        from tools.stability.ou3_theorem.matrix_certificates import add, identity, ldlt, matmul, transpose
        g, field, dt = F("9.80665"), F(75), F(4, 125)
        # Literal skew row groups; cross blocks are retained in the audit.
        c = [[0, -g, 0], [g, 0, 0], [0, 0, 0],
             [0, 0, 0], [0, 0, -field], [0, field, 0]]
        o = [v + [F(0)] * 3 for v in c] + [v + [dt * x for x in v] for v in c]
        gram = matmul(transpose(o), o)
        _, pivots = ldlt(add(gram, identity(6), -F(18, 127)**2))
        self.assertTrue(all(x > 0 for x in pivots))

    def test_independent_exact_matrix_audit_does_not_promote_source(self):
        row = certificate()
        self.assertTrue(all(F(x) > 0 for x in row["supplied_example_exact_gram_residual_pivots"]))
        self.assertFalse(row["physical_span_supplies_uniform_row_and_defect_budget"])
        self.assertFalse(row["uniform_historical_AG_readout_action"])
        self.assertFalse(row["rho0_certified"])


if __name__ == "__main__":
    unittest.main()
