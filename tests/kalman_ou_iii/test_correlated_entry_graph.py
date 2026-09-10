"""Exact tests for the fresh-entry/source graph, not sampled stability tests."""
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools" / "stability"))
import ou3_p4_correlated_entry_graph as graph


class CorrelatedEntryGraphTest(unittest.TestCase):
    def test_all_polynomial_coefficients(self):
        evidence = graph.certificate()
        self.assertTrue(evidence["identity_certificate_pass"])
        self.assertEqual(evidence["gain_polynomial_coefficients_checked"], 36)
        self.assertEqual(evidence["prediction_polynomial_coefficients_checked"], 144)

    def test_nonzero_true_S_is_not_nonzero_innovation_at_fresh_entry(self):
        # A 300 m*s true primitive is illustrative algebra, not a BRMM admission.
        root = [Q(0)] * 12
        root[6] = Q(300)
        joint = graph.multiply(graph.fresh_entry_lift(), graph.matrix([[x] for x in root]))
        self.assertEqual(joint[6][0], Q(300))  # actual e_S is NOT zero
        self.assertEqual(joint[18][0], Q(300))
        residual = graph.multiply(graph.innovation_rows(), joint)
        self.assertEqual(residual, graph.matrix([[0], [0], [0]]))

    def test_detaching_true_source_creates_spurious_residual(self):
        detached = [[0] for _ in range(24)]
        detached[6][0] = 300
        residual = graph.multiply(graph.innovation_rows(), graph.matrix(detached))
        self.assertEqual(residual, graph.matrix([[300], [0], [0]]))

    def test_zero_entry_residual_does_not_delete_error_storage(self):
        storage = [[0]*24 for _ in range(24)]
        for i in range(12):
            storage[i][i] = 1
        pulled = graph.restrict_joint_quadratic(graph.matrix(storage), graph.fresh_entry_lift())
        self.assertEqual(pulled, graph.identity(12))
        c = graph.innovation_rows()
        residual_storage = graph.multiply(graph.transpose(c), c)
        pulled_residual = graph.restrict_joint_quadratic(residual_storage, graph.fresh_entry_lift())
        self.assertEqual(pulled_residual, graph.matrix([[0]*12 for _ in range(12)]))

    def test_graph_must_not_be_reimposed_after_other_measurements(self):
        # An intervening sensor correction changes the estimate, not the truth.
        joint = [[0] for _ in range(24)]
        joint[6][0] = Q(-1, 2)
        residual = graph.multiply(graph.innovation_rows(), graph.matrix(joint))
        self.assertEqual(residual[0][0], Q(-1, 2))

    def test_finite_measurement_sign(self):
        # x_hat=2, x_true=5 -> e=3, r=-2; K=1/4 -> e_plus=7/2.
        k = graph.matrix([[Q(1, 4)]])
        update = graph.finite_linear_measurement_map(k, (0,))
        out = graph.multiply(update, graph.matrix([[3], [5]]))
        self.assertEqual(out, graph.matrix([[Q(7, 2)], [5]]))

    def test_same_driver_cancellation_and_prediction(self):
        f = graph.matrix([[1, Q(1, 200)], [0, 1]])
        g = graph.fresh_entry_lift(2)
        self.assertEqual(graph.multiply(graph.joint_prediction_map(f), g), graph.multiply(g, f))
        self.assertEqual(graph.multiply(graph.innovation_rows(2, (0,)), g), graph.matrix([[0, 0]]))

    def test_floating_point_inputs_cannot_masquerade_as_exact_certificate(self):
        with self.assertRaises(TypeError):
            graph.matrix([[0.1]])
        with self.assertRaises(ValueError):
            graph.innovation_rows(12, (6, 6))
        with self.assertRaises(ValueError):
            graph.restrict_joint_quadratic(graph.matrix([[1, 2], [0, 1]]), graph.identity(2))


if __name__ == "__main__":
    unittest.main()
