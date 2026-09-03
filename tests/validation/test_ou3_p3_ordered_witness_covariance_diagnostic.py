import math
import unittest

import ou3_p3_ordered_witness_covariance_diagnostic as D


class OrderedWitnessCovarianceDiagnosticTests(unittest.TestCase):
    def test_period_change_matches_shipping_fmod_semantics(self):
        self.assertAlmostEqual(D._set_period(0.17, 0.05), math.fmod(0.17, 0.05))

    def test_periodic_due_waits_then_fires(self):
        due, elapsed = D._due(0.005, 0.015, 0.0)
        self.assertFalse(due)
        self.assertAlmostEqual(elapsed, 0.005)
        due, elapsed = D._due(0.005, 0.015, elapsed)
        self.assertFalse(due)
        due, elapsed = D._due(0.005, 0.015, elapsed)
        self.assertTrue(due)
        self.assertGreaterEqual(elapsed, 0.0)
        self.assertLess(elapsed, 0.015)

    def test_extension_prefers_exact_self_edge_and_ends_on_target_phase(self):
        rt = {
            "nodes": [{}, {}],
            "gaps": [13, 14],
            "labelled_successors": [
                [[0], [1]],
                [[1], [0]],
            ],
        }
        witness = [{
            "source": 0,
            "successor": 0,
            "gap_samples": 13,
            "cumulative_samples": 13,
        }]
        segs = D.extend_witness_to_target(witness, rt, 30)
        self.assertEqual([row[0] for row in segs], [0, 0, 0])
        self.assertEqual([row[1] for row in segs], [13, 13, 4])
        self.assertEqual(segs[-1][2], 13)
        self.assertFalse(segs[-1][4])

    def test_synthetic_tuple_reconstructs_process_maxima(self):
        summary = {
            "sigma_squared_upper": 36.0,
            "q_c_upper": 72.0,
            "S_measurement_variance_upper": 16.0,
            "pseudo_update_cadence_s": [0.005, 0.25],
        }
        p = D._synthetic(summary)
        self.assertAlmostEqual(p["sigma"] ** 2, 36.0)
        self.assertAlmostEqual(2.0 * p["sigma"] ** 2 / p["tau"], 72.0)
        self.assertEqual(p["period"], 0.25)
        self.assertEqual(p["Rstd"], 4.0)
        self.assertIsNone(p["node"])


if __name__ == "__main__":
    unittest.main()
