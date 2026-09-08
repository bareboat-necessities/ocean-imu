from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests/kalman_ou_iii"))
import ou3_p4_connected_motion as C  # noqa: E402
import ou3_p4_trace_overlay as T  # noqa: E402


class ConnectedMotionTests(unittest.TestCase):
    def test_failure_reporter_binds_each_event_and_mode(self):
        first, second = C.ParityChecks(), C.ParityChecks()
        first.row = {"index": 5, "stage": "prediction"}
        first("mismatch", [1.], [0.])
        first.row = {"index": 6, "stage": "reset"}
        first("mismatch", [2.], [0.])
        second.row = {"index": 100, "stage": "projection"}
        second("mismatch", [3.], [0.])
        self.assertEqual([(r["index"], r["stage"]) for r in first.failures],
                         [(5, "prediction"), (6, "reset")])
        self.assertEqual(len(second.failures), 1)
        self.assertEqual(second.failures[0]["index"], 100)

    def row(self, kind):
        return {"kind": kind, "R_hat": np.eye(3).ravel().tolist(),
                "R_true": np.eye(3).ravel().tolist(), "x_hat": [0.]*21,
                "linear_true": [0.]*12, "true_bias": [0.]*3,
                "measured": [0.]*3, "temperature_mean": [0.]*3,
                "mag_reference": [20., 0., 40.], "gravity": 9.80665,
                "P": np.eye(21).ravel().tolist()}

    def test_S_physical_forcing_is_not_the_RS_corrective_gain(self):
        row = self.row("S_zero")
        row["linear_true"][6] = 1.
        error = np.zeros(21)
        error[12] = .25
        residual, source, h = C.residual_graph(row, error)
        np.testing.assert_array_equal(source, [-1, 0, 0])
        np.testing.assert_array_equal(residual, [-.75, 0, 0])
        # Both the source offset and the actual R_S-dependent gain matter.
        p = np.eye(21)
        p[6, 12] = p[12, 6] = .5
        results = []
        for rs in (1., 4.):
            gain = p@h.T@np.linalg.inv(h@p@h.T+rs*np.eye(3))
            after = C.correct(error, gain, residual, [0.]*3, .4)
            results.append(after[6])
            self.assertGreater(after[6], 0.)
            omitted_source = C.correct(error, gain, h@error, [0.]*3, .4)
            self.assertLess(omitted_source[6], 0.)
            self.assertNotEqual(after[6], error[6])
        self.assertGreater(results[0], results[1])

    def test_accelerometer_source_and_bias_share_full_nonlinear_residual(self):
        row = self.row("accelerometer")
        error = np.zeros(21)
        error[:3] = [.2, -.1, .05]
        error[15:18] = [.1, -.2, .3]
        error[18:21] = [.02, -.03, .01]
        row["R_true"] = C.rotation(error[:3]).ravel().tolist()
        row["linear_true"][9:12] = error[15:18].tolist()
        row["true_bias"] = error[18:21].tolist()
        source = np.array([.001, .002, -.003])
        row["measured"] = (C.rotation(error[:3])@(error[15:18]-[0, 0, row["gravity"]])
                           + error[18:21]+source).tolist()
        residual, defect, h = C.residual_graph(row, error)
        np.testing.assert_allclose(defect, source, atol=2e-15)
        np.testing.assert_allclose(residual, np.asarray(row["measured"])+[0, 0, row["gravity"]], atol=2e-15)
        np.testing.assert_array_equal(h[:, 18:21], np.eye(3))

    def test_magnetic_reference_defect_survives_zero_noise(self):
        row = self.row("magnetometer")
        row["measured"] = [19., 2., 40.]
        residual, source, _ = C.residual_graph(row, np.zeros(21))
        np.testing.assert_array_equal(source, [-1., 2., 0.])
        np.testing.assert_array_equal(residual, source)

    def test_finite_rotation_and_projection_keep_all_21_rows(self):
        error = np.zeros(21)
        error[:3] = [.3, -.1, .2]
        gain = np.zeros((21, 3))
        gain[:3] = .2*np.eye(3)
        gain[6:9] = .1*np.eye(3)
        gain[18:21] = np.eye(3)
        residual = np.array([.7, -.2, .1])
        after = C.correct(error, gain, residual, [0.]*3, .4)
        np.testing.assert_allclose(C.rotation(after[:3]),
                                   C.rotation(error[:3])@C.deployed_rotation(.2*residual).T, atol=3e-16)
        np.testing.assert_allclose(after[6:9], -.1*residual, atol=1e-16)
        self.assertAlmostEqual(np.linalg.norm(after[18:21]), .4)

    def test_performance_storage_does_not_discard_active_covariance_cross_terms(self):
        row = self.row("S_zero")
        row["x_hat"][6] = 1.
        row["x_hat"][18] = .2
        p = np.eye(21)
        p[6, 18] = p[18, 6] = .5
        row["P"] = p.ravel().tolist()
        self.assertAlmostEqual(C.motion_energy(row), 4/3)
        self.assertAlmostEqual(float(C.exact_motion_energy(row)), 4/3)
        self.assertEqual(row["x_hat"][18], .2)

    def test_small_angle_branch_is_shipping_polynomial_not_linear_surrogate(self):
        d = np.array([.008, 0., 0.])
        t2, t4 = np.dot(d, d), np.dot(d, d)**2
        expected = C.rotation(2*(.5-t2/48+t4/3840)*d/(1-t2/8+t4/384))
        np.testing.assert_array_equal(C.deployed_rotation(d), expected)
        self.assertGreater(np.max(np.abs(C.deployed_rotation(d)-C.rotation(d))), 1e-9)

    def test_trace_overlay_only_inserts_uniquely_identified_read_callbacks(self):
        original = (ROOT / "src/kalman_ou_iii/Kalman3D_Wave_OU_III.h").read_text()
        traced, edits = T.instrument(original)
        self.assertEqual(len(edits), 10)
        for edit in edits:
            self.assertEqual(traced.count(edit), 1)
            traced = traced.replace(edit, "")
        self.assertEqual(traced, original)
        with self.assertRaises(ValueError):
            T.instrument(original.replace("    project_acc_bias_();", ""))

    def test_missing_complete_source_word_cannot_reach_motion_gate(self):
        with self.assertRaises(ValueError):
            C.audit({"dt": .005}, [], [])


if __name__ == "__main__":
    unittest.main()
