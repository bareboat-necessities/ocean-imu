from pathlib import Path
import json
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests/kalman_ou_iii"))
import ou3_p4_source_endpoint as E  # noqa: E402


class SourceEndpointTests(unittest.TestCase):
    def test_forced_endpoint_cannot_be_promoted_to_homogeneous(self):
        self.assertEqual(E.decision(source_matches=True, common_premises=True,
                                    forcing_zero=False),
                         "ACTUAL_FORCED_ENDPOINT_NOT_HOMOGENEOUS_P4_TEST")

    def test_failed_source_and_common_premises_remain_blockers(self):
        self.assertEqual(E.decision(source_matches=False, common_premises=True,
                                    forcing_zero=True), "SOURCE_ATTACHMENT_REJECTED")
        self.assertEqual(E.decision(source_matches=True, common_premises=False,
                                    forcing_zero=True),
                         "FIXED_WORD_OUTSIDE_CHECKED_NORMAL_LIVE_PREMISES")

    def test_zero_forcing_does_not_supply_missing_subevent_graph(self):
        self.assertEqual(E.decision(source_matches=True, common_premises=True,
                                    forcing_zero=True),
                         "SUBEVENT_GRAPH_AND_HOMOGENEOUS_ATTACHMENT_STILL_REQUIRED")

    def test_source_antiderivative_and_orbital_velocity_share_root(self):
        root = {"atoms": [[1., 1/9.8, .05, .4, .6, .8]]}
        x, rotation = E.stokes(root, .3)
        eps = 1e-5
        xp, _ = E.stokes(root, .3+eps)
        xm, _ = E.stokes(root, .3-eps)
        derivative = (xp-xm)/(2*eps)
        np.testing.assert_allclose(derivative[6:9], x[3:6], atol=1e-10)
        np.testing.assert_allclose(derivative[3:6], x[:3], atol=1e-10)
        np.testing.assert_allclose(derivative[:3], x[9:12], atol=1e-10)
        np.testing.assert_allclose(rotation@rotation.T, np.eye(3), atol=1e-14)

    def test_full_21_state_storage_keeps_bias_and_cross_terms(self):
        p = np.eye(21)
        p[6, 18] = p[18, 6] = .25
        x = np.zeros(21)
        x[6], x[18] = 2, .1
        row = {"R_true": np.eye(3).ravel().tolist(),
               "R_hat": np.eye(3).ravel().tolist(),
               "x_hat": x.tolist(), "linear_true": np.zeros(12).tolist(),
               "true_bias": [0, 0, 0], "P": p.ravel().tolist()}
        e, covariance = E.error_and_covariance(row, 21)
        expected = float(e @ np.linalg.solve(covariance, e))
        self.assertAlmostEqual(float(E.high_precision_energy(row, 21)), expected, places=13)
        self.assertNotAlmostEqual(expected, float(e@e), places=5)

    def test_complete_report_is_native_json_without_silent_promotion(self):
        root = {"atoms": [[1., 1/9.8, 0., .4, .6, .8]]*128,
                "gravity": 9.8, "branch": "STOKES_WAVE_FOLLOWING"}
        estimate = np.zeros(21)
        estimate[6] = 1
        rows = []
        for mode in ("H18", "A21"):
            for index in range(601):
                rows.append({"word": mode, "event": "root" if index == 0 else "imu",
                             "index": index, "source_time": index*.005, "wind_heel": 0.,
                             "R_true": np.eye(3).ravel().tolist(),
                             "R_hat": np.eye(3).ravel().tolist(),
                             "x_hat": estimate.tolist(), "P": np.eye(21).ravel().tolist(),
                             "linear_true": [0.]*12, "prediction_forcing": [0.]*12,
                             "physical_acceleration": [0., 0., 0.], "physical_gyro": [0., 0., 0.],
                             "true_bias": [0., 0., 0.], "tau_b": 5000., "live": True,
                             "active": mode == "A21", "mag_lock": True, "mag_refined": True,
                             "acc_accepted": True, "s_due": False})
        report = json.loads(json.dumps(E.audit(root, rows), allow_nan=False))
        self.assertFalse(report["P4_promoted"])
        for mode in ("H18", "A21"):
            self.assertIsNone(report["modes"][mode]["homogeneous_P4_endpoint_ratio"])
            self.assertEqual(report["modes"][mode]["actual_forced_endpoint_ratio"], "1")


if __name__ == "__main__":
    unittest.main()
