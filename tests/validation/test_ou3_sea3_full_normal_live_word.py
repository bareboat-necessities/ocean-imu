from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_full_normal_live_word as mod  # noqa: E402
from ou3_interval import Interval, matrix_identity, matrix_point  # noqa: E402


class Sea3LiteralFullNormalLiveWordTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = mod.build()
        cls.failures = mod.validate(cls.d)

    def test_shipping_order_is_method_body_bound(self):
        self.assertEqual(self.failures, [])
        p = self.d["shipping_event_order_parity"]
        self.assertTrue(all(p.values()), p)
        self.assertTrue(p["commit_precedes_prediction"])
        self.assertTrue(p["pending_aw_floor_applied_inside_prediction_before_S_service"])
        self.assertTrue(p["prediction_precedes_accelerometer"])
        self.assertTrue(p["accelerometer_precedes_tuner_candidate_update"])
        self.assertTrue(p["magnetometer_is_separate_external_call"])

    def test_literal_word_keeps_every_required_event_family(self):
        self.assertGreaterEqual(self.d["imu_samples_upper"], 600)
        self.assertGreaterEqual(self.d["guaranteed_S_updates_lower_over_word"], 4)
        self.assertGreater(self.d["guaranteed_aw_floor_applications_lower_over_word"], 0)
        self.assertTrue(self.d["every_valid_imu_sample_requires_prediction"])
        self.assertTrue(self.d["every_valid_imu_sample_requires_accelerometer_Joseph"])
        self.assertTrue(self.d["S_scheduler_is_executed_not_replaced_by_selected_four"])
        self.assertTrue(self.d["magnetometer_is_asynchronous_external_event_family"])
        self.assertFalse(self.d["hardware_magnetometer_ODR_used_as_PE_recurrence"])

    def test_accelerometer_packs_full_cross_block_jacobian(self):
        for mode in ("H", "A"):
            H = mod.H_accelerometer(
                mode,
                [Interval.point(1.0), Interval.point(2.0), Interval.point(3.0)],
                matrix_identity(3),
            )
            self.assertEqual(len(H), 3)
            self.assertEqual(len(H[0]), mod.state_dimension(mode))
            # J_aw = I is retained.
            for i in range(3):
                self.assertEqual(H[i][mod.OFF_AW + i], Interval.point(1.0))
            # J_theta=-skew(f) is retained and nonzero.
            self.assertNotEqual(H[0][mod.OFF_TH + 1], Interval.point(0.0))
            if mode == "A":
                for i in range(3):
                    self.assertEqual(H[i][mod.OFF_BA + i], Interval.point(1.0))

    def test_one_literal_sample_executes_shipping_kalman_order(self):
        mode = "H"
        n = mod.state_dimension(mode)
        P0 = matrix_point([[2.0 if i == j else 0.0 for j in range(n)] for i in range(n)])
        w = mod.initialize_word(mode, P0)
        Faa = matrix_identity(6)
        Qaa = matrix_point([[0.01 if i == j else 0.0 for j in range(6)] for i in range(6)])
        Fll = matrix_identity(12)
        Qll = matrix_point([[0.01 if i == j else 0.0 for j in range(12)] for i in range(12)])
        F, Q = mod.pack_prediction(mode, Faa, Qaa, Fll, Qll)
        mod.apply_imu_sample(
            w,
            F=F,
            Q=Q,
            f_cog_body=[Interval.point(0.0), Interval.point(0.0), Interval.point(-9.80665)],
            R_wb=matrix_identity(3),
            Racc=mod.diagonal_R([0.2, 0.2, 0.2]),
            due_S=True,
            rs_std_xyz=[Interval.point(0.72), Interval.point(0.72), Interval.point(1.0)],
            Delta_aw=matrix_point([[0.001 if i == j else 0.0 for j in range(3)] for i in range(3)]),
        )
        self.assertEqual(w.event_log, ["prediction", "aw_floor", "S_zero", "accelerometer"])
        self.assertTrue(mod.BACKEND.decomposition_identity_enclosed(w.riccati))

    def test_reduced_routes_cannot_promote(self):
        r = self.d["no_reduced_promotion_routes"]
        self.assertTrue(r)
        self.assertTrue(all(v is False for v in r.values()), r)
        self.assertFalse(self.d["SOURCE_REACHABLE_EVENT_FAMILY_MATERIALIZED"])
        self.assertFalse(self.d["FULL_H18_WORD_EXECUTED"])
        self.assertFalse(self.d["FULL_A21_WORD_EXECUTED"])
        self.assertFalse(self.d["FULL_H18_A21_LDLT_CLOSED"])
        self.assertFalse(self.d["P3_CANONICAL_PASS"])


if __name__ == "__main__":
    unittest.main()
