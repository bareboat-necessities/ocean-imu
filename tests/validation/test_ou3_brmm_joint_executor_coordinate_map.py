import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools" / "stability"))

from ou3_interval import Interval
import ou3_brmm_joint_executor_coordinate_map as M


class JointExecutorCoordinateMapTests(unittest.TestCase):
    def test_status_closes_map_not_provider_or_p4(self):
        d = M.build()
        self.assertEqual(M.validate(d), [])
        self.assertTrue(d["joint_source_output_map_closed"])
        self.assertTrue(d["raw_gyro_and_corrected_rate_distinct"])
        self.assertTrue(d["measurement_coordinates_and_nominal_geometry_distinct"])
        self.assertFalse(d["complete_601_sample_provider_materialized_here"])
        self.assertFalse(d["sensor_forcing_hard_bound_closed_here"])
        self.assertFalse(d["P4_PASS"])
        self.assertEqual(d["P3_delta"], 1e-18)

    def test_corrected_rate_is_raw_minus_estimated_bias(self):
        I = Interval.point
        got = M.corrected_gyro_rate_body(
            (I(1.0), I(2.0), I(3.0)),
            (I(0.1), I(-0.2), I(0.3)),
        )
        self.assertEqual([(x.lo, x.hi) for x in got], [(0.9,0.9),(2.2,2.2),(2.7,2.7)])

    def test_world_down_gravity_sign_and_world_to_body_map(self):
        I = Interval.point
        R = [[I(1.0),I(0.0),I(0.0)], [I(0.0),I(1.0),I(0.0)], [I(0.0),I(0.0),I(1.0)]]
        got = M.nominal_cog_specific_force_body(R, (I(0.0),I(0.0),I(0.0)), I(9.80665))
        self.assertAlmostEqual(got[0].lo, 0.0)
        self.assertAlmostEqual(got[1].lo, 0.0)
        self.assertAlmostEqual(got[2].lo, -9.80665)

    def test_truth_and_nominal_geometry_are_not_collapsed(self):
        I = Interval.point
        Z = (I(0.0),I(0.0),I(0.0))
        R = [[I(1.0),I(0.0),I(0.0)], [I(0.0),I(1.0),I(0.0)], [I(0.0),I(0.0),I(1.0)]]
        s = M.materialize_sample(
            R_true_wb=R, R_hat_wb=R,
            a_true_world=(I(1.0),I(0.0),I(0.0)),
            a_w_hat_world=(I(2.0),I(0.0),I(0.0)),
            omega_true_bprime=Z,
            true_gyro_bias_bprime=Z,
            estimated_gyro_bias_bprime=Z,
            true_accel_bias_bprime=Z,
            gyro_forcing_bprime=Z,
            accel_forcing_bprime=Z,
            gravity_mps2=I(9.80665),
        )
        self.assertEqual(s.specific_force_body[0], I(1.0))
        self.assertEqual(s.f_cog_body[0], I(2.0))

    def test_sensor_forcing_stays_in_measurement_not_nominal_geometry(self):
        I = Interval.point
        Z = (I(0.0),I(0.0),I(0.0))
        R = [[I(1.0),I(0.0),I(0.0)], [I(0.0),I(1.0),I(0.0)], [I(0.0),I(0.0),I(1.0)]]
        s = M.materialize_sample(
            R_true_wb=R, R_hat_wb=R,
            a_true_world=Z, a_w_hat_world=Z, omega_true_bprime=Z,
            true_gyro_bias_bprime=Z, estimated_gyro_bias_bprime=Z,
            true_accel_bias_bprime=Z, gyro_forcing_bprime=Z,
            accel_forcing_bprime=(I(1.0),I(0.0),I(0.0)),
            gravity_mps2=I(9.80665),
        )
        self.assertEqual(s.specific_force_body[0], I(1.0))
        self.assertEqual(s.f_cog_body[0], I(0.0))

    def test_false_shortcuts_and_promotion_rejected(self):
        d = M.build()
        d["independent_raw_and_corrected_gyro_boxes_allowed"] = True
        d["P4_PASS"] = True
        f = M.validate(d)
        self.assertIn("independent_raw_and_corrected_gyro_boxes_allowed not false", f)
        self.assertIn("P4_PASS not false", f)


if __name__ == "__main__":
    unittest.main()
