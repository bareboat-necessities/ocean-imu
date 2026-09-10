#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_p5_live_entry_reachability as REACH


class LiveEntryReachabilityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = REACH.build()

    def test_validates_and_promotes_nothing(self):
        self.assertEqual(REACH.validate(self.d), [])
        self.assertFalse(self.d["P4_promoted_here"])
        self.assertFalse(self.d["P5_MAY_START"])
        self.assertFalse(self.d["P5_capture_assumed"])
        self.assertFalse(self.d["reachability_argument_is_replay_or_fit"])
        self.assertFalse(self.d["inventing_primitive_values_here"])

    def test_every_shipping_parity_string_is_present(self):
        self.assertTrue(self.d["shipping_source_parity_pass"])
        for name, ok in self.d["shipping_source_parity"].items():
            self.assertTrue(ok, name)

    def test_linear_block_is_not_propagated_before_go_live(self):
        self.assertFalse(self.d["mekf_linear_block_propagated_before_go_live"])
        self.assertTrue(self.d["shipping_source_parity"]["mekf_drive_calls_are_inside_the_guard"])

    def test_guard_detector_rejects_an_unguarded_drive_call(self):
        good = "if (drive_mekf) {\n mekf_->time_update(gyro, dt);\n mekf_->measurement_update_acc_only(acc_in, tempC);\n}\n"
        self.assertTrue(REACH._guarded_mekf_calls(good))
        moved_out = "if (drive_mekf) {\n}\nmekf_->time_update(gyro, dt);\nmekf_->measurement_update_acc_only(acc_in, tempC);\n"
        self.assertFalse(REACH._guarded_mekf_calls(moved_out))
        missing = "if (drive_mekf) {\n mekf_->time_update(gyro, dt);\n}\n"
        self.assertFalse(REACH._guarded_mekf_calls(missing))

    def test_entry_relation_is_fully_correlated(self):
        rel = self.d["deployed_live_entry_relation"]
        self.assertFalse(rel["coordinates_are_independent"])
        self.assertTrue(rel["integral_is_the_running_integral_of_the_same_position_history"])
        self.assertEqual(rel["cross_covariance_attitude_to_linear_at_entry"], 0.0)
        self.assertEqual(rel["first_S_zero_event_attitude_gain"], 0.0)

    def test_every_non_attitude_coordinate_is_pinned_to_negated_truth(self):
        pinned = self.d["pinned_entry_coordinates"]
        self.assertEqual(set(pinned), set(REACH.PINNED_COORDINATES))
        for name, row in pinned.items():
            with self.subTest(coordinate=name):
                self.assertTrue(row["equals"].startswith("-"))
                self.assertTrue(row["equals"].endswith("_true(T)"))

    def test_open_radii_are_exactly_the_uninstantiated_BRMM_primitives(self):
        self.assertEqual(sorted(self.d["uninstantiated_BRMM_primitives"]), ["P_m", "S_m", "V_m"])
        self.assertFalse(self.d["BRMM_motion_primitives_instantiated"])
        self.assertEqual(set(self.d["entry_radii_still_open"]),
                         {"velocity", "position", "integral_displacement"})

    def test_bias_entry_radius_is_the_true_bias_envelope_not_the_projection_radius(self):
        row = self.d["pinned_entry_coordinates"]["accelerometer_bias"]
        self.assertTrue(row["instantiated"])
        declared = self.d["declared_independent_handoff_radii"]["accelerometer_bias"]
        self.assertLess(float(row["value"]), declared)
        self.assertTrue(self.d["deployed_pin_is_tighter_than_declared"]["accelerometer_bias"])

    def test_whole_file_state_writer_audit_is_closed(self):
        allowed = set(self.d["allowed_unguarded_state_writers"])
        self.assertEqual(allowed, set(REACH.ALLOWED_UNGUARDED_STATE_WRITERS))
        found = set(self.d["mekf_state_writers_outside_the_drive_guard"])
        self.assertTrue(found <= allowed, sorted(found - allowed))
        # Each allow-listed writer carries its own recorded argument.
        self.assertTrue(self.d["pre_live_magnetometer_cannot_move_the_linear_block"])
        self.assertTrue(self.d["pre_live_accel_relock_writes_only_the_attitude_head"])

    def test_an_unargued_unguarded_writer_fails_validation(self):
        d = dict(self.d)
        d["mekf_state_writers_outside_the_drive_guard"] = sorted(
            set(self.d["mekf_state_writers_outside_the_drive_guard"]) | {"time_update"})
        failures = REACH.validate(d)
        self.assertTrue(any("unargued MEKF state writer" in x for x in failures), failures)

    def test_writer_scan_finds_a_call_moved_out_of_the_guard(self):
        guarded = ("if (drive_mekf) {\n mekf_->time_update(gyro, dt);\n}\n")
        self.assertEqual(REACH._unguarded_state_writers(guarded), set())
        moved = ("if (drive_mekf) {\n}\nmekf_->time_update(gyro, dt);\n")
        self.assertEqual(REACH._unguarded_state_writers(moved), {"time_update"})

    def test_capture_timing_and_attitude_radius_remain_open(self):
        self.assertFalse(self.d["go_live_timing_established_here"])
        self.assertFalse(self.d["proxy_attitude_radius_established_here"])


if __name__ == "__main__":
    unittest.main()
