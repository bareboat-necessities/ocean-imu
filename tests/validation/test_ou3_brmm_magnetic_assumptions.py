from pathlib import Path
import sys
import unittest
from fractions import Fraction as F

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "tools" / "stability"))

import ou3_brmm_magnetic_source_envelope as ENV  # noqa: E402
import ou3_brmm_magnetic_call_schedule as SCHED  # noqa: E402
import ou3_brmm_magnetic_startup_yaw_capture as CAP  # noqa: E402


class MagneticSourceEnvelopeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = ENV.build()

    def test_named_class_matches_the_declared_numbers(self):
        self.assertEqual(ENV.validate(self.d), [])
        self.assertEqual(self.d["assumption_id"], "MAG-BMM150-DET-v1")
        e = self.d["envelope"]
        self.assertEqual(e["world_field_norm_lower_uT"], 20.0)
        self.assertEqual(e["world_field_norm_upper_uT"], 75.0)
        self.assertEqual(e["world_field_horizontal_lower_uT"], 15.0)
        self.assertEqual(e["hard_iron_norm_upper_uT"], 5.0)
        self.assertEqual(e["residual_norm_upper_uT"], 2.0)

    def test_declared_PE_magnetic_bounds_become_consequences(self):
        self.assertEqual(self.d["measured_body_field_norm_lower_uT"], 13.0)
        self.assertEqual(self.d["measured_body_field_norm_upper_uT"], 82.0)
        self.assertTrue(self.d["declared_PE_magnetic_floor_is_now_derived"])
        self.assertTrue(self.d["declared_PE_magnetic_ceiling_is_now_derived"])
        self.assertTrue(self.d["shipping_mag_norm_guard_cleared_unconditionally"])

    def test_sine_separation_is_not_silently_derived(self):
        self.assertFalse(self.d["magnetic_vector_sine_separation_derivable_from_this_class"])
        self.assertTrue(self.d["sine_separation_remains_declared_PE_hypothesis"])

    def test_norm_ratio_gate_is_not_forced_by_the_declared_class(self):
        g = self.d["norm_ratio_gate"]
        self.assertAlmostEqual(g["shipping_ratio"], 0.35)
        self.assertFalse(g["universal_sample_admission_forced"])
        # ratio*Bmin/(2+ratio) = 0.35*20/2.35
        self.assertAlmostEqual(
            g["universal_admission_perturbation_budget_uT"], 0.35 * 20.0 / 2.35, places=12
        )

    def test_derived_tilt_requirements_are_exact_and_ordered(self):
        r = self.d["derived_tilt_requirements"]
        # min_horizontal_fraction: E <= (H - f*Bmax)/(1+f) = 75/7 uT, hence
        # sin(delta/2) <= (75/7 - 7)/150 = 13/525.
        self.assertEqual(
            F(r["max_sin_half_tilt_for_horizontal_fraction_gate"]).limit_denominator(10**6),
            F(13, 525),
        )
        # north non-vanishing: E < H, hence sin(delta/2) < (15-7)/150 = 4/75.
        self.assertEqual(
            F(r["max_sin_half_tilt_for_north_capture"]).limit_denominator(10**6), F(4, 75)
        )
        self.assertEqual(r["binding_requirement"], "min_horizontal_fraction")
        self.assertLess(
            r["max_tilt_deg_for_horizontal_fraction_gate"], r["max_tilt_deg_for_north_capture"]
        )

    def test_declared_startup_direction_error_reproduces_17_over_30(self):
        a = self.d["at_declared_startup_direction_error"]
        self.assertEqual(a["mean_perturbation_upper_uT"], 8.5)
        self.assertEqual(a["sin_yaw_error_upper_exact"], "17/30")
        self.assertTrue(a["north_nonvanishing"])
        self.assertTrue(a["horizontal_fraction_gate_satisfied"])

    def test_supply_is_not_silently_discharged(self):
        self.assertFalse(self.d["declared_startup_direction_error_is_the_accumulation_frame_error"])
        self.assertTrue(self.d["accumulation_frame_error_supply_is_open_obligation"])
        self.assertFalse(self.d["P4_promoted_here"])
        self.assertFalse(self.d["P5_promoted_here"])

    def test_derive_rejects_out_of_range_tilt_supply(self):
        c = ENV.shipping_constants()
        with self.assertRaises(ValueError):
            ENV.derive(F(-1, 100), c)
        with self.assertRaises(ValueError):
            ENV.derive(F(2), c)


class MagneticCallScheduleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = SCHED.build()

    def test_release_time_is_finite_and_forces_the_strict_guard(self):
        self.assertEqual(SCHED.validate(self.d), [])
        self.assertEqual(self.d["assumption_id"], "MAG-CALL-SCHEDULE-v1")
        r = self.d["release"]
        self.assertEqual(r["unlock_count"], 250)
        # 0.04 + 249*0.04 = 10.0 s from Live; 9.96 s from the first call.
        self.assertAlmostEqual(r["live_to_unlock_upper_s"], 10.0, places=12)
        self.assertAlmostEqual(r["elapsed_first_to_unlock_upper_s"], 9.96, places=12)
        self.assertTrue(r["strict_one_second_guard_forced"])
        self.assertTrue(self.d["H18_A21_RELEASE_TIME_CLOSED"])

    def test_schedule_is_only_25_hz_and_separate_from_the_value_class(self):
        s = self.d["schedule"]
        self.assertAlmostEqual(s["effective_rate_hz"], 25.0)
        self.assertTrue(s["is_deployment_timing_assumption_not_sensor_guarantee"])
        self.assertFalse(self.d["magnetic_value_class_used_here"])

    def test_counter_parity_and_open_downstream_obligations(self):
        p = self.d["shipping_parity"]
        self.assertTrue(p["counter_increment_is_unconditional"])
        self.assertTrue(p["release_requires_live_stage"])
        self.assertTrue(p["external_hold_gates_enable"])
        self.assertFalse(self.d["innovation_acceptance_used_for_counter"])
        self.assertFalse(self.d["eventual_A21_under_arbitrary_external_acc_bias_hold_claimed"])
        self.assertFalse(self.d["H18_A21_covariance_transport_closed_here"])

    def test_release_reachability_rejects_nonpositive_count(self):
        with self.assertRaises(ValueError):
            SCHED.release_reachability(0)


class StartupYawCaptureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = CAP.build()

    def test_declared_supply_lands_inside_the_declared_entrance_set(self):
        self.assertEqual(CAP.validate(self.d), [])
        b = self.d["declared_supply_branch"]
        self.assertEqual(b["sin_yaw_error_upper_exact"], "17/30")
        self.assertAlmostEqual(b["yaw_rad_upper"], 0.61, places=12)
        self.assertAlmostEqual(b["full_attitude_rad_upper"], 0.63, places=12)
        self.assertLess(b["full_attitude_deg_upper"], self.d["declared_entrance_full_attitude_deg"])
        self.assertTrue(self.d["DECLARED_SUPPLY_ENTRANCE_CLOSED"])
        self.assertFalse(b["sqrtN_statistical_reduction_used"])

    def test_taylor_bound_dominates_the_magnetic_sine_bound(self):
        lower = CAP._taylor_sin_lower(F(61, 100))
        self.assertGreater(lower, F(17, 30))
        with self.assertRaises(ValueError):
            CAP._taylor_sin_lower(F(3, 2))

    def test_certified_supply_branch_stays_fail_closed_with_a_reported_gap(self):
        c = self.d["certified_supply_branch"]
        self.assertFalse(c["invariant_closed_at_padded_envelope"])
        self.assertFalse(c["magnetic_capture_forced_from_certified_supply"])
        self.assertFalse(c["certified_tilt_alone_inside_declared_entrance_set"])
        self.assertFalse(c["perfect_yaw_gauge_would_repair_entrance"])
        self.assertGreater(c["gate_shortfall_factor"], 1.0)
        self.assertGreater(c["entrance_shortfall_factor"], 1.0)
        self.assertTrue(c["bound_is_level_set_radius_not_accuracy_bound"])
        self.assertFalse(self.d["CERTIFIED_SUPPLY_ENTRANCE_CLOSED"])
        self.assertFalse(self.d["accumulation_frame_tilt_supply_discharged"])
        self.assertFalse(self.d["timeout_branch_yaw_gauge_forced"])

    def test_scaling_argument_rules_out_metric_refinement_alone(self):
        s = self.d["forcing_scaling_argument"]
        self.assertTrue(s["is_scaling_argument_not_certificate"])
        self.assertTrue(s["metric_refinement_alone_cannot_reach_requirement"])
        self.assertGreater(s["shortfall_factor"], 1.0)
        # theta_ss = -m/s from the observer's own integral channel, plus the
        # 0.1*xi primitive term carried by z = x - B*xi.
        self.assertAlmostEqual(
            s["metric_free_tilt_floor_rad"],
            s["declared_mean_direction_chord_norm_upper"] / s["endpoint_sector_s_lower"]
            + 0.1 * s["declared_direction_primitive_norm_upper_s"],
            places=12,
        )

    def test_three_qualitatively_different_alternatives_are_recorded(self):
        self.assertGreaterEqual(len(self.d["alternatives"]), 3)
        self.assertEqual(
            self.d["limiting_quantity"],
            "private-observer accumulation tilt-frame error bound",
        )

    def test_shipping_handoff_parity_is_observed_not_assumed(self):
        p = self.d["shipping_parity"]
        self.assertTrue(p["seed_composes_proxy_tilt_with_gauge_yaw"])
        self.assertTrue(p["timeout_branch_needs_only_aligned_branch"])
        self.assertTrue(p["quality_branch_needs_north_ready"])
        self.assertTrue(p["ungauged_seed_uses_free_yaw_sigma"])
        self.assertTrue(self.d["declared_entrance_ungauged_branch_is_tilt_only"])


class NoReplayOrRetuningTests(unittest.TestCase):
    def test_new_modules_do_not_replay_or_move_frozen_gates(self):
        for name in (
            "ou3_brmm_magnetic_source_envelope",
            "ou3_brmm_magnetic_call_schedule",
            "ou3_brmm_magnetic_startup_yaw_capture",
        ):
            text = (ROOT / "tools" / "stability" / f"{name}.py").read_text()
            for forbidden in ("ou3_exact_replay", "np.quantile", "observed_min", "replay_min", "1e-19"):
                self.assertNotIn(forbidden, text, name)


if __name__ == "__main__":
    unittest.main()
