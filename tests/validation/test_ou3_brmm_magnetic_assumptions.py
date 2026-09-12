from pathlib import Path
import math
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

    def test_excursion_parameter_covers_vessel_heading(self):
        r = self.d["derived_tilt_requirements"]
        # The accumulation frame strips the estimator's yaw, not the vessel's,
        # so a heading change during the window smears the mean.  The parameter
        # is the total excursion and the heading supply is recorded as missing.
        self.assertEqual(
            r["parameter"], "TOTAL_ACCUMULATION_FRAME_ROTATION_EXCURSION_TILT_PLUS_HEADING"
        )
        self.assertTrue(self.d["accumulation_frame_retains_vessel_heading"])
        self.assertFalse(self.d["heading_excursion_bounded_by_this_class"])
        self.assertFalse(self.d["startup_heading_excursion_declared_in_operating_domain"])
        self.assertTrue(self.d["declared_startup_direction_error_carries_no_heading_content"])
        a = self.d["at_declared_startup_direction_error"]
        self.assertFalse(a["supply_covers_heading_excursion"])
        self.assertTrue(a["row_is_conditional_on_a_total_excursion_supply"])
        self.assertTrue(ENV.derive(F(0), ENV.shipping_constants())["excursion_includes_vessel_heading"])

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

    def test_release_time_is_derived_by_a_case_split(self):
        self.assertEqual(SCHED.validate(self.d), [])
        self.assertEqual(self.d["assumption_id"], "MAG-CALL-SCHEDULE-v1")
        r = self.d["release"]
        self.assertEqual(r["unlock_count"], 250)
        # The gap bound is an UPPER bound, so `elapsed > guard` cannot be
        # inferred from it: 250 calls 1 ms apart clear the count at 0.249 s with
        # the shipping `> 1.0f` predicate still false.  Both conditions are
        # monotone once true, so the release fires by the later of the two.
        self.assertAlmostEqual(r["count_condition_satisfied_by_upper_s"], 10.0, places=12)
        self.assertAlmostEqual(r["guard_condition_satisfied_by_upper_s"], 1.08, places=12)
        self.assertAlmostEqual(
            r["live_to_unlock_upper_s"],
            max(r["count_condition_satisfied_by_upper_s"], r["guard_condition_satisfied_by_upper_s"]),
            places=12,
        )
        self.assertEqual(r["binding_condition"], "count")
        self.assertTrue(r["release_derived_by_case_split_not_upper_bound_comparison"])
        self.assertFalse(r["lower_bound_on_call_spacing_assumed"])

    def test_release_is_conditional_on_north_lock(self):
        # The only call site of the counter-owning inner method sits behind
        # `if (mag_ref_set_ && stage_ == Stage::Live)`, so a host call schedule
        # alone does not advance the counter.
        self.assertTrue(self.d["shipping_parity"]["counter_call_is_behind_north_lock_gate"])
        self.assertTrue(self.d["H18_A21_RELEASE_TIME_CLOSED_AFTER_NORTH_LOCK"])
        self.assertFalse(self.d["H18_A21_RELEASE_TIME_CLOSED_FROM_LIVE"])
        self.assertFalse(self.d["counter_advances_without_north_lock"])
        self.assertFalse(self.d["ungauged_timeout_path_reaches_the_release"])
        self.assertTrue(self.d["north_lock_is_an_unmet_prerequisite_on_this_route"])
        self.assertTrue(self.d["release"]["measured_from_north_lock_not_from_Live"])

    def test_release_reachability_rejects_nonpositive_gap(self):
        with self.assertRaises(ValueError):
            SCHED.release_reachability(250, gap=F(0))

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

    def test_magnetic_algebra_closes_but_the_entrance_claim_is_retracted(self):
        self.assertEqual(CAP.validate(self.d), [])
        b = self.d["declared_supply_branch"]
        self.assertEqual(b["sin_yaw_error_upper_exact"], "17/30")
        self.assertAlmostEqual(b["yaw_rad_upper"], 0.61, places=12)
        self.assertAlmostEqual(b["full_attitude_rad_upper"], 0.63, places=12)
        self.assertLess(b["full_attitude_deg_upper"], self.d["declared_entrance_full_attitude_deg"])
        self.assertFalse(b["sqrtN_statistical_reduction_used"])
        # The algebra holds; the supply does not.  The declared quantity is a
        # gravity-direction error and the perturbation bound needs the total
        # accumulation-frame excursion, heading included.
        self.assertTrue(self.d["magnetic_algebra_closed"])
        self.assertTrue(b["inside_declared_entrance_set_if_supply_covered_heading"])
        self.assertTrue(b["supply_is_gravity_direction_error_only"])
        self.assertFalse(b["supply_covers_total_excursion"])
        self.assertFalse(self.d["DECLARED_SUPPLY_ENTRANCE_CLOSED"])
        self.assertFalse(self.d["total_excursion_supply_exists"])
        self.assertFalse(self.d["heading_excursion_charged_as_tilt"])

    def test_taylor_bound_dominates_the_magnetic_sine_bound(self):
        lower = CAP._taylor_sin_lower(F(61, 100))
        self.assertGreater(lower, F(17, 30))
        with self.assertRaises(ValueError):
            CAP._taylor_sin_lower(F(3, 2))

    def test_certified_supply_branch_stays_fail_closed_with_a_reported_gap(self):
        c = self.d["certified_supply_branch"]
        # The two-phase Mahony certificate now closes, so the supply exists -
        # it is simply far too weak, which is a different fact and is reported
        # as such rather than as a missing certificate.
        self.assertTrue(c["invariant_closed_at_padded_envelope"])
        self.assertEqual(c["invariant_validation_failures"], [])
        self.assertTrue(c["two_phase_certificate_consumed"])
        self.assertFalse(c["magnetic_capture_forced_from_certified_supply"])
        self.assertFalse(c["magnetic_capture_forced_from_ultimate_supply"])
        self.assertFalse(c["certified_tilt_alone_inside_declared_entrance_set"])
        self.assertFalse(c["perfect_yaw_gauge_would_repair_entrance"])
        self.assertGreater(c["gate_shortfall_factor"], 1.0)
        self.assertGreater(c["ultimate_gate_shortfall_factor"], 1.0)
        self.assertGreater(c["ultimate_north_capture_shortfall_factor"], 1.0)
        self.assertGreater(c["entrance_shortfall_factor"], 1.0)
        # The accumulation frame can start well before the deployed horizon, so
        # the sound supply is the all-time bound, not the asymptotic one.
        self.assertTrue(c["accumulation_frame_may_start_before_the_horizon"])
        self.assertTrue(c["sound_supply_for_accumulation_frame_is_all_time_bound"])
        self.assertLess(c["certified_ultimate_tilt_deg_upper"], c["certified_tilt_deg_upper"])
        self.assertFalse(self.d["CERTIFIED_SUPPLY_ENTRANCE_CLOSED"])
        self.assertFalse(self.d["accumulation_frame_tilt_supply_discharged"])
        self.assertFalse(self.d["timeout_branch_yaw_gauge_forced"])

    def test_formulation_floor_rules_out_the_whole_route_not_just_a_metric(self):
        s = self.d["declared_formulation_tilt_floor"]
        self.assertTrue(s["is_formulation_lower_bound_not_physical_claim"])
        self.assertTrue(s["frozen_sector_member_s_equals_one"])
        self.assertTrue(s["no_metric_or_level_can_reach_requirement"])
        self.assertGreater(s["shortfall_factor"], 1.0)
        # The floor is the in-band primitive-channel peak plus an admitted DC
        # mean chord; both are superposed on the linear z-dynamics.
        self.assertTrue(s["peak_frequency_inside_declared_wave_band"])
        lo, hi = s["declared_wave_band_hz"]
        self.assertLessEqual(lo, s["peak_frequency_hz"])
        self.assertLessEqual(s["peak_frequency_hz"], hi)
        self.assertAlmostEqual(
            s["declared_formulation_tilt_floor_rad"],
            s["primitive_channel_tilt_floor_rad"] + s["mean_chord_tilt_floor_rad"],
            places=12,
        )
        # It must also dominate the weaker north-capture requirement.
        self.assertGreater(
            math.degrees(s["declared_formulation_tilt_floor_rad"]),
            self.d["certified_supply_branch"]["required_tilt_deg_for_north_capture"],
        )
        self.assertFalse(self.d["private_observer_can_supply_the_declared_gates"])
        self.assertTrue(self.d["quadratic_metric_class_ruled_out"])

    def test_required_envelope_narrowing_is_quantified_at_the_floor(self):
        n = self.d["required_envelope_narrowing_at_the_floor"]
        self.assertAlmostEqual(
            n["supplied_tilt_rad"],
            self.d["declared_formulation_tilt_floor"]["declared_formulation_tilt_floor_rad"],
            places=12,
        )
        self.assertEqual(n["declared_combined_perturbation_uT"], 7.0)
        self.assertAlmostEqual(n["shipping_min_horizontal_fraction"], 0.05)
        rows = n["rows"]
        self.assertGreaterEqual(len(rows), 3)
        # A narrower total-field ceiling needs a lower horizontal floor, and
        # every declared row must stay inside its own total field.
        for row in rows:
            self.assertTrue(row["feasible_within_total_field"], row)
            self.assertGreater(
                row["required_horizontal_floor_uT"], n["declared_envelope_row"]["declared_horizontal_floor_uT"]
            )
        ceilings = [r["world_field_norm_upper_uT"] for r in rows]
        floors = [r["required_horizontal_floor_uT"] for r in rows]
        self.assertEqual(ceilings, sorted(ceilings))
        self.assertEqual(floors, sorted(floors))

    def test_rational_sqrt_lower_bound_is_sound(self):
        from fractions import Fraction as Frac

        for x in (Frac(2), Frac(1, 3), Frac(157, 1000), Frac(0)):
            lo = CAP._sqrt_lower(x)
            self.assertLessEqual(lo * lo, x, x)               # sound
            step = Frac(1, 10**28)
            self.assertGreater((lo + step) * (lo + step), x, x)  # and tight
        with self.assertRaises(ValueError):
            CAP._sqrt_lower(Frac(-1))

    def test_qualitatively_different_alternatives_are_recorded(self):
        self.assertGreaterEqual(len(self.d["alternatives"]), 4)
        joined = " ".join(self.d["alternatives"])
        self.assertIn("IQC", joined)
        self.assertIn("heading-excursion", joined)
        self.assertEqual(
            self.d["limiting_quantity"],
            "declared gravity-direction forcing pair (mean chord, primitive) and the commissioned magnetic band",
        )

    def test_shipping_handoff_parity_compares_complete_expressions(self):
        p = self.d["shipping_parity"]
        self.assertTrue(p["parity_compares_complete_expressions_not_tokens"])
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
