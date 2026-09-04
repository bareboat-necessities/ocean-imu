import copy
import importlib.util
import math
from pathlib import Path
import struct
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
spec = importlib.util.spec_from_file_location(
    "ou3_source_domain_contract", ROOT / "tools" / "ou3_source_domain_contract.py"
)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

import ou3_sea3_directional_p2_ha_feasibility as sea3  # noqa: E402
import ou3_sea3_finite_horizon_concentration as sea3_concentration  # noqa: E402
import ou3_sea3_finite_window_response_admission as sea3_window  # noqa: E402
import ou3_sea3_p1_compatibility as sea3_p1  # noqa: E402
import ou3_sea3_physical_admissibility as sea3_phys  # noqa: E402
import ou3_sea3_wave_period_spectral_identity as sea3_period_identity  # noqa: E402


def f32(value):
    return struct.unpack(">f", struct.pack(">f", float(value)))[0]


class SourceDomainContractTests(unittest.TestCase):
    def test_contract_uses_shipping_clamps_and_keeps_theorem_unpromoted(self):
        d = mod.build(mod.DEFAULT_HEADER)
        self.assertEqual(d["schema"], 3)
        self.assertTrue(d["source_generated_not_trajectory_fit"])
        self.assertTrue(d["source_complete_parameter_domain"])
        self.assertFalse(d["validated_arithmetic"])
        self.assertFalse(d["outward_rounded"])
        self.assertEqual(d["implementation_scalar_semantics"]["type"], "IEEE754_BINARY32")
        self.assertEqual(d["continuous_parameters"]["tau_aw_s"], [f32(0.02), f32(12.0)])
        self.assertEqual(d["continuous_parameters"]["sigma_aw_mps2"], [f32(0.05), f32(4.0)])
        self.assertEqual(set(d["discrete_source_branches"]["mode"]), {"H", "A"})

    def test_constexpr_arithmetic_rounds_as_binary32_after_each_operation(self):
        text = """
        constexpr float A = 0.1f;
        constexpr float B = A + A;
        constexpr float C = B + A;
        constexpr float DT = 1.0f / 200.0f;
        """
        a = f32(0.1)
        b = f32(a + a)
        c = f32(b + a)
        dt = f32(f32(1.0) / f32(200.0))
        self.assertEqual(mod.parse_const(text, "A"), a)
        self.assertEqual(mod.parse_const(text, "B"), b)
        self.assertEqual(mod.parse_const(text, "C"), c)
        self.assertEqual(mod.parse_const(text, "DT"), dt)
        self.assertNotEqual(dt, 0.005)

    def test_validated_parameter_box_outwardly_contains_every_source_endpoint(self):
        d = mod.build(mod.DEFAULT_HEADER)
        box = d["validated_parameter_box"]
        self.assertTrue(box["validated_arithmetic"])
        self.assertTrue(box["outward_rounded"])
        self.assertEqual(box["theorem_promotion"], "NOT_ESTABLISHED")
        self.assertFalse(box["continuous_word_enclosed"])
        self.assertFalse(box["nonlinear_word_enclosed"])
        for name, source_bounds in d["continuous_parameters"].items():
            lo, hi = box["continuous_parameters"][name]
            self.assertLess(lo, source_bounds[0])
            self.assertGreater(hi, source_bounds[1])
            self.assertEqual(lo, math.nextafter(source_bounds[0], -math.inf))
            self.assertEqual(hi, math.nextafter(source_bounds[1], math.inf))
        for name, source_value in d["timing_constants_s"].items():
            lo, hi = box["timing_constants_s"][name]
            self.assertLess(lo, source_value)
            self.assertGreater(hi, source_value)
            self.assertEqual(lo, math.nextafter(source_value, -math.inf))
            self.assertEqual(hi, math.nextafter(source_value, math.inf))

    def test_configured_runtime_sampling_assumption_is_explicit_and_source_bound(self):
        d = mod.build(mod.DEFAULT_HEADER)
        runtime = d["configured_runtime_assumption"]
        expected = f32(f32(1.0) / f32(200.0))
        self.assertEqual(runtime["qualification"], "CONFIGURED_VALIDATION_RUNTIME_ASSUMPTION")
        self.assertEqual(runtime["sample_period_contract"], "FIXED_SOURCE_NOMINAL")
        self.assertEqual(runtime["imu_dt_s"], expected)
        self.assertFalse(runtime["api_enforces_this_bound"])
        lo, hi = runtime["imu_dt_outward_interval_s"]
        self.assertLess(lo, expected)
        self.assertGreater(hi, expected)
        self.assertEqual(d["validated_parameter_box"]["configured_runtime"], runtime)

    def test_contract_names_every_hybrid_transition_required_for_deployment(self):
        d = mod.build(mod.DEFAULT_HEADER)
        self.assertEqual(
            set(d["hybrid_obligations"]),
            {
                "startup_handoff", "held_to_active", "magnetic_lock",
                "magnetic_regauge_refinement", "tilt_reset", "tilt_relock",
                "cooldown_reentry", "periodic_aw_covariance_sync",
            },
        )
        self.assertEqual(d["periodic_aw_covariance_sync_proof"]["required_mode"], "PSD_NONEXPANSIVE")

    def test_sea3_physical_height_period_coupling_is_part_of_the_source_contract(self):
        d = sea3_phys.build()
        self.assertEqual(sea3_phys.validate(d), [])
        self.assertTrue(
            d["three_partition_contract"]["independent_H_r_and_T_p_rectangular_extrema_forbidden"]
        )
        self.assertTrue(
            d["three_partition_contract"]["independent_three_partition_H_maxima_forbidden"]
        )
        self.assertEqual(d["repository_total_Hs_upper_m"], 8.5)
        self.assertFalse(d["left_language_inclusion_closed"])

    def test_sea3_cartesian_sea_x_rao_domain_is_rejected_before_p1(self):
        d = sea3_p1.build()
        self.assertEqual(sea3_p1.validate(d), [])
        self.assertTrue(d["cartesian_product_refuted_by_analytical_witness"])
        self.assertTrue(d["coupled_SEA3_domain_required"])
        self.assertFalse(d["independent_cartesian_sea_x_RAO_domain_is_P1_sound"])
        w = d["witness"]
        self.assertTrue(w["PM_is_JONSWAP_gamma_1_boundary"])
        self.assertTrue(w["witness_is_inside_declared_JONSWAP_gamma_interval_1_to_7"])
        self.assertTrue(w["RAO_parameter_bounds_finite"])
        self.assertTrue(w["witness_is_inside_declared_RAO_parameter_ranges"])
        self.assertTrue(w["all_nondyadic_witness_constants_outward_enclosed"])
        self.assertEqual(w["x_interval"], [1.0, 9.0])
        self.assertGreater(
            w["validated_acceleration_mean_square_lower_m2_s4"],
            w["P1_cap_squared_upper_m2_s4"],
        )
        self.assertGreater(w["validated_acceleration_RMS_lower_mps2"], 4.5)
        self.assertTrue(
            d["coupled_domain_contract"]["finite_window_deterministic_response_certificate_required"]
        )
        self.assertFalse(d["finite_window_realization_certificate_closed"])
        self.assertFalse(d["L_actual_sea_subset_Lhat_SEA3_closed"])

    def test_sea3_finite_window_admission_interface_is_exact_and_nonpromoting(self):
        d = sea3_window.build_contract()
        self.assertEqual(sea3_window.validate(d), [])
        self.assertTrue(d["finite_window_admission_predicate_defined"])
        self.assertFalse(d["gaussian_spectrum_alone_can_admit_window"])
        self.assertFalse(d["RMS_or_PSD_moment_alone_can_admit_window"])
        self.assertFalse(d["finite_window_realization_producer_closed"])
        self.assertFalse(d["L_actual_sea_subset_Lhat_SEA3_closed"])
        self.assertEqual(
            d["normal_live_caps"]["non_gravitational_cog_acceleration_norm_upper_mps2"],
            4.0,
        )
        self.assertEqual(d["normal_live_caps"]["body_rate_norm_upper_deg_s"], 30.0)

    def test_sea3_finite_window_admission_rejects_psd_only_and_accepts_validated_pathwise_bounds(self):
        contract = sea3_window.build_contract()
        base = {
            "validated_arithmetic": True,
            "outward_rounded": True,
            "post_rao_response_enclosed": True,
            "all_valid_imu_samples_covered": True,
            "trajectory_replay_used": False,
            "gaussian_spectrum_only": False,
            "rms_or_psd_only": False,
            "response_parameter_box_sha256": contract["response_parameter_box_sha256"],
            "window_samples": 200,
            "post_rao_cog_acceleration_norm_upper_mps2": 3.9,
            "body_rate_norm_upper_deg_s": 29.0,
        }
        accepted = sea3_window.evaluate_window(base)
        self.assertTrue(accepted["window_admitted_to_Lhat_SEA3"])
        self.assertFalse(accepted["global_left_inclusion_promoted_by_this_decision"])

        psd_only = copy.deepcopy(base)
        psd_only["rms_or_psd_only"] = True
        rejected = sea3_window.evaluate_window(psd_only)
        self.assertFalse(rejected["window_admitted_to_Lhat_SEA3"])
        self.assertTrue(any("RMS/PSD" in x for x in rejected["validation_failures"]))

        over_cap = copy.deepcopy(base)
        over_cap["post_rao_cog_acceleration_norm_upper_mps2"] = math.nextafter(4.0, math.inf)
        rejected = sea3_window.evaluate_window(over_cap)
        self.assertFalse(rejected["window_admitted_to_Lhat_SEA3"])
        self.assertTrue(any("exceeds Normal-Live P1 cap" in x for x in rejected["validation_failures"]))

    def test_sea3_gaussian_concentration_kernel_respects_finite_horizon_budget(self):
        one_second = sea3_concentration.build(200)
        self.assertEqual(sea3_concentration.validate(one_second), [])
        self.assertFalse(one_second["temporal_independence_required"])
        self.assertFalse(one_second["cross_axis_independence_required"])
        self.assertLessEqual(
            one_second["combined_failure_probability_upper"],
            one_second["finite_horizon_failure_probability_budget"],
        )
        self.assertGreater(
            one_second["acceleration"]["required_trace_covariance_upper_m2_s4"],
            0.0,
        )

        twenty_minutes = sea3_concentration.build(240000)
        self.assertEqual(sea3_concentration.validate(twenty_minutes), [])
        self.assertLess(
            twenty_minutes["acceleration"]["required_trace_covariance_upper_m2_s4"],
            one_second["acceleration"]["required_trace_covariance_upper_m2_s4"],
        )
        self.assertFalse(twenty_minutes["deterministic_left_inclusion_closed"])
        self.assertFalse(twenty_minutes["infinite_horizon_Gaussian_hard_bound_claimed"])

    def test_sea3_gaussian_concentration_candidate_is_fail_closed(self):
        d = sea3_concentration.build(200)
        digest = d["response_parameter_box_sha256"]
        acc_threshold = d["acceleration"]["required_trace_covariance_upper_m2_s4"]
        rate_threshold = d["body_rate"]["required_trace_covariance_upper_deg2_s2"]

        good = sea3_concentration.evaluate_covariance_candidate(
            d,
            acceleration_trace_covariance_upper_m2_s4=0.5 * acc_threshold,
            body_rate_trace_covariance_upper_deg2_s2=0.5 * rate_threshold,
            validated_covariance_trace_enclosures=True,
            response_parameter_box_sha256=digest,
        )
        self.assertTrue(good["finite_horizon_good_event_candidate_pass"])
        self.assertFalse(good["deterministic_left_inclusion_promoted"])

        unvalidated = sea3_concentration.evaluate_covariance_candidate(
            d,
            acceleration_trace_covariance_upper_m2_s4=0.5 * acc_threshold,
            body_rate_trace_covariance_upper_deg2_s2=0.5 * rate_threshold,
            validated_covariance_trace_enclosures=False,
            response_parameter_box_sha256=digest,
        )
        self.assertFalse(unvalidated["finite_horizon_good_event_candidate_pass"])

        over = sea3_concentration.evaluate_covariance_candidate(
            d,
            acceleration_trace_covariance_upper_m2_s4=math.nextafter(acc_threshold, math.inf),
            body_rate_trace_covariance_upper_deg2_s2=0.5 * rate_threshold,
            validated_covariance_trace_enclosures=True,
            response_parameter_box_sha256=digest,
        )
        self.assertFalse(over["finite_horizon_good_event_candidate_pass"])

    def test_wave_period_leak_subtraction_is_exact_for_admissible_steady_spectra(self):
        d = sea3_period_identity.build()
        self.assertEqual(sea3_period_identity.validate(d), [])
        self.assertEqual(set(d["source_parity"]), set(sea3_period_identity.SOURCE_PARITY_KEYS))
        self.assertTrue(all(d["source_parity"].values()))
        ident = d["continuous_time_steady_state_identity"]
        self.assertTrue(ident["input_spectrum_nonnegative"])
        self.assertTrue(ident["weighted_denominator_finite_and_strictly_positive_required"])
        self.assertTrue(ident["weighted_second_moment_finite_required"])
        self.assertTrue(ident["holds_for_any_input_spectrum_satisfying_these_preconditions"])
        self.assertFalse(ident["narrow_band_approximation"])
        self.assertFalse(ident["single_sinusoid_approximation"])
        self.assertFalse(d["single_frequency_assumption_used"])
        self.assertFalse(d["single_RAO_used"])
        self.assertFalse(d["finite_RAO_grid_used"])
        self.assertFalse(d["promotion"]["SEA0_full_certificate_promoted"])
        self.assertFalse(d["promotion"]["P2_pruning_promoted"])
        self.assertFalse(d["promotion"]["finite_EWMA_transient_enclosed"])
        self.assertFalse(d["promotion"]["discrete_estimator_identified_with_continuous_steady_state"])

    def test_wave_period_identity_validator_rejects_overclaims(self):
        d = sea3_period_identity.build()
        bad = copy.deepcopy(d)
        bad["single_RAO_used"] = True
        bad["promotion"]["discrete_estimator_identified_with_continuous_steady_state"] = True
        bad["continuous_time_steady_state_identity"]["weighted_denominator_finite_and_strictly_positive_required"] = False
        failures = sea3_period_identity.validate(bad)
        self.assertTrue(any("single_RAO_used" in x for x in failures))
        self.assertTrue(any("discrete_estimator_identified" in x for x in failures))
        self.assertTrue(any("weighted_denominator" in x for x in failures))

    def test_sea3_rao_family_right_inclusion_is_part_of_the_canonical_source_contract(self):
        """`ou3-proof` executes this file, so the RAO bridge has no side workflow."""
        d = sea3.build_inclusion()
        self.assertEqual(sea3.validate(d), [])

        r = d["response_enclosure"]
        self.assertTrue(r["single_worst_envelope_proves_entire_parameter_box_by_monotonicity"])
        self.assertFalse(r["single_nominal_RAO_used"])
        self.assertFalse(r["finite_RAO_grid_used"])
        self.assertTrue(r["uniform_moment_theorem"]["unbanded_acceleration_moment_finite"])
        self.assertGreater(r["acceleration_moment_tightening_vs_flat_6Hz_corner_lower"], 624.0)

        p = d["p2_inclusion"]
        self.assertEqual(p["SEA3_TO_P2_INCLUSION_CERTIFICATE"], "PASS")
        self.assertTrue(p["Lhat_SEA3_subset_L_current_source"])
        self.assertFalse(p["single_RAO_selected_for_inclusion"])
        self.assertFalse(p["P2_pruned_by_SEA3"])
        self.assertEqual(p["P2_physical_source_states"], 800)
        self.assertEqual(p["RAO_parameter_box_consumed"], r["rao_envelope_parameter_box"])

    def test_sea3_validator_rejects_mismatched_consumed_rao_box(self):
        d = sea3.build_inclusion()
        bad = copy.deepcopy(d)
        # build_inclusion() intentionally reuses the response-box object in memory.
        # Detach the consumed copy here to emulate a serialized/stale candidate
        # artifact, which is the mismatch the validator is required to reject.
        bad["p2_inclusion"]["RAO_parameter_box_consumed"] = copy.deepcopy(
            bad["p2_inclusion"]["RAO_parameter_box_consumed"]
        )
        bad["p2_inclusion"]["RAO_parameter_box_consumed"]["peak_translation_gain"] = [0.0, 3.0]
        failures = sea3.validate(bad)
        self.assertIn("P2 inclusion consumed a different RAO parameter box", failures)


if __name__ == "__main__":
    unittest.main()
