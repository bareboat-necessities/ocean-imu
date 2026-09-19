from __future__ import annotations
import math,sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2]; sys.path.insert(0,str(ROOT))
from tools.stability.ou3_theorem.tail_stability import (
    A21TailPremises, RefinementPremises, ShippingScheduleBounds,
    active_bias_homogeneous_ratio, active_tilt_bias_minor,
    additive_supply_practical_radius, a21_entry_bound,
    attitude_gyro_quotient_uniformly_observable,
    attitude_information_from_marine_diversity, bias_projection_inactive_margin,
    compact_uniform_information_exists, coupled_information_floor,
    detectable_tail_margin, explicit_small_gain_radius, finite_bridge_bound,
    gravity_tilt_sensitivity_floor, information_contraction_ratio,
    linear_kalman_tail_exponentially_stable, local_nonlinear_radius_exists,
    marine_magnetic_diversity_window_min, marine_magnetic_vector_diversity_floor,
    neutral_quotient_uniform_mu_exists, nonlinear_margin_from_information,
    nonlinear_small_gain_exists_from_linear, nonlinear_tail_ratio,
    normalized_translation_information_floor, conservative_float32_kernel_error,
    bias_release_error_bound, projection_sector_retained_error_bound,
    vector_attitude_information_floor, rotation_integral_singular_floor,
    two_epoch_attitude_gyro_floor, explicit_neutral_information_floor,
    rho_from_explicit_mu, psd_service_schur_floor,
    aggregate_repeated_row_information_floor, aggregate_translation_information_floor,
    certified_rho_and_margin, vibration_inflated_accel_std_ceiling,
    root_attitude_information_floor, attitude_gyro_information_floor_from_windows,
    full_neutral_information_floor, covariance_normalized_information_floor,
    rho_from_covariance_normalized_mu, classical_riccati_covariance_bounds,
    practical_radius_bound, projection_sector_is_dissipative, proof_route_status,
    pseudo_decay_exponent_per_gap, release_can_guarantee_projection_inactive,
    small_gain_budget,
    refinement_completion_bound, refinement_gate_margins,
    refinement_sample_gate_uniform, refinement_tilt_limit_rad,
    release_bound_from_service, shipping_covariance_hard_events_retain_compactness,
    smooth_a21_remainder_vanishes_locally, tilt_gyro_bias_quotient_minor,
    time_varying_translation_minor_floor, translational_observability_determinant,
    translational_observability_nonsingular, uniform_controllability_exists,
)

class TailStabilityTests(unittest.TestCase):
    def test_h18_is_finite_bridge(self):
        self.assertAlmostEqual(finite_bridge_bound(2,3,1,.5),3.5)
        self.assertAlmostEqual(finite_bridge_bound(2,2,2,1),11)
    def test_service_makes_count_guard_finite_after_reference_release(self):
        self.assertEqual(release_bound_from_service(250,1,1,12),263)
    def test_nonlinear_small_gain(self):
        p=A21TailPremises(.81,.05,2,.5,3)
        self.assertAlmostEqual(nonlinear_tail_ratio(p),.9025)
        self.assertTrue(math.isfinite(practical_radius_bound(p,.1)))
    def test_noncontractive_tail_refused(self):
        p=A21TailPremises(.99,.01,1,1,1)
        self.assertGreaterEqual(nonlinear_tail_ratio(p),1)
        with self.assertRaises(ValueError): practical_radius_bound(p,.1)
    def test_route_fail_closed(self):
        r=proof_route_status(reference_refinement_finite=True,h18_bridge_retained=True,
            a21_linear_uniform=True,a21_nonlinear_bound=False,a21_prefix_retained=True)
        self.assertEqual(r["h18_role"],"finite bridge only"); self.assertFalse(r["route_closed"])
    def test_invalid_inputs(self):
        with self.assertRaises(ValueError): finite_bridge_bound(1,-1,1,0)
        with self.assertRaises(ValueError): release_bound_from_service(-1,1,1,1)
        with self.assertRaises(ValueError): A21TailPremises(-1,0,1,1,1)



class TranslationalObservabilityTests(unittest.TestCase):
    def test_exact_constant_parameter_ou_minor_is_nonzero(self):
        for dt in (.004,.005,.006):
            for tau in (.25,1.0,10.0):
                self.assertTrue(translational_observability_nonsingular(dt,tau))
                self.assertLess(translational_observability_determinant(dt,tau),0)

    def test_invalid_ou_minor_inputs_refused(self):
        with self.assertRaises(ValueError): translational_observability_determinant(0,1)
        with self.assertRaises(ValueError): translational_observability_determinant(.005,0)

class TimeVaryingShippingObservabilityTests(unittest.TestCase):
    def schedule(self):
        return ShippingScheduleBounds(.004,.006,.02,12.0,.015/1.1,.005,.15)

    def test_tau_scaled_scheduler_prevents_decay_collapse(self):
        s=self.schedule()
        self.assertAlmostEqual(pseudo_decay_exponent_per_gap(s),.55,places=12)
        self.assertGreater(math.exp(-3*pseudo_decay_exponent_per_gap(s)),.19)

    def test_time_varying_translation_minor_has_uniform_positive_floor(self):
        floor=time_varying_translation_minor_floor(self.schedule())
        self.assertGreater(floor,0)
        self.assertGreater(floor,7.8e-16)

    def test_gravity_tilt_floor_is_uniform(self):
        self.assertAlmostEqual(gravity_tilt_sensitivity_floor(9.80665,8.8),1.00665)

    def test_cross_coupling_must_fit_block_floors(self):
        self.assertGreater(coupled_information_floor(2.0,3.0,1.0),0)
        self.assertLessEqual(coupled_information_floor(1.0,1.0,1.0),0)

class QuantitativeCertificateTests(unittest.TestCase):
    def test_normalized_translation_floor_is_constructive(self):
        mu=normalized_translation_information_floor(
            word_s=16.0,event_gap_max_s=.156,v_scale=5.5,p_scale=8.1,
            S_scale=1100.0,S_noise_std_max=100.0)
        self.assertGreater(mu,.00204)
        self.assertLess(mu,.00205)

    def test_sparse_two_epoch_bound_is_diagnostic_only(self):
        racc=math.hypot(.2,.75*.3)
        j=vector_attitude_information_floor(
            vertical_specific_force_floor=9.80665-8.8,
            horizontal_specific_force_ceiling=8.8,
            accel_noise_std_ceiling=racc,heading_information_floor=1.0)
        b=rotation_integral_singular_floor(1.0,.6108652381980153)
        mu_ag=two_epoch_attitude_gyro_floor(j,.02,b)
        self.assertGreater(mu_ag,0.0)
        # This is deliberately not composed into the theorem's full mu_N:
        # magnetic service is a transported 2-D window Gramian, not an
        # instantaneous pure-heading row.

    def test_information_adds_over_guaranteed_events(self):
        self.assertEqual(aggregate_repeated_row_information_floor(.25,16),4.0)
        self.assertEqual(aggregate_translation_information_floor(.002,3),.002)
        c=certified_rho_and_margin(.002)
        self.assertLess(c["rho0"],.999)
        self.assertGreater(c["nonlinear_norm_margin"],.0009)

    def test_long_word_closes_fixed_coordinate_neutral_floor(self):
        racc=vibration_inflated_accel_std_ceiling(
            nominal_std=.2,vibration_gain=.75,
            detector_residual_rms_ceiling=.3)
        self.assertLess(racc,.302)
        j=root_attitude_information_floor(
            vertical_specific_force_floor=9.80665-8.8,
            horizontal_specific_force_ceiling=8.8,
            accel_std_ceiling=racc,magnetic_heading_floor=1.0)
        self.assertGreater(j,.0129)
        b=rotation_integral_singular_floor(1.0,.6108652381980153)
        self.assertGreater(b,.984)
        mu_ag=attitude_gyro_information_floor_from_windows(
            root_attitude_floor=j,service_windows=2048,
            gyro_bias_coordinate_scale=.02,
            one_window_rotation_integral_floor=b)
        self.assertGreater(mu_ag,.0021)
        mu_t=normalized_translation_information_floor(
            word_s=16.0,event_gap_max_s=.156,v_scale=5.5,p_scale=8.1,
            S_scale=1100.0,S_noise_std_max=100.0)
        mu=full_neutral_information_floor(
            translation_floor=mu_t,attitude_gyro_floor=mu_ag)
        self.assertGreater(mu,.00204)
        self.assertLess(mu,.00205)

    def test_corrected_classical_riccati_bounds_are_fail_closed(self):
        b=classical_riccati_covariance_bounds(
            observability_lower=.1,observability_upper=2.0,
            controllability_lower=.2,controllability_upper=3.0,
            horizon_steps=10)
        self.assertGreater(b["p_min"],0)
        self.assertGreater(b["p_max"],b["p_min"])

    def test_rho_requires_covariance_normalization(self):
        mu=covariance_normalized_information_floor(.002, .5)
        self.assertAlmostEqual(mu,.001)
        self.assertAlmostEqual(rho_from_covariance_normalized_mu(mu),1/1.001)

    def test_transported_service_uses_schur_not_instantaneous_row(self):
        mu=psd_service_schur_floor(
            controlled_block_floor=.1,controlled_block_ceiling=10.0,
            nuisance_block_ceiling=10.0,service_floor=1.0)
        self.assertGreater(mu,0.0)
        self.assertLess(mu,.1)

    def test_float32_kernel_gamma_bound(self):
        self.assertLess(conservative_float32_kernel_error(100,1.0),6e-6)

    def test_release_cannot_start_in_old_point_one_five_bias_ball(self):
        e=bias_release_error_bound(.22516660498395405,0.0)
        self.assertGreater(e,.225)
        self.assertLess(e,.226)
        self.assertAlmostEqual(projection_sector_retained_error_bound(
            .22516660498395405,.4),.6251666049839541)

class ProjectionReleaseTests(unittest.TestCase):
    def test_universal_release_is_not_projection_inactive(self):
        self.assertFalse(release_can_guarantee_projection_inactive(
            physical_bias_bound=.22516660498395405,projection_radius=.4,
            pre_release_estimate_bound=0.0))
        self.assertTrue(projection_sector_is_dissipative(
            truth_bound=.22516660498395405,projection_radius=.4))

    def test_small_gain_budget_is_explicit(self):
        b=small_gain_budget(.81,.05)
        self.assertAlmostEqual(b["available_eta"],.1)
        self.assertAlmostEqual(b["slack"],.05)
        self.assertTrue(b["closed"])

class RiccatiAndNonlinearClosureTests(unittest.TestCase):
    def test_shipping_positive_noise_envelope_is_uniformly_controllable(self):
        self.assertTrue(uniform_controllability_exists(
            dt_min_s=.004,tau_min_s=.02,tau_max_s=12.0,aw_sigma_floor=.05,
            gyro_noise_floor=.00157,gyro_bias_rw_floor=1e-5,
            accel_bias_drive_floor=math.sqrt(2.5e-7)))

    def test_shipping_covariance_events_retain_compactness(self):
        self.assertTrue(shipping_covariance_hard_events_retain_compactness(
            aw_stationary_std_floor=.05,aw_stationary_std_ceiling=4.0,
            release_bias_variance_ceiling=1.0))

    def test_linear_tail_waits_for_hard_event_retention(self):
        self.assertFalse(linear_kalman_tail_exponentially_stable(
            uniform_detectability=True,uniform_controllability=True,
            covariance_hard_events_retained=False))
        self.assertTrue(linear_kalman_tail_exponentially_stable(
            uniform_detectability=True,uniform_controllability=True,
            covariance_hard_events_retained=True))

    def test_constructive_radius_and_additive_supply_formulae(self):
        r=explicit_small_gain_radius(.81,2.0,.1)
        self.assertAlmostEqual(r,.05)
        self.assertTrue(math.isfinite(additive_supply_practical_radius(.81,.05,.5,1e-4)))

    def test_same_history_exogenous_schedule_leaves_smooth_local_remainder(self):
        smooth=smooth_a21_remainder_vanishes_locally(
            projection_sector_certified=True,tuner_exogenous_same_history=True,
            magnetic_reference_exogenous_same_history=True,
            finite_operation_domain=True)
        self.assertTrue(smooth)
        self.assertTrue(nonlinear_small_gain_exists_from_linear(.99,smooth))

class MarineMagneticDiversityTests(unittest.TestCase):
    def test_bounded_velocity_forces_vector_diversity(self):
        threshold=marine_magnetic_diversity_window_min(9.80665,15.0,75.0,5.5)
        self.assertGreater(threshold,5.6); self.assertLess(threshold,5.7)
        floor=marine_magnetic_vector_diversity_floor(8.0,9.80665,15.0,75.0,5.5)
        self.assertAlmostEqual(floor,43.97475)
        self.assertTrue(attitude_information_from_marine_diversity(
            diversity_floor=floor,magnetic_service_floor=1.0,
            magnetic_service_window_s=1.0,angular_rate_ceiling=.6108652382))

    def test_too_short_window_does_not_claim_diversity(self):
        self.assertLess(marine_magnetic_vector_diversity_floor(
            5.0,9.80665,15.0,75.0,5.5),0.0)

class NeutralQuotientDetectabilityTests(unittest.TestCase):
    def test_tilt_gyro_bias_quotient_minor_is_uniform(self):
        self.assertAlmostEqual(abs(tilt_gyro_bias_quotient_minor(1.00665,.004)),
                               1.00665**2*.004)

    def test_attitude_gyro_quotient_is_uniformly_observable(self):
        self.assertTrue(attitude_gyro_quotient_uniformly_observable(
            gravity_floor=1.00665,magnetic_floor=1.0,
            sample_spacing_floor=.004,angular_rate_ceiling=.6108652382))

    def test_neutral_quotient_has_uniform_mu_when_all_components_are_strict(self):
        self.assertTrue(neutral_quotient_uniform_mu_exists(
            translation_minor_floor=1e-15,gravity_tilt_floor=1.00665,
            tilt_gyro_minor_abs_floor=1e-3,magnetic_information_floor=1.0,
            strict_guard_margin=1e-6))

    def test_active_bias_is_strictly_stable_without_observability(self):
        q=active_bias_homogeneous_ratio(16.0,5000.0)
        self.assertLess(q,.994)
        self.assertGreater(detectable_tail_margin(1.0,16.0,5000.0),.003)

class HybridInformationClosureTests(unittest.TestCase):
    def test_active_tilt_bias_stationary_minor_is_strict(self):
        d=active_tilt_bias_minor(1.00665,.5,5000.0)
        self.assertLess(d,0)
        self.assertGreater(abs(d),4e-9)

    def test_compactness_needs_every_strict_margin(self):
        self.assertTrue(compact_uniform_information_exists(
            translation_minor_floor=1e-15,gravity_tilt_floor=1.0,
            magnetic_information_floor=1.0,active_bias_minor_abs_floor=1e-9,
            strict_branch_guard_margin=1e-6))
        self.assertFalse(compact_uniform_information_exists(
            translation_minor_floor=1e-15,gravity_tilt_floor=1.0,
            magnetic_information_floor=1.0,active_bias_minor_abs_floor=1e-9,
            strict_branch_guard_margin=0.0))

    def test_positive_information_has_a_smooth_local_radius(self):
        self.assertTrue(local_nonlinear_radius_exists(.01,0.0))

class A21InformationTests(unittest.TestCase):
    def test_information_floor_gives_strict_linear_ratio(self):
        self.assertAlmostEqual(information_contraction_ratio(1.0),.5)
        self.assertAlmostEqual(nonlinear_margin_from_information(1.0),1-math.sqrt(.5))

    def test_projection_inactive_on_inner_bias_domain(self):
        margin=bias_projection_inactive_margin(.4,.22516660498395405)
        self.assertGreater(margin,.174)
        self.assertLess(margin,.175)

    def test_zero_information_refused(self):
        with self.assertRaises(ValueError): information_contraction_ratio(0)

class RefinementTests(unittest.TestCase):
    def deployed_like(self, **kw):
        d=dict(min_samples=128,min_window_s=30.0,service_window_s=1.0,true_field_norm_lower=20.0,
               true_field_norm_upper=75.0,measurement_residual_norm=2.0,
               true_horizontal_lower=15.0,tilt_error_upper_rad=math.radians(2.0),
               max_norm_ratio_from_mean=.35,min_horizontal_fraction=.05)
        d.update(kw); return RefinementPremises(**d)

    def test_physical_residual_closes_norm_gate(self):
        p=self.deployed_like()
        m=refinement_gate_margins(p)
        self.assertLess(m["norm_ratio_upper"],.35)
        self.assertGreater(m["horizontal_fraction_lower"],.05)
        self.assertTrue(refinement_sample_gate_uniform(p))
        self.assertGreater(math.degrees(refinement_tilt_limit_rad(p)),6.9)
        self.assertLess(math.degrees(refinement_tilt_limit_rad(p)),7.1)

    def test_service_closes_count_and_window_gates(self):
        p=self.deployed_like()
        # Conservative theorem bound: one usable event per one-second service
        # window. Real shipping cadence is much faster but is not needed.
        self.assertEqual(refinement_completion_bound(90.0,p),218.0)

    def test_capture_to_a21_entry_is_finite(self):
        p=self.deployed_like()
        # refinement: 128 s worst case; then at most 250 service windows + 1 s guard
        self.assertEqual(a21_entry_bound(80.0,90.0,p,250,1.0),469.0)

    def test_large_tilt_error_refuses_horizontal_gate(self):
        p=self.deployed_like(tilt_error_upper_rad=math.radians(15.0))
        self.assertFalse(refinement_sample_gate_uniform(p))
        with self.assertRaises(ValueError): refinement_completion_bound(90.0,p)

    def test_residual_too_large_refuses_norm_gate(self):
        p=self.deployed_like(measurement_residual_norm=4.0)
        self.assertFalse(refinement_sample_gate_uniform(p))

if __name__=="__main__": unittest.main()
