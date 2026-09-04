from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_directional_response_moments as response  # noqa: E402
import ou3_sea3_spectral_moment_bridge as bridge  # noqa: E402


class Sea3DirectionalResponseMomentsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = response.build()

    def test_subcertificate_is_replay_free_and_non_promoting(self) -> None:
        payload = self.payload
        self.assertEqual(response.validate(payload), [])
        self.assertEqual(
            payload["schema_version"],
            "OU3_SEA3_DIRECTIONAL_RESPONSE_MOMENTS_V1",
        )
        self.assertFalse(payload["trajectory_replay_used"])
        self.assertFalse(payload["filter_changed"])
        self.assertFalse(payload["declared_operating_domain_shrunk"])
        self.assertFalse(payload["SEA0_full_certificate_promoted"])
        for stage in ("P2", "P3", "P4", "P5"):
            self.assertFalse(payload[f"{stage}_promoted_from_this_artifact"])
        self.assertFalse(
            payload["sea3_source_language_included_in_frozen_p2_contract"]
        )

    def test_sea_domain_is_not_narrower_than_the_declared_deployment_envelope(
        self,
    ) -> None:
        domain = json.loads(
            (ROOT / "tools/ou3_proof_operating_domain.json").read_text(
                encoding="utf-8"
            )
        )
        declared_hs = domain["initial_filter_entrance"]["position"][
            "significant_wave_height_Hs_upper_m"
        ]
        self.assertGreaterEqual(response.HS_MAX_M, declared_hs)
        self.assertFalse(self.payload["declared_operating_domain_shrunk"])

    def test_sea_family_matches_the_frozen_spectral_bridge(self) -> None:
        # One JONSWAP convention across every SEA0 artifact.
        self.assertEqual(response.M_MAX, bridge.M_MAX)
        self.assertEqual(response.GAMMA_MIN, bridge.GAMMA_MIN)
        self.assertEqual(response.GAMMA_MAX, bridge.GAMMA_MAX)
        elevation = response.partition_elevation_spectrum(4.0, 10.0, 1.0)
        # m0 = H^2/16 = 1 m^2 by construction.
        self.assertAlmostEqual(response._integrate(elevation), 1.0, places=9)
        surface_ratio = (
            response.surface_zero_crossing_period(elevation) / 10.0
        )
        self.assertAlmostEqual(
            surface_ratio,
            ((5.0 / 4.0) * math.pi) ** -0.25,
            places=9,
        )

    def test_directional_spreading_is_normalized_and_axis_selective(self) -> None:
        head = response.directional_gram(0.0, 8.0)
        beam = response.directional_gram(math.radians(90.0), 8.0)
        self.assertAlmostEqual(head[0][0], 1.0, places=9)
        self.assertAlmostEqual(beam[0][0], 1.0, places=9)
        # A head sea drives surge, a beam sea drives sway.
        self.assertGreater(head[1][1], head[2][2])
        self.assertGreater(beam[2][2], beam[1][1])
        # cos^2 spreading at s = 1 has closed-form Gram entries.
        broad = response.directional_gram(0.0, 1.0)
        self.assertAlmostEqual(broad[1][1], 0.5, places=9)
        self.assertAlmostEqual(broad[0][1], 2.0 / math.pi, places=9)

    def test_cross_axis_coupling_survives_the_response_enclosure(self) -> None:
        moments = self.payload["matrix_moments"]
        self.assertGreater(moments["worst_normalized_offdiagonal"], 0.5)
        self.assertFalse(moments["per_axis_scalarization_valid"])
        for row in moments["example_acceleration_matrix_moment_m2s4"]:
            for value in row:
                self.assertTrue(math.isfinite(value))

    def test_leak_inversion_identity_is_exact_not_narrow_band(self) -> None:
        # omega_est^2 = var_v/var_eta - lambda^2 must equal the mu-weighted
        # mean square frequency for a broadband, strongly multimodal input.
        totals = [0.0] * len(response.OMEGA_NODES)
        for height, peak_period, gamma in (
            (3.0, 4.0, 1.0),
            (2.0, 13.0, 7.0),
        ):
            spectrum = response.partition_elevation_spectrum(
                height, peak_period, gamma
            )
            totals = [a + b for a, b in zip(totals, spectrum)]

        source = response.up_specific_force_spectrum(totals)
        proxy_variance = response._integrate(
            [
                w * s
                for w, s in zip(response.PROXY_ELEVATION_WEIGHT, source)
            ]
        )
        velocity_variance = response.deployed_velocity_proxy_variance(totals)
        inverted = (
            velocity_variance / proxy_variance
            - response.LAMBDA_LEAK * response.LAMBDA_LEAK
        )
        from_identity = (
            2.0 * math.pi / response.deployed_period_s(totals)
        ) ** 2
        self.assertAlmostEqual(inverted / from_identity, 1.0, places=9)

    def test_mixture_period_stays_inside_single_partition_extremes(self) -> None:
        convexity = self.payload["mixture_convexity_check"]
        self.assertLessEqual(convexity["worst_relative_excess"], 1e-9)
        for case in convexity["cases"]:
            self.assertTrue(case["inside_component_extremes"])
        self.assertTrue(
            self.payload["deployed_period_channel"][
                "extremes_valid_for_three_partition_class"
            ]
        )

    def test_deployed_period_is_not_the_surface_period(self) -> None:
        period = self.payload["deployed_period_channel"]
        bias_lo, bias_hi = period["deployed_over_surface_Tz"]
        self.assertLess(bias_lo, 1.0)
        self.assertGreater(bias_hi, 1.0)
        deployed_lo, deployed_hi = period["deployed_Tz_over_Tp"]
        surface_lo, surface_hi = bridge.build()["gamma_continuum_screen"][
            "surface_elevation_Tz_over_Tp_outer"
        ]
        # The response and the leak widen the ratio on both sides, so the
        # surface screen may not be substituted for the deployed one.
        self.assertLess(deployed_lo, surface_lo)
        self.assertGreater(deployed_hi, surface_hi)

    def test_induced_tuner_frequency_stays_in_the_committed_channel(self) -> None:
        period = self.payload["deployed_period_channel"]
        low, high = period["induced_tuner_frequency_hz"]
        self.assertGreaterEqual(low, response.MIN_TUNE_FREQ_HZ)
        self.assertLessEqual(high, response.MAX_TUNE_FREQ_HZ)
        self.assertTrue(period["inside_committed_tuning_channel"])

    def test_sampling_alone_does_not_band_limit_the_sea(self) -> None:
        alias = self.payload["alias_obligation"]
        self.assertFalse(alias["sampling_alone_band_limits_the_sea"])
        self.assertTrue(alias["response_rolloff_is_mandatory"])
        self.assertTrue(alias["flat_response_folded_power_is_divergent"])
        self.assertTrue(alias["flat_response_exceeds_sigma_clamp_within_decades"])
        # The declared roll-off must reduce the fold to a negligible quantity.
        self.assertTrue(
            math.isfinite(alias["folded_power_with_declared_rolloff_m2s4"])
        )
        self.assertLess(
            alias["folded_sigma_with_declared_rolloff_ms2"],
            0.01 * response.SIGMA_AW_CLAMP_MS2,
        )

    def test_flat_response_alias_power_grows_without_bound(self) -> None:
        constant = self.payload["alias_obligation"][
            "asymptotic_acceleration_constant_m2s3"
        ]
        three = response.unbounded_response_alias_power(3.0, constant)
        six = response.unbounded_response_alias_power(6.0, constant)
        self.assertAlmostEqual(six / three, 2.0, places=9)

    def test_sigma_clamp_contains_the_source_coordinate_by_saturation(self) -> None:
        sigma = self.payload["deployed_sigma_channel"]
        self.assertTrue(sigma["per_partition_and_mixture_steepness_both_imposed"])
        self.assertGreaterEqual(
            sigma["worst_admissible_sigma_a_ms2"],
            sigma["single_partition_max_sigma_a_ms2"],
        )
        # sigma_aw cannot leave its interval whatever the sea does.
        self.assertLessEqual(
            sigma["worst_admissible_sigma_aw_ms2"], response.SIGMA_AW_CLAMP_MS2
        )
        self.assertTrue(sigma["sigma_aw_coordinate_stays_inside_declared_interval"])
        self.assertAlmostEqual(
            sigma["sigma_a_saturating_the_clamp_ms2"],
            response.SIGMA_AW_CLAMP_MS2 / response.SIGMA_COEFF,
            places=12,
        )

    def test_clamp_rail_is_reachable_only_off_the_settled_reference(self) -> None:
        sigma = self.payload["deployed_sigma_channel"]
        rail = sigma["sigma_a_saturating_the_clamp_ms2"]
        self.assertTrue(sigma["clamp_rail_reachable_from_declared_sea_domain"])
        self.assertGreaterEqual(sigma["worst_admissible_sigma_a_ms2"], rail)
        self.assertTrue(sigma["rail_makes_sea_to_source_map_non_injective"])
        # Neither frequency-source mode at rest reaches the rail on its own.
        self.assertFalse(sigma["clamp_rail_reachable_at_settled_band_reference"])
        self.assertLess(sigma["self_consistent_max_sigma_a_ms2"], rail)
        self.assertFalse(sigma["clamp_rail_reachable_at_startup_prior"])
        self.assertLess(sigma["startup_prior_max_sigma_a_ms2"], rail)
        # Reachability stays expressed as a response/domain quantity.
        self.assertLess(
            sigma["rho_max_below_which_the_rail_is_unreachable"], response.RHO_MAX
        )
        self.assertLess(
            sigma["steepness_below_which_the_rail_is_unreachable"],
            response.SIGNIFICANT_STEEPNESS_MAX,
        )

    def test_band_corners_follow_the_shipping_bandpass_clamps(self) -> None:
        # Well inside the clamps the corners are exact ratio multiples.
        low, high = response.sigma_band_corners_hz(0.4)
        self.assertAlmostEqual(low, 0.2, places=12)
        self.assertAlmostEqual(high, 1.6, places=12)
        # The 6 Hz absolute ceiling is unreachable inside the committed
        # channel: the high corner tops out at 4 * 1.2 Hz.
        _, high = response.sigma_band_corners_hz(response.MAX_TUNE_FREQ_HZ)
        self.assertAlmostEqual(
            high,
            response.SIGMA_BAND_HIGH_RATIO * response.MAX_TUNE_FREQ_HZ,
            places=12,
        )
        self.assertLess(high, response.SIGMA_BAND_MAX_HZ)
        self.assertFalse(
            self.payload["deployed_sigma_channel"][
                "absolute_band_ceiling_reachable_in_committed_channel"
            ]
        )
        # The absolute floor wins at low reference frequencies.
        low, _ = response.sigma_band_corners_hz(0.01)
        self.assertAlmostEqual(low, response.SIGMA_BAND_MIN_HZ, places=12)

    def test_committed_artifact_is_current(self) -> None:
        artifact = json.loads(
            (ROOT / "tools/ou3_sea3_directional_response_moments.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(artifact, self.payload)


if __name__ == "__main__":
    unittest.main()
