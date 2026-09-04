from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_wave_period_frontend as frontend  # noqa: E402


class Sea3WavePeriodFrontendTest(unittest.TestCase):
    def test_frontend_identity_and_startup_semantics_match_shipping_source(self) -> None:
        payload = frontend.build(ROOT)
        failures = frontend.validate(payload)
        self.assertEqual(failures, [])
        self.assertEqual(
            payload["schema_version"],
            "OU3_SEA3_WAVE_PERIOD_FRONTEND_V2",
        )
        self.assertTrue(all(payload["source_parity"].values()))
        self.assertTrue(all(payload["operating_domain_parity"].values()))
        exact = payload["exact_transfer_identity"]
        self.assertTrue(exact["two_shared_high_pass_stages_cancel_from_ratio"])
        self.assertTrue(exact["source_leak_square_subtraction_used"])

    def test_validated_warping_is_small_but_not_promoted_to_tuner_period(self) -> None:
        payload = frontend.build(ROOT)
        intervals = payload["validated_intervals"]
        omega_lo, omega_hi = intervals["omega_hat_over_omega"]
        period_lo, period_hi = intervals["period_hat_over_period"]

        self.assertGreater(omega_lo, 0.999)
        self.assertLessEqual(omega_hi, 1.001)
        self.assertGreater(period_lo, 0.999)
        self.assertLess(period_hi, 1.001)
        self.assertTrue(intervals["validated_transcendentals_used"])
        self.assertFalse(
            intervals["ordinary_libm_transcendentals_used_for_enclosure"]
        )
        self.assertFalse(payload["SEA0_full_certificate_promoted"])
        self.assertFalse(payload["P2_promoted_from_this_artifact"])
        self.assertFalse(
            payload["interpretation"]
            ["surface_Tz_or_sinusoid_period_may_replace_tuner_Tz"]
        )

    def test_prior_and_estimator_takeover_are_separate_source_modes(self) -> None:
        payload = frontend.build(ROOT)
        startup = payload["startup_source_language"]
        settle_lo, settle_hi = payload["validated_intervals"][
            "wave_period_integrator_settle_lower_bound_s"
        ]

        self.assertFalse(startup["tuner_ready_requires_wave_period_estimator_ready"])
        self.assertTrue(startup["live_entry_may_precede_wave_period_estimator_first_valid_period"])
        self.assertFalse(startup["wave_period_takeover_waits_for_isReady"])
        self.assertTrue(startup["tuner_consumes_previous_sample_wave_period_state"])
        self.assertTrue(startup["current_sample_wave_period_update_occurs_after_tuner_update"])
        self.assertTrue(
            startup[
                "first_newly_finite_wave_period_can_affect_tuner_no_earlier_than_next_valid_sample"
            ]
        )
        self.assertTrue(
            startup[
                "first_valid_tuner_update_can_satisfy_debiased_variance_ready_threshold"
            ]
        )
        self.assertGreater(settle_lo, 47.0)
        self.assertLess(settle_hi, 49.0)

    def test_remaining_estimator_obligations_are_explicit(self) -> None:
        payload = frontend.build(ROOT)
        interpretation = payload["interpretation"]
        self.assertFalse(interpretation["discrete_frontend_warping_is_current_limiter"])
        self.assertTrue(interpretation["startup_prior_to_estimator_takeover_is_now_source_certified"])
        self.assertTrue(interpretation["multimodal_response_moment_enclosure_still_required"])
        self.assertTrue(interpretation["finite_EW_moment_transient_still_required"])
        self.assertTrue(interpretation["canonical_log_period_EMA_still_required"])
        self.assertTrue(interpretation["target_float_libm_rounding_still_required"])


if __name__ == "__main__":
    unittest.main()
