from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_wave_period_frontend as frontend  # noqa: E402


class Sea3WavePeriodFrontendTest(unittest.TestCase):
    def test_frontend_identity_matches_shipping_source(self) -> None:
        payload = frontend.build(ROOT)
        failures = frontend.validate(payload)
        self.assertEqual(failures, [])
        self.assertTrue(all(payload["source_parity"].values()))
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

    def test_remaining_estimator_obligations_are_explicit(self) -> None:
        payload = frontend.build(ROOT)
        interpretation = payload["interpretation"]
        self.assertFalse(interpretation["discrete_frontend_warping_is_current_limiter"])
        self.assertTrue(interpretation["multimodal_response_moment_enclosure_still_required"])
        self.assertTrue(interpretation["finite_EW_moment_transient_still_required"])
        self.assertTrue(interpretation["canonical_log_period_EMA_still_required"])
        self.assertTrue(interpretation["target_float_libm_rounding_still_required"])


if __name__ == "__main__":
    unittest.main()
