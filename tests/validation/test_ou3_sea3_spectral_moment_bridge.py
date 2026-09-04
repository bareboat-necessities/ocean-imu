from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_sea3_spectral_moment_bridge as bridge  # noqa: E402


class Sea3SpectralMomentBridgeTest(unittest.TestCase):
    def test_surface_bridge_is_replay_free_and_non_promoting(self) -> None:
        payload = bridge.build()
        failures = bridge.validate(payload)
        self.assertEqual(failures, [])
        self.assertEqual(
            payload["schema_version"],
            "OU3_SEA3_SPECTRAL_MOMENT_BRIDGE_V1",
        )
        self.assertFalse(payload["trajectory_replay_used"])
        self.assertFalse(payload["SEA0_full_certificate_promoted"])
        self.assertFalse(payload["gamma_continuum_screen"]["promotion_use"])
        self.assertFalse(
            payload["tuner_bridge_contract"]
            ["surface_elevation_Tz_may_be_substituted_for_tuner_Tz"]
        )
        self.assertTrue(
            payload["tuner_bridge_contract"]["directional_vessel_IMU_RAO_required"]
        )

    def test_pm_exact_ratio_is_inside_gamma_continuum_screen(self) -> None:
        payload = bridge.build()
        lo, hi = payload["gamma_continuum_screen"][
            "surface_elevation_Tz_over_Tp_outer"
        ]
        exact = ((5.0 / 4.0) * math.pi) ** (-0.25)
        self.assertAlmostEqual(
            payload["analytical_lemmas"]["pm_gamma_1_exact_Tz_over_Tp"],
            exact,
            places=11,
        )
        self.assertLess(lo, exact)
        self.assertLess(exact, hi)
        self.assertGreater(lo, 0.70)
        self.assertLess(hi, 0.84)

    def test_multimodal_zero_crossing_period_uses_energy_inverse_square(self) -> None:
        value = bridge.mixture_zero_crossing_period(
            [3.0, 4.0, 0.0],
            [5.0, 10.0, 7.0],
        )
        expected = math.sqrt(25.0 / (9.0 / 25.0 + 16.0 / 100.0))
        self.assertAlmostEqual(value, expected, places=12)
        self.assertGreaterEqual(value, 5.0)
        self.assertLessEqual(value, 10.0)
        self.assertNotAlmostEqual(value, 0.36 * 5.0 + 0.64 * 10.0, places=6)

    def test_unbanded_acceleration_shortcut_is_explicitly_forbidden(self) -> None:
        payload = bridge.build()
        analytical = payload["analytical_lemmas"]
        self.assertFalse(analytical["unbanded_surface_acceleration_variance_finite"])
        self.assertTrue(
            analytical["band_or_directional_response_required_for_acceleration_moments"]
        )

    def test_committed_bridge_artifact_is_current(self) -> None:
        payload = bridge.build()
        failures = bridge.validate(payload)
        payload["validation_pass"] = not failures
        payload["validation_failures"] = failures
        artifact = json.loads(
            (ROOT / "tools/ou3_sea3_spectral_moment_bridge.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(artifact, payload)


if __name__ == "__main__":
    unittest.main()
