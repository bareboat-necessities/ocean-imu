"""Late contract for the analytical OU-III robustness coupling.

The archived robustness bundle predates SpectralMSE and remains byte-for-byte
restatable.  New study runs, however, must couple sigma/tau perturbations to the
deployed analytical law rather than the historical cubic relation.
"""
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))

import ou_robustness as robustness  # noqa: E402
import ou_validation as validation  # noqa: E402
import test_ou_robustness as robustness_tests  # noqa: E402
import test_ou_robustness_bounds as bound_tests  # noqa: E402


def _coupled_sweeps_follow_deployed_spectral_mse(self):
    baseline = validation.TuningPoint(tau_s=1.2, sigma_a_mps2=0.8, RS_ms=3.5)

    coupled_sigma = robustness.scaled_tuning_point(
        baseline, "sigma_aw_rs", 0.5
    )
    self.assertAlmostEqual(coupled_sigma.tau_s, 1.2)
    self.assertAlmostEqual(coupled_sigma.sigma_a_mps2, 0.4)
    self.assertAlmostEqual(coupled_sigma.RS_ms, 3.5 * 0.5 ** (6.0 / 7.0))

    coupled_tau = robustness.scaled_tuning_point(
        baseline, "tau_rs", 1.5
    )
    expected_ratio = robustness.spectral_mse_tau_ratio(1.2, 1.8)
    self.assertAlmostEqual(coupled_tau.tau_s, 1.8)
    self.assertAlmostEqual(coupled_tau.sigma_a_mps2, 0.8)
    self.assertAlmostEqual(coupled_tau.RS_ms, 3.5 * expected_ratio)

    # This nominal point is away from cadence clamps, so the exact ratio must
    # collapse to the cadence-normalized exponent 41/14.
    self.assertAlmostEqual(expected_ratio, 1.5 ** (41.0 / 14.0), places=6)
    self.assertEqual(
        robustness.SENSITIVITY_PARAMETERS,
        robustness.OFAT_PARAMETERS + robustness.COUPLED_PARAMETERS,
    )


def _sigma_coupled_half_scale_stays_above_floor(self):
    point = robustness.scaled_tuning_point(bound_tests._baseline(), "sigma_aw_rs", 0.5)
    expected = bound_tests._baseline().RS_ms * 0.5 ** (6.0 / 7.0)
    self.assertGreater(point.RS_ms, bound_tests.R_S_FLOOR)
    self.assertAlmostEqual(point.RS_ms, expected)
    self.assertAlmostEqual(point.sigma_a_mps2, 0.4)
    robustness.validate_tuning_point(point)


def _tau_coupled_half_scale_realizes_floor(self):
    point = robustness.scaled_tuning_point(bound_tests._baseline(), "tau_rs", 0.5)
    self.assertAlmostEqual(point.tau_s, 0.6)
    self.assertAlmostEqual(point.RS_ms, bound_tests.R_S_FLOOR)
    robustness.validate_tuning_point(point)


def _tau_floor_is_only_clipping(self):
    baseline = bound_tests._baseline()
    unclipped = baseline.RS_ms * robustness.spectral_mse_tau_ratio(
        baseline.tau_s, baseline.tau_s * 0.5
    )
    self.assertLess(unclipped, bound_tests.R_S_FLOOR)


robustness_tests.RobustnessDesignTests.test_coupled_sweeps_follow_the_deployed_regularization_law = (
    _coupled_sweeps_follow_deployed_spectral_mse
)
bound_tests.RobustnessBoundTests.test_sigma_coupled_half_scale_stays_above_the_r_s_floor = (
    _sigma_coupled_half_scale_stays_above_floor
)
bound_tests.RobustnessBoundTests.test_tau_coupled_half_scale_realizes_r_s_floor = (
    _tau_coupled_half_scale_realizes_floor
)
bound_tests.RobustnessBoundTests.test_the_floor_is_the_only_thing_clipping_that_request = (
    _tau_floor_is_only_clipping
)


class AnalyticalRobustnessCouplingContract(unittest.TestCase):
    def test_active_driver_states_current_spectral_mse_powers(self):
        text = (REPO_ROOT / "tools" / "ou_robustness.py").read_text(encoding="utf-8")
        self.assertIn("6.0 / 7.0", text)
        self.assertIn("24.0 / 7.0", text)
        self.assertIn("41/14", text)
        self.assertNotIn("scale**3", text)
        self.assertNotIn("sigma_aw * tau^3", text)

    def test_exact_tau_ratio_handles_cadence_clamps(self):
        old_tau = 0.30
        new_tau = 0.15
        exact = robustness.spectral_mse_tau_ratio(old_tau, new_tau)
        # Both points hit the 5-ms cadence floor, so no -1/2 tau exponent is
        # contributed by cadence and the ratio is the base 24/7 power.
        self.assertAlmostEqual(exact, (new_tau / old_tau) ** (24.0 / 7.0), places=6)


if __name__ == "__main__":
    unittest.main()
