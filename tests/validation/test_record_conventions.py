"""Conventions of the shipped wave records, asserted rather than assumed.

The underway records retain incident propagation-to azimuth in filenames.
Their orbital correlations include signed encounter frequency and changing
heading, so stationary-hull phase matching cannot validate propagation sign.
Check the delivered total-motion kinematics and declared track instead.

The second contract here is the Python mirror of the deployed fixed-tuning
operating point, pinned against the shipping headers.
"""

import sys
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))

SAMPLE_RATE_HZ = 200.0
ANALYSIS_SECONDS = 400.0


def _records():
    for directory in (Path(__file__).resolve().parent, REPO_ROOT / "tests" / "kalman_ou_iii"):
        found = sorted(directory.glob("wave_data_jonswap_*.csv")) + sorted(directory.glob("wave_data_pmstokes_*.csv"))
        if found:
            return found
    return []


def _azimuth_from_name(path: Path) -> float:
    for token in path.stem.split("_"):
        if token.startswith("A"):
            return float(token[1:])
    raise ValueError(f"no azimuth in {path.name}")


class GeneratorMotionConventionTests(unittest.TestCase):
    def test_underway_records_have_consistent_total_motion_and_heading(self):
        from underway_reference import track
        records = _records()
        if not records:
            self.skipTest("simulation records not fetched")
        for path in records:
            with self.subTest(record=path.name):
                with path.open() as stream:
                    columns = stream.readline().strip().split(',')
                rows = int(SAMPLE_RATE_HZ * ANALYSIS_SECONDS)
                data = np.loadtxt(path, delimiter=',', skiprows=1, max_rows=rows)
                times = data[:, columns.index('time')]
                position = data[:, [columns.index('disp_' + a) for a in 'xyz']]
                velocity = data[:, [columns.index('vel_' + a) for a in 'xyz']]
                acceleration = data[:, [columns.index('acc_' + a) for a in 'xyz']]
                # These are independent numerical derivatives of the raw total
                # history, including forward travel and steering acceleration.
                dp = np.gradient(position, 1 / SAMPLE_RATE_HZ, axis=0)
                dv = np.gradient(velocity, 1 / SAMPLE_RATE_HZ, axis=0)
                # Upstream exports total travel through float32: its north
                # position ULP grows with distance. Bound derivative error
                # by that quantization, independently of estimator accuracy.
                ulp = np.abs(np.spacing(position.astype(np.float32))).astype(float)
                bound = .0001 + 2 * ulp / (2 / SAMPLE_RATE_HZ)
                self.assertTrue(np.all(np.abs(dp[2:-2] - velocity[2:-2]) <= bound[2:-2]))
                np.testing.assert_allclose(dv[2:-2], acceleration[2:-2], atol=.003, rtol=.003)
                cruise, _, _, heading = track(times)
                np.testing.assert_allclose(data[:, columns.index('yaw_deg')], heading, atol=2e-6)
                # Total displacement includes kilometre-scale travel; its
                # wave residual remains bounded rather than accumulating it.
                self.assertGreater(position[-1, 0] - position[0, 0], 900)
                self.assertLess(np.max(np.abs(position - cruise)), 20)


class DeployedLawMirrorTests(unittest.TestCase):
    """The fixed-tuning modes derive their frozen operating point in Python.

    That derivation mirrors the C++ adaptation law.  When the two drift apart
    every fixed mode is scored at a point the deployed filter would never
    choose, silently: nothing fails, the numbers are just wrong.  This pins the
    mirror to the header.
    """

    HEADER = REPO_ROOT / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"
    HEADER_OU_II = REPO_ROOT / "src" / "kalman_ou_ii" / "SeaStateFusionFilter_OU_II.h"

    def _header_value(self, pattern: str) -> float:
        return self._value_from(self.HEADER, pattern)

    def _value_from(self, header, pattern: str) -> float:
        import re

        from wrapper_sources import wrapper_source

        # The filter together with the shared kalman_common headers it is
        # built from: shared defaults are defined there once.
        family = "ou2" if header == self.HEADER_OU_II else "ou3"
        text = wrapper_source(family)
        match = re.search(pattern, text)
        self.assertIsNotNone(match, f"{pattern} not found in {header.name}")
        return self._resolve(text, match.group(1), header)

    def _resolve(self, text: str, token: str, header) -> float:
        """Read a float, following one level of named constant.

        A default written as a named constant says why it has the value it has,
        which is worth more than a bare literal -- but it means the mirror test
        cannot simply float() what it matched.  Resolve the name against the
        same header rather than forcing the constant back into a literal.
        """
        from wrapper_sources import resolved_value

        try:
            return float(token)
        except ValueError:
            pass
        named = resolved_value(text, token.split("::")[-1])
        self.assertIsNotNone(
            named, f"{token} is not a float constant defined in {header.name}"
        )
        return named

    def test_rs_coefficient_matches_the_filter_default(self):
        import ou_validation as validation

        self.assertAlmostEqual(
            validation.OU_III_RS_COEFF,
            self._header_value(r"float\s+R_S_coeff_\s*=\s*([A-Za-z0-9._]+?)f?\s*;"),
            places=6,
        )

    def test_mse_coefficient_and_sigma_coeff_match_the_filter_defaults(self):
        """The deployed law is SpectralMSE, so its C_J is the one that matters.

        c_sigma is mirrored too because the law divides it back out to recover
        the physical band RMS; a drift in either would move every fixed-tuning
        operating point.
        """
        import ou_validation as validation

        self.assertAlmostEqual(
            validation.OU_III_RS_MSE_COEFF,
            self._header_value(r"R_S_MSE_COEFF_DEFAULT\s*=\s*([A-Za-z0-9._]+?)f?\s*;"),
            places=6,
        )
        self.assertAlmostEqual(
            validation.OU_III_SIGMA_COEFF,
            self._header_value(r"float\s+sigma_coeff_\s*=\s*([A-Za-z0-9._]+?)f?\s*;"),
            places=6,
        )

    def test_rs_acceleration_scale_matches_the_filter_default(self):
        """The OU-III base schedule scales with sqrt(R_a), not with sigma_aw.

        The mirror carries sqrt(R_a) as its own constant, so it can drift from
        the header exactly the way the coefficient can.  The header states the
        reduced density r_a = R_a * h; recover sqrt(R_a) from it.
        """
        import ou_validation as validation

        sigma_a = self._header_value(
            r"R_S_ACCEL_NOISE_DENSITY_DEFAULT\s*=\s*([0-9.]+)f\s*\*"
        )
        self.assertAlmostEqual(
            validation.OU_III_RS_ACCEL_SIGMA_MPS2, sigma_a, places=6
        )

    def test_rs_bounds_match_the_filter_clamps(self):
        import ou_validation as validation

        self.assertAlmostEqual(
            validation.OU_III_RS_BOUNDS_MS[0],
            self._header_value(r"MIN_R_S\s*=\s*([0-9.]+)f"),
            places=6,
        )
        self.assertAlmostEqual(
            validation.OU_III_RS_BOUNDS_MS[1],
            self._header_value(r"MAX_R_S\s*=\s*([0-9.]+)f"),
            places=6,
        )

    def test_pseudo_update_cadence_bounds_match_the_filter_clamps(self):
        """T_S sits under a square root in both laws, so its clamps are load bearing.

        The mirror clips the pseudo-update period itself, so a ceiling that
        drifts from the header silently moves r_S on every tau long enough to
        reach it -- above 11 s for OU-III.  The two families carry the ceiling
        independently, so pin both.
        """
        import ou_validation as validation

        for bounds, header in (
            (validation.OU_III_PSEUDO_PERIOD_BOUNDS_S, self.HEADER),
            (validation.OU_II_PSEUDO_PERIOD_BOUNDS_S, self.HEADER_OU_II),
        ):
            with self.subTest(header=header.name):
                # The floor is the IMU schedule itself, written as a rate.
                rate = self._value_from(
                    header, r"NOMINAL_IMU_DT_S\s*=\s*1\.0f\s*/\s*([0-9.]+)f"
                )
                self.assertAlmostEqual(bounds[0], 1.0 / rate, places=9)
                self.assertAlmostEqual(
                    bounds[1],
                    self._value_from(
                        header,
                        r"PSEUDO_UPDATE_PERIOD_MAX_S_DEFAULT\s*=\s*([0-9.]+)f",
                    ),
                    places=6,
                )

    def test_ou_ii_coefficients_match_the_filter_defaults(self):
        """OU-II's mirror has the same failure mode and had no test.

        Its two coefficients are the ones a re-fit moves, so pin both.
        """
        import ou_validation as validation

        self.assertAlmostEqual(
            validation.OU_II_RP0_COEFF,
            self._value_from(self.HEADER_OU_II,
                             r"float\s+R_p0_coeff_\s*=\s*([0-9.]+)f"),
            places=6,
        )
        self.assertAlmostEqual(
            validation.OU_II_RV0_COEFF,
            self._value_from(self.HEADER_OU_II,
                             r"float\s+R_v0_coeff_\s*=\s*([0-9.]+)f"),
            places=6,
        )

    def test_ou_ii_mse_law_constants_match_the_filter_defaults(self):
        """The deployed OU-II law is PhysicalMSE, so its constants matter most.

        C_P, the channel ratio C_P/C_V and c_sigma all move every fixed-tuning
        operating point; the law divides c_sigma back out to recover the
        physical band RMS, so a drift in it is as damaging as a drift in C_P.
        """
        import ou_validation as validation

        for value, pattern in (
            (validation.OU_II_PSEUDO_MSE_COEFF,
             r"R_PSEUDO_MSE_COEFF_DEFAULT\s*=\s*([0-9.]+)f"),
            (validation.OU_II_PSEUDO_MSE_RATIO,
             r"R_PSEUDO_MSE_RATIO_DEFAULT\s*=\s*([0-9.]+)f"),
            (validation.OU_II_SIGMA_COEFF,
             r"float\s+sigma_coeff_\s*=\s*([0-9.]+)f"),
        ):
            self.assertAlmostEqual(
                value, self._value_from(self.HEADER_OU_II, pattern), places=6)

    def test_ou_ii_mse_law_acceleration_density_matches_the_filter(self):
        """q_eff = 2 r_a is mirrored as its own constant and can drift alone.

        The header states it as the square of the acceleration noise floor, so
        the mirror is checked against that constant rather than a literal: the
        two are the same physical number and a re-characterization of the
        platform has to move them together.
        """
        import ou_validation as validation

        sigma_a = self._value_from(
            self.HEADER_OU_II,
            r"R_PSEUDO_ACCEL_NOISE_DENSITY_DEFAULT\s*=\s*\n?\s*"
            r"([A-Za-z0-9._]+?)f?\s*\*",
        )
        self.assertAlmostEqual(
            validation.OU_II_PSEUDO_QEFF, 2.0 * sigma_a ** 2 / 200.0, places=12
        )
        self.assertAlmostEqual(
            sigma_a,
            self._value_from(self.HEADER_OU_II,
                             r"\bACC_NOISE_FLOOR_SIGMA_DEFAULT\s*=\s*([A-Za-z0-9_:.]+?)f?\s*;"),
            places=6,
        )

    def test_ou_ii_bounds_match_the_filter_clamps(self):
        import ou_validation as validation

        for value, pattern in (
            (validation.OU_II_RP0_BOUNDS_M[0], r"MIN_R_p0_std\s*=\s*([0-9.]+)f"),
            (validation.OU_II_RP0_BOUNDS_M[1], r"MAX_R_p0_std\s*=\s*([0-9.]+)f"),
            (validation.OU_II_RV0_BOUNDS_MPS[0], r"MIN_R_v0_std\s*=\s*([0-9.]+)f"),
            (validation.OU_II_RV0_BOUNDS_MPS[1], r"MAX_R_v0_std\s*=\s*([0-9.]+)f"),
        ):
            self.assertAlmostEqual(
                value, self._value_from(self.HEADER_OU_II, pattern), places=6)

    def test_robustness_bounds_match_the_filter_clamps(self):
        import ou_robustness as robustness

        self.assertAlmostEqual(
            robustness.TAU_BOUNDS_S[1],
            self._header_value(r"MAX_TAU_S\s*=\s*([0-9.]+)f"),
            places=6,
        )
        self.assertAlmostEqual(
            robustness.R_S_BOUNDS_MS[1],
            self._header_value(r"MAX_R_S\s*=\s*([0-9.]+)f"),
            places=6,
        )


if __name__ == "__main__":
    unittest.main()
