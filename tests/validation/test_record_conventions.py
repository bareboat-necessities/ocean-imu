"""Conventions of the shipped wave records, asserted rather than assumed.

The generator azimuth in a record name is the direction the waves come *from*,
so travel-sense scoring has to compare against ``azimuth + 180``.  That was
never written down and it mattered.  The test needs a record and skips when the
simulation data has not been fetched.

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
        found = sorted(directory.glob("wave_data_jonswap_*.csv"))
        if found:
            return found
    return []


def _azimuth_from_name(path: Path) -> float:
    for token in path.stem.split("_"):
        if token.startswith("A"):
            return float(token[1:])
    raise ValueError(f"no azimuth in {path.name}")


class GeneratorAzimuthConventionTests(unittest.TestCase):
    def test_azimuth_is_the_direction_waves_come_from(self):
        records = _records()
        if not records:
            self.skipTest("simulation records not fetched")

        for path in records:
            with self.subTest(record=path.name):
                columns = path.read_text(encoding="utf-8").split("\n", 1)[0].strip().split(",")
                index = {name: position for position, name in enumerate(columns)}
                rows = int(SAMPLE_RATE_HZ * ANALYSIS_SECONDS)
                data = np.loadtxt(
                    path,
                    delimiter=",",
                    skiprows=1,
                    max_rows=rows,
                    usecols=[index["disp_x"], index["disp_y"], index["disp_z"]],
                )
                horizontal = data[:, :2] - data[:, :2].mean(axis=0)
                vertical = data[:, 2] - data[:, 2].mean()

                # Principal axis of the horizontal orbital displacement.
                _, vectors = np.linalg.eigh(np.cov(horizontal.T))
                axis = vectors[:, -1]

                # For a deep-water wave travelling along +axis with z up, the
                # along-axis displacement and the vertical rate are in
                # antiphase, so the correlation picks the propagation end.
                vertical_rate = np.gradient(vertical, 1.0 / SAMPLE_RATE_HZ)
                along = horizontal @ axis
                if float(np.mean(along * vertical_rate)) > 0.0:
                    axis = -axis

                travel_deg = np.degrees(np.arctan2(axis[1], axis[0])) % 360.0
                expected_deg = (_azimuth_from_name(path) + 180.0) % 360.0
                error = (travel_deg - expected_deg + 180.0) % 360.0 - 180.0
                self.assertLess(
                    abs(error),
                    20.0,
                    f"{path.name}: propagation-to direction {travel_deg:.1f} deg is not "
                    f"azimuth + 180 = {expected_deg:.1f} deg",
                )


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

        text = header.read_text(encoding="utf-8")
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
        import re

        try:
            return float(token)
        except ValueError:
            pass
        named = re.search(
            rf"\b{re.escape(token)}\s*=\s*([0-9.eE+-]+)f?\s*;", text
        )
        self.assertIsNotNone(
            named, f"{token} is not a float constant defined in {header.name}"
        )
        return float(named.group(1))

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
                    header, r"FREQ_SMOOTHER_DT\s*=\s*1\.0f\s*/\s*([0-9.]+)f"
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
                             r"ACC_NOISE_FLOOR_SIGMA_DEFAULT\s*=\s*([0-9.]+)f"),
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
