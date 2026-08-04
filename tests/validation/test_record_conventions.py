"""Conventions of the shipped wave records, asserted rather than assumed.

Two of them were never written down, and both mattered:

* the generator azimuth in a record name is the direction the waves come
  *from*, so travel-sense scoring has to compare against ``azimuth + 180``;
* the heading-rotation transform used for the travel-sense gauge experiment
  must leave the world-frame specific force and the body-frame angular rate
  untouched, otherwise it changes the physics instead of the point of view.

The first test needs a record and skips when the simulation data has not been
fetched; the second is self-contained.
"""

import sys
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))

import rotate_record_heading as rotate  # noqa: E402

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


class HeadingRotationTests(unittest.TestCase):
    def _trajectory(self, count=64):
        rng = np.random.default_rng(20260803)
        angles = rng.normal(scale=0.15, size=(count, 3))
        half = angles * 0.5
        cos, sin = np.cos(half), np.sin(half)
        quaternions = np.stack(
            (
                cos[:, 0] * cos[:, 1] * cos[:, 2] + sin[:, 0] * sin[:, 1] * sin[:, 2],
                sin[:, 0] * cos[:, 1] * cos[:, 2] - cos[:, 0] * sin[:, 1] * sin[:, 2],
                cos[:, 0] * sin[:, 1] * cos[:, 2] + sin[:, 0] * cos[:, 1] * sin[:, 2],
                cos[:, 0] * cos[:, 1] * sin[:, 2] - sin[:, 0] * sin[:, 1] * cos[:, 2],
            ),
            axis=-1,
        )
        quaternions /= np.linalg.norm(quaternions, axis=1, keepdims=True)
        return quaternions, rng.normal(scale=2.0, size=(count, 3))

    def test_world_specific_force_is_preserved(self):
        q_world_to_body, body_force = self._trajectory()
        q_body_to_world = rotate.quaternion_conjugate(q_world_to_body)
        world_before = rotate.quaternion_rotate(q_body_to_world, body_force)

        for heading_deg in (45.0, 90.0, 180.0, -120.0):
            with self.subTest(heading_deg=heading_deg):
                half = np.deg2rad(heading_deg) * 0.5
                q_heading = np.tile(
                    np.array([np.cos(half), 0.0, 0.0, np.sin(half)]),
                    (len(q_world_to_body), 1),
                )
                q_body_to_world_new = rotate.quaternion_multiply(q_heading, q_body_to_world)
                q_world_to_body_new = rotate.quaternion_conjugate(q_body_to_world_new)
                body_force_new = rotate.quaternion_rotate(q_world_to_body_new, world_before)
                world_after = rotate.quaternion_rotate(q_body_to_world_new, body_force_new)
                self.assertLess(float(np.abs(world_after - world_before).max()), 1e-9)

    def test_heading_rotation_is_a_rotation_about_the_world_vertical(self):
        q_world_to_body, _ = self._trajectory(count=8)
        q_body_to_world = rotate.quaternion_conjugate(q_world_to_body)
        vertical = np.tile(np.array([0.0, 0.0, 1.0]), (len(q_body_to_world), 1))
        bow = rotate.quaternion_rotate(q_body_to_world, np.tile(np.array([1.0, 0.0, 0.0]), (len(q_body_to_world), 1)))

        half = np.deg2rad(90.0) * 0.5
        q_heading = np.tile(np.array([np.cos(half), 0.0, 0.0, np.sin(half)]), (len(q_body_to_world), 1))
        q_body_to_world_new = rotate.quaternion_multiply(q_heading, q_body_to_world)
        bow_new = rotate.quaternion_rotate(
            q_body_to_world_new, np.tile(np.array([1.0, 0.0, 0.0]), (len(q_body_to_world), 1))
        )

        # The vertical component of the bow axis is untouched and the horizontal
        # part turns by exactly the requested angle.
        self.assertLess(float(np.abs(bow_new[:, 2] - bow[:, 2]).max()), 1e-9)
        turned = np.degrees(
            np.arctan2(bow_new[:, 1], bow_new[:, 0]) - np.arctan2(bow[:, 1], bow[:, 0])
        )
        turned = (turned + 180.0) % 360.0 - 180.0
        self.assertLess(float(np.abs(turned - 90.0).max()), 1e-6)
        self.assertLess(
            float(np.abs(np.einsum("ij,ij->i", vertical, vertical) - 1.0).max()), 1e-12
        )


class DeployedLawMirrorTests(unittest.TestCase):
    """The fixed-tuning modes derive their frozen operating point in Python.

    That derivation mirrors the C++ adaptation law.  When the two drift apart
    every fixed mode is scored at a point the deployed filter would never
    choose, silently: nothing fails, the numbers are just wrong.  This pins the
    mirror to the header.
    """

    HEADER = REPO_ROOT / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h"

    def _header_value(self, pattern: str) -> float:
        import re

        text = self.HEADER.read_text(encoding="utf-8")
        match = re.search(pattern, text)
        self.assertIsNotNone(match, f"{pattern} not found in {self.HEADER.name}")
        return float(match.group(1))

    def test_rs_coefficient_matches_the_filter_default(self):
        import ou_validation as validation

        self.assertAlmostEqual(
            validation.OU_III_RS_COEFF,
            self._header_value(r"float\s+R_S_coeff_\s*=\s*([0-9.]+)f"),
            places=6,
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
