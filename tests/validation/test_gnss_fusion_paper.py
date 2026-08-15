"""Publication/source contract for the GNSS-aided marine inertial fusion paper."""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
PAPER = DOC / "kalman-gnss-fusion.tex"
OU2 = REPO_ROOT / "src" / "kalman_ou_ii" / "Kalman3D_Wave_OU_II.h"
OU3 = REPO_ROOT / "src" / "kalman_ou_iii" / "Kalman3D_Wave_OU_III.h"
TFG = REPO_ROOT / "src" / "kalman_tfg" / "Kalman3D_Wave_TFG.h"
NLO = REPO_ROOT / "src" / "nlo" / "TimeVaryingGainNLO.h"


def compact(text: str) -> str:
    return " ".join(text.split())


class GnssFusionPaperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = PAPER.read_text(encoding="utf-8")
        cls.flat = compact(cls.text)
        cls.ou2 = OU2.read_text(encoding="utf-8")
        cls.ou3 = OU3.read_text(encoding="utf-8")
        cls.tfg = TFG.read_text(encoding="utf-8")
        cls.nlo = NLO.read_text(encoding="utf-8")

    def test_standalone_scope_and_bibliography(self):
        self.assertIn(
            "GNSS-Aided Fusion in Marine Inertial Navigation Filters", self.flat
        )
        self.assertIn("No new GNSS performance result is claimed here", self.flat)
        self.assertIn(r"\bibliography{w3d,w3d-baselines}", self.text)
        self.assertNotIn(r"\begin{table*}", self.text)
        self.assertNotIn(r"\begin{figure*}", self.text)

    def test_trajectory_wave_split_is_explicit(self):
        for label in (
            r"\label{eq:sum-position}",
            r"\label{eq:sum-velocity}",
            r"\label{eq:sum-acceleration}",
        ):
            self.assertIn(label, self.text)
        self.assertIn(
            "absolute navigation and wave motion are separate latent components",
            self.flat.lower(),
        )
        self.assertIn(
            "must never be applied to $(\\vct r_n,\\vct u_n)$", self.flat
        )
        self.assertIn(
            "not the integral of absolute position", self.flat
        )

    def test_both_sensor_lever_arms_are_modeled(self):
        for token in (
            r"\vct \ell_I^{\mathcal B}",
            r"\vct \ell_G^{\mathcal B}",
            r"\label{eq:imu-rot-acc}",
            r"\label{eq:gnss-position}",
            r"\label{eq:gnss-velocity}",
        ):
            self.assertIn(token, self.text)
        self.assertIn(
            r"\vct\alpha^{\mathcal B}\times\vct\ell_I^{\mathcal B}", self.text
        )
        self.assertIn(
            r"\vct\omega^{\mathcal B}\times", self.text
        )
        self.assertIn("IMU-to-CoG offset modifies measured specific force", self.flat)
        self.assertIn("GNSS-antenna-to-CoG offset", self.flat)

    def test_ou2_ou3_source_really_has_the_documented_imu_lever_arm(self):
        required_api = (
            "set_imu_lever_arm_body",
            "clear_imu_lever_arm",
            "set_alpha_smoothing_tau",
        )
        for src in (self.ou2, self.ou3):
            for token in required_api:
                self.assertIn(token, src)
            self.assertIn("alpha_bprime.cross(r_imu_bprime)", src)
            self.assertIn(
                "omega_bprime.cross(omega_bprime.cross(r_imu_bprime))", src
            )
            self.assertIn("deheel_vector_(r_imu_wrt_cog_body_phys_)", src)

    def test_ou2_ou3_existing_position_velocity_updates_are_not_misrepresented(self):
        for src in (self.ou2, self.ou3):
            self.assertIn("measurement_update_position_pseudo", src)
            self.assertIn("measurement_update_velocity_pseudo", src)
        self.assertIn("a direct call with absolute GNSS position would still be wrong", self.flat)
        self.assertIn("those methods currently target the wave blocks", self.flat)

    def test_family_specific_sections_exist(self):
        for heading in (
            r"\section{Application to OU--II}",
            r"\section{Application to OU--III}",
            r"\section{Application to the Two-Frame Lie-Group EKF}",
            r"\section{Application to the Time-Varying-Gain NLO}",
        ):
            self.assertIn(heading, self.text)

    def test_tfg_claim_matches_current_source(self):
        self.assertIn("TwoFrameGroup<T, 4", self.tfg)
        self.assertNotIn("measurement_update_position_pseudo", self.tfg)
        self.assertIn("does not currently expose the OU-style generic", self.flat)
        self.assertIn("seven-column two-frame state", self.flat)
        self.assertIn("full retraction/reset covariance transport", self.flat)

    def test_nlo_claim_matches_current_source(self):
        self.assertIn("struct GnssXY", self.nlo)
        self.assertIn("Vec2 position_n_m", self.nlo)
        self.assertIn("R rms_xy_m", self.nlo)
        self.assertIn("usesVirtualHorizontal_", self.nlo)
        self.assertNotIn("velocity_n_mps", self.nlo)
        self.assertIn("does not yet provide GNSS velocity", self.flat)
        self.assertIn("sample age, freshness timeout", self.flat)

    def test_asynchronous_outage_and_reacquisition_are_first_class(self):
        for phrase in (
            "fuse a GNSS sample exactly once at its measurement time",
            "state-history ring buffer",
            "AIDED:",
            "COAST:",
            "REACQUIRE:",
            "velocity-first reacquisition",
            "preserving $\\vct p_w$",
        ):
            self.assertIn(phrase, self.flat)

    def test_required_studies_are_not_placeholder_prose(self):
        for heading in (
            r"\subsection{Reference-point and lever-arm studies}",
            r"\subsection{Absolute-navigation studies}",
            r"\subsection{Wave-plus-trajectory studies}",
            r"\subsection{GNSS-rate, latency, and asynchrony studies}",
            r"\subsection{Outage and reacquisition studies}",
            r"\subsection{GNSS quality and vertical-aiding studies}",
            r"\subsection{Architecture ablations}",
        ):
            self.assertIn(heading, self.text)
        for token in ("1, 5, 10, and 20~Hz", "1, 5, 30, 120, and 600~s", "250~ms"):
            self.assertIn(token, self.text)

    def test_results_section_remains_a_reporting_contract(self):
        self.assertIn(r"\section{Results Reporting Contract}", self.text)
        self.assertIn(
            "presently contains no GNSS-aided performance table", self.flat
        )
        self.assertIn("navigation-to-wave and wave-to-navigation leakage", self.flat)
        self.assertRegex(
            self.flat,
            re.compile(r"result that reports only total antenna or CoG position error"),
        )


if __name__ == "__main__":
    unittest.main()
