"""Publication contract for the wind-heel / wave-roll concept paper."""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
PAPER = REPO_ROOT / "doc" / "kalman_ou_iii" / "kalman-wind-heel.tex"


def compact(text: str) -> str:
    return " ".join(text.split())


class WindHeelPaperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = PAPER.read_text(encoding="utf-8")
        cls.flat = compact(cls.text)

    def test_title_and_prospective_scope(self):
        self.assertIn(
            "Separation of Steady Wind Heel and Wave-Induced Roll Using Marine IMU Sensor Fusion",
            self.flat,
        )
        self.assertIn("concept paper", self.flat.lower())
        self.assertIn("No performance result is claimed in this draft", self.flat)
        self.assertIn("No experimental or simulation results are claimed in this draft", self.flat)
        self.assertIn(r"\bibliography{w3d}", self.text)

    def test_requested_eighteen_section_structure_is_present(self):
        sections = re.findall(r"\\section\{([^}]*)\}", self.text)
        expected = [
            "Introduction",
            "Sailboat Heel and Wave-Roll Dynamics",
            "Coordinate Frames and Measurements",
            "Steady-Heel State/Model",
            "Dynamic Wave-Attitude Model",
            "Apparent-Wind Aiding",
            "Heel-Frame Transformation",
            "Observability and Ambiguity with Horizontal Acceleration",
            "Adaptive Separation of Mean Heel and Wave Roll",
            "Maneuver/Tack Detection",
            "Validation Methodology",
            "Constant-Wind Cases",
            "Gusts and Changing Wind",
            "Tacks and Gybes",
            "Irregular Directional Seas",
            "Results",
            "Discussion",
            "Conclusion",
        ]
        self.assertEqual(expected, sections)

    def test_frame_transform_preserves_dynamic_roll_and_physical_attitude(self):
        for token in (
            r"\vct v_{B'}=\mat R_x(-\phi_h)\vct v_B",
            r"\mat R_{WB}=\mat R_{WB'}\mat R_x(-\phi_h)",
            r"\mat R_{WB'}^{+}=\mat R_{WB'}^{-}\mat R_x(\Delta\phi_h)",
        ):
            self.assertIn(token, self.text)
        self.assertIn("The dynamic wave roll is retained inside", self.flat)

    def test_observability_ambiguity_is_explicit(self):
        self.assertIn(r"\label{eq:acc-linearized}", self.text)
        self.assertIn(r"\label{eq:gramian}", self.text)
        self.assertIn("horizontal acceleration", self.flat.lower())
        self.assertIn("accelerometer bias", self.flat.lower())
        self.assertIn("magnetometer and a gyroscope improve physical-attitude estimation", self.flat)
        self.assertIn("they do not fully solve the decomposition problem", self.flat)

    def test_wind_is_a_soft_aid_not_a_direct_heel_sensor(self):
        self.assertIn(r"\label{eq:wind-pseudo}", self.text)
        self.assertIn("wind as a broad, adaptive prior", self.flat)
        self.assertIn("not as a hard pseudo-measurement", self.flat)
        self.assertIn("masthead lever arm", self.flat.lower())
        self.assertIn("sensor latency", self.flat.lower())

    def test_validation_program_contains_required_controls_and_metrics(self):
        for phrase in (
            "no heel separation",
            "fixed-cutoff low-pass subtraction of roll",
            "slow heel state without wind aiding",
            "wind-aided adaptive heel state",
            "oracle heel transform",
            "deliberately over-fast heel state",
            "paired seeds",
            "wave-band leakage ratio",
            "observability diagnostics",
        ):
            self.assertIn(phrase, self.flat)

    def test_all_requested_study_classes_have_negative_or_stress_controls(self):
        for label in (
            r"\label{sec:constant-wind}",
            r"\label{sec:gusts}",
            r"\label{sec:tack-gybe-study}",
            r"\label{sec:directional-seas}",
        ):
            self.assertIn(label, self.text)
        self.assertIn("constant transverse acceleration with no wind change", self.flat)
        self.assertIn("large wave-roll excursion with no maneuver", self.flat)
        self.assertIn("Long-period swell is a required edge case", self.flat)

    def test_results_section_is_a_reporting_contract_not_fake_evidence(self):
        results = self.text.split(r"\section{Results}", 1)[1].split(r"\section{Discussion}", 1)[0]
        flat = compact(results)
        self.assertIn("No experimental or simulation results are claimed in this draft", flat)
        self.assertIn("reporting contract", flat)
        self.assertIn("primary table for constant-wind", flat)
        self.assertIn("observability condition-number map", flat)

    def test_two_column_layout_avoids_full_width_floats(self):
        self.assertNotIn(r"\begin{table*}", self.text)
        self.assertNotIn(r"\begin{figure*}", self.text)


if __name__ == "__main__":
    unittest.main()
