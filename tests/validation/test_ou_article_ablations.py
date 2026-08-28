#!/usr/bin/env python3

import csv
import math
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "doc" / "kalman_ou_iii"
MISMATCH_RESULTS = ROOT / "reports" / "results" / "model_mismatch_ablation"
MISMATCH_SUMMARY = MISMATCH_RESULTS / "model_mismatch_summary.csv"
MISMATCH_FIGURES = (
    "ou_model_mismatch_floor.svg",
    "ou_model_mismatch_scaling.svg",
)
ENGINE_RESULTS = ROOT / "reports" / "results" / "engine_noise_degradation"
ENGINE_SUMMARY = ENGINE_RESULTS / "engine_noise_summary.csv"
ENGINE_FIGURES = (
    "ou_engine_noise_speed.svg",
    "ou_engine_noise_mechanism.svg",
)


class OuArticleAblationContractTests(unittest.TestCase):
    def test_post_results_section_keeps_current_investigations_only(self):
        text = (DOC / "w3d-post-results-investigations.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertIn(r"\input{w3d-roundtrip-transition-ablation.tex-part}", text)
        self.assertIn(r"\input{w3d-model-mismatch-ablation.tex-part}", text)
        self.assertIn(r"\input{w3d-engine-noise-degradation.tex-part}", text)
        self.assertIn(r"\input{w3d-full-state-h2.tex-part}", text)
        self.assertNotIn(
            r"\input{w3d-adaptation-coefficient-investigation.tex-part}", text
        )

    def test_roundtrip_protocol_is_low_high_same_low(self):
        text = (DOC / "w3d-roundtrip-transition-ablation.tex-part").read_text(
            encoding="utf-8"
        )
        for token in ("400", "520", "800", "920"):
            self.assertIn(token, text)
        self.assertIn("same phase-randomized", text)
        self.assertIn("low-start", text)
        self.assertIn("low-return", text)
        self.assertIn("rise", text)
        self.assertIn("fall", text)
        self.assertIn("No historical", text)
        self.assertIn("w3d-roundtrip-transition-scores-generated.tex-part", text)

    def test_only_bidirectional_transition_renders(self):
        baseline = (DOC / "w3d-baseline-comparison.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("ou_validation_transition", baseline)
        self.assertNotIn(r"fig:ou_transition", baseline)
        self.assertIn(r"Sec.~\ref{sec:roundtrip-transition-ablation}", baseline)

        transition = (DOC / "w3d-roundtrip-transition-ablation.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertIn(
            r"\includesvg[width=\columnwidth,inkscapelatex=false]{ou_rs_roundtrip_transition}",
            transition,
        )
        self.assertNotIn(r"Fig.~\ref{fig:ou_transition}", transition)
        self.assertFalse((DOC / "ou_validation_transition.svg").exists())

    def test_roundtrip_figure_is_mirrored_from_generated_evidence(self):
        generated = (
            ROOT / "reports" / "results" / "ou_rs_law"
            / "ou_rs_roundtrip_transition.svg"
        )
        mirrored = DOC / "ou_rs_roundtrip_transition.svg"
        self.assertTrue(generated.exists(), generated)
        self.assertEqual(generated.read_bytes(), mirrored.read_bytes())

    def test_roundtrip_generator_exports_full_segment_scores(self):
        tool = (ROOT / "tools" / "ou_roundtrip_transition.py").read_text(
            encoding="utf-8"
        )
        for token in (
            "ou_rs_roundtrip_scores.csv",
            "ou_rs_roundtrip_scores.tex",
            '"rise"',
            '"fall"',
            "disp_z_pct_refrms",
            "disp_3d_rms_m",
        ):
            self.assertIn(token, tool)

    def test_model_mismatch_table_matches_committed_summary(self):
        text = (DOC / "w3d-model-mismatch-ablation.tex-part").read_text(
            encoding="utf-8"
        )
        normalized = " ".join(text.split())
        self.assertIn(r"\texttt{--no-noise}", text)
        self.assertIn("not pure plant-model mismatch", text)

        with MISMATCH_SUMMARY.open(encoding="utf-8", newline="") as stream:
            rows = {row["family"]: row for row in csv.DictReader(stream)}

        tex_names = {"OU-II": "OU--II", "OU-III": "OU--III", "TFG": "TFG"}
        for family in ("OU-II", "OU-III", "TFG"):
            row = rows[family]
            expected = " ".join(
                (
                    tex_names[family], "&",
                    f"{float(row['disp_x_rms_m']):.3f}", "&",
                    f"{float(row['disp_y_rms_m']):.3f}", "&",
                    f"{float(row['disp_z_rms_m']):.3f}", "&",
                    f"{float(row['disp_3d_rms_m']):.3f}", "&",
                    f"{float(row['disp_z_pct_refrms']):.2f}",
                )
            )
            self.assertIn(expected, normalized)

    def test_model_mismatch_figures_are_mirrored_from_generated_evidence(self):
        for name in MISMATCH_FIGURES:
            with self.subTest(name=name):
                generated = MISMATCH_RESULTS / name
                mirrored = DOC / name
                self.assertTrue(generated.exists(), generated)
                self.assertEqual(generated.read_bytes(), mirrored.read_bytes(), name)

    def test_model_mismatch_section_renders_both_figures(self):
        text = (DOC / "w3d-model-mismatch-ablation.tex-part").read_text(
            encoding="utf-8"
        )
        for name in MISMATCH_FIGURES:
            stem = name.removesuffix(".svg")
            with self.subTest(name=name):
                self.assertIn(rf"\IfFileExists{{{name}}}", text)
                self.assertIn(
                    rf"\includesvg[width=\columnwidth,inkscapelatex=false]{{{stem}}}",
                    text,
                )
        self.assertIn(r"\label{fig:model-mismatch-floor}", text)
        self.assertIn(r"\label{fig:model-mismatch-scaling}", text)
        self.assertIn(r"Fig.~\ref{fig:model-mismatch-scaling}", text)

    def test_model_mismatch_generator_publishes_the_article_figures(self):
        tool = (ROOT / "tools" / "model_mismatch_ablation.py").read_text(
            encoding="utf-8"
        )
        for token in (
            "write_floor_plot",
            "write_scaling_plot",
            "ou_model_mismatch_floor.svg",
            "ou_model_mismatch_scaling.svg",
            'rcParams["svg.hashsalt"]',
            'metadata={"Date": None}',
        ):
            self.assertIn(token, tool)


class OuArticleEngineNoiseContractTests(unittest.TestCase):
    """The engine-noise section must agree with the committed study."""

    def summary(self) -> dict[tuple[str, str], dict[str, str]]:
        with ENGINE_SUMMARY.open(encoding="utf-8", newline="") as stream:
            return {
                (row["family"], row["label"]): row for row in csv.DictReader(stream)
            }

    def test_cruise_table_matches_committed_summary(self):
        text = (DOC / "w3d-engine-noise-degradation.tex-part").read_text(
            encoding="utf-8"
        )
        normalized = " ".join(text.split())
        rows = self.summary()

        tex_names = {"OU-II": "OU--II", "OU-III": "OU--III", "TFG": "TFG"}
        for family in ("OU-II", "OU-III", "TFG"):
            for label, engine in (("engine off", "off"), ("speed 2400 rpm", "on")):
                row = rows[(family, label)]
                offset = math.sqrt(
                    sum(
                        float(row[name]) ** 2
                        for name in ("disp_x_mean_m", "disp_y_mean_m", "disp_z_mean_m")
                    )
                )
                expected = " ".join(
                    (
                        tex_names[family], "&", engine, "&",
                        f"{float(row['disp_3d_rms_m']):.3f}", "&",
                        f"{offset:.3f}", "&",
                        f"{float(row['pitch_rms_deg']):.3f}", "&",
                        f"{float(row['yaw_rms_deg']):.2f}",
                    )
                )
                with self.subTest(family=family, engine=engine):
                    self.assertIn(expected, normalized)

    def test_matched_power_control_is_the_stated_comparison(self):
        """The central claim is power, not fold placement; check both spreads."""

        rows = self.summary()
        for family in ("OU-II", "OU-III"):
            as_is = [
                float(rows[(family, label)]["disp_3d_rms_m"])
                for label in (
                    "bandwidth 20 Hz", "bandwidth 40 Hz",
                    "speed 2400 rpm", "bandwidth 160 Hz",
                )
            ]
            matched = [
                float(rows[(family, label)]["disp_3d_rms_m"])
                for label in (
                    "matched 20 Hz", "matched 40 Hz",
                    "speed 2400 rpm", "matched 160 Hz",
                )
            ]
            with self.subTest(family=family):
                # Equalizing recorded power must collapse the bandwidth spread.
                self.assertGreater(max(as_is) / min(as_is), 5.0)
                self.assertLess(max(matched) / min(matched), 3.0)

        # And every matched cell really is at one recorded level.
        recorded = {
            round(float(rows[(family, label)]["recorded_rms_mps2"]), 4)
            for family in ("OU-II", "OU-III", "TFG")
            for label in (
                "matched 20 Hz", "matched 40 Hz",
                "speed 2400 rpm", "matched 160 Hz",
            )
        }
        self.assertEqual(len(recorded), 1, recorded)

    def test_gyro_path_ablation_leaves_the_result_unchanged(self):
        """The section attributes the whole effect to the accelerometer."""

        rows = self.summary()
        for family in ("OU-II", "OU-III", "TFG"):
            full = float(rows[(family, "speed 2400 rpm")]["disp_3d_rms_m"])
            accel_only = float(rows[(family, "accelerometer only")]["disp_3d_rms_m"])
            with self.subTest(family=family):
                self.assertAlmostEqual(accel_only / full, 1.0, delta=0.01)

    def test_error_is_dominated_by_the_static_offset(self):
        rows = self.summary()
        for family in ("OU-II", "OU-III", "TFG"):
            row = rows[(family, "speed 2400 rpm")]
            offset = math.sqrt(
                sum(
                    float(row[name]) ** 2
                    for name in ("disp_x_mean_m", "disp_y_mean_m", "disp_z_mean_m")
                )
            )
            with self.subTest(family=family):
                self.assertGreater(offset / float(row["disp_3d_rms_m"]), 0.8)

    def test_figures_are_mirrored_from_generated_evidence(self):
        for name in ENGINE_FIGURES:
            with self.subTest(name=name):
                generated = ENGINE_RESULTS / name
                mirrored = DOC / name
                self.assertTrue(generated.exists(), generated)
                self.assertEqual(generated.read_bytes(), mirrored.read_bytes(), name)

    def test_section_renders_both_figures(self):
        text = (DOC / "w3d-engine-noise-degradation.tex-part").read_text(
            encoding="utf-8"
        )
        for name in ENGINE_FIGURES:
            stem = name.removesuffix(".svg")
            with self.subTest(name=name):
                self.assertIn(rf"\IfFileExists{{{name}}}", text)
                self.assertIn(
                    rf"\includesvg[width=\columnwidth,inkscapelatex=false]{{{stem}}}",
                    text,
                )
        self.assertIn(r"\label{fig:engine-noise-speed}", text)
        self.assertIn(r"\label{fig:engine-noise-mechanism}", text)
        self.assertIn("sensor-path", text)

    def test_generator_publishes_the_article_figures(self):
        tool = (ROOT / "tools" / "engine_noise_degradation.py").read_text(
            encoding="utf-8"
        )
        for token in (
            "write_speed_plot",
            "write_mechanism_plot",
            "ou_engine_noise_speed.svg",
            "ou_engine_noise_mechanism.svg",
            'rcParams["svg.hashsalt"]',
            'metadata={"Date": None}',
            "VALIDATION_METRICS_MEANS",
        ):
            self.assertIn(token, tool)

    def test_engine_model_is_off_by_default(self):
        """Every historical realization depends on this staying opt-in."""

        source = (ROOT / "src" / "util" / "W3dSimCommon.cpp").read_text(
            encoding="utf-8"
        )
        self.assertIn('w3d_engine_float_from_env("W3D_ENGINE_RPM", cfg.rpm)', source)
        self.assertIn("return std::nullopt;", source)

        header = (ROOT / "src" / "util" / "W3dSimCommon.h").read_text(encoding="utf-8")
        self.assertIn("float rpm = 0.0f;", header)

    def test_means_stay_off_the_fingerprinted_metrics_line(self):
        """Adding these columns to VALIDATION_METRICS would rewrite evidence."""

        source = (ROOT / "src" / "util" / "W3dSimCommon.cpp").read_text(
            encoding="utf-8"
        )
        metrics_line = source.split('<< record_tag\n')[1].split('<< "\\n";')[0]
        for name in ("pitch_mean_deg", "disp_z_mean_m", "yaw_mean_deg"):
            self.assertNotIn(name, metrics_line)

        parser = (ROOT / "tools" / "ou_validation.py").read_text(encoding="utf-8")
        for name in ("pitch_mean_deg", "disp_z_mean_m", "yaw_mean_deg"):
            self.assertNotIn(name, parser)


if __name__ == "__main__":
    unittest.main()
