#!/usr/bin/env python3

import csv
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


class OuArticleAblationContractTests(unittest.TestCase):
    def test_post_results_section_keeps_current_investigations_only(self):
        text = (DOC / "w3d-post-results-investigations.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertIn(r"\input{w3d-roundtrip-transition-ablation.tex-part}", text)
        self.assertIn(r"\input{w3d-model-mismatch-ablation.tex-part}", text)
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


if __name__ == "__main__":
    unittest.main()
