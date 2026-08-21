#!/usr/bin/env python3

import csv
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "doc" / "kalman_ou_iii"
MISMATCH_SUMMARY = (
    ROOT / "reports" / "results" / "model_mismatch_ablation"
    / "model_mismatch_summary.csv"
)


class OuArticleAblationContractTests(unittest.TestCase):
    def test_post_results_section_includes_roundtrip_and_model_mismatch(self):
        text = (DOC / "w3d-post-results-investigations.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertIn(r"\input{w3d-roundtrip-transition-ablation.tex-part}", text)
        self.assertIn(r"\input{w3d-model-mismatch-ablation.tex-part}", text)

    def test_roundtrip_protocol_is_low_high_same_low(self):
        text = (DOC / "w3d-roundtrip-transition-ablation.tex-part").read_text(
            encoding="utf-8"
        )
        for token in ("400", "520", "800", "920"):
            self.assertIn(token, text)
        self.assertIn("same phase-randomized", text)
        self.assertIn("low-start", text)
        self.assertIn("low-return", text)
        self.assertIn("high-to-low", text)

    def test_each_transition_protocol_renders_its_own_diagnostic(self):
        # Both protocols are shown the same way -- sea state, wave shapes,
        # vertical error and the applied tuning traces -- so the round-trip
        # figure cannot regress to an abstract plot of the blend weight while
        # the one-way figure keeps its time-domain evidence.
        baseline = (DOC / "w3d-baseline-comparison.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertIn(
            r"\includesvg[width=\columnwidth,inkscapelatex=false]{ou_validation_transition}",
            baseline,
        )
        self.assertIn(r"Sec.~\ref{sec:roundtrip-transition-ablation}", baseline)

        ablation = (DOC / "w3d-roundtrip-transition-ablation.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertIn(
            r"\includesvg[width=\columnwidth,inkscapelatex=false]{ou_rs_roundtrip_transition}",
            ablation,
        )
        self.assertNotIn("high-sea blend weight", ablation)
        self.assertIn(r"Fig.~\ref{fig:ou_transition}", ablation)

    def test_roundtrip_figure_is_mirrored_from_the_generated_evidence(self):
        generated = (
            ROOT / "reports" / "results" / "ou_rs_law"
            / "ou_rs_roundtrip_transition.svg"
        )
        mirrored = DOC / "ou_rs_roundtrip_transition.svg"
        self.assertTrue(generated.exists(), generated)
        self.assertEqual(generated.read_bytes(), mirrored.read_bytes())

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


if __name__ == "__main__":
    unittest.main()
