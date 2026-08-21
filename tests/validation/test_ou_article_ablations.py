#!/usr/bin/env python3

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DOC = ROOT / "doc" / "kalman_ou_iii"


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

    def test_article_no_longer_renders_legacy_one_way_transition_svg(self):
        text = (DOC / "w3d-baseline-comparison.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertNotIn(r"\includesvg[width=\columnwidth,inkscapelatex=false]{ou_validation_transition}", text)
        self.assertIn(r"Sec.~\ref{sec:roundtrip-transition-ablation}", text)

    def test_model_mismatch_table_matches_committed_summary(self):
        text = (DOC / "w3d-model-mismatch-ablation.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertIn(r"\texttt{--no-noise}", text)
        self.assertIn("OU--III & 0.246 & 0.174 & 0.140 & 0.332 & 11.87", text)
        self.assertIn("OU--II  & 0.333 & 0.224 & 0.240 & 0.467 & 20.32", text)
        self.assertIn("TFG     & 0.542 & 0.199 & 0.155 & 0.598 & 13.15", text)
        self.assertIn("not pure plant-model mismatch", text)


if __name__ == "__main__":
    unittest.main()
