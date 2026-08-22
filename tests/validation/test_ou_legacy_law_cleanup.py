"""Publication and study contract for retiring legacy OU adaptation laws."""

import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"

OU2_FORBIDDEN = ("PseudoAdaptationLaw::Empirical",)
OU3_FORBIDDEN = (
    "RSAdaptationLaw::Cubic",
    "RSAdaptationLaw::StrongRiccati",
    "RSAdaptationLaw::PosteriorRiccati",
)


class OULegacyLawCleanupTests(unittest.TestCase):
    def test_ou2_tests_use_only_physical_mse(self):
        for path in (REPO_ROOT / "tests" / "kalman_ou_ii").glob("*-test.cpp"):
            text = path.read_text(encoding="utf-8")
            for token in OU2_FORBIDDEN:
                self.assertNotIn(token, text, path.name)
        self.assertFalse((REPO_ROOT / "tools" / "ou2_pseudo_law_compare.py").exists())
        self.assertFalse((REPO_ROOT / "tools" / "ou2_pseudo_mse_scale_sweep.py").exists())

    def test_ou3_tests_use_only_spectral_mse(self):
        for path in (REPO_ROOT / "tests" / "kalman_ou_iii").glob("*-test.cpp"):
            text = path.read_text(encoding="utf-8")
            for token in OU3_FORBIDDEN:
                self.assertNotIn(token, text, path.name)

    def test_arduino_sketches_never_select_retired_laws(self):
        for path in REPO_ROOT.rglob("*.ino"):
            text = path.read_text(encoding="utf-8", errors="ignore")
            for token in OU2_FORBIDDEN + OU3_FORBIDDEN:
                self.assertNotIn(token, text, str(path.relative_to(REPO_ROOT)))
            self.assertNotIn("OU_II_PSEUDO_LAW", text, str(path.relative_to(REPO_ROOT)))
            self.assertNotIn("OU_III_RS_LAW", text, str(path.relative_to(REPO_ROOT)))

    def test_full_validation_has_no_legacy_law_selector(self):
        texts = [
            (REPO_ROOT / ".github" / "workflows" / "ou-validation.yml").read_text(encoding="utf-8"),
            (REPO_ROOT / "tools" / "ou_validation.py").read_text(encoding="utf-8"),
            (REPO_ROOT / "tools" / "ou_robustness.py").read_text(encoding="utf-8"),
        ]
        joined = "\n".join(texts)
        for token in OU2_FORBIDDEN + OU3_FORBIDDEN:
            self.assertNotIn(token, joined)
        # Validation may exercise channel freezes, but it must never request a
        # different law through either simulator's compatibility selector.
        self.assertNotIn("OU_II_PSEUDO_LAW", joined)
        self.assertNotIn("OU_III_RS_LAW", joined)
        self.assertNotIn("OU_III_RS_SIGMA_EXP", joined)

    def test_publication_tree_has_no_legacy_numerical_law_studies(self):
        publication_files = [
            DOC / "w3d-ou-validation-results-publication.tex-part",
            DOC / "w3d-ou-robustness.tex-part",
            DOC / "w3d-ou-robustness-results-generated.tex-part",
            DOC / "w3d-ou-robustness-macros-generated.tex-part",
            DOC / "w3d-roundtrip-transition-ablation.tex-part",
        ]
        publication = "\n".join(
            p.read_text(encoding="utf-8") for p in publication_files if p.exists()
        )
        for retired in (
            r"0.35\,\sigma_{aw}\tau",
            "legacy cubic",
            "historical coupled",
            "OU--III tuning sensitivity",
            r"\label{tab:ou_robustness_sensitivity}",
            r"\OURobustnessWorstParameter",
        ):
            self.assertNotIn(retired, publication)
        self.assertFalse((DOC / "ou_robustness_sensitivity.svg").exists())
        self.assertFalse((DOC / "ou_validation_transition.svg").exists())

    def test_only_roundtrip_transition_is_publishable(self):
        sync = (REPO_ROOT / "tools" / "ou_publication_sync.py").read_text(encoding="utf-8")
        baseline = (DOC / "w3d-baseline-comparison.tex-part").read_text(encoding="utf-8")
        roundtrip = (DOC / "w3d-roundtrip-transition-ablation.tex-part").read_text(encoding="utf-8")
        self.assertIn("curate_validation_for_article", sync)
        self.assertIn("tab:ou_transition_segments", sync)
        self.assertNotIn("fig:ou_transition", baseline)
        self.assertNotIn("ou_validation_transition.svg", baseline)
        self.assertIn("ou_rs_roundtrip_transition.svg", roundtrip)
        self.assertIn("low--high--low", roundtrip)
        self.assertIn("rise", roundtrip)
        self.assertIn("fall", roundtrip)

    def test_current_analytical_derivations_remain_publication_laws(self):
        ou3 = (DOC / "w3d-adaptation-deployed-law.tex-part").read_text(encoding="utf-8")
        self.assertIn(r"\label{eq:adapt-default-rs-law}", ou3)
        self.assertIn("not an adaptation law", ou3)
        ou2 = (REPO_ROOT / "src" / "kalman_ou_ii" / "SeaStateFusionFilter_OU_II.h").read_text(encoding="utf-8")
        ou3src = (REPO_ROOT / "src" / "kalman_ou_iii" / "SeaStateFusionFilter_OU_III.h").read_text(encoding="utf-8")
        self.assertIn("PseudoAdaptationLaw::PhysicalMSE", ou2)
        self.assertIn("RSAdaptationLaw::SpectralMSE", ou3src)

    def test_roundtrip_instrument_has_no_law_comparison_mode(self):
        tool = (REPO_ROOT / "tools" / "ou_roundtrip_transition.py").read_text(encoding="utf-8")
        self.assertIn("does not sweep or compare regularizer laws", tool)
        for token in ("--exponents", "OU_III_RS_SIGMA_EXP", "LAW_CUBIC"):
            self.assertNotIn(token, tool)


if __name__ == "__main__":
    unittest.main()
