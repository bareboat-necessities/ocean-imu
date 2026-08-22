"""Publication contract for retiring the OU-III legacy regularizer-law studies."""

import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"


class OUIIILegacyLawCleanupTests(unittest.TestCase):
    def test_main_sources_do_not_restore_legacy_numerical_comparisons(self):
        sources = {
            name: (DOC / name).read_text(encoding="utf-8")
            for name in (
                "w3d-reduced-mse-envelope.tex-part",
                "w3d-reduced-mse-validation.tex-part",
                "w3d-adaptation-deployed-law.tex-part",
                "w3d-conclusion-summary.tex-part",
                "w3d-post-results-investigations.tex-part",
            )
        }
        joined = "\n".join(sources.values())
        for retired in (
            r"\label{tab:adapt-reduced-mse-envelope}",
            r"\label{eq:adapt-reduced-vs-applied-ratio}",
            r"\label{eq:adapt-reduced-direct-comparison}",
            r"\label{tab:adapt-mse-paired}",
            r"\label{eq:deployed-effective-law}",
            "Legacy dimensional cubic schedule",
            "Coefficient-Law Investigation and Transition Tradeoff",
            "Sensor-floor plus sea-dependent leakage",
            "supported low-cost cubic",
        ):
            self.assertNotIn(retired, joined)

    def test_publication_tree_has_no_legacy_numerical_law_studies(self):
        validation = (DOC / "w3d-ou-validation-results-generated.tex-part").read_text(
            encoding="utf-8"
        )
        robustness = (DOC / "w3d-ou-robustness.tex-part").read_text(
            encoding="utf-8"
        )
        robustness_results = (
            DOC / "w3d-ou-robustness-results-generated.tex-part"
        ).read_text(encoding="utf-8")
        robustness_macros = (
            DOC / "w3d-ou-robustness-macros-generated.tex-part"
        ).read_text(encoding="utf-8")
        publication = "\n".join(
            (validation, robustness, robustness_results, robustness_macros)
        )

        for retired in (
            r"0.35\,\sigma_{aw}\tau",
            r"\tau,r_S\propto\tau^3",
            "legacy cubic",
            "historical coupled",
            "OU--III tuning sensitivity",
            r"\label{tab:ou_robustness_sensitivity}",
            r"\OURobustnessWorstParameter",
            r"\OURobustnessTauSpan",
            r"\OURobustnessSigmaSpan",
            r"\OURobustnessRSSpan",
            r"\OURobustnessCoupledSigmaSpan",
            r"\OURobustnessCoupledTauSpan",
        ):
            self.assertNotIn(retired, publication)

        self.assertIn("SpectralMSE", validation)
        self.assertIn(r"\label{tab:ou_mc_channels}", validation)
        self.assertIn(r"\label{tab:ou_robustness_stress}", robustness_results)
        self.assertFalse((DOC / "ou_robustness_sensitivity.svg").exists())

    def test_current_spectral_mse_derivation_remains_the_publication_law(self):
        envelope = (DOC / "w3d-reduced-mse-envelope.tex-part").read_text(
            encoding="utf-8"
        )
        validation = (DOC / "w3d-reduced-mse-validation.tex-part").read_text(
            encoding="utf-8"
        )
        deployed = (DOC / "w3d-adaptation-deployed-law.tex-part").read_text(
            encoding="utf-8"
        )
        for token in (
            r"q_{\rm eff}^{1/14}m_{-4}^{3/7}T_S^{-1/2}",
            r"\sigma_a^{6/7}\tau^{41/14}",
        ):
            self.assertIn(token, envelope)
        self.assertIn(r"\boxed{C_J\simeq0.0538.}", validation)
        self.assertIn(r"\label{eq:adapt-default-rs-law}", deployed)

    def test_sigma_tau_cubed_survives_only_as_dimensional_context(self):
        deployed = (DOC / "w3d-adaptation-deployed-law.tex-part").read_text(
            encoding="utf-8"
        )
        self.assertIn(r"r_{S,\mathrm{dim}}", deployed)
        self.assertIn(r"\sigma_{aw}\tau^3", deployed)
        self.assertIn("coordinate-scale statement only", deployed)
        self.assertIn("not an adaptation law", deployed)

    def test_retired_publication_fragments_and_study_drivers_are_absent(self):
        for relative in (
            "doc/kalman_ou_iii/w3d-adaptation-coefficient-investigation.tex-part",
            "doc/kalman_ou_iii/w3d-rs-design-interpretation.tex-part",
            "docs/ou-iii-rs-amplitude-retune.md",
            "tools/ou_rs_amplitude_retune_compare.py",
            "tools/ou_rs_amplitude_retune_sweep.py",
            "tools/ou_rs_law_ablation.py",
            "tools/ou_rs_spectral_mse_sweep.py",
            "tools/rs_law_bias_variance.py",
        ):
            self.assertFalse((REPO_ROOT / relative).exists(), relative)

    def test_roundtrip_instrument_has_no_law_comparison_mode(self):
        tool = (REPO_ROOT / "tools" / "ou_roundtrip_transition.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("does not sweep or compare regularizer laws", tool)
        for token in ("--exponents", "OU_III_RS_SIGMA_EXP", "LAW_CUBIC"):
            self.assertNotIn(token, tool)


if __name__ == "__main__":
    unittest.main()
