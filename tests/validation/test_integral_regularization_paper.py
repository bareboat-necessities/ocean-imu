"""Publication contract for the standalone integral-state regularization paper."""

import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
PAPER = DOC / "adaptive-integral-state-regularization.tex"


def compact(text: str) -> str:
    return " ".join(text.split())


class IntegralRegularizationPaperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = PAPER.read_text(encoding="utf-8")
        cls.flat = compact(cls.text)

    def test_scope_is_current_law_design(self):
        self.assertIn(
            "Sea-State Scaling and Reduced Physical-MSE Design",
            self.flat,
        )
        self.assertIn("regularizer rather than a physical sensor", self.flat)
        self.assertIn("Historical tuning-law comparisons are intentionally outside", self.flat)

    def test_similarity_theorem_is_dimensional_context_only(self):
        self.assertIn(r"\label{thm:similarity}", self.text)
        self.assertIn(r"\boxed{\sigma_a\tau^3}", self.text)
        self.assertIn("coordinate-scale statement", self.flat)
        self.assertIn("not an adaptation law", self.flat)

    def test_reduced_mse_theorem_and_current_scaling_are_present(self):
        for token in (
            r"\label{thm:reduced-mse}",
            r"\frac{3q_{\rm eff}}{2\omega_R^3}",
            r"4m_{-4}\omega_R^4",
            r"\frac{9}{32}\frac{q_{\rm eff}}{m_{-4}}",
            r"q_{\rm eff}^{1/14}m_{-4}^{3/7}T_S^{-1/2}",
            r"q_{\rm eff}^{4/7}m_{-4}^{3/7}",
            r"\sigma_a^{6/7}\tau^{41/14}",
        ):
            self.assertIn(token, self.text)
        self.assertIn("optimum of $J_{\\rm red}$, not of the complete MEKF", self.flat)

    def test_exact_spectral_evaluation_sets_current_coefficient(self):
        self.assertIn(r"\label{tab:exact-mse}", self.text)
        self.assertIn(r"\label{eq:deployed-spectral-mse}", self.text)
        self.assertIn(r"C_J=0.0538", self.text)
        self.assertIn("reproduce $(6/7,41/14)$", self.flat)

    def test_analysis_limit_is_specific(self):
        self.assertIn("Where the Analytical Logic Stops", self.flat)
        self.assertIn("design covariance", self.flat)
        self.assertIn("no random $n_S$ is generated", self.flat)
        self.assertIn(
            "filter-internal Riccati covariance and the actual estimation-error covariance",
            self.flat,
        )
        self.assertIn("Full-estimator experiments remain the validation authority", self.flat)

    def test_legacy_numerical_studies_do_not_return(self):
        for token in (
            "legacy cubic comparator",
            "Paired ten-seed comparison",
            r"\label{tab:p-ablation}",
            r"\label{tab:channel-ablation}",
            "Amplitude-free law",
            "sensor-floor plus sea-dependent leakage",
            r"r_S\propto\sigma_a\tau^{5/2}",
        ):
            self.assertNotIn(token, self.text)

    def test_horizontal_diagnostic_and_mechanism_are_preserved(self):
        self.assertIn(r"\label{eq:q-horizontal-ratios}", self.text)
        self.assertIn("65.4", self.text)
        self.assertIn("102.3", self.text)
        self.assertIn(r"q_i^{1/14}", self.text)
        self.assertIn("six times the logarithmic sensitivity", self.flat)
        self.assertIn("reduce $q_{\\rm eff,H}$ itself", self.flat)

    def test_horizontal_claim_is_direction_limited_and_isotropic(self):
        self.assertIn(
            "Sparse directional simulations therefore cannot support a universal numerical $x/y$ optimum",
            self.flat,
        )
        self.assertIn("Isotropic integral regularization therefore performs adequately", self.flat)

    def test_two_column_artifacts_stay_single_column(self):
        self.assertNotIn(r"\begin{table*}", self.text)
        self.assertNotIn(r"\begin{figure*}", self.text)


if __name__ == "__main__":
    unittest.main()
