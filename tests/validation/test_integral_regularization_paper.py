"""Publication contract for the standalone integral-state regularization paper."""

import re
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
PAPER = DOC / "adaptive-integral-state-regularization.tex"


def compact(text: str) -> str:
    return " ".join(text.split())


def prose(text: str) -> str:
    flat = compact(text)
    return re.sub(r"\\(?:emph|textbf|textit)\{([^{}]*)\}", r"\1", flat)


class IntegralRegularizationPaperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = PAPER.read_text(encoding="utf-8")
        cls.flat = compact(cls.text)
        cls.prose = prose(cls.text)

    def test_scope(self):
        self.assertIn(
            "Sea-State Scaling, Reduced-Order Analysis, and Empirical Tuning",
            self.flat,
        )
        self.assertIn("regularizer rather than a physical sensor", self.prose)

    def test_similarity_theorem(self):
        self.assertIn(r"\label{thm:similarity}", self.text)
        self.assertIn(r"\boxed{\sigma_a\tau^3}", self.text)
        self.assertIn("coordinate-scale statement", self.prose)

    def test_reduced_mse_theorem_is_restored(self):
        for token in (
            r"\label{thm:reduced-mse}",
            r"\frac{3q_{\rm eff}}{2\omega_R^3}",
            r"4m_{-4}\omega_R^4",
            r"\frac{9}{32}\frac{q_{\rm eff}}{m_{-4}}",
            r"q_{\rm eff}^{1/14}m_{-4}^{3/7}T_S^{-1/2}",
            r"q_{\rm eff}^{4/7}m_{-4}^{3/7}",
        ):
            self.assertIn(token, self.text)
        self.assertIn("optimum of $J_{\\rm red}$, not of the complete MEKF", self.flat)

    def test_analysis_limit_is_specific(self):
        self.assertIn("Where the analytical logic stops", self.flat)
        self.assertIn("design covariance", self.prose)
        self.assertIn("no random $n_S$ is generated", self.flat)
        self.assertIn(
            "filter-internal Riccati covariance and the actual estimation-error covariance",
            self.flat,
        )
        self.assertIn("full-estimator experiments remain the tuning authority", self.flat)

    def test_deployed_law_remains_empirical(self):
        self.assertIn(r"r_S\propto\sigma_a\tau^{5/2}", self.text)
        self.assertIn(
            "empirical design choices validated in the complete estimator", self.flat
        )
        self.assertNotIn("deployed law is analytically optimal", self.prose.lower())

    def test_amplitude_ablation(self):
        self.assertIn(r"\label{tab:p-ablation}", self.text)
        self.assertIn(
            r"$\mathbf{p=1}$ & \textbf{4.747} & --- & \textbf{0.651}",
            self.text,
        )
        self.assertIn(r"p_{\rm red}=\frac{m+12}{14}", self.text)
        self.assertIn("used directly as empirical support", self.flat)

    def test_channel_ablation(self):
        self.assertIn(r"\label{tab:channel-ablation}", self.text)
        self.assertIn(r"Transition    & 4.70 & 4.91 & 8.66 & 10.31", self.text)
        self.assertIn("adaptive integral anchor is", self.flat)

    def test_awkward_sensitivity_table_stays_removed(self):
        self.assertNotIn(r"\label{tab:sensitivity}", self.text)

    def test_horizontal_diagnostic_and_mechanism_are_preserved(self):
        self.assertIn(r"\label{eq:q-horizontal-ratios}", self.text)
        self.assertIn("65.4", self.text)
        self.assertIn("102.3", self.text)
        self.assertIn(r"q_{\rm eff}^{1/14}", self.text)
        self.assertIn("eight times the sensitivity", self.flat)
        self.assertIn("reduce $q_{\\rm eff,H}$ itself", self.flat)

    def test_horizontal_claim_is_direction_limited_and_isotropic(self):
        self.assertIn(
            "sparse directional simulations cannot support a universal numerical $x/y$ optimum",
            self.flat,
        )
        self.assertIn("isotropic integral regularization therefore performs adequately", self.flat)
        self.assertIn("full MEKF experiments decide the practical choice", self.flat)

    def test_legacy_axis_optimum_does_not_return(self):
        for token in (
            r"\rho_{xy}",
            r"\label{eq:axis-opt}",
            r"\label{tab:axis-opt}",
            r"\label{tab:aniso-arms}",
            r"\label{tab:ssigma}",
            "1.138",
            "0.861",
        ):
            self.assertNotIn(token, self.text)

    def test_two_column(self):
        self.assertNotIn(r"\begin{table*}", self.text)
        self.assertNotIn(r"\begin{figure*}", self.text)


if __name__ == "__main__":
    unittest.main()
