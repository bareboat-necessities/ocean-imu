"""Publication contract for the standalone integral-state regularization paper."""

import re
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
PAPER = DOC / "adaptive-integral-state-regularization.tex"


def compact(text: str) -> str:
    """Normalize source wrapping without changing LaTeX tokens."""
    return " ".join(text.split())


def prose(text: str) -> str:
    """Normalize wrapping and simple emphasis markup for prose assertions."""
    flat = compact(text)
    return re.sub(r"\\(?:emph|textbf|textit)\{([^{}]*)\}", r"\1", flat)


class IntegralRegularizationPaperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = PAPER.read_text(encoding="utf-8")
        cls.flat = compact(cls.text)
        cls.prose = prose(cls.text)

    def test_paper_has_standalone_publication_scope(self):
        self.assertIn(
            "Adaptive Integral-State Regularization for Drift-Bounded Inertial Motion Estimation",
            self.flat,
        )
        self.assertIn("Sea-State Scaling and Empirical Tuning", self.flat)
        self.assertIn("integral-state regularization", self.prose.lower())
        self.assertIn("regularizer rather than a physical sensor", self.prose)
        self.assertIn(r"\bibliography{w3d}", self.text)

    def test_similarity_theorem_is_scaling_not_optimality_claim(self):
        self.assertIn(r"\label{thm:similarity}", self.text)
        self.assertIn(r"\boxed{\sigma_a\tau^3}", self.text)
        self.assertIn("coordinate-scale statement", self.prose)
        self.assertIn("does not prove that this is the optimal Kalman tuning", self.prose)

    def test_reduced_tradeoff_is_explicitly_illustrative(self):
        self.assertIn(r"\label{prop:reduced-tradeoff}", self.text)
        self.assertIn(r"\frac{3q_{\rm eff}}{2\omega_R^3}", self.text)
        self.assertIn(r"4m_{-4}\omega_R^4", self.text)
        self.assertIn(r"q_{\rm eff}^{1/14}m_{-4}^{3/7}T_S^{-1/2}", self.text)
        self.assertIn("not asserted to be the optimum of the complete MEKF", self.prose)
        self.assertIn("Where the analytical logic stops", self.flat)

    def test_deployed_law_is_empirical_not_theorem_derived(self):
        self.assertIn(r"r_S\propto\sigma_a\tau^{5/2}", self.text)
        self.assertIn("The coefficient and the decision to retain this schedule", self.flat)
        self.assertIn("are empirical design choices", self.flat)
        self.assertIn("full-estimator evidence", self.flat)
        self.assertNotIn(r"\label{thm:mse-opt}", self.text)
        self.assertNotIn(r"r_S^\star", self.text)

    def test_amplitude_exponent_ablation_matches_committed_study(self):
        self.assertIn(r"\label{tab:p-ablation}", self.text)
        for row in (
            r"$p=0$    & 5.047 & $+0.300\pm0.071$ & 0.712",
            r"$p=0.5$  & 4.823 & $+0.076\pm0.025$ & 0.670",
            r"$\mathbf{p=1}$ & \textbf{4.747} & --- & \textbf{0.651}",
            r"$p=1.25$ & 4.779 & $+0.032\pm0.012$ & 0.653",
            r"$p=1.5$  & 4.817 & $+0.070\pm0.016$ & 0.660",
        ):
            self.assertIn(row, self.text)
        self.assertNotIn(r"p=\frac{m+12}{14}", self.text)
        self.assertIn("used directly as empirical support", self.flat)

    def test_channel_ablation_keeps_rs_as_the_dominant_adaptive_channel(self):
        self.assertIn(r"\label{tab:channel-ablation}", self.text)
        for row in (
            r"$H_s=0.27$ m & 5.13 & 5.13 & 15.82 & 15.72",
            r"$H_s=4.00$ m & 4.85 & 4.84 & 7.05 & 7.02",
            r"$H_s=8.50$ m & 4.03 & 4.02 & 11.28 & 11.14",
            r"Transition    & 4.66 & 4.65 & 8.82 & 8.90",
        ):
            self.assertIn(row, self.text)
        self.assertIn("adaptive integral anchor is", self.flat)

    def test_awkward_local_sensitivity_table_is_removed(self):
        self.assertNotIn(r"\label{tab:sensitivity}", self.text)
        self.assertIn("A separate local sensitivity sweep reaches the same mechanistic conclusion", self.flat)

    def test_horizontal_claim_is_isotropic_and_direction_limited(self):
        self.assertIn(r"\section{Three-Axis Regularization}", self.text)
        self.assertIn("limited set of wave headings", self.flat)
        self.assertIn("separate universal $x$ and $y$ regularizer coefficients", self.flat)
        self.assertIn("isotropic integral regularization performed adequately", self.flat)
        self.assertIn("same integral-regularization scale on all three axes", self.flat)

    def test_legacy_horizontal_tuning_and_axis_optimum_are_removed(self):
        for token in (
            r"\rho_{xy}",
            r"\label{eq:axis-opt}",
            r"\label{tab:axis-opt}",
            r"\label{tab:aniso-arms}",
            r"\label{tab:ssigma}",
            "1.138",
            "0.861",
            "historical horizontal factor",
        ):
            self.assertNotIn(token, self.text)

    def test_abstract_contains_no_reduced_optimum_numbers(self):
        abstract = self.text.split(r"\begin{abstract}", 1)[1].split(r"\end{abstract}", 1)[0]
        for token in (
            r"q_{\rm eff}^{1/14}",
            r"m_{-4}^{3/7}",
            "1.138",
            "0.861",
            "two orders of magnitude",
        ):
            self.assertNotIn(token, abstract)
        self.assertIn("justified empirically", abstract)

    def test_two_column_paper_does_not_use_full_width_tables(self):
        self.assertNotIn(r"\begin{table*}", self.text)
        self.assertNotIn(r"\begin{figure*}", self.text)


if __name__ == "__main__":
    unittest.main()
