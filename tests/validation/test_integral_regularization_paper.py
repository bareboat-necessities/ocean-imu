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


class IntegralRegularizationPaperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = PAPER.read_text(encoding="utf-8")
        cls.flat = compact(cls.text)

    def test_paper_has_standalone_publication_scope(self):
        self.assertIn(
            "Adaptive Integral-State Regularization for Drift-Bounded Inertial Motion Estimation",
            self.flat,
        )
        self.assertIn("integral-state regularization", self.flat.lower())
        self.assertIn("regularizer, not as a physical sensor", self.flat)
        self.assertIn(r"\bibliography{w3d}", self.text)

    def test_bias_variance_optimum_is_explicit(self):
        for token in (
            r"\label{thm:mse-opt}",
            r"\frac{3q_{\rm eff}}{2\omega_R^3}",
            r"4m_{-4}\omega_R^4",
            r"\frac{9}{32}\frac{q_{\rm eff}}{m_{-4}}",
            r"q_{\rm eff}^{1/14}m_{-4}^{3/7}T_S^{-1/2}",
            r"q_{\rm eff}^{4/7}m_{-4}^{3/7}",
        ):
            self.assertIn(token, self.text)

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
        self.assertIn(r"p=\frac{m+12}{14}", self.text)

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

    def test_horizontal_vertical_optimum_is_measured_not_assumed(self):
        self.assertIn(r"\label{eq:axis-opt}", self.text)
        self.assertIn(
            r"\left(\frac{q_{{\rm eff},i}}{q_{{\rm eff},j}}\right)^{1/14}",
            self.text,
        )
        self.assertIn(
            r"\left(\frac{m_{-4,i}}{m_{-4,j}}\right)^{3/7}", self.text
        )
        self.assertRegex(
            self.flat,
            re.compile(
                r"Geometric mean\s*&\s*65\.4\s*&\s*102\.3\s*&\s*0\.674\s*&\s*0\.326\s*&\s*\\textbf\{1\.138\}\s*&\s*\\textbf\{0\.861\}"
            ),
        )
        self.assertIn("two orders of magnitude", self.flat)
        self.assertIn("do not attempt to compensate", self.flat)

    def test_anisotropy_ablation_preserves_counterintuitive_result(self):
        self.assertIn(r"\label{tab:aniso-arms}", self.text)
        self.assertIn(r"$(1.87,1.87)$ & +2.3 & +27.6 & -0.2 & +13.7", self.text)
        self.assertIn(r"$(1.00,1.00)$ & -4.2 & -1.1 & +0.1 & -2.6", self.text)
        self.assertIn(r"\label{tab:ssigma}", self.text)
        self.assertIn(
            r"\textbf{1.00} & \textbf{-2.11} & \textbf{-0.54} & +0.16 & \textbf{-1.28}",
            self.text,
        )

    def test_historical_validation_scope_is_not_hidden(self):
        self.assertIn("generated before the final horizontal OU-prior factor", self.flat)
        self.assertIn("changed from 1.87 to the currently deployed isotropic value 1.0", self.flat)
        self.assertIn("A future full regeneration should collapse that historical distinction", self.flat)

    def test_two_column_paper_does_not_use_full_width_tables(self):
        self.assertNotIn(r"\begin{table*}", self.text)
        self.assertNotIn(r"\begin{figure*}", self.text)
        # Wide numeric tables, where needed, are explicitly fit to one IEEE column.
        self.assertGreaterEqual(self.text.count(r"\resizebox{\columnwidth}{!}"), 2)


if __name__ == "__main__":
    unittest.main()
