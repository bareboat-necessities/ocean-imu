from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC_RESULTS = REPO_ROOT / "doc" / "kalman_ou_iii" / "w3d-results.tex-part"
PLOT_SCRIPT = REPO_ROOT / "plots" / "kalman_ou_iii" / "kalman_ou_iii-plots.py"
DRAW_SCRIPT = REPO_ROOT / "plots" / "kalman_ou_iii" / "draw_plots.sh"


class OUPublicationLayoutContractTests(unittest.TestCase):
    def test_attitude_tuning_and_bias_figures_are_vertical_single_column(self):
        text = DOC_RESULTS.read_text(encoding="utf-8")

        attitude = text.split(r"\subsection{Attitude and adaptation behavior}", 1)[1]
        attitude = attitude.split(r"\subsection{Bias tracking}", 1)[0]
        self.assertIn(r"\begin{figure}[t]", attitude)
        self.assertNotIn(r"\begin{figure*}", attitude)
        self.assertEqual(attitude.count(r"\begin{minipage}[t]{\columnwidth}"), 2)
        self.assertIn(r"\vspace{0.5em}", attitude)

        bias = text.split(r"\subsection{Bias tracking}", 1)[1]
        bias = bias.split(r"\subsection{Comparative results and adaptation ablations}", 1)[0]
        self.assertIn(r"\begin{figure}[t]", bias)
        self.assertNotIn(r"\begin{figure*}", bias)
        self.assertEqual(bias.count(r"\begin{minipage}[t]{\columnwidth}"), 2)
        self.assertIn(r"\vspace{0.5em}", bias)

    def test_plot_dataset_headers_are_opt_in(self):
        source = PLOT_SCRIPT.read_text(encoding="utf-8")
        self.assertIn('"--show-titles"', source)
        self.assertIn('action="store_true"', source)
        self.assertIn("SHOW_TITLES = args.show_titles", source)
        self.assertIn("if SHOW_TITLES:\n        fig.suptitle(title)", source)
        self.assertNotIn("default=True", source)

    def test_axes_table_publication_view_is_compacted_without_touching_evidence(self):
        source = DRAW_SCRIPT.read_text(encoding="utf-8")
        self.assertIn("def compact_axes_table(text: str) -> str:", source)
        self.assertIn(r'block.replace(spacing, r"\setlength{\tabcolsep}{1.6pt}"', source)
        self.assertIn(r're.subn(r"\$H_s=([0-9.]+)\$ m", r"\1 m", block)', source)
        self.assertIn(r'r"$\1{\\pm}\2$"', source)
        self.assertIn("text = compact_axes_table(text)", source)
        self.assertIn("w3d-ou-validation-results-publication.tex-part", source)

    def test_publication_view_uses_compact_validation_table_captions(self):
        source = DRAW_SCRIPT.read_text(encoding="utf-8")
        self.assertIn("def replace_table_caption(text: str, label: str, caption: str) -> str:", source)
        for caption in (
            r"Ten-seed paired OU-family comparison over the final \SI{900}{s}.",
            r"Ten-seed adaptive displacement RMS by axis over the final \SI{900}{s}.",
            r"OU--III adaptation-channel ablation over the final \SI{900}{s}.",
            r"Controlled transition scored by interval (Adaptive mode).",
            r"Ten-seed paired OU-family comparison on PM--Stokes seas.",
        ):
            self.assertIn(caption, source)


if __name__ == "__main__":
    unittest.main()
