"""Publication contract for generated figures in kalman-wave-dir.tex."""

import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DOC = REPO_ROOT / "doc" / "kalman_ou_iii"
PLOTS = REPO_ROOT / "plots" / "kalman_ou_iii"


class WaveDirectionChartContractTests(unittest.TestCase):
    def test_direction_article_uses_generated_single_column_pgfs(self):
        article = (DOC / "kalman-wave-dir.tex").read_text(encoding="utf-8")
        results = (DOC / "w3d-wave-direction-results.tex-part").read_text(
            encoding="utf-8"
        )
        frequency = (DOC / "w3d-frequency-tracking.tex-part").read_text(
            encoding="utf-8"
        )
        charts = (DOC / "w3d-wave-direction-charts.tex-part").read_text(
            encoding="utf-8"
        )
        plotter = (PLOTS / "kalman_ou_iii-plots.py").read_text(encoding="utf-8")

        self.assertIn(r"\usepackage{graphicx}", article)
        self.assertIn(r"\usepackage{pgf}", article)
        self.assertIn(r"\usepackage{pgfplots}", article)
        self.assertIn(r"\providecommand{\mathdefault}[1]{#1}", article)
        self.assertIn(r"\input{w3d-wave-direction-charts.tex-part}", results)

        # The existing plot pipeline creates these two PGF families before the
        # LaTeX step; the publication must consume them rather than duplicate
        # plotting logic in the manuscript.
        self.assertIn('finalize_plot(fig, outbase, "_dir")', plotter)
        self.assertIn('finalize_plot(fig, outbase, "_tuner")', plotter)
        for name in (
            "w3d_ou3_jonswap_medium_dir.pgf",
            "w3d_ou3_jonswap_low_dir.pgf",
            "w3d_ou3_pmstokes_medium_dir.pgf",
        ):
            self.assertIn(name, charts)
        self.assertIn("w3d_ou3_jonswap_medium_tuner.pgf", frequency)

        # IEEE conference output is two-column.  Publication plots stay inside
        # one column; no figure* or text-width scaling may creep back in.
        figure_sources = charts + frequency
        self.assertNotIn(r"\begin{figure*}", figure_sources)
        self.assertNotIn(r"\textwidth", figure_sources)
        self.assertEqual(figure_sources.count(r"\resizebox{\columnwidth}{!}"), 4)


if __name__ == "__main__":
    unittest.main()
