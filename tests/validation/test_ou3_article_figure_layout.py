#!/usr/bin/env python3
"""The five requested article figures are horizontal, single-column panels."""
from pathlib import Path
import csv
import re
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import ou3_lever_arm_study as lever
import ou3_engine_noise_mitigation as guard
from matplotlib.figure import Figure


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


class ArticleHorizontalPanelTests(unittest.TestCase):
    def test_actual_panels_share_one_row_and_fit_the_column(self):
        summaries = read_rows(lever.DEFAULT_OUTPUT_DIR / "lever_arm_summary.csv")
        runs = read_rows(lever.DEFAULT_OUTPUT_DIR / "lever_arm_runs.csv")
        guard_paths = list((ROOT / "reports/results").rglob("guard_summary.csv"))
        self.assertEqual(len(guard_paths), 1)
        guards = read_rows(guard_paths[0])
        d = lever.worst_unmodeled_axis(summaries, 0.3, "disp_3d_ratio_to_baseline")
        t = lever.worst_unmodeled_axis(summaries, 0.3, "tilt_ratio_to_baseline")
        cases = [
            (lever.write_ratio_plot, (summaries, "disp_3d_ratio_to_baseline", "Pooled 3-D RMS / CG baseline"), 3),
            (lever.write_ratio_plot, (summaries, "tilt_ratio_to_baseline", "Pooled max tilt RMS / CG baseline"), 3),
            (lever.write_mechanism_plot, (summaries,), 2),
            (lever.write_sea_state_plot, (runs, d, t, 0.3), 2),
            (guard.write_plot, (guards,), 2),
        ]
        with tempfile.TemporaryDirectory() as directory:
            for function, args, count in cases:
                with self.subTest(function=function.__name__, metric=args[1] if function == lever.write_ratio_plot else ""):
                    figures = []
                    with patch.object(Figure, "savefig", autospec=True,
                                      side_effect=lambda fig, *a, **kw: figures.append(fig)):
                        function(Path(directory) / "figure.svg", *args)
                    self.assertEqual(len(figures), 1)
                    fig = figures[0]
                    self.assertAlmostEqual(fig.get_figwidth(), 3.5)
                    self.assertEqual(len(fig.axes), count)
                    boxes = [axis.get_position() for axis in fig.axes]
                    for box in boxes:
                        self.assertAlmostEqual(box.y0, boxes[0].y0)
                        self.assertAlmostEqual(box.y1, boxes[0].y1)
                    for left, right in zip(boxes, boxes[1:]):
                        self.assertLess(left.x1, right.x0)
                    fig.canvas.draw()
                    renderer = fig.canvas.get_renderer()
                    # Tight extents include rendered labels, not out-of-view locator ticks.
                    artists = [*fig.axes, *fig.texts, *fig.legends]
                    for artist in artists:
                        box = artist.get_tightbbox(renderer)
                        if box is None:
                            continue
                        self.assertGreaterEqual(box.x0, -1)
                        self.assertGreaterEqual(box.y0, -1)
                        self.assertLessEqual(box.x1, fig.bbox.width + 1)
                        self.assertLessEqual(box.y1, fig.bbox.height + 1)


    def test_article_keeps_single_column_floats(self):
        doc = ROOT / "doc/kalman_ou_iii"
        sources = (doc / "w3d-imu-lever-arm-study.tex-part").read_text() + (doc / "w3d-engine-noise-degradation.tex-part").read_text()
        names = ["ou3_lever_arm_penalty", "ou3_lever_arm_tilt", "ou3_lever_arm_mechanism", "ou3_lever_arm_sea_state", "ou_engine_noise_guard"]
        blocks = re.findall(r"\\begin\{figure\}.*?\\end\{figure\}", sources, re.S)
        for name in names:
            selected = [block for block in blocks if "{" + name + "}" in block]
            self.assertEqual(len(selected), 1, name)
            self.assertIn(r"width=\columnwidth", selected[0])


if __name__ == "__main__":
    unittest.main()
