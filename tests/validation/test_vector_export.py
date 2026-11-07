"""Rendering budgets must not change the scientific data behind a figure."""
import ast
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "plots"))
from vector_export import MAX_PLOT_POINTS, envelope_indices, save_figure


class VectorExportTests(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def test_extrema_endpoints_and_gaps(self):
        t = np.arange(100003)
        y = np.zeros(len(t))
        y[12345], y[54321], y[75000] = 100, -100, np.nan
        i = envelope_indices((t, y))
        self.assertLessEqual(len(i), MAX_PLOT_POINTS)
        self.assertTrue(np.all(np.diff(i) > 0))
        for index in (0, 12345, 54321, 75000, len(t) - 1):
            self.assertIn(index, i)
        np.testing.assert_array_equal(envelope_indices((t[:10], y[:10])), np.arange(10))

    def test_svg_and_pgf_share_geometry_and_restore_data(self):
        fig, ax = plt.subplots()
        x = np.arange(20003)
        y = np.sin(x)
        line, = ax.plot(x, y)
        band = ax.fill_between(x, y - 1, y + 1)
        path = band.get_paths()[0]
        original_vertices = path.vertices
        original_codes = path.codes
        rendered = []

        def capture(*args, **kwargs):
            self.assertLessEqual(len(line.get_xdata()), MAX_PLOT_POINTS)
            self.assertLessEqual(len(path.vertices), MAX_PLOT_POINTS)
            rendered.append(line.get_xydata().copy())

        with patch.object(fig, "savefig", side_effect=capture):
            save_figure(fig, "chart.svg")
            save_figure(fig, "chart.pgf")
        np.testing.assert_array_equal(rendered[0], rendered[1])
        np.testing.assert_array_equal(line.get_ydata(), y)
        self.assertIs(path.vertices, original_vertices)
        self.assertIs(path.codes, original_codes)
        with patch.object(fig, "savefig", side_effect=RuntimeError("backend failed")):
            with self.assertRaises(RuntimeError):
                save_figure(fig, "chart.pgf")
        np.testing.assert_array_equal(line.get_ydata(), y)
        self.assertIs(path.vertices, original_vertices)

    def test_overlaid_filters_share_the_axes_budget(self):
        fig, ax = plt.subplots()
        t = np.arange(10001)
        for offset in range(4):
            ax.plot(t, np.sin(t) + offset)

        def capture(*args, **kwargs):
            self.assertLessEqual(sum(len(line.get_xdata()) for line in ax.lines), MAX_PLOT_POINTS)

        with patch.object(fig, "savefig", side_effect=capture):
            save_figure(fig, "filters.pgf")
        self.assertTrue(all(len(line.get_xdata()) == len(t) for line in ax.lines))

    def test_plot_generators_use_the_shared_exporter(self):
        root = Path(__file__).resolve().parents[2]
        for directory in (root / "plots", root / "tools"):
            for source in directory.rglob("*.py"):
                if source.name == "vector_export.py":
                    continue
                tree = ast.parse(source.read_text())
                direct_exports = [node.lineno for node in ast.walk(tree)
                                  if isinstance(node, ast.Call)
                                  and isinstance(node.func, ast.Attribute)
                                  and node.func.attr == "savefig"]
                self.assertFalse(direct_exports, f"{source}: bypasses bounded export at {direct_exports}")

    def test_trajectory_coordinates_stay_aligned(self):
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        t = np.arange(10001)
        line, = ax.plot(t, t * 2, t * 3)

        def capture(*args, **kwargs):
            x, y, z = line.get_data_3d()
            self.assertLessEqual(len(x), MAX_PLOT_POINTS)
            np.testing.assert_array_equal(y, x * 2)
            np.testing.assert_array_equal(z, x * 3)

        with patch.object(fig, "savefig", side_effect=capture):
            save_figure(fig, "trajectory.svg")
        np.testing.assert_array_equal(line.get_data_3d()[0], t)


if __name__ == "__main__":
    unittest.main()
