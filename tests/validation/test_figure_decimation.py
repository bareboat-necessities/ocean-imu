import importlib.util
from pathlib import Path
import re
import shutil
import tempfile
import unittest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("figure_decimation", REPO_ROOT / "plots" / "figure_decimation.py")
assert SPEC is not None and SPEC.loader is not None
figure_decimation = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(figure_decimation)


def _trace(n=240_000):
    t = np.arange(n) / 200.0
    y = np.sin(2 * np.pi * 0.3 * t) + 0.01 * np.cos(2 * np.pi * 37.0 * t)
    y[123_457] = 5.0          # an isolated spike must survive
    return t, y


class FigureDecimationTests(unittest.TestCase):
    def test_time_series_keeps_budget_endpoints_and_extremes(self):
        t, y = _trace()
        fig, ax = plt.subplots()
        (line,) = ax.plot(t, y)
        self.assertEqual(figure_decimation.decimate_figure(fig, 2000), 1)
        x2, y2 = line.get_xdata(), line.get_ydata()
        self.assertLessEqual(len(x2), 2002)
        self.assertEqual(x2[0], t[0])
        self.assertEqual(x2[-1], t[-1])
        self.assertTrue(np.all(np.diff(x2) > 0))
        self.assertEqual(y2.max(), y.max())
        self.assertEqual(y2.min(), y.min())
        plt.close(fig)

    def test_zoomed_window_keeps_resolution_and_limits(self):
        t, y = _trace()
        fig, ax = plt.subplots()
        (line,) = ax.plot(t, y)
        ax.set_xlim(600.0, 660.0)                   # 12 000 samples in view
        ylim = ax.get_ylim()
        self.assertEqual(figure_decimation.decimate_figure(fig, 2000), 1)
        x2, y2 = line.get_xdata(), line.get_ydata()
        in_view = (x2 >= 600.0) & (x2 <= 660.0)
        self.assertGreaterEqual(in_view.sum(), 1900)
        self.assertLessEqual(len(x2), 2020)
        self.assertEqual((x2[0], x2[-1]), (t[0], t[-1]))
        self.assertEqual((y2.min(), y2.max()), (y.min(), y.max()))   # spike outside the view
        self.assertEqual(ax.get_ylim(), ylim)
        self.assertEqual(ax.get_xlim(), (600.0, 660.0))
        plt.close(fig)

    def test_short_lines_markers_and_trajectories(self):
        fig, ax = plt.subplots()
        (short,) = ax.plot(np.arange(100), np.arange(100))
        (markers,) = ax.plot(np.arange(10_000), np.arange(10_000), "o")
        theta = np.linspace(0, 40 * np.pi, 50_000)
        (spiral,) = ax.plot(np.cos(theta) * theta, np.sin(theta) * theta)
        figure_decimation.decimate_figure(fig, 1000)
        self.assertEqual(len(short.get_xdata()), 100)
        self.assertEqual(len(markers.get_xdata()), 10_000)
        self.assertLessEqual(len(spiral.get_xdata()), 1001)
        self.assertEqual(spiral.get_xdata()[-1], np.cos(theta[-1]) * theta[-1])
        plt.close(fig)

    def test_install_decimates_saved_svg_and_pgf(self):
        figure_decimation.install(matplotlib)
        figure_decimation.install(matplotlib)          # idempotent
        t, y = _trace()
        tex = next((e for e in ("lualatex", "pdflatex", "xelatex") if shutil.which(e)), None)
        formats = ("svg", "pgf") if tex else ("svg",)
        with tempfile.TemporaryDirectory() as tmp, matplotlib.rc_context({"pgf.texsystem": tex or "xelatex"}):
            for ext in formats:
                fig, ax = plt.subplots()
                ax.plot(t, y)
                out = Path(tmp) / f"trace.{ext}"
                fig.savefig(out)
                plt.close(fig)
                text = out.read_text(encoding="utf-8")
                # Each vertex is one coordinate pair in either format.
                if ext == "svg":
                    path = max(re.findall(r'd="([^"]+)"', text), key=len)
                    vertices = len(re.findall(r"[ML]", path))
                else:
                    vertices = text.count("\\pgfpathlineto")
                self.assertLess(vertices, 2100, ext)


if __name__ == "__main__":
    unittest.main()
