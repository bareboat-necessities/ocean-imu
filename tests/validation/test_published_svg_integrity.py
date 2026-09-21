#!/usr/bin/env python3
"""Published SVG figures are complete documents, not truncated fragments.

A figure that only gets part of the way to disk still satisfies the
byte-for-byte mirror contract, because the evidence copy and the article copy
are truncated together.  Nothing downstream complains either: the article
guards every include with ``\\IfFileExists``, the file does exist, and Inkscape
converts the fragment without an error, so the figure reaches the PDF as a
blank rectangle.  The OU-III vibration-level figure shipped that way, with only
a stub of one axis rendered.  These tests read the committed bytes.
"""
from pathlib import Path
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ElementTree
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
import engine_noise_degradation as engine
from matplotlib.figure import Figure

SVG_NAMESPACE = "http://www.w3.org/2000/svg"
PUBLISHED_ROOTS = ("doc", "img", "reports")
ENGINE_RESULTS = ROOT / "reports" / "results" / "engine_noise_degradation"
ENGINE_SUMMARY = ENGINE_RESULTS / "engine_noise_summary.csv"


def published_svgs() -> list[Path]:
    found: list[Path] = []
    for name in PUBLISHED_ROOTS:
        found.extend(sorted((ROOT / name).rglob("*.svg")))
    return found


def read_summary() -> list[dict[str, str]]:
    import csv

    with ENGINE_SUMMARY.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def rendered_panel_count(writer, *args) -> int:
    """Panels the generator draws, taken from the figure it would have saved."""

    figures: list[Figure] = []
    with tempfile.TemporaryDirectory() as directory:
        with patch.object(
            Figure,
            "savefig",
            autospec=True,
            side_effect=lambda fig, *a, _figures=figures, **kw: _figures.append(fig),
        ):
            writer(Path(directory) / "figure.svg", *args)
    if len(figures) != 1:
        raise AssertionError(f"{writer.__name__} saved {len(figures)} figures")
    return len(figures[0].axes)


class PublishedSvgIntegrityTests(unittest.TestCase):
    def test_every_published_svg_parses_and_closes(self):
        paths = published_svgs()
        # A guard that silently matches nothing is worse than no guard.
        self.assertGreater(len(paths), 20)
        for path in paths:
            with self.subTest(path=str(path.relative_to(ROOT))):
                payload = path.read_bytes()
                self.assertTrue(
                    payload.rstrip().endswith(b"</svg>"),
                    f"{path} does not end with a closing svg tag",
                )
                try:
                    root = ElementTree.fromstring(payload)
                except ElementTree.ParseError as error:
                    self.fail(f"{path} is not well-formed XML: {error}")
                self.assertEqual(root.tag, f"{{{SVG_NAMESPACE}}}svg", path)
                for attribute in ("width", "height"):
                    value = root.get(attribute)
                    self.assertIsNotNone(value, f"{path} declares no {attribute}")
                    self.assertNotEqual(value.strip(), "", path)

    def test_engine_noise_figures_keep_every_panel(self):
        summaries = read_summary()
        cases = (
            (engine.SPEED_PLOT_NAME, engine.write_speed_plot),
            (engine.MECHANISM_PLOT_NAME, engine.write_mechanism_plot),
        )
        for name, writer in cases:
            with self.subTest(name=name):
                expected = rendered_panel_count(writer, summaries)
                # Matplotlib labels each Axes group `axes_1`, `axes_2`, ...
                for path in (ENGINE_RESULTS / name, ROOT / "doc/kalman_ou_iii" / name):
                    text = path.read_text(encoding="utf-8")
                    # A failure here prints counts, never the figure source.
                    drawn = sum(
                        f'<g id="axes_{index}">' in text
                        for index in range(1, expected + 2)
                    )
                    self.assertEqual(
                        drawn,
                        expected,
                        f"{path.relative_to(ROOT)} carries {drawn} of the "
                        f"{expected} panels {writer.__name__} draws",
                    )


if __name__ == "__main__":
    unittest.main()
