import hashlib
import importlib.util
from pathlib import Path
import tempfile
import time
import unittest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO_ROOT / "plots" / "svg_determinism.py"
SPEC = importlib.util.spec_from_file_location("svg_determinism", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
svg_determinism = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(svg_determinism)


class SvgReproducibilityTests(unittest.TestCase):
    def _render(self, path: Path) -> str:
        with matplotlib.rc_context():
            svg_determinism.configure_svg(matplotlib)
            figure, axis = plt.subplots(figsize=(4.0, 3.0))
            axis.plot(
                [0.0, 1.0, 2.0],
                [0.0, 1.0, 0.25],
                marker="o",
                label="deterministic input",
            )
            axis.set_xlabel("time")
            axis.set_ylabel("value")
            axis.set_title("SVG determinism sentinel")
            axis.legend()
            figure.tight_layout()
            svg_determinism.save_svg(figure, path, bbox_inches="tight")
            plt.close(figure)

        payload = path.read_bytes()
        self.assertNotIn(b"<dc:date>", payload)
        return hashlib.sha256(payload).hexdigest()

    def test_same_input_renders_to_identical_sha256_svg(self):
        with tempfile.TemporaryDirectory() as tmp:
            first = Path(tmp) / "first.svg"
            second = Path(tmp) / "second.svg"

            first_sha = self._render(first)
            # Ensure a default Matplotlib creation timestamp would have moved.
            time.sleep(0.05)
            second_sha = self._render(second)

            self.assertEqual(first_sha, second_sha)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_readme_svg_producers_use_shared_deterministic_export(self):
        spectrum = (
            REPO_ROOT / "plots" / "spectrum" / "spectrum-plots.py"
        ).read_text(encoding="utf-8")
        wrapper = (
            REPO_ROOT / "plots" / "kalman_ou_iii" / "draw_plots.sh"
        ).read_text(encoding="utf-8")

        for text in (spectrum, wrapper):
            self.assertIn(
                "from svg_determinism import configure_svg, save_svg", text
            )
            self.assertIn("configure_svg(mpl)", text)
            self.assertIn("save_svg(fig,", text)


if __name__ == "__main__":
    unittest.main()
