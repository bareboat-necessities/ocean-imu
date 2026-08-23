from pathlib import Path
import re
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
README = REPO_ROOT / "README.md"
GALLERY = REPO_ROOT / "reports" / "results" / "readme"
PROVENANCE = REPO_ROOT / "reports" / "readme-results-provenance.md"
PUBLISH_WORKFLOW = REPO_ROOT / ".github" / "workflows" / "readme-results-publish.yml"

EXPECTED_GALLERY = (
    "w3d_ou3_pmstokes_medium.svg",
    "w3d_ou3_pmstokes_medium_zkin.svg",
    "w3d_ou3_pmstokes_medium_xykin.svg",
    "w3d_ou3_pmstokes_medium_acc_bias.svg",
    "w3d_ou3_pmstokes_medium_gyro_bias.svg",
    "spectrum_pmstokes_medium_3d.svg",
)


class ReadmeResultGalleryTests(unittest.TestCase):
    def test_expected_ou3_gallery_is_committed_as_svg(self):
        for name in EXPECTED_GALLERY:
            with self.subTest(name=name):
                path = GALLERY / name
                self.assertTrue(path.is_file(), f"missing README result SVG: {path}")
                self.assertEqual(path.suffix, ".svg")
                header = path.read_bytes()[:4096].lower()
                self.assertIn(b"<svg", header, f"not an SVG document: {path}")

    def test_readme_references_every_expected_ou3_svg(self):
        text = README.read_text(encoding="utf-8")
        for name in EXPECTED_GALLERY:
            with self.subTest(name=name):
                expected = f"./reports/results/readme/{name}"
                self.assertIn(expected, text)

    def test_every_local_readme_image_resolves_in_repository(self):
        text = README.read_text(encoding="utf-8")
        sources = re.findall(r'<img\b[^>]*\bsrc=["\']([^"\']+)["\']', text, flags=re.I)
        self.assertTrue(sources, "README contains no <img> sources to validate")

        for source in sources:
            with self.subTest(source=source):
                if source.startswith(("http://", "https://", "data:")):
                    continue
                relative = source.split("?", 1)[0].split("#", 1)[0]
                path = REPO_ROOT / relative.removeprefix("./")
                self.assertTrue(path.is_file(), f"README image target does not exist: {source}")

    def test_run_specific_provenance_is_outside_scientific_results(self):
        self.assertFalse(
            (GALLERY / "PROVENANCE.md").exists(),
            "build-run provenance must not live under reports/results",
        )
        self.assertTrue(PROVENANCE.is_file())

        workflow = PUBLISH_WORKFLOW.read_text(encoding="utf-8")
        self.assertIn("reports/readme-results-provenance.md", workflow)
        self.assertNotIn("reports/ou_evidence_fingerprint.json", workflow)
        self.assertNotIn("ou_replay_fingerprint.py", workflow)
        self.assertNotIn("sim-data-files.zip", workflow)


if __name__ == "__main__":
    unittest.main()
