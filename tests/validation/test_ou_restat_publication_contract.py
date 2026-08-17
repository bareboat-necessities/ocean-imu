from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
GENERATOR = REPO_ROOT / "tools" / "ou_validation.py"
PUBLICATION = (
    REPO_ROOT
    / "reports"
    / "results"
    / "ou_validation"
    / "ou_validation_publication.tex"
)

RETIRED_CLAIM = "The vertical gain of OU--III is paid for in the horizontal channels."
CURRENT_CAPTION_END = r"positive $\Delta$ favors OU--II.}"


class OUValidationRestatementPublicationContractTests(unittest.TestCase):
    def test_generator_and_committed_publication_share_current_axis_caption(self):
        generator = GENERATOR.read_text(encoding="utf-8")
        publication = PUBLICATION.read_text(encoding="utf-8")

        self.assertNotIn(RETIRED_CLAIM, generator)
        self.assertNotIn(RETIRED_CLAIM, publication)
        self.assertIn(CURRENT_CAPTION_END, generator)
        self.assertIn(CURRENT_CAPTION_END, publication)


if __name__ == "__main__":
    unittest.main()
