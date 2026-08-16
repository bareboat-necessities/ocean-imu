"""Stable title constraint requested for the main OU-III publication."""

from __future__ import annotations

import re
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
MAIN = REPO_ROOT / "doc/kalman_ou_iii/kalman_ou-w3d.tex"


class OUTitleContractTests(unittest.TestCase):
    def test_main_ou_iii_title_retains_navigation(self):
        text = MAIN.read_text(encoding="utf-8")
        match = re.search(r"\\title\{([^}]*)\}", text, re.IGNORECASE | re.DOTALL)
        self.assertIsNotNone(match, "OU-III manuscript has no \\title{...}")
        self.assertIn("navigation", match.group(1).lower())


if __name__ == "__main__":
    unittest.main()
