#!/usr/bin/env python3
"""Pin the implemented OU-III startup proxy to the theorem constants."""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
VERTICAL = (ROOT / "src/tuner/VerticalAccelComplementary.h").read_text(encoding="utf-8")
THEOREM = (ROOT / "doc/kalman_ou_iii/w3d-semiglobal-stability.tex-part").read_text(encoding="utf-8")
STARTUP_DOC = (ROOT / "docs/ou-iii-startup-init.md").read_text(encoding="utf-8")


class OU3StartupTheoremParityTests(unittest.TestCase):
    def test_default_proxy_gains_match_theorem(self):
        self.assertRegex(
            VERTICAL,
            r"VerticalAccelComplementary\(float\s+two_kp\s*=\s*0\.2f,\s*"
            r"float\s+two_ki\s*=\s*0\.02f",
        )
        flat = re.sub(r"\s+", " ", THEOREM)
        self.assertIn(r"2k_P=0.2", flat)
        self.assertIn(r"2k_I=0.02", flat)

    def test_startup_study_and_implementation_use_same_integral_gain(self):
        flat = re.sub(r"\s+", " ", STARTUP_DOC)
        self.assertIn("two_kp = 0.2", flat)
        self.assertIn("two_ki = 0.02", flat)
        self.assertNotIn("float two_ki = 0.0f", VERTICAL)


if __name__ == "__main__":
    unittest.main()
