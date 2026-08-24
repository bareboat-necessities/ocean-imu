#!/usr/bin/env python3
"""Pin the implemented OU-III startup proxy to the theorem constants."""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WRAP = (ROOT / "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h").read_text(encoding="utf-8")
THEOREM = (ROOT / "doc/kalman_ou_iii/w3d-semiglobal-stability.tex-part").read_text(encoding="utf-8")
STARTUP_DOC = (ROOT / "docs/ou-iii-startup-init.md").read_text(encoding="utf-8")


class OU3StartupTheoremParityTests(unittest.TestCase):
    def test_wrapper_binds_proxy_gains_to_theorem_constants(self):
        self.assertRegex(
            WRAP,
            r"STARTUP_PROXY_TWO_KP_DEFAULT\s*=\s*0\.2f",
        )
        self.assertRegex(
            WRAP,
            r"STARTUP_PROXY_TWO_KI_DEFAULT\s*=\s*0\.02f",
        )
        self.assertIn(
            "VerticalAccelComplementary      vertical_accel_comp_{\n"
            "        STARTUP_PROXY_TWO_KP_DEFAULT,\n"
            "        STARTUP_PROXY_TWO_KI_DEFAULT};",
            WRAP,
        )

        flat = re.sub(r"\s+", " ", THEOREM)
        self.assertIn(r"2k_P=0.2", flat)
        self.assertIn(r"2k_I=0.02", flat)

    def test_startup_study_matches_the_bound_wrapper(self):
        flat = re.sub(r"\s+", " ", STARTUP_DOC)
        self.assertIn("two_kp = 0.2", flat)
        self.assertIn("two_ki = 0.02", flat)


if __name__ == "__main__":
    unittest.main()
