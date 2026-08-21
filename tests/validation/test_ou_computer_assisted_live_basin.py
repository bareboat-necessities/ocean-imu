"""Contract tests for the computer-assisted OU-III Live-basin proof."""

from __future__ import annotations

import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "ou_live_basin_interval_proof.py"
DOC = ROOT / "doc" / "kalman_ou_iii" / "w3d-computer-assisted-live-basin.tex-part"
ARTICLE = ROOT / "doc" / "kalman_ou_iii" / "kalman_ou-w3d.tex"


class ComputerAssistedLiveBasinProofTests(unittest.TestCase):
    def run_proof(self):
        cp = subprocess.run(
            [sys.executable, str(TOOL), "--repo-root", str(ROOT), "--json"],
            check=True,
            text=True,
            capture_output=True,
        )
        return {k: Decimal(v) for k, v in json.loads(cp.stdout).items()}

    def test_verified_uniform_margin_is_strictly_positive(self):
        v = self.run_proof()
        self.assertGreater(v["qH_lower"], 0)
        self.assertGreater(v["s18_lower"], 0)
        self.assertGreater(v["deltaH_lower"], 0)
        self.assertGreater(v["one_minus_chi_lower"], 0)

    def test_frozen_schedule_is_only_a_reported_strengthening(self):
        v = self.run_proof()
        self.assertGreater(v["q_frozen_horizon_lower"], v["q_step_lower"])
        self.assertGreater(v["delta_frozen_lower"], v["deltaH_lower"])
        # The uniform proof intentionally remains extremely conservative.
        self.assertLess(v["deltaH_lower"], Decimal("1e-90"))

    def test_article_states_scope_and_no_simulation_dependency(self):
        text = DOC.read_text(encoding="utf-8")
        for marker in (
            "Computer-Assisted Analytical Closure",
            r"\label{eq:iss-ca-delta-num}",
            r"\label{eq:iss-ca-chi-margin}",
            "not a useful numerical handoff radius",
            "no reference-sea simulation is needed",
            "production tuner has no positive lower clamp",
            "frozen-schedule comparison",
        ):
            self.assertIn(marker, text)

    def test_publication_wires_both_constructive_certificate_parts(self):
        text = ARTICLE.read_text(encoding="utf-8")
        self.assertIn(r"\input{w3d-live-basin-certificate.tex-part}", text)
        self.assertIn(r"\input{w3d-computer-assisted-live-basin.tex-part}", text)


if __name__ == "__main__":
    unittest.main()
