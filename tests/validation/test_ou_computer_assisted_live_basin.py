"""Contract tests for the internal OU-III Live-basin interval proof utility."""

from __future__ import annotations

import json
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "ou_live_basin_interval_proof.py"
ARTICLE = ROOT / "doc" / "kalman_ou_iii" / "kalman_ou-w3d.tex"
BROAD_BOX_TEX = ROOT / "doc" / "kalman_ou_iii" / "w3d-computer-assisted-live-basin.tex-part"


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
        self.assertGreater(v["alpha_controllability_lower"], 0)
        self.assertGreater(v["qH_lower"], 0)
        self.assertGreater(v["s18_lower"], 0)
        self.assertGreater(v["deltaH_lower"], 0)
        self.assertGreater(v["one_minus_chi_lower"], 0)
        self.assertGreater(v["riccati_cert_radius_lower"], 0)

    def test_finite_horizon_controllability_replaces_h7_primary_floor(self):
        v = self.run_proof()
        self.assertEqual(v["sigma_aw_cert_min"], Decimal("0.05"))
        self.assertGreater(
            v["alpha_controllability_lower"],
            Decimal("1e10") * v["q_step_lower"],
        )
        self.assertGreater(
            v["alpha_controllability_lower"],
            v["q_frozen_horizon_lower"],
        )
        self.assertGreater(v["deltaH_lower"], v["delta_frozen_lower"])
        self.assertGreater(v["qba_step_lower"], v["q_step_lower"])

    def test_detectability_and_lifted_nonlinear_closures_are_computed(self):
        v = self.run_proof()
        self.assertGreater(v["measurement_dissipation_lower"], 0)
        self.assertGreater(v["detectability_db_lower"], 0)
        self.assertGreater(v["detectability_G_lower"], 0)
        self.assertGreater(v["delta_detectability_lower"], 0)
        self.assertGreater(v["cxi_upper"], 0)
        self.assertGreater(v["rxi_lower"], 0)
        self.assertGreater(v["CH_upper"], 0)
        self.assertGreater(v["riccati_cert_radius_lower"], 0)
        # This utility remains only an internal worst-case diagnostic.  The
        # publication now uses an exact finite-step regional proof instead.
        self.assertLess(v["riccati_cert_radius_lower"], Decimal("1e-100"))

    def test_publication_excludes_broad_box_certificate(self):
        text = ARTICLE.read_text(encoding="utf-8")
        self.assertIn(r"\input{w3d-live-basin-certificate.tex-part}", text)
        self.assertNotIn(r"\input{w3d-computer-assisted-live-basin.tex-part}", text)
        self.assertFalse(BROAD_BOX_TEX.exists())


if __name__ == "__main__":
    unittest.main()
