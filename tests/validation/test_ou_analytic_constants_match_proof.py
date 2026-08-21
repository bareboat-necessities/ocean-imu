"""Every analytical constant printed in the article must be the one the proof
program actually verifies.

The computer-assisted subsection quotes about fifteen numbers.  They are the
output of `tools/ou_live_basin_interval_proof.py`, but nothing was checking
that: when the directed arithmetic behind
`beta_process_measurement_upper` was tightened, six downstream constants in the
document kept their old values and the article claimed margins three times
stronger than its own program verified.  The conclusions did not change -- the
margins are strictly positive either way, and the broad-box radius stays
operationally useless -- but a proof document that overstates its own verified
bound is a defect regardless of whether the overstatement matters.

This test compares the printed value against the program's own formatter, so a
future tightening either updates the article or fails here.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
import unittest
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOL = ROOT / "tools" / "ou_live_basin_interval_proof.py"
DOC = ROOT / "doc" / "kalman_ou_iii" / "w3d-computer-assisted-live-basin.tex-part"

# Printed constant -> proof-program key.  The pattern is anchored on the key's
# surrounding LaTeX so a number that moves to a different claim is caught too.
PRINTED = {
    "alpha_controllability_lower": r"\\alpha_c\\ge([0-9.]+)\\times10\^\{(-?\d+)\}",
    "q_step_lower": r"q_\{\\rm step\}\\ge([0-9.]+)\\times10\^\{(-?\d+)\}",
    "q_frozen_horizon_lower": r"gives\s*\n?\$?([0-9.]+)\\times10\^\{(-?\d+)\}\$?\s*\n?and is therefore weaker",
    "sL_lower": r"\\sigma_\{\\min\}\(\\mathcal O_L\)\s*\n?\s*\\ge([0-9.]+)\\times10\^\{(-?\d+)\}",
    "sA_lower": r"\\sigma_\{\\min\}\(\\mathcal O_A\)\s*\n?\s*\\ge([0-9.]+)\\times10\^\{(-?\d+)\}",
    "s18_lower": r"\\sigma_\{\\min\}\(\\mathcal O_\{18\}\)\s*\n?\s*\\ge([0-9.]+)\\times10\^\{(-?\d+)\}",
    "pbar_upper": r"\\overline p\\le([0-9.]+)\\times10\^\{(-?\d+)\}",
    "beta_process_measurement_upper": r"\\beta_c\\le([0-9.]+)\\times10\^\{(-?\d+)\}",
    "qH_lower": r"q_H\\ge([0-9.]+)\\times10\^\{(-?\d+)\}",
    "deltaH_lower": r"\\delta_H\s*\n?\s*\\ge([0-9.]+)\\times10\^\{(-?\d+)\}",
    "one_minus_chi_lower": r"1-\\chi_H\s*\n?\s*\\ge([0-9.]+)\\times10\^\{(-?\d+)\}",
    "delta_detectability_lower": r"\\delta_\{\\rm det\}\s*\n?\s*\\ge([0-9.]+)\\times10\^\{(-?\d+)\}",
    "cxi_upper": r"c_\\xi\\le([0-9.]+)\\times10\^\{(-?\d+)\}",
    "rxi_lower": r"r_\\xi\\ge([0-9.]+)\\times10\^\{(-?\d+)\}",
    "riccati_cert_radius_lower": r"R_\{\\rm cert\}\s*\n?\s*\\ge([0-9.]+)\\times10\^\{(-?\d+)\}",
}


def proof_values() -> dict[str, Decimal]:
    cp = subprocess.run(
        [sys.executable, str(TOOL), "--repo-root", str(ROOT), "--json"],
        check=True,
        text=True,
        capture_output=True,
    )
    return {k: Decimal(v) for k, v in json.loads(cp.stdout).items()}


def as_printed(x: Decimal) -> str:
    """The proof program's own reporting format."""
    return f"{x:.12E}"


class AnalyticConstantsMatchProofTests(unittest.TestCase):
    def setUp(self):
        self.values = proof_values()
        self.text = DOC.read_text(encoding="utf-8")

    def test_every_printed_constant_is_the_verified_one(self):
        for key, pattern in PRINTED.items():
            with self.subTest(constant=key):
                m = re.search(pattern, self.text)
                self.assertIsNotNone(
                    m, f"{key} is no longer printed where this test looks for it"
                )
                printed = Decimal(f"{m.group(1)}E{m.group(2)}")
                verified = self.values[key]
                self.assertEqual(
                    as_printed(printed),
                    as_printed(verified),
                    f"{key}: article prints {printed:.12E}, proof verifies "
                    f"{verified:.12E}",
                )

    def test_a_lower_bound_is_never_printed_above_what_is_verified(self):
        """The direction that matters: overstating a margin is unsound."""
        for key in (k for k in PRINTED if k.endswith("_lower")):
            with self.subTest(constant=key):
                m = re.search(PRINTED[key], self.text)
                printed = Decimal(f"{m.group(1)}E{m.group(2)}")
                # Allow only last-digit rounding of the shared formatter.
                verified = self.values[key]
                self.assertLessEqual(
                    printed,
                    verified * Decimal("1.000000000001"),
                    f"{key} is advertised above the verified bound",
                )


if __name__ == "__main__":
    unittest.main()
