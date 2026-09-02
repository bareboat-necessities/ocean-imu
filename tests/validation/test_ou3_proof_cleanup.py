#!/usr/bin/env python3
"""Contracts for the flattened OU-III proof-tool tree."""
from __future__ import annotations

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"


class Ou3ProofCleanupTest(unittest.TestCase):
    def test_no_retired_p5_aliases_or_imports_remain(self):
        offenders = []
        for path in sorted(TOOLS.glob("ou3_*.py")):
            text = path.read_text(encoding="utf-8")
            if path.name.startswith("ou3_p5_"):
                offenders.append(f"retired filename: {path.name}")
            if "ou3_p5_" in text:
                offenders.append(f"retired reference: {path.name}")
        self.assertEqual([], offenders, "\n".join(offenders))


if __name__ == "__main__":
    unittest.main()
