#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p3_tau_decay_budget as D


class TauDecayBudgetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = D.build(window_samples=(100,))

    def test_source_contract_is_fail_closed(self):
        self.assertEqual(D.validate(self.payload), [])
        self.assertFalse(self.payload["P3_PROMOTED"])
        self.assertEqual(self.payload["clock_phase_gap_alphabet_samples"], list(range(13, 27)))

    def test_long_tau_endpoint_beats_global_lambda_max_bound(self):
        row = self.payload["windows"]["100"][9]
        H = row["window_s_nominal"]
        global_bound = H / (1.0 / 3.0)
        self.assertGreater(global_bound, 0.0)
        self.assertLess(row["decay_exponent_upper"], global_bound)

    def test_short_tau_endpoint_contains_frozen_hold_branch(self):
        row = self.payload["windows"]["100"][0]
        self.assertGreaterEqual(
            row["decay_exponent_upper"], row["frozen_clock_decay_exponent_upper"]
        )
        self.assertTrue(math.isfinite(row["decay_exponent_upper"]))


if __name__ == "__main__":
    unittest.main()
