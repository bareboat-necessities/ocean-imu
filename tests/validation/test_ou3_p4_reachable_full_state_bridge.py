#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tools"))

import ou3_p4_reachable_full_state_bridge as BRIDGE


class P4ReachableFullStateBridgeTests(unittest.TestCase):
    def test_cross_block_budget_is_schur_threshold_and_fail_closed(self):
        translation = {
            "P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS": "PASS",
            "modes": {
                "H": {"complete_word_translation_margin_lower": 9.0e-24},
                "A": {"complete_word_translation_margin_lower": 16.0e-24},
            },
        }
        bottleneck = {
            "modes": {
                "H": {"existing_direct_nontranslation_margin_lower": 4.0e-28},
                "A": {"existing_direct_nontranslation_margin_lower": 9.0e-28},
            }
        }
        path = {
            "path_graph_ready": True,
            "partition": {"states": 800},
            "recurrent_states": 800,
            "strongly_connected_components": 1,
            "old_worst_corner_has_internal_recurrent_cycle": True,
        }
        d = BRIDGE.build(translation, bottleneck, path)
        self.assertEqual([], BRIDGE.validate(d))
        self.assertEqual("NOT_ESTABLISHED", d["P4_USABLE_CERTIFICATE_STATUS"])
        self.assertAlmostEqual(6.0e-26, d["modes"]["H"]["normalized_cross_block_spectral_norm_budget_upper_open"])
        self.assertAlmostEqual(1.2e-25, d["modes"]["A"]["normalized_cross_block_spectral_norm_budget_upper_open"])
        for mode in ("H", "A"):
            self.assertFalse(d["modes"][mode]["cross_block_bound_validated"])
            self.assertFalse(d["modes"][mode]["full_state_linear_certificate_established"])

    def test_nonpass_translation_is_rejected(self):
        d = BRIDGE.build(
            {"P4_COMPLETE_TRANSLATION_WORST_CELL_STATUS": "NOT_ESTABLISHED", "modes": {}},
            {"modes": {}},
            {"path_graph_ready": True, "partition": {"states": 1}, "recurrent_states": 1},
        )
        f = BRIDGE.validate(d)
        self.assertTrue(any("translation input is not PASS" in x for x in f))


if __name__ == "__main__":
    unittest.main()
