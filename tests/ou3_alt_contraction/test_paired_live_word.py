"""ALT paired Live-word probe remains non-promoting and same-history."""
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT / "tools/stability")]
from tools.stability.ou3_alt_contraction import paired_live_word as probe


class PairedLiveWordTests(unittest.TestCase):
    def test_paired_measurement_graphs_close_without_promotion(self):
        d = probe.build()
        self.assertEqual(probe.validate(d), [])
        self.assertTrue(d["paired_r_N_S_q_materialized"])
        self.assertTrue(d["nonzero_nominal_residual_exercised"])
        self.assertTrue(d["endogenous_N_or_S_increment_exercised"])
        self.assertTrue(d["rank3_product_port_finite_increment_identities_hold"])
        self.assertTrue(
            d["outward_product_graphs_close_for_all_materialized_paired_cells"])
        self.assertEqual(len(d["event_reports"]), 6)
        self.assertEqual(
            {(e["mode"], e["kind"]) for e in d["event_reports"]},
            {(m, k) for m in ("H", "A")
             for k in ("S_zero", "accelerometer", "magnetometer")},
        )
        for event in d["event_reports"]:
            graph = event["interval_product_graph"]
            self.assertTrue(graph["hard_product_graph_for_this_paired_cell"])
            self.assertTrue(graph["solve_graph_zero_enclosed"])
            self.assertTrue(graph["correction_graph_zero_enclosed"])
        self.assertFalse(d["source_uniform_interval_product_graph_closed"])
        self.assertFalse(
            d["ALT_ACTUAL_SOURCE_UNIFORM_FINITE_INCREMENT_WORD_ATTACHED"])
        self.assertFalse(d["ALT_LIVE_PASS"])
        self.assertFalse(d["P4_promoted"])


if __name__ == "__main__":
    unittest.main()
