from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_h18_source_node_screen as SCREEN


class H18SourceNodeScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = SCREEN.build(source_node_index=0, samples=2, cell_limit=1)

    def test_exact_P2_node_is_bound_into_shared_H18_engine(self):
        self.assertEqual(SCREEN.validate(self.d), [])
        self.assertTrue(self.d["exact_P2_source_node_cell_used"])
        self.assertTrue(self.d["shared_H18_differential_operations_used"])
        self.assertEqual(self.d["P2_source_node_count_available"], 800)
        self.assertEqual(self.d["P2_source_node_index"], 0)
        self.assertEqual(self.d["P2_source_node"]["index"], 0)

    def test_single_node_screen_stays_fail_closed(self):
        self.assertFalse(self.d["all_P2_source_nodes_checked"])
        self.assertFalse(self.d["actual_per_node_Sigma_KF_whitening_used"])
        self.assertFalse(self.d["P4_USABLE_CERTIFICATE_PROMOTED"])
        self.assertIn("all reachable g->h edges", self.d["next_obligation"])

    def test_source_node_preserves_raw_vs_filter_sigma_distinction(self):
        n = self.d["P2_source_node"]
        self.assertLess(n["sigma_tuner_raw_mps2"][0], self.d["P2_source_node"]["sigma_filter_committed_mps2"][0])
        self.assertGreaterEqual(n["sigma_filter_committed_mps2"][0], 0.05)


if __name__ == "__main__":
    unittest.main()
