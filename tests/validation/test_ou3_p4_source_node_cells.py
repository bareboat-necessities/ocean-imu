from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_source_node_cells as NODES


class P4SourceNodeCellsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = NODES.build()

    def test_exact_P2_partition_and_order(self):
        self.assertEqual(NODES.validate(self.d), [])
        self.assertEqual(
            self.d["partition"],
            {"tau": 10, "sigma_tuner_raw": 8, "R_S": 10, "states": 800},
        )
        self.assertEqual(self.d["nodes"][0]["index"], 0)
        self.assertEqual(self.d["nodes"][-1]["index"], 799)
        for i in (0, 9, 10, 79, 80, 799):
            n = self.d["nodes"][i]
            self.assertEqual(i, ((n["tau_index"] * 8) + n["sigma_raw_index"]) * 10 + n["R_S_index"])

    def test_raw_tuner_sigma_floor_is_not_confused_with_filter_floor(self):
        self.assertLess(
            self.d["raw_tuner_sigma_partition_lower"],
            self.d["filter_sigma_floor_mps2"],
        )
        self.assertGreater(self.d["raw_tuner_sigma_subfloor_node_count"], 0)
        self.assertGreater(self.d["filter_sigma_floor_intersecting_node_count"], 0)
        floor = self.d["filter_sigma_floor_mps2"]
        for n in self.d["nodes"]:
            raw_lo, raw_hi = n["sigma_tuner_raw_mps2"]
            filt_lo, filt_hi = n["sigma_filter_committed_mps2"]
            for target in (max(floor, raw_lo), max(floor, raw_hi)):
                self.assertLessEqual(filt_lo, target)
                self.assertGreaterEqual(filt_hi, target)

    def test_h18_source_cell_preserves_node_coordinates_and_tau_cadence_coupling(self):
        for i in (0, 399, 799):
            n = self.d["nodes"][i]
            c = NODES.h18_source_cell(i, self.d)
            self.assertEqual(c["source_node_index"], i)
            self.assertEqual(c["tau_s"].as_list(), n["tau_s"])
            self.assertEqual(c["sigma_aw_mps2"].as_list(), n["sigma_filter_committed_mps2"])
            self.assertEqual(c["R_S_filter_std"].as_list(), n["R_S_filter_std"])
            self.assertEqual(c["pseudo_period_s"].as_list(), n["pseudo_update_period_s"])

    def test_materialization_is_only_a_metric_attachment_prerequisite(self):
        self.assertFalse(self.d["source_graph_rebuilt_or_pruned_here"])
        self.assertFalse(self.d["P4_metric_attached_here"])
        self.assertFalse(self.d["P4_USABLE_CERTIFICATE_PROMOTED"])
        self.assertIn("actual source-correlated H/A covariance-information metric", self.d["next_obligation"])


if __name__ == "__main__":
    unittest.main()
