from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_p4_h18_shared_word_screen as SHARED


class H18SharedWordScreenTests(unittest.TestCase):
    def test_short_screen_uses_shared_operations_and_stays_fail_closed(self):
        d = SHARED.build(samples=2, cell_limit=1)
        self.assertEqual(SHARED.validate(d), [])
        self.assertTrue(d["shared_H18_differential_operations_used"])
        self.assertEqual(d["shared_operation_module"], "ou3_p4_h18_differential_operations")
        self.assertEqual(d["dimension"], 18)
        self.assertFalse(d["actual_per_node_Sigma_KF_whitening_used"])
        self.assertFalse(d["source_graph_all_reachable_edges_checked"])
        self.assertFalse(d["P4_USABLE_CERTIFICATE_PROMOTED"])


if __name__ == "__main__":
    unittest.main()
