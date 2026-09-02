#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from ou3_interval import Interval
import ou3_p3_source_uniform_translation_covariance as U
import ou3_p4_source_node_cells as NODES
import ou3_source_reachable_matrix_p3 as BASE


class SourceUniformTranslationCovarianceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = U.build()

    def test_source_uniform_contract_passes(self):
        self.assertEqual(U.validate(self.payload), [])
        self.assertTrue(self.payload["time_varying_source_parameters_covered_by_pointwise_extrema"])
        self.assertFalse(self.payload["P3_PROMOTED"])

    def test_global_upper_dominates_diagnostic_endpoint_cells(self):
        global_upper = self.payload["Sigma_translation_diagonal_upper"]
        sched = BASE.source_schedule()
        nodes = NODES.build()
        for index in (0, 729):
            node = NODES.node(index, nodes)
            tau = Interval(*map(float, node["tau_s"]))
            sigma = Interval(*map(float, node["sigma_filter_committed_mps2"]))
            rs = Interval(*map(float, node["R_S_filter_std"]))
            local, _timing = BASE.translation_upper(tau, sigma, rs, 1.0, sched)
            for g, x in zip(global_upper, local):
                self.assertGreaterEqual(float(g), float(x))


if __name__ == "__main__":
    unittest.main()
