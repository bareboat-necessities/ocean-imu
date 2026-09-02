#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

import ou3_p2_correlation_path_memory as CORR
import ou3_p2_p3_correlation_interface as I


class P2P3CorrelationInterfaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = I.build()

    def test_source_timing_is_kept_but_endpoint_only_p2_is_rejected(self):
        self.assertEqual(I.validate(self.d), [])
        self.assertTrue(self.d["P2_source_timing_mathematics_retained"])
        self.assertTrue(self.d["P2_source_timing_certificate_pass"])
        self.assertTrue(self.d["physical_800_state_partition_retained"])
        self.assertFalse(self.d["endpoint_only_800_state_quotient_sufficient_for_P3"])

    def test_versioned_path_memory_makes_refined_p2_ready(self):
        self.assertTrue(self.d["P2_CORRELATION_INTERFACE_READY"])
        self.assertTrue(self.d["P2_READY_FOR_CANONICAL_P3"])
        self.assertEqual(self.d["correlation_interface_version"], CORR.INTERFACE_VERSION)
        self.assertEqual(self.d["P3_required_correlation_interface_version"], CORR.INTERFACE_VERSION)
        self.assertFalse(self.d["P3_PROMOTED_HERE"])

    def test_same_history_contract_forbids_cartesian_extrema(self):
        self.assertTrue(self.d["independent_cartesian_tau_sigma_RS_extrema_forbidden"])
        self.assertTrue(self.d["P3_must_use_one_common_source_history_for_all_bounds"])
        contract = self.d["correlation_consumer_contract"]
        self.assertTrue(contract["correlated_quantities_must_come_from_same_segment_node"])
        self.assertEqual(contract["independent_tau_sigma_R_S_extremization_before_propagation"], "FORBIDDEN")
        self.assertEqual(contract["global_800_ancestor_hull_as_P3_covariance_information_input"], "FORBIDDEN")

    def test_external_metadata_cannot_replace_repository_certificate(self):
        with self.assertRaises(ValueError):
            I.build(candidate={"qualification": "fake"})


if __name__ == "__main__":
    unittest.main()
