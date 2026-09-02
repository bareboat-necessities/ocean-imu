#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / "tools"
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

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

    def test_current_interface_is_deliberately_not_ready_without_correlation_certificate(self):
        self.assertFalse(self.d["correlation_candidate"]["provided"])
        self.assertFalse(self.d["P2_CORRELATION_INTERFACE_READY"])
        self.assertFalse(self.d["P2_READY_FOR_CANONICAL_P3"])
        self.assertFalse(self.d["P3_PROMOTED_HERE"])

    def test_same_history_contract_forbids_cartesian_extrema(self):
        self.assertTrue(self.d["independent_cartesian_tau_sigma_RS_extrema_forbidden"])
        self.assertTrue(self.d["P3_must_use_one_common_source_history_for_all_bounds"])
        required = set(self.d["required_correlated_quantities"])
        for name in (
            "tau_applied", "sigma_aw_applied", "R_S_applied",
            "process_excitation_lower", "covariance_upper", "measurement_R_S_lower",
        ):
            self.assertIn(name, required)

    def test_metadata_only_candidate_cannot_pass(self):
        fake = {
            "qualification": I.CANDIDATE_QUALIFICATION,
            "source_only": True,
            "physical_source_states": 800,
            "stage_boundary_pair_states": self.d["stage_boundary_pair_states"],
        }
        x = I.build(candidate=fake)
        self.assertFalse(x["P2_CORRELATION_INTERFACE_READY"])
        self.assertFalse(x["P2_READY_FOR_CANONICAL_P3"])
        self.assertTrue(x["correlation_candidate"]["reasons"])


if __name__ == "__main__":
    unittest.main()
