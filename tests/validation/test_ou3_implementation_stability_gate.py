import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_implementation_stability_gate as GATE
import ou3_source_domain_contract as SOURCE


class Ou3ImplementationStabilityGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = SOURCE.build(SOURCE.DEFAULT_HEADER.resolve())
        cls.out = GATE.compose({}, source, {})

    def test_live_only_or_missing_downstream_proof_cannot_report_stable(self):
        out = self.out
        self.assertTrue(out["implementation_manifest_pass"])
        self.assertTrue(out["startup_certificate_pass"])
        self.assertTrue(out["source_complete_word_language_pass"])
        self.assertFalse(out["downstream_deployment_theorem_pass"])
        self.assertEqual(out["implementation_stability_certificate"], "FAIL")

    def test_current_staged_P5_bridge_blocks_final_promotion(self):
        out = self.out
        self.assertTrue(out["generic_deployment_capture_is_not_P5"])
        self.assertTrue(out["P5_local_capture_identifier_is_not_completion_certificate"])
        self.assertFalse(out["P5_finite_startup_capture_pass"])
        self.assertEqual(
            out["P5_first_obstruction"],
            "GRAVITY_QUOTIENT_NONLINEAR_SOURCE_WORD_NOT_CERTIFIED",
        )
        self.assertEqual(out["P5_identification"]["P5_OBSTRUCTION_IDENTIFIED"], "PASS")
        bridge = out["P5_bridge"]
        self.assertEqual(bridge["P5_OUTER_H_BRIDGE_CERTIFICATE"], "NOT_ESTABLISHED")
        self.assertEqual(
            bridge["gauged_full_heading_first_failure"],
            "COMPLETE_SOURCE_WORD_EXACT_CAYLEY_INFORMATION_TRANSPORT_NOT_CERTIFIED",
        )
        self.assertEqual(
            bridge["raw_V_R_large_angle_sector_audit"]["status"],
            "DISPROVED_ON_DECLARED_SOURCE_FAMILY",
        )
        self.assertEqual(bridge["finite_angle_information_geometry"]["status"], "PASS")
        self.assertEqual(
            bridge["yaw_only_quotient_audit"]["obstruction_identified"],
            "PASS",
        )
        self.assertEqual(bridge["detectable_gravity_quotient"]["status"], "PASS")
        self.assertIsNone(bridge["N_H_words"])
        self.assertTrue(any("P5 finite startup-to-inner-funnel capture not established" in x for x in out["failures"]))
        self.assertEqual(out["implementation_stability_certificate"], "FAIL")

    def test_final_gate_names_existing_replay_as_evidence_not_theorem_derivation(self):
        role = self.out["performance_and_replay_role"]
        self.assertIn("regression/falsification", role)
        self.assertIn("not used to derive theorem bounds", role)


if __name__ == "__main__":
    unittest.main()
