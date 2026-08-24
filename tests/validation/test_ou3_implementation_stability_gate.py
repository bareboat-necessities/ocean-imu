import sys
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))

import ou3_implementation_stability_gate as GATE
import ou3_source_domain_contract as SOURCE


class Ou3ImplementationStabilityGateTests(unittest.TestCase):
    def test_live_only_or_missing_downstream_proof_cannot_report_stable(self):
        # Empty downstream proof artifacts must fail, even though the source
        # manifest, startup certificate, and declared word language are rebuilt
        # and may themselves pass.
        source = SOURCE.build(SOURCE.DEFAULT_HEADER.resolve())
        out = GATE.compose({}, source, {})
        self.assertTrue(out["implementation_manifest_pass"])
        self.assertTrue(out["startup_certificate_pass"])
        self.assertTrue(out["source_complete_word_language_pass"])
        self.assertFalse(out["downstream_deployment_theorem_pass"])
        self.assertEqual(out["implementation_stability_certificate"], "FAIL")

    def test_final_gate_names_existing_replay_as_evidence_not_theorem_derivation(self):
        source = SOURCE.build(SOURCE.DEFAULT_HEADER.resolve())
        out = GATE.compose({}, source, {})
        role = out["performance_and_replay_role"]
        self.assertIn("regression/falsification", role)
        self.assertIn("not used to derive theorem bounds", role)


if __name__ == "__main__":
    unittest.main()
