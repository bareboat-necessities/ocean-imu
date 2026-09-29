"""Proof-record refreshes must validate the same final tree as publication."""
from fnmatch import fnmatchcase
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / ".github/workflows/ou-full-evidence-branch.yml"


class BranchEvidenceTriggerTests(unittest.TestCase):
    def test_proof_refresh_starts_full_validation(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        push = workflow.split("  push:\n", 1)[1].split("\nconcurrency:", 1)[0]
        paths = re.findall(r'^      - "([^"\n]+)"$', push, re.MULTILINE)
        self.assertTrue(paths)
        for path in ("reports/results/ou3_stability/provenance.json",
                     "reports/results/ou3_stability/world-frame-source-feasibility.json"):
            with self.subTest(path=path):
                self.assertTrue(any(fnmatchcase(path, pattern) for pattern in paths))
        # The publishing bot must not recursively regenerate its own output.
        for path in ("reports/results/ou_validation/ou_validation_manifest.json",
                     "reports/results/ou_robustness/ou_robustness_manifest.json",
                     "reports/results/tfg_comparison/tfg_comparison_manifest.json",
                     "reports/ou_evidence_fingerprint.json"):
            with self.subTest(path=path):
                self.assertFalse(any(fnmatchcase(path, pattern) for pattern in paths))
        self.assertIn("branches-ignore: [main]", push)
        self.assertIn("cancel-in-progress: true", workflow)
        self.assertIn("validation_mode: full", workflow)


if __name__ == "__main__":
    unittest.main()
