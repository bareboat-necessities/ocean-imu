"""Run the fail-closed proof provenance gate in the shared evidence suite."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tools.stability.ou3_theorem import build_evidence  # noqa: E402


class ProofProvenanceTests(unittest.TestCase):
    def test_committed_proof_sources_and_certificates_match(self):
        report = build_evidence.validate()
        self.assertEqual(report["failures"], [])
        self.assertTrue(report["validation_pass"])

    def test_stale_document_hash_fails_without_promoting_theorem(self):
        manifest = json.loads(build_evidence.PROVENANCE.read_text(encoding="utf-8"))
        document = "docs/ou3-construction-mean-action.md"
        row = next(
            item for item in manifest["operation_lemma_sources"]
            if item["path"] == document
        )
        row["git_blob_sha"] = "0" * 40
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "provenance.json"
            path.write_text(json.dumps(manifest), encoding="utf-8")
            with mock.patch.object(build_evidence, "PROVENANCE", path):
                report = build_evidence.validate()
        self.assertFalse(report["validation_pass"])
        self.assertFalse(report["theorem_closed"])
        self.assertTrue(any(
            failure.startswith(f"source provenance changed: {document} ")
            for failure in report["failures"]
        ))


if __name__ == "__main__":
    unittest.main()
