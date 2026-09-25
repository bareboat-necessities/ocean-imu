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

from tools.stability.ou3_theorem import build_evidence, construction_mean_action  # noqa: E402


class ProofProvenanceTests(unittest.TestCase):
    def test_committed_proof_sources_and_certificates_match(self):
        report = build_evidence.validate()
        self.assertEqual(report["failures"], [])
        self.assertTrue(report["validation_pass"])

    def test_stale_bound_hash_fails_without_promoting_theorem(self):
        for group, source in (
            ("authoritative_shipping_sources", "src/kalman_ou_iii/Kalman3D_Wave_OU_III.h"),
            ("operation_lemma_sources", "src/kalman_ou_common/KalmanOUCoreMath.h"),
            ("operation_lemma_sources", "docs/ou3-construction-mean-action.md"),
            ("operation_lemma_sources", "reports/results/ou3_stability/construction-mean-action.json"),
        ):
            with self.subTest(source=source):
                manifest = json.loads(build_evidence.PROVENANCE.read_text(encoding="utf-8"))
                row = next(item for item in manifest[group] if item["path"] == source)
                row["git_blob_sha"] = "0" * 40
                with tempfile.TemporaryDirectory() as temporary:
                    path = Path(temporary) / "provenance.json"
                    path.write_text(json.dumps(manifest), encoding="utf-8")
                    with mock.patch.object(build_evidence, "PROVENANCE", path):
                        report = build_evidence.validate()
                self.assertFalse(report["validation_pass"])
                self.assertFalse(report["theorem_closed"])
                self.assertEqual(len(report["failures"]), 1)
                self.assertTrue(report["failures"][0].startswith(
                    f"source provenance changed: {source} "
                ))

    def test_stale_instrumented_header_fails_without_promoting_theorem(self):
        with mock.patch.object(
            construction_mean_action, "instrument", return_value="stale instrumented header\n"
        ):
            report = build_evidence.validate()
        self.assertFalse(report["validation_pass"])
        self.assertFalse(report["theorem_closed"])
        self.assertEqual(report["failures"], ["construction mean header fingerprint changed"])


if __name__ == "__main__":
    unittest.main()
