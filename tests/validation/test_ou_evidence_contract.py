"""Publication/provenance contract for committed OU statistical evidence."""

from __future__ import annotations

import copy
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))

import ou_evidence_contract as contract  # noqa: E402
import ou_evidence_provenance as provenance  # noqa: E402


class OUEvidenceContractTests(unittest.TestCase):
    def test_committed_bundles_match_immutable_replay_provenance(self):
        for study in ("validation", "robustness"):
            self.assertEqual(contract.check(study), [], study)

    def test_statistical_rows_do_not_export_legacy_regression_gate_fields(self):
        for study, directory, csv_name, json_name in (
            ("validation", REPO_ROOT / "reports/results/ou_validation", "ou_validation_raw.csv", "ou_validation.json"),
            ("robustness", REPO_ROOT / "reports/results/ou_robustness", "ou_robustness_raw.csv", "ou_robustness.json"),
        ):
            with (directory / csv_name).open(newline="", encoding="utf-8") as handle:
                fields = csv.DictReader(handle).fieldnames or []
            self.assertNotIn("quality_gate_pass", fields, study)
            self.assertNotIn("simulator_return_code", fields, study)
            bundle = json.loads((directory / json_name).read_text(encoding="utf-8"))
            for row in bundle["raw_runs"]:
                self.assertNotIn("quality_gate_pass", row, study)
                self.assertNotIn("simulator_return_code", row, study)
            self.assertFalse(bundle["protocol"]["simulator_regression_gates_exported"])
            self.assertIn("all completed replays", bundle["protocol"]["replay_inclusion_rule"])

    def test_manifest_separates_replay_from_restatement(self):
        manifest = json.loads((REPO_ROOT / "reports/results/ou_validation/ou_validation_manifest.json").read_text(encoding="utf-8"))
        self.assertEqual(manifest["schema_version"], 2)
        replay = manifest["replay_provenance"]
        self.assertRegex(replay["git_commit"], r"^[0-9a-f]{40}$")
        self.assertRegex(replay["raw_rows_sha256"], r"^[0-9a-f]{64}$")
        implementation = replay["implementation_files"]
        required = {
            "tests/kalman_ou_ii/kalman_ou_ii-sim.cpp",
            "tests/kalman_ou_iii/kalman_ou_iii-sim.cpp",
            "tests/kalman_ou_iii/Makefile",
            "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
            "src/kalman_ou_iii/Kalman3D_Wave_OU_III.h",
            "src/util/W3dSimCommon.cpp",
        }
        self.assertTrue(required.issubset(implementation), sorted(required - set(implementation)))
        for record in implementation.values():
            self.assertRegex(record["sha256"], r"^[0-9a-f]{64}$")
            self.assertGreater(record["size_bytes"], 0)
        self.assertEqual(manifest["git_commit"], replay["git_commit"])
        if manifest.get("restatement"):
            self.assertIn("analysis_pipeline_files", manifest["restatement"])
            self.assertIn("source_bundle_sha256", manifest["restatement"])

    def test_gate_count_macros_are_absent_from_publication_inputs(self):
        for path in (
            REPO_ROOT / "reports/results/ou_validation/ou_validation_macros.tex",
            REPO_ROOT / "doc/kalman_ou_iii/w3d-ou-validation-macros-generated.tex-part",
        ):
            text = path.read_text(encoding="utf-8")
            self.assertNotIn("OUValidationGatePasses", text)
            self.assertNotIn("OUValidationGateFailures", text)


class ProvenanceAdversarialTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        (self.root / "src/kalman_ou_iii").mkdir(parents=True)
        (self.root / "src/util").mkdir(parents=True)
        (self.root / "tests/kalman_ou_iii").mkdir(parents=True)
        (self.root / "tools").mkdir(parents=True)
        self.out = self.root / "reports/results/ou_validation"
        self.out.mkdir(parents=True)

        self.header = self.root / "src/kalman_ou_iii/Kalman3D_Wave_OU_III.h"
        self.header.write_text("#pragma once\ninline int model_version(){return 1;}\n", encoding="utf-8")
        self.sim = self.root / "tests/kalman_ou_iii/kalman_ou_iii-sim.cpp"
        self.sim.write_text('#include "kalman_ou_iii/Kalman3D_Wave_OU_III.h"\nint main(){return model_version()-1;}\n', encoding="utf-8")
        (self.root / "src/util/W3dSimCommon.cpp").write_text("// fixture\n", encoding="utf-8")
        (self.root / "tests/kalman_ou_iii/Makefile").write_text("CXXFLAGS=-O3 -std=c++20\n", encoding="utf-8")
        self.analysis = self.root / "tools/ou_validation.py"
        self.analysis.write_text("# analysis v1\n", encoding="utf-8")

        self.raw = self.out / "ou_validation_raw.csv"
        self.raw.write_text("scenario,value\nA,1.0\n", encoding="utf-8")
        self.bundle = self.out / "ou_validation.json"
        self.bundle.write_text(json.dumps({
            "protocol": {
                "simulator_regression_gates_exported": False,
                "replay_inclusion_rule": "all completed replays with machine-readable metrics",
            },
            "raw_runs": [{"scenario": "A", "value": 1.0}],
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        self.manifest = self.out / "ou_validation_manifest.json"

        self.cfg = {
            "validation": {
                "dir": self.out,
                "raw_csv": "ou_validation_raw.csv",
                "json": "ou_validation.json",
                "manifest": "ou_validation_manifest.json",
                "macro": None,
                "mirrored_macro": None,
                "cpp_roots": (self.sim, self.root / "src/util/W3dSimCommon.cpp"),
                "build_files": (self.root / "tests/kalman_ou_iii/Makefile",),
                "analysis_files": (self.analysis,),
            }
        }
        self.patches = (
            mock.patch.object(provenance, "REPO_ROOT", self.root),
            mock.patch.object(provenance, "STUDIES", self.cfg),
            mock.patch.object(contract, "STUDIES", self.cfg),
        )
        for patch in self.patches:
            patch.start()
        replay = {
            "git_commit": "0" * 40,
            "implementation_files": provenance.implementation_records("validation"),
            "input_files": {},
            "raw_rows_sha256": provenance.sha256_file(self.raw),
            "raw_rows_size_bytes": self.raw.stat().st_size,
            "raw_rows_schema": "normalized",
            "environment": {},
        }
        payload = {
            "schema_version": 2,
            "replay_provenance": replay,
            "git_commit": replay["git_commit"],
            "implementation_files": {name: provenance.legacy_record(record) for name, record in replay["implementation_files"].items()},
            "source_files": {},
            "result_files": {
                self.raw.name: provenance.legacy_record(provenance.file_record(self.raw)),
                self.bundle.name: provenance.legacy_record(provenance.file_record(self.bundle)),
            },
        }
        self.manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    def tearDown(self):
        for patch in reversed(self.patches):
            patch.stop()
        self.temp.cleanup()

    def _manifest_payload(self):
        return json.loads(self.manifest.read_text(encoding="utf-8"))

    def test_safe_restatement_keeps_replay_provenance_and_records_new_analysis(self):
        before = copy.deepcopy(self._manifest_payload()["replay_provenance"])
        self.analysis.write_text("# analysis v2\n", encoding="utf-8")
        context = provenance.begin_restatement("validation", self.bundle)
        with mock.patch.object(provenance, "git_output", return_value="1" * 40):
            after_manifest = provenance.finalize_restatement_manifest("validation", self._manifest_payload(), context)
        self.assertEqual(after_manifest["replay_provenance"], before)
        self.assertEqual(after_manifest["git_commit"], before["git_commit"])
        self.assertEqual(after_manifest["restatement"]["git_commit"], "1" * 40)
        self.assertEqual(
            after_manifest["restatement"]["analysis_pipeline_files"]["tools/ou_validation.py"],
            provenance.file_record(self.analysis),
        )

    def test_estimator_change_requires_full_replay(self):
        before = self.manifest.read_bytes()
        self.header.write_text("#pragma once\ninline int model_version(){return 2;}\n", encoding="utf-8")
        with self.assertRaisesRegex(provenance.ProvenanceError, "full replay is required"):
            provenance.begin_restatement("validation", self.bundle)
        self.assertEqual(self.manifest.read_bytes(), before)
        self.assertNotEqual(contract._auto(("validation",)), 0)

    def test_simulator_change_requires_full_replay(self):
        self.sim.write_text(self.sim.read_text(encoding="utf-8") + "// changed\n", encoding="utf-8")
        with self.assertRaisesRegex(provenance.ProvenanceError, "full replay is required"):
            provenance.begin_restatement("validation", self.bundle)

    def test_restatement_never_replaces_replay_commit(self):
        context = provenance.begin_restatement("validation", self.bundle)
        old = context.replay_provenance["git_commit"]
        with mock.patch.object(provenance, "git_output", return_value="f" * 40):
            manifest = provenance.finalize_restatement_manifest("validation", self._manifest_payload(), context)
        self.assertEqual(manifest["replay_provenance"]["git_commit"], old)
        self.assertEqual(manifest["git_commit"], old)
        self.assertEqual(manifest["restatement"]["git_commit"], "f" * 40)

    def test_raw_rows_identity_is_immutable(self):
        self.raw.write_text("scenario,value\nA,9.0\n", encoding="utf-8")
        with self.assertRaisesRegex(provenance.ProvenanceError, "raw replay rows"):
            provenance.begin_restatement("validation", self.bundle)

    def test_archive_mode_needs_no_git_metadata(self):
        with mock.patch.object(provenance, "git_available", return_value=False):
            self.assertEqual(contract.check("validation"), [])


if __name__ == "__main__":
    unittest.main()
