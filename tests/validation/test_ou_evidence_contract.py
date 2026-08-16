"""Publication/provenance contract for committed OU statistical evidence."""

import csv
import json
import sys
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools"))

import ou_evidence_contract as contract  # noqa: E402


class OUEvidenceContractTests(unittest.TestCase):
    def test_committed_bundles_match_implementation_and_pipeline_hashes(self):
        for study in ("validation", "robustness"):
            self.assertEqual(contract.check(study), [], study)

    def test_statistical_rows_do_not_export_legacy_regression_gate_fields(self):
        for study, directory, csv_name, json_name in (
            (
                "validation",
                REPO_ROOT / "reports" / "results" / "ou_validation",
                "ou_validation_raw.csv",
                "ou_validation.json",
            ),
            (
                "robustness",
                REPO_ROOT / "reports" / "results" / "ou_robustness",
                "ou_robustness_raw.csv",
                "ou_robustness.json",
            ),
        ):
            with (directory / csv_name).open(newline="", encoding="utf-8") as handle:
                fields = csv.DictReader(handle).fieldnames or []
            self.assertNotIn("quality_gate_pass", fields, study)
            self.assertNotIn("simulator_return_code", fields, study)

            bundle = json.loads((directory / json_name).read_text(encoding="utf-8"))
            for row in bundle["raw_runs"]:
                self.assertNotIn("quality_gate_pass", row, study)
                self.assertNotIn("simulator_return_code", row, study)
            protocol = bundle["protocol"]
            self.assertFalse(protocol["simulator_regression_gates_exported"])
            self.assertIn("all completed replays", protocol["replay_inclusion_rule"])

    def test_manifest_covers_filter_implementation_not_only_wave_inputs(self):
        manifest = json.loads(
            (
                REPO_ROOT
                / "reports"
                / "results"
                / "ou_validation"
                / "ou_validation_manifest.json"
            ).read_text(encoding="utf-8")
        )
        implementation = manifest["implementation_files"]
        required = {
            "tests/kalman_ou_ii/kalman_ou_ii-sim.cpp",
            "tests/kalman_ou_iii/kalman_ou_iii-sim.cpp",
            "src/kalman_ou_iii/SeaStateFusionFilter_OU_III.h",
            "src/kalman_ou_iii/Kalman3D_Wave_OU_III.h",
            "src/util/W3dSimCommon.cpp",
        }
        self.assertTrue(required.issubset(implementation), sorted(required - set(implementation)))
        for record in implementation.values():
            self.assertRegex(record["sha256"], r"^[0-9a-f]{64}$")
            self.assertGreater(record["bytes"], 0)

        analysis = manifest["analysis_pipeline_files"]
        self.assertIn("tools/ou_validation.py", analysis)

    def test_gate_count_macros_are_absent_from_publication_inputs(self):
        for path in (
            REPO_ROOT / "reports" / "results" / "ou_validation" / "ou_validation_macros.tex",
            REPO_ROOT / "doc" / "kalman_ou_iii" / "w3d-ou-validation-macros-generated.tex-part",
        ):
            text = path.read_text(encoding="utf-8")
            self.assertNotIn("OUValidationGatePasses", text)
            self.assertNotIn("OUValidationGateFailures", text)


if __name__ == "__main__":
    unittest.main()
