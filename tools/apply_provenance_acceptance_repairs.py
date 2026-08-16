#!/usr/bin/env python3
from __future__ import annotations

import ast
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def replace_method(text: str, class_name: str, method_name: str, replacement: str) -> str:
    tree = ast.parse(text)
    target = None
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == class_name:
            for child in node.body:
                if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == method_name:
                    target = child
                    break
    if target is None:
        raise SystemExit(f"missing {class_name}.{method_name}")
    lines = text.splitlines(keepends=True)
    start = target.lineno - 1
    end = target.end_lineno
    if not replacement.endswith("\n"):
        replacement += "\n"
    lines[start:end] = [replacement]
    return "".join(lines)


# Restore the complete validation test file first; a previous broad textual
# replacement accidentally deleted unrelated fixtures. Then change only the
# specific method bodies that formerly depended on import-order monkeypatches.
validation_path = ROOT / "tests/validation/test_ou_validation.py"
baseline = subprocess.check_output(
    ["git", "show", "origin/main:tests/validation/test_ou_validation.py"],
    cwd=ROOT,
    text=True,
)
text = baseline
for name, helper in (
    ("test_abstract_reports_committed_stationary_aggregate", "_abstract_reports_committed_stationary_aggregate"),
    ("test_fixed_reference_and_transition_limits_are_stated", "_fixed_reference_and_transition_limits_are_stated"),
    ("test_inference_is_qualified_rather_than_asserted", "_inference_is_qualified_rather_than_asserted"),
    ("test_three_dimensional_and_channel_results_are_reported", "_three_dimensional_and_channel_results_are_reported"),
    ("test_transition_and_secondary_ensembles_are_rescored", "_transition_and_secondary_ensembles_are_rescored"),
    ("test_baseline_fairness_thresholds_and_hardware_limits_are_recorded", "_baseline_fairness_thresholds_and_hardware_limits_are_recorded"),
    ("test_the_contribution_is_framed_at_the_width_of_the_evidence", "_contribution_is_framed_at_the_width_of_the_evidence"),
):
    text = replace_method(
        text,
        "CommittedFullResultsTests" if name == "test_abstract_reports_committed_stationary_aggregate" else "ManuscriptMethodologyTests",
        name,
        f"    def {name}(self):\n        from test_zz_publication_contract import {helper}\n        {helper}(self)\n",
    )
for name, helper in (
    ("test_the_deterministic_table_is_generated_not_typed", "_deterministic_table_remains_generated_evidence"),
    ("test_singer_relationship_and_contribution_wording", "_singer_relationship_and_contribution_boundary"),
    ("test_tracker_and_heel_material_are_scoped_as_unevaluated_appendices", "_unevaluated_tracker_and_heel_material_stays_out_of_main_claims"),
):
    text = replace_method(
        text,
        "ManuscriptMethodologyTests",
        name,
        f"    def {name}(self):\n        from test_reorg_compat_overrides import {helper}\n        {helper}(self)\n",
    )
text = replace_method(
    text,
    "RestatTests",
    "test_a_restat_cannot_manufacture_an_ensemble",
    '''    def test_a_restat_cannot_manufacture_an_ensemble(self):
        with tempfile.TemporaryDirectory() as tmp:
            empty = Path(tmp) / "empty.json"
            empty.write_text(json.dumps({"raw_runs": []}), encoding="utf-8")
            with contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaises(validation.evidence_provenance.ProvenanceError):
                    validation.restat_bundle(empty, Path(tmp), 100, 1)
''',
)
validation_path.write_text(text, encoding="utf-8")

# Publication helper modules are ordinary semantic libraries now. Remove the
# obsolete test that asserted monkeypatch identity; vanilla discovery itself is
# the contract that the direct delegates above work.
pub_path = ROOT / "tests/validation/test_zz_publication_contract.py"
pub_text = pub_path.read_text(encoding="utf-8")
tree = ast.parse(pub_text)
node = next((n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "PublicationOverrideInstallationTests"), None)
if node is None:
    raise SystemExit("PublicationOverrideInstallationTests not found")
lines = pub_text.splitlines(keepends=True)
del lines[node.lineno - 1:node.end_lineno]
pub_path.write_text("".join(lines), encoding="utf-8")

# Robustness tests now assert the actual schema-v2 namespaces instead of the
# deleted moved-source/re-pin mechanism.
robust_path = ROOT / "tests/validation/test_ou_robustness.py"
rtext = robust_path.read_text(encoding="utf-8")
rtext = replace_method(
    rtext,
    "RestatTests",
    "test_a_restat_keeps_covering_the_files_it_did_not_rewrite",
    '''    def test_a_restat_keeps_covering_the_files_it_did_not_rewrite(self):
        """Result inventory stays complete while code provenance has explicit namespaces."""
        manifest = json.loads(
            (self.RESULTS / "ou_robustness_manifest.json").read_text(encoding="utf-8")
        )
        covered = set(manifest["result_files"])
        for name in ("ou_robustness_sensitivity.svg", "ou_robustness_stress.svg"):
            self.assertIn(name, covered, name)
        replay_impl = set(manifest["replay_provenance"]["implementation_files"])
        self.assertIn("tests/kalman_ou_iii/kalman_ou_iii-sim.cpp", replay_impl)
        self.assertIn("src/kalman_ou_iii/Kalman3D_Wave_OU_III.h", replay_impl)
        analysis = set(manifest["restatement"]["analysis_pipeline_files"])
        self.assertIn("tools/ou_validation.py", analysis)
        self.assertIn("tools/ou_robustness.py", analysis)
''',
)
rtext = replace_method(
    rtext,
    "RestatTests",
    "test_a_restat_records_the_sources_that_moved_under_it",
    '''    def test_a_restat_records_the_sources_that_moved_under_it(self):
        """The obsolete warn-and-repin path must not survive schema v2."""
        self.assertFalse(hasattr(robustness, "_restated_source_files"))
        manifest = json.loads(
            (self.RESULTS / "ou_robustness_manifest.json").read_text(encoding="utf-8")
        )
        self.assertNotIn("sources_moved_since_rows", manifest)
        self.assertIn("replay_provenance", manifest)
        self.assertIn("restatement", manifest)
''',
)
rtext = replace_method(
    rtext,
    "RestatTests",
    "test_a_restat_cannot_manufacture_an_ensemble",
    '''    def test_a_restat_cannot_manufacture_an_ensemble(self):
        with tempfile.TemporaryDirectory() as tmp:
            empty = Path(tmp) / "empty.json"
            empty.write_text(json.dumps({"raw_runs": []}), encoding="utf-8")
            with contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaises(robustness.evidence_provenance.ProvenanceError):
                    robustness.restat_bundle(empty, Path(tmp), 100, 1)
''',
)
old_block = re.compile(
    r"        # Source pins are provenance, and the two kinds of source drift for\n.*?(?=    def test_manuscript_copies_match_generated_evidence)",
    re.DOTALL,
)
new_block = '''        # Schema v2 separates replay-producing code and versioned inputs from
        # later analysis code. Replay dependencies are enforced byte-for-byte;
        # analysis provenance belongs to the restatement block.
        replay = manifest["replay_provenance"]
        implementation = replay["implementation_files"]
        self.assertIn("tests/kalman_ou_iii/kalman_ou_iii-sim.cpp", implementation)
        self.assertIn("src/kalman_ou_iii/Kalman3D_Wave_OU_III.h", implementation)
        self.assertIn("tests/kalman_ou_iii/Makefile", implementation)
        for name, metadata in implementation.items():
            path = REPO_ROOT / name
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(), metadata["sha256"], name
            )
            self.assertEqual(path.stat().st_size, metadata["size_bytes"], name)
        for name, metadata in replay["input_files"].items():
            path = REPO_ROOT / name
            if not path.is_file():
                continue
            self.assertEqual(
                hashlib.sha256(path.read_bytes()).hexdigest(), metadata["sha256"], name
            )
            self.assertEqual(path.stat().st_size, metadata["size_bytes"], name)
        analysis = manifest["restatement"]["analysis_pipeline_files"]
        self.assertIn("tools/ou_validation.py", analysis)
        self.assertIn("tools/ou_robustness.py", analysis)

'''
rtext, count = old_block.subn(new_block, rtext, count=1)
if count != 1:
    raise SystemExit("robustness legacy source-pin block not found")
robust_path.write_text(rtext, encoding="utf-8")

# Canonical raw CSV is immutable and hash-pinned. JSON `raw_runs` historically
# omits fields whose CSV cell is blank, so reconcile semantically: JSON may omit
# only blank canonical cells and may not add or contradict any CSV field.
prov_path = ROOT / "tools/ou_evidence_provenance.py"
ptext = prov_path.read_text(encoding="utf-8")
start = ptext.index("def _assert_bundle_rows_match_raw_csv(")
end = ptext.index("\n\ndef begin_restatement", start)
replacement = '''def _scalar_matches_csv(value: Any, csv_value: str) -> bool:
    if value is None:
        return csv_value == ""
    if isinstance(value, bool):
        return csv_value.lower() in ({"true", "1"} if value else {"false", "0"})
    if isinstance(value, (int, float)) and csv_value != "":
        try:
            lhs = float(value)
            rhs = float(csv_value)
            if lhs != lhs and rhs != rhs:  # NaN
                return True
            return lhs == rhs
        except ValueError:
            pass
    return str(value) == csv_value


def _assert_bundle_rows_match_raw_csv(source_bundle: Path, raw_csv: Path) -> None:
    bundle = json.loads(source_bundle.read_text(encoding="utf-8"))
    rows = bundle.get("raw_runs", [])
    with raw_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        csv_rows = list(reader)
    if len(rows) != len(csv_rows):
        raise ProvenanceError("bundle raw_runs count differs from canonical raw replay CSV")
    expected_fields = set(fields)
    for index, (row, csv_row) in enumerate(zip(rows, csv_rows)):
        extra = set(row) - expected_fields
        if extra:
            raise ProvenanceError(
                f"bundle raw_runs has fields outside canonical raw CSV at row {index}: {sorted(extra)}"
            )
        for field in fields:
            if field not in row:
                if csv_row[field] != "":
                    raise ProvenanceError(
                        f"bundle raw_runs omits non-empty canonical field at row {index}, field {field}"
                    )
                continue
            if not _scalar_matches_csv(row[field], csv_row[field]):
                raise ProvenanceError(
                    f"bundle raw_runs contradicts canonical raw replay CSV at row {index}, field {field}"
                )
'''
ptext = ptext[:start] + replacement + ptext[end:]
prov_path.write_text(ptext, encoding="utf-8")

print("applied acceptance repairs")
