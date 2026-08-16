#!/usr/bin/env python3
from pathlib import Path

path = Path(__file__).resolve().parents[1] / "tests/validation/test_ou_robustness.py"
text = path.read_text(encoding="utf-8")
old = '''    def test_a_restat_records_the_sources_that_moved_under_it(self):
        """The obsolete warn-and-repin path must not survive schema v2."""
        self.assertFalse(hasattr(robustness, "_restated_source_files"))
        manifest = json.loads(
            (self.RESULTS / "ou_robustness_manifest.json").read_text(encoding="utf-8")
        )
        self.assertNotIn("sources_moved_since_rows", manifest)
        self.assertIn("replay_provenance", manifest)
        self.assertIn("restatement", manifest)
'''
new = '''    def test_committed_manifest_has_no_obsolete_moved_source_state(self):
        """Schema v2 has one immutable replay model; restatement exists only after one occurs."""
        self.assertFalse(hasattr(robustness, "_restated_source_files"))
        manifest = json.loads(
            (self.RESULTS / "ou_robustness_manifest.json").read_text(encoding="utf-8")
        )
        self.assertNotIn("sources_moved_since_rows", manifest)
        replay = manifest["replay_provenance"]
        self.assertEqual(manifest["git_commit"], replay["git_commit"])
        if "restatement" in manifest:
            self.assertIn("analysis_pipeline_files", manifest["restatement"])
        else:
            # A freshly generated bundle has no reason to manufacture a
            # restatement event merely to satisfy the schema.
            self.assertIn("analysis_pipeline_files", manifest)
'''
if old not in text:
    raise SystemExit("stale robustness provenance test anchor not found")
path.write_text(text.replace(old, new, 1), encoding="utf-8")
print("fixed fresh-generation/restatement test semantics")
