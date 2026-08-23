from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW = REPO_ROOT / ".github" / "workflows" / "ou-validation.yml"


class WorkflowContractTests(unittest.TestCase):
    def test_full_regeneration_validates_before_commit(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        check_marker = "- name: Check the manuscript against the regenerated evidence"
        commit_marker = "- name: Commit the regenerated evidence"
        check = workflow.index(check_marker)
        commit = workflow.index(commit_marker)

        self.assertLess(check, commit)
        gate = workflow[check:commit]
        self.assertIn("make -C tests/validation test", gate)
        self.assertNotIn("continue-on-error", gate)

    def test_validation_publication_is_aligned_before_bundle_upload(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        combine = workflow.index("- name: Combine paired validation shards")
        align = workflow.index(
            "- name: Align validation publication wording with the article"
        )
        upload = workflow.index("- name: Upload regenerated bundle")

        self.assertLess(combine, align)
        self.assertLess(align, upload)
        stage = workflow[align:upload]
        self.assertIn("tools/ou_publication_sync.py", stage)
        self.assertIn("reports/results/ou_validation", stage)

    def test_full_replay_is_gated_by_replay_and_results_fingerprints(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        fingerprint = workflow.index("  fingerprint:")
        regenerate = workflow.index("  regenerate:")
        self.assertLess(fingerprint, regenerate)

        gate = workflow[fingerprint:regenerate]
        self.assertIn("tools/ou_replay_fingerprint.py", gate)
        self.assertIn("sim-data-files.zip", gate)
        self.assertIn("reports/ou_evidence_fingerprint.json", gate)
        self.assertIn("replay_required=false", gate)
        self.assertIn("replay_required=true", gate)
        self.assertIn("complete results tree", gate)
        self.assertIn("make -C tests/validation test", gate)

        regen_header = workflow[regenerate:workflow.index("    runs-on:", regenerate)]
        self.assertIn("needs: fingerprint", regen_header)
        self.assertIn(
            "needs.fingerprint.outputs.replay_required == 'true'", regen_header
        )

    def test_every_full_replay_uses_the_fingerprinted_zip_bytes(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        full = workflow[workflow.index("  fingerprint:"):]
        self.assertIn("- name: Preserve fingerprinted simulation archive", full)
        self.assertIn("name: ou-fingerprinted-simulation-data", full)
        self.assertGreaterEqual(
            full.count("name: ou-fingerprinted-simulation-data"),
            4,
        )
        regenerate = full[full.index("  regenerate:"):]
        self.assertNotIn("gh release download v1.1.3", regenerate)

    def test_regenerated_evidence_hashes_the_final_validated_tree(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        check = workflow.index("- name: Check the manuscript against the regenerated evidence")
        record = workflow.index("- name: Record replay and results fingerprints")
        verify = workflow.index("- name: Verify final evidence fingerprints")
        commit = workflow.index("- name: Commit the regenerated evidence")
        self.assertLess(check, record)
        self.assertLess(record, verify)
        self.assertLess(verify, commit)

        stage = workflow[record:commit]
        self.assertIn("tools/ou_replay_fingerprint.py", stage)
        self.assertIn("--write reports/ou_evidence_fingerprint.json", stage)
        self.assertIn("--check reports/ou_evidence_fingerprint.json", stage)

        commit_stage = workflow[commit:]
        self.assertIn("reports/ou_evidence_fingerprint.json", commit_stage)
        self.assertIn("reports/results/ou_validation", commit_stage)
        self.assertIn("reports/results/ou_robustness", commit_stage)

    def test_results_tree_changes_are_in_workflow_trigger(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        self.assertIn('- "reports/results/**"', workflow)

    def test_push_retry_revalidates_after_rebase(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        start = workflow.index("for attempt in 1 2 3 4; do")
        end = workflow.index("\n          done", start)
        loop = workflow[start:end]

        rebase = loop.index("git pull --rebase")
        fingerprint = loop.index("tools/ou_replay_fingerprint.py")
        validate = loop.index("make -C tests/validation test")
        self.assertLess(rebase, fingerprint)
        self.assertLess(fingerprint, validate)

    def test_failure_message_cannot_run_after_a_push(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        self.assertNotIn("The regenerated bundle is committed, but", workflow)
        self.assertNotIn("Fail if the manuscript no longer matches", workflow)


if __name__ == "__main__":
    unittest.main()
