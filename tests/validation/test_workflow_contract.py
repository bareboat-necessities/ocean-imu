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

    def test_full_replay_is_gated_by_broad_content_fingerprint(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        fingerprint = workflow.index("  fingerprint:")
        regenerate = workflow.index("  regenerate:")
        self.assertLess(fingerprint, regenerate)

        gate = workflow[fingerprint:regenerate]
        self.assertIn("tools/ou_replay_fingerprint.py", gate)
        self.assertIn("sim-data-files.zip", gate)
        self.assertIn("reports/results/ou_replay_fingerprint.json", gate)
        self.assertIn("replay_required=false", gate)
        self.assertIn("replay_required=true", gate)
        self.assertIn("make -C tests/validation test", gate)

        regen_header = workflow[regenerate:workflow.index("    runs-on:", regenerate)]
        self.assertIn("needs: fingerprint", regen_header)
        self.assertIn(
            "needs.fingerprint.outputs.replay_required == 'true'", regen_header
        )

    def test_regenerated_evidence_records_the_replay_fingerprint(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        record = workflow.index("- name: Record replay fingerprint")
        check = workflow.index("- name: Check the manuscript against the regenerated evidence")
        commit = workflow.index("- name: Commit the regenerated evidence")
        self.assertLess(record, check)
        self.assertLess(check, commit)

        stage = workflow[record:commit]
        self.assertIn("tools/ou_replay_fingerprint.py", stage)
        self.assertIn("--write reports/results/ou_replay_fingerprint.json", stage)

        commit_stage = workflow[commit:]
        self.assertIn("reports/results/ou_replay_fingerprint.json", commit_stage)

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
