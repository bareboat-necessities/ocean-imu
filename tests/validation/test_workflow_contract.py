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

    def test_push_retry_revalidates_after_rebase(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        start = workflow.index("for attempt in 1 2 3 4; do")
        end = workflow.index("\n          done", start)
        loop = workflow[start:end]

        rebase = loop.index("git pull --rebase")
        validate = loop.index("make -C tests/validation test")
        self.assertLess(rebase, validate)

    def test_failure_message_cannot_run_after_a_push(self):
        workflow = WORKFLOW.read_text(encoding="utf-8")
        self.assertNotIn("The regenerated bundle is committed, but", workflow)
        self.assertNotIn("Fail if the manuscript no longer matches", workflow)


if __name__ == "__main__":
    unittest.main()
