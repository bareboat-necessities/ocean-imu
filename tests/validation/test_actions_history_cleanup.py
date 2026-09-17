import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / ".github/scripts/actions_history_cleanup.py"
WORKFLOW = ROOT / ".github/workflows/actions-history-cleanup.yml"


class TestActionsHistoryCleanupSource(unittest.TestCase):
    def test_workflow_uses_sharded_cleanup_script(self):
        workflow = WORKFLOW.read_text()
        source = SCRIPT.read_text()
        self.assertIn('params["created"] = created', source)
        self.assertIn("if count <= 1000:", source)
        self.assertIn("span > dt.timedelta(days=1)", source)
        self.assertIn("actions_history_cleanup.py", workflow)
        self.assertNotIn("visible in capped search", workflow)

    def test_cancelled_cleanup_includes_pull_request_runs(self):
        source = SCRIPT.read_text()
        # The 11k+ backlog is predominantly superseded PR checks. Passing
        # exclude_pull_requests=true silently omitted that population.
        self.assertNotIn("exclude_pull_requests", source)
        self.assertIn('enumerate_status(conclusion)', source)
        self.assertIn('(\"cancelled\", \"skipped\")', source)

    def test_large_backlog_is_deleted_in_bounded_batches(self):
        workflow = WORKFLOW.read_text()
        source = SCRIPT.read_text()
        self.assertIn('DELETE_BATCH_SIZE", "250"', source)
        self.assertIn("for offset in range(0, len(ids), BATCH):", source)
        self.assertIn("MAX_DELETIONS: ${{ inputs.max_deletions || '20000' }}", workflow)
        self.assertIn("DELETE_BATCH_SIZE: ${{ inputs.delete_batch_size || '250' }}", workflow)
        self.assertIn("timeout-minutes: 180", workflow)

    def test_cleanup_fix_runs_on_its_own_pr_change(self):
        workflow = WORKFLOW.read_text()
        self.assertIn("pull_request:", workflow)
        self.assertIn("'.github/scripts/actions_history_cleanup.py'", workflow)
        self.assertIn("actions: write", workflow)


if __name__ == "__main__":
    unittest.main()
