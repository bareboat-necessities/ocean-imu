import unittest
from pathlib import Path

P=Path(__file__).resolve().parents[2]/'.github/scripts/actions_history_cleanup.py'

class TestActionsHistoryCleanupSource(unittest.TestCase):
    def test_workflow_uses_sharded_cleanup_script(self):
        workflow=(Path(__file__).resolve().parents[2]/'.github/workflows/actions-history-cleanup.yml').read_text()
        source=P.read_text()
        self.assertIn('created_at',source)
        self.assertIn("if n<=1000:return page_all",source)
        self.assertIn("span>dt.timedelta(days=1)",source)
        self.assertIn('actions_history_cleanup.py',workflow)
        self.assertNotIn('visible in capped search',workflow)

if __name__=='__main__': unittest.main()
