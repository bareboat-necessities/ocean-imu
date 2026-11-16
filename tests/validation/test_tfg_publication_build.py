"""Article builds reuse valid evidence and fail closed on replay/check errors."""
import os
from pathlib import Path
import subprocess
import tempfile
import textwrap
import unittest

ROOT = Path(__file__).resolve().parents[2]


class TFGPublicationBuildTests(unittest.TestCase):
    def run_step(self, case, initial=0, replay=0, final=0):
        workflow = (ROOT / '.github/workflows/build.yml').read_text()
        step = workflow.split('      - name: Prepare case-specific full TFG publication evidence\n', 1)[1]
        script = textwrap.dedent(step.split('        run: |\n', 1)[1].split('\n      - name:', 1)[0])
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stub = root / 'python3'
            stub.write_text('''#!/bin/bash
printf '%s %s\\n' "$OCEAN_IMU_SIM_CASE" "$*" >> "$CALL_LOG"
case "$*" in
  *--publish*) exit "$FINAL_STATUS" ;;
  *--mode*) exit "$REPLAY_STATUS" ;;
  *) exit "$INITIAL_STATUS" ;;
esac
''')
            stub.chmod(0o755)
            env = {**os.environ, 'PATH': tmp + ':' + os.environ['PATH'],
                   'OCEAN_IMU_SIM_CASE': case, 'CALL_LOG': str(root / 'calls'),
                   'INITIAL_STATUS': str(initial), 'REPLAY_STATUS': str(replay),
                   'FINAL_STATUS': str(final)}
            result = subprocess.run(['bash', '-e', '-o', 'pipefail', '-c', script], env=env)
            return result.returncode, (root / 'calls').read_text().splitlines()

    def test_current_bundle_is_reused_for_each_case(self):
        for case in ('following', 'head'):
            status, calls = self.run_step(case)
            self.assertEqual(status, 0)
            self.assertEqual(len(calls), 2)
            self.assertTrue(all(c.startswith(case + ' ') for c in calls))
            self.assertTrue(all('--check --require-ci' in c for c in calls))

    def test_invalid_bundle_requires_full_replay_then_check(self):
        status, calls = self.run_step('head', initial=1)
        self.assertEqual(status, 0)
        self.assertEqual(len(calls), 3)
        self.assertIn('--mode full --jobs 2 --data-dir plots/kalman_ou_iii', calls[1])
        self.assertIn('--check --require-ci --publish', calls[2])

    def test_failed_replay_prevents_publication(self):
        status, calls = self.run_step('head', initial=1, replay=7)
        self.assertEqual(status, 7)
        self.assertEqual(len(calls), 2)
        self.assertNotIn('--publish', ' '.join(calls))

    def test_final_provenance_failure_fails_build(self):
        status, calls = self.run_step('following', initial=1, final=9)
        self.assertEqual(status, 9)
        self.assertEqual(len(calls), 3)
