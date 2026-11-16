"""Study diagnostics survive output-directory replacement and launch failures."""
import os
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
from simulation_studies import run_commands


class StudyDiagnosticsTests(unittest.TestCase):
    def test_successful_study_can_replace_its_output_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)/'study'
            out.mkdir()
            (out/'stale.txt').write_text('previous results')
            code = ('import pathlib,shutil,sys; p=pathlib.Path(sys.argv[1]); '
                    'print("before replacement",flush=True); shutil.rmtree(p); p.mkdir(); '
                    '(p/"metrics.csv").write_text("rms\\n1\\n"); print("after replacement")')
            result = run_commands([[sys.executable, '-c', code, str(out)]], os.environ.copy(), out)
            self.assertEqual(result.returncode, 0)
            self.assertEqual((out/'metrics.csv').read_text(), 'rms\n1\n')
            self.assertEqual((out/'study.log').read_text(), 'before replacement\nafter replacement\n')
            self.assertEqual((out.parent/'study.log').read_bytes(), (out/'study.log').read_bytes())
            self.assertFalse((out/'stale.txt').exists())

    def test_failure_preserves_log_and_does_not_launch_postprocessing(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)/'study'
            out.mkdir()
            code = ('import shutil,sys; shutil.rmtree(sys.argv[1]); '
                    'print("failed study",file=sys.stderr); sys.exit(3)')
            commands = [[sys.executable, '-c', code, str(out)],
                        [sys.executable, '-c', 'raise RuntimeError("must not run")']]
            result = run_commands(commands, os.environ.copy(), out)
            self.assertEqual(result.returncode, 3)
            self.assertEqual((out/'study.log').read_text(), 'failed study\n')

    def test_launch_failure_is_recorded(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)/'study'
            result = run_commands([[str(out/'missing-executable')]], os.environ.copy(), out)
            self.assertEqual(result.returncode, 127)
            self.assertIn('Unable to launch', (out/'study.log').read_text())


if __name__ == '__main__':
    unittest.main()
