"""Study diagnostics survive output-directory replacement and launch failures."""
import os
from fnmatch import fnmatchcase
from pathlib import Path
import re
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
from simulation_studies import prepare_output, run_commands


class StudyTriggerTests(unittest.TestCase):
    def test_shared_runtime_and_report_dependencies_regenerate_all_studies(self):
        root = Path(__file__).resolve().parents[2]
        workflow = (root/'.github/workflows/simulation-studies.yml').read_text()
        push = workflow.split('  push:\n', 1)[1].split('\npermissions:', 1)[0]
        paths = re.findall(r'^      - "([^"\n]+)"$', push, re.MULTILINE)
        self.assertTrue(paths)
        for path in (
            'src/tuner/GyroMagCalibrationImpl.h',
            'src/ahrs/FrameConversions.h',
            'src/kalman_common/MarineTurnModel.h',
            'src/lie/TangentBiasGroup.h',
            'tests/kalman_tfg_shared_modes/kalman_tfg_shared_modes-sim.cpp',
            'tools/ou_validation.py',
            'tools/ou_robustness.py',
            'plots/vector_export.py',
            'Makefile',
        ):
            with self.subTest(path=path):
                self.assertTrue(any(fnmatchcase(path, pattern) for pattern in paths))
        # Publishing completed studies must not restart all expensive replays.
        for case in ('following', 'head'):
            for path in (f'reports/results/sim_v1.4.0/{case}/studies/ou_validation/case-receipt.json',
                         f'reports/results/sim_v1.4.0/{case}/studies/ou_validation/ou_validation.json'):
                with self.subTest(path=path):
                    self.assertFalse(any(fnmatchcase(path, pattern) for pattern in paths))


class StudyDiagnosticsTests(unittest.TestCase):
    def test_protocol_bound_replays_discard_committed_results_and_cache(self):
        for study, protocol in (
                ('tfg_shared_modes_dataset', 'dataset_protocol.json'),
                ('carma_shared_modes_validation', 'protocol.json')):
            with self.subTest(study=study), tempfile.TemporaryDirectory() as tmp:
                out = Path(tmp)/'study'
                cache = out/'dataset_cache'
                cache.mkdir(parents=True)
                (out/protocol).write_text('{"binaries":"old"}\n')
                (cache/'partial.json').write_text('{"stale":true}\n')
                prepare_output(study, out)
                self.assertTrue(out.is_dir())
                self.assertEqual(list(out.iterdir()), [])

    def test_other_study_output_is_not_reset_by_dataset_preparation(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)/'study'
            out.mkdir()
            marker = out/'metrics.json'
            marker.write_text('{}\n')
            prepare_output('ou_validation', out)
            self.assertTrue(marker.exists())

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
