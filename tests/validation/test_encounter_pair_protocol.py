"""Paired amplitude experiments must request and score the same complete record."""
import os
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
import encounter_scaling_study as study


class EncounterPairProtocolTests(unittest.TestCase):
    def job(self, previous=False):
        return dict(family='csm', previous=previous, previous_binary=None,
                    engine_rpm=None, dataset='moving', seed='7', noise='normal',
                    binary=Path('/test/replay'), input=Path('/test/wave.csv'))

    def output(self, **changes):
        fields = dict(input='wave.csv', samples='180000', start_s='300', window_s='900')
        fields.update({key: '0.1' for key in study.METRICS + (
            'disp_3d_rms_m', 'accel_bias_3d_rms_mps2', 'gyro_bias_3d_rms_radps')})
        fields.update(changes)
        return 'VALIDATION_METRICS ' + ' '.join(f'{k}={v}' for k, v in fields.items()) + '\n'

    def test_same_stream_explicit_record_and_only_previous_disables_mapping(self):
        environments = []
        with patch.dict(os.environ, {'W3D_SEED': '42', 'W3D_CARMA_BAND_NORMALIZATION': '0'}):
            for previous in (True, False):
                with patch.object(study.subprocess, 'run', return_value=subprocess.CompletedProcess(
                        [], 0, self.output(), '')) as run, patch('builtins.print'):
                    row = study.run_one(self.job(previous))
                self.assertEqual(run.call_args.args[0], ['/test/replay', '--input', 'wave.csv'])
                environments.append(run.call_args.kwargs['env'])
                self.assertEqual(row['arm'], 'previous' if previous else 'band_normalized')
        self.assertEqual(environments[0].pop('W3D_CARMA_BAND_NORMALIZATION'), '0')
        self.assertNotIn('W3D_CARMA_BAND_NORMALIZATION', environments[1])
        self.assertEqual(environments[0], environments[1])
        self.assertEqual(environments[1]['W3D_SEED'], '7')

    def test_mismatched_input_incomplete_window_and_nonfinite_scores_fail_closed(self):
        for changes in ({'input': 'other.csv'}, {'samples': '179999'},
                        {'start_s': '0'}, {'disp_z_rms_m': 'nan'}):
            with self.subTest(changes=changes), patch.object(
                    study.subprocess, 'run', return_value=subprocess.CompletedProcess(
                        [], 0, self.output(**changes), '')), self.assertRaises(ValueError):
                study.run_one(self.job())

    def test_failed_gates_remain_visible_in_diagnostic_rows(self):
        error = 'ERROR: yaw RMS above limit'
        with patch.object(study.subprocess, 'run', return_value=subprocess.CompletedProcess(
                [], 1, self.output(), error)), patch('builtins.print'):
            row = study.run_one(self.job())
        self.assertEqual(row['exit_code'], 1)
        self.assertEqual(row['gate_failures'], [error])


if __name__ == '__main__':
    unittest.main()
