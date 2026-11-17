"""Independently check the physical histories used for steering attribution."""
from pathlib import Path
import importlib.util
import sys
import unittest
import tempfile
import hashlib
from unittest.mock import patch

import numpy as np
from scipy.spatial.transform import Rotation

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tools'))
from heading_kinematics_diagnosis import control, harmonic_amplitude, indices, cached_trace_valid
import heading_turn_validation
from underway_reference import track


def record(case):
    t = np.arange(12000)/200
    p, v, a, yaw = track(t, case)
    wave = .08*np.sin(2*t)
    p[:, 2] = wave
    v[:, 2] = .16*np.cos(2*t)
    a[:, 2] = -.32*np.sin(2*t)
    roll, pitch = .03*np.sin(t), .02*np.cos(t)
    rd, pd = .03*np.cos(t), -.02*np.sin(t)
    yd = np.deg2rad(3)*2*np.pi/20*np.cos(2*np.pi*t/20)
    rotation = Rotation.from_euler('ZYX', np.column_stack((np.deg2rad(yaw), pitch, roll))).inv()
    gyro = np.column_stack((rd-yd*np.sin(pitch),
                            pd*np.cos(roll)+yd*np.sin(roll)*np.cos(pitch),
                            -pd*np.sin(roll)+yd*np.cos(roll)*np.cos(pitch)))
    q = rotation.as_quat()[:, [3, 0, 1, 2]]
    columns = ['time']+[prefix+x for prefix in ('disp_', 'vel_', 'acc_', 'acc_b', 'gyro_') for x in 'xyz']
    columns += ['roll_deg', 'pitch_deg', 'yaw_deg']+['q_wb_zu_'+x for x in 'wxyz']
    data = np.column_stack((t, p, v, a, rotation.apply(a+[0, 0, 9.80665]), gyro,
                            np.rad2deg(roll), np.rad2deg(pitch), yaw, q))
    return columns, data


class HeadingControlTests(unittest.TestCase):
    def test_each_counterfactual_is_a_complete_physical_history(self):
        for case in ('following', 'head'):
            columns, original = record(case)
            original_copy = original.copy()
            for fixed in (False, True):
                with self.subTest(case=case, fixed_yaw=fixed):
                    data = control(columns, original, case, fixed)
                    np.testing.assert_array_equal(original, original_copy)
                    p, v, a = [data[:, indices(columns, prefix)] for prefix in ('disp_', 'vel_', 'acc_')]
                    np.testing.assert_allclose(p[:, :2], 0, atol=1e-13)
                    np.testing.assert_allclose(v[:, :2], 0, atol=1e-13)
                    np.testing.assert_allclose(a[:, :2], 0, atol=1e-13)
                    np.testing.assert_array_equal(p[:, 2], original[:, columns.index('disp_z')])
                    np.testing.assert_allclose(np.gradient(p, .005, axis=0)[2:-2], v[2:-2], atol=3e-6)
                    np.testing.assert_allclose(np.gradient(v, .005, axis=0)[2:-2], a[2:-2], atol=6e-6)
                    q = data[:, [columns.index('q_wb_zu_'+x) for x in 'xyzw']]
                    rotation = Rotation.from_quat(q)
                    np.testing.assert_allclose(rotation.apply(a+[0, 0, 9.80665]),
                                               data[:, indices(columns, 'acc_b')], atol=1e-12)
                    # A central relative rotation measures instantaneous body rate.
                    gyro = (rotation[:-2]*rotation[2:].inv()).as_rotvec()/.01
                    np.testing.assert_allclose(gyro, data[1:-1, indices(columns, 'gyro_')], atol=2e-6)
                    if not fixed:
                        for prefix in ('gyro_', 'q_wb_zu_'):
                            selected = [i for i, name in enumerate(columns) if name.startswith(prefix)]
                            np.testing.assert_array_equal(data[:, selected], original[:, selected])

    def test_receipt_cannot_reuse_a_missing_or_changed_trace(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'w3d-test.csv'
            content = b'time,disp_est_y\n0,1\n'
            receipt = {'trace_sha256': hashlib.sha256(content).hexdigest()}
            self.assertFalse(cached_trace_valid(Path(folder), receipt))
            path.write_bytes(content)
            self.assertTrue(cached_trace_valid(Path(folder), receipt))
            path.write_bytes(content + b'1,2\n')
            self.assertFalse(cached_trace_valid(Path(folder), receipt))

    def test_steering_signature(self):
        t = np.arange(12000)/200
        p, _, a, _ = track(t)
        self.assertAlmostEqual(harmonic_amplitude(t, p[:, 1]), .428557, places=5)
        self.assertAlmostEqual(np.max(np.abs(a[:, 1])), .04231145, places=6)

    def test_comparison_rejects_refreshed_baseline_metrics(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            content = b'{"records": {}}\n'
            expected = {}
            for case in ('following', 'head'):
                relative = Path('reports/results/heading_kinematics/pre_fix')/case/'measurements.json'
                path = root/relative
                path.parent.mkdir(parents=True)
                path.write_bytes(content)
                expected[case] = hashlib.sha256(content).hexdigest()
            with (patch.object(heading_turn_validation, 'ROOT', root),
                  patch.object(heading_turn_validation, 'BASELINE_SHA256', expected)):
                self.assertEqual(set(heading_turn_validation.load_baseline()), {'following', 'head'})
                path.write_text(path.read_text()+'\n')
                with self.assertRaisesRegex(ValueError, 'baseline changed'):
                    heading_turn_validation.load_baseline()

    def test_harmonic_baseline_cannot_be_replaced_by_fresh_runs(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            path = root/'reports/results/heading_kinematics/measurements.json'
            path.parent.mkdir(parents=True)
            content = b'{"runs": []}\n'
            path.write_bytes(content)
            with (patch.object(heading_turn_validation, 'ROOT', root),
                  patch.object(heading_turn_validation, 'ATTRIBUTION_SHA256',
                               {'measurements': hashlib.sha256(content).hexdigest()})):
                self.assertEqual(heading_turn_validation.load_attribution('measurements'), {'runs': []})
                path.write_bytes(content+b'\n')
                with self.assertRaisesRegex(ValueError, 'attribution baseline changed'):
                    heading_turn_validation.load_attribution('measurements')

    def test_tfg_bias_chart_distinguishes_estimates_from_injected_truth(self):
        import pandas as pd
        path = Path(__file__).resolve().parents[2]/'plots/kalman_tfg/kalman_tfg-plots.py'
        spec = importlib.util.spec_from_file_location('tfg_plots', path)
        plots = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(plots)
        frame = pd.DataFrame({'time': [0., 1., 2.]})
        for i, axis in enumerate('xyz'):
            frame['acc_bias_'+axis] = [i, i+.1, i+.2]
            frame['acc_bias_est_'+axis] = [i+1, i+1.2, i+1.3]
        with patch.object(plots, 'finalize') as saved:
            plots.plot_bias(frame, 'unused', 'test')
        fig = saved.call_args.args[0]
        self.assertEqual(len(fig.axes), 3)
        for ax, axis in zip(fig.axes, 'xyz'):
            self.assertEqual([line.get_label() for line in ax.lines], ['true', 'estimate'])
            np.testing.assert_array_equal(ax.lines[0].get_ydata(), frame['acc_bias_'+axis])
            np.testing.assert_array_equal(ax.lines[1].get_ydata(), frame['acc_bias_est_'+axis])
            self.assertIn('Body', ax.get_ylabel())
        plots.plt.close(fig)


if __name__ == '__main__':
    unittest.main()
