"""Counterfactual construction and scoring of tools/moving_vessel_diagnosis.py."""
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'tools'))
import moving_vessel_diagnosis as mvd  # noqa: E402


def synthetic_record(n=4000, period=3.0):
    """Self-consistent wave record: heave/surge sinusoids, roll/pitch rocking."""
    t = mvd.sample_times(n)
    w = 2 * np.pi / period
    frame = pd.DataFrame({'time': np.round(t, 2)})
    for axis, amp, phase in (('x', 0.05, 0.3), ('y', 0.02, 1.1), ('z', 0.10, 0.0)):
        frame[f'disp_{axis}'] = amp * np.sin(w * t + phase)
        frame[f'vel_{axis}'] = amp * w * np.cos(w * t + phase)
        frame[f'acc_{axis}'] = -amp * w * w * np.sin(w * t + phase)
    roll, pitch = 0.02 * np.sin(w * t), 0.015 * np.cos(w * t)
    # q_wb is world -> body: conjugate of body -> world Rx(roll) then Ry(pitch).
    qx = np.column_stack((np.cos(roll / 2), np.sin(roll / 2), 0 * t, 0 * t))
    qy = np.column_stack((np.cos(pitch / 2), 0 * t, np.sin(pitch / 2), 0 * t))
    q_bw = mvd.qmul(qy, qx)
    q_wb = q_bw * np.array([1, -1, -1, -1])
    dq = np.gradient(q_wb, mvd.DT, axis=0)
    gyro = 2 * mvd.qmul(q_wb, dq * np.array([1, -1, -1, -1]))[:, 1:]
    mvd.set_body_fields(frame, q_wb, gyro)
    return frame


def body_consistency(frame):
    residual = mvd.body_residuals(frame)
    return residual['specific_force_residual_mps2'], residual['angular_rate_residual_radps']


class TrackTest(unittest.TestCase):
    def test_primitive_and_derivatives_match(self):
        t = mvd.sample_times(24000 * 2)
        p, v, a, heading, rate = mvd.track(t)
        self.assertLess(np.abs(np.gradient(p, mvd.DT, axis=0) - v)[5:-5].max(), 1e-6)
        self.assertLess(np.abs(np.gradient(v, mvd.DT, axis=0) - a)[5:-5].max(), 1e-7)
        np.testing.assert_allclose(np.linalg.norm(v, axis=1), mvd.SPEED_MPS)
        self.assertAlmostEqual(np.degrees(np.abs(heading).max()), 3.0, places=6)
        np.testing.assert_allclose(np.gradient(heading, mvd.DT)[5:-5], rate[5:-5], atol=1e-9)

    def test_rounded_csv_time_misplaces_track_by_speed_times_rounding(self):
        n = 240000
        frame = pd.DataFrame({'time': [float(f'{x:.6g}') for x in mvd.sample_times(n)]})
        artifact = mvd.time_stamp_artifact(frame)
        # Past 1000 s, every other sample is rounded by 5 ms.
        expected = mvd.SPEED_MPS * 0.005 * np.sqrt(0.5 * 200 / 900)
        self.assertAlmostEqual(artifact, expected, delta=0.03 * expected)


class CounterfactualTest(unittest.TestCase):
    def test_retrack_is_consistent_and_invertible(self):
        base = synthetic_record()
        moving = mvd.retrack(base, +1)
        for frame in (base, moving):
            f_err, w_err = body_consistency(frame)
            self.assertLess(f_err, 1e-9)
            self.assertLess(w_err, 1e-6)
        np.testing.assert_allclose(moving['yaw_deg'], np.degrees(mvd.track(mvd.sample_times(len(base)))[3]),
                                   atol=1e-6)
        back = mvd.retrack(moving, -1)
        columns = [c for c in base.columns if c != 'time']
        np.testing.assert_allclose(back[columns].to_numpy(), base[columns].to_numpy(), atol=1e-9)

    def test_stretch_keeps_heave_and_scales_derivatives(self):
        base = synthetic_record(n=8000)
        slow = mvd.stretch(base, 1.5)
        f_err, w_err = body_consistency(slow)
        self.assertLess(f_err, 1e-9)
        self.assertLess(w_err, 1e-5)
        tail = slice(-4000, None)
        self.assertAlmostEqual(slow['disp_z'][tail].abs().max(), 0.10, delta=1e-3)
        self.assertAlmostEqual(slow['acc_z'][tail].abs().max(), 0.10 * (2 * np.pi / 4.5) ** 2, delta=1e-4)


class StillnessReplicaTest(unittest.TestCase):
    def test_threshold_is_absolute_acceleration(self):
        t = mvd.sample_times(mvd.WINDOW_SAMPLES + 2000)
        w = 2 * np.pi / 5.0
        # Same heave amplitude; doubling the period divides acceleration by 4.
        fast = mvd.stillness(0.3 * w * w * np.sin(w * t))
        slow = mvd.stillness(0.3 * (w / 2) ** 2 * np.sin(w / 2 * t))
        self.assertLess(fast['still_fraction'], 0.05)
        self.assertGreater(slow['still_fraction'], 0.95)


if __name__ == '__main__':
    unittest.main()
