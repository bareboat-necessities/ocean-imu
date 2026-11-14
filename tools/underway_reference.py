"""Declared v1.3.0 5 kn track; truth-side decomposition, never filter aid."""
import numpy as np
from scipy.special import jv

SPEED_MPS = 5 * 1852 / 3600
AMPLITUDE_RAD = np.deg2rad(3)
FREQUENCY_RADPS = 2 * np.pi / 120


def track(times):
    times = np.asarray(times, dtype=float)
    phase = FREQUENCY_RADPS * times
    x = jv(0, AMPLITUDE_RAD) * times
    y = np.zeros_like(times)
    for n in range(1, 13):
        coefficient = 2 * jv(n, AMPLITUDE_RAD) / (n * FREQUENCY_RADPS)
        if n % 2 == 0:
            x = x + coefficient * np.sin(n * phase)
        else:
            y = y + coefficient * (1 - np.cos(n * phase))
    heading = AMPLITUDE_RAD * np.sin(phase)
    rate = AMPLITUDE_RAD * FREQUENCY_RADPS * np.cos(phase)
    zero = np.zeros_like(times)
    position = SPEED_MPS * np.column_stack((x, y, zero))
    velocity = SPEED_MPS * np.column_stack((np.cos(heading), np.sin(heading), zero))
    acceleration = SPEED_MPS * np.column_stack((-rate * np.sin(heading), rate * np.cos(heading), zero))
    return position, velocity, acceleration, np.rad2deg(heading)


def change_track(columns, data, sign):
    """Remove/add one deterministic p/v/a/yaw chain for surrogate generation."""
    result = data.copy()
    position, velocity, acceleration, heading = track(result[:, columns.index('time')])
    for prefix, values in zip(('disp', 'vel', 'acc'), (position, velocity, acceleration)):
        indices = [columns.index(prefix + '_' + a) for a in 'xyz']
        result[:, indices] += sign * values
    result[:, columns.index('yaw_deg')] += sign * heading
    return result
