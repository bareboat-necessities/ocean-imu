#!/usr/bin/env python3
"""Plot retained sigma time series; the reported statistics use full-rate rows."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

root = Path(__file__).resolve().parent
for folder, name, configs in [
    ('diagnostics', 'sigma-averaging', [('baseline_cap35', '4 periods · 35 s cap'),
                                       ('K16_cap180', '16 periods · 180 s cap'),
                                       ('K32_cap180', '32 periods · 180 s cap')]),
    ('diagnostics-decay', 'sigma-decay', [('baseline_cap35', '4 periods · 1 s decay'),
                                        ('decay5_only', '4 periods · 5 s decay'),
                                        ('K16_decay5', '16 periods · 5 s decay')]),
]:
    frame = pd.read_csv(root / folder / 'series-1hz.csv.gz')
    frame = frame[(frame.family == 'OU-II') & frame.input.str.contains('jonswap_H0.270')
                  & (frame.time >= 300)]
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True, layout='constrained')
    for config, label in configs:
        rows = frame[frame.config == config]
        axes[0].plot(rows.time, rows.accel_var_tuner, label=label, linewidth=1)
        axes[1].plot(rows.time, rows.sigma_a_applied, label=label, linewidth=1)
    axes[0].set_ylabel('Variance estimate (m²/s⁴)')
    axes[1].set_ylabel('Applied sigma (m/s²)')
    axes[1].set_xlabel('Time (s)')
    axes[0].legend(loc='upper right', ncols=3, fontsize=8)
    for axis in axes:
        axis.grid(alpha=.2)
    fig.suptitle('OU-II · low-wave JONSWAP RAO record · default sensor draw')
    fig.savefig(root / f'{name}.png', dpi=150)
    fig.savefig(root / f'{name}.svg')
