#!/usr/bin/env python3
"""Render the shared-mode TFG degradation tables and figure from the scored CSVs."""
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sim_dataset import CASES, ROOT, REFERENCE_CSV_SHA256_BY_CASE, digest
from tfg_shared_modes_degradation import conditions

sys.path.insert(0, str(ROOT / 'plots'))
from svg_determinism import configure_svg, save_svg
from vector_export import save_figure

OUT = ROOT / 'reports/results/tfg_shared_modes/degradation'
ARTICLE = ROOT / 'doc/kalman_tfg_shared_modes'
AXES = {'x-forward': 'Forward', 'y-starboard': 'Starboard', 'z-up': 'Above'}
SEGMENTS = ('low_start', 'rise', 'high_recover', 'high', 'fall', 'low_recover', 'low_return')
METRICS = ['disp_z_rms_m', 'disp_3d_rms_m', 'yaw_rms_deg']


def pooled(values):
    return np.sqrt(np.mean(values ** 2))


def read_case(case):
    directory = OUT / case
    manifest = json.loads((directory / 'manifest.json').read_text())
    for name, expected in manifest['files'].items():
        if digest(directory / name) != expected:
            raise ValueError(f'{case}/{name} differs from its manifest')
    rows = pd.read_csv(directory / 'steady_seas.csv')
    engine = rows[rows['mode'].isin(['baseline', 'engine'])]
    lever = rows[rows['mode'].isin(['unmodeled', 'wrapper', 'exact'])]
    ids = ['record', 'mode', 'rpm', 'axis', 'distance_m']
    expected = {(record, mode, rpm, axis, distance)
                for record in REFERENCE_CSV_SHA256_BY_CASE[case]
                for mode, rpm, axis, distance, _ in conditions()}
    if (rows.duplicated(ids).any() or set(rows[ids].itertuples(index=False, name=None)) != expected
            or manifest['dataset']['case'] != case):
        raise ValueError(f'Incomplete paired design for {case}')
    engines = engine.groupby('rpm')[METRICS].agg(pooled)
    installs = lever.groupby(['distance_m', 'axis', 'mode'])[METRICS].agg(pooled)
    segments = pd.read_csv(directory / 'roundtrip.csv').set_index('segment')
    if tuple(segments.index) != SEGMENTS:
        raise ValueError(f'Incomplete round trip for {case}')
    return engines, installs / engines.loc[0], segments


def table(name, caption, headings, rows):
    columns = 'l' + 'r' * (len(headings) - 1)
    lines = [r'\begin{table}[!ht]\centering\small',
             rf'\caption{{{caption}}}\label{{tab:tcsm-{name}}}',
             r'\setlength{\tabcolsep}{3pt}',
             rf'\begin{{tabular}}{{{columns}}}\toprule',
             ' & '.join(headings) + r' \\ \midrule']
    lines.extend(' & '.join(row) + r' \\' for row in rows)
    lines.extend([r'\bottomrule\end{tabular}\end{table}', ''])
    return '\n'.join(lines)


def main():
    data = {case: read_case(case) for case in CASES}
    f_engine, f_lever, f_segments = data['following']
    head_engine, head_lever, head_segments = data['head']
    phase_names = ['Low', 'Rise', 'Early high', 'High', 'Fall', 'Early low', 'Return']
    engine_rows = [[str(rpm)] + [f'{v:.3f}' for engines in (f_engine, head_engine)
                                 for v in engines.loc[rpm, METRICS]] for rpm in f_engine.index]
    transition_rows = []
    for phase, name in zip(SEGMENTS, phase_names):
        segment = f_segments.loc[phase]
        transition_rows.append([name, f'{segment.start_s:.0f}--{segment.stop_s:.0f}'] +
                               [f'{segments.loc[phase, m]:.3f}' for segments in (f_segments, head_segments)
                                for m in METRICS[:2]])
    lever_rows = []
    for distance, axis in f_lever.index.droplevel('mode').unique():
        lever_rows.append([AXES[axis], f'{distance:.1f}'] +
                          [f'{levers.loc[(distance, axis, mode), "disp_3d_rms_m"]:.3f}'
                           for levers in (f_lever, head_lever) for mode in ('unmodeled', 'wrapper')])
    engine_tex = table('engine', 'Engine sweep. Pooled heave and 3-D displacement RMS (m), and yaw RMS '
                       '(degrees). F and H denote following and head seas.', ['rpm', '$z_F$', '3-D$_F$',
                       'Yaw$_F$', '$z_H$', '3-D$_H$', 'Yaw$_H$'], engine_rows)
    transition_tex = table('roundtrip', 'Continuous round trip. Heave and 3-D displacement RMS (m) '
                           'within each scoring interval.', ['Phase', 'Time (s)', '$z_F$', '3-D$_F$',
                           '$z_H$', '3-D$_H$'], transition_rows)
    lever_tex = table('mount', 'Pooled 3-D displacement RMS divided by the at-CoG RMS. Raw denotes '
                      'an uncorrected IMU and gyro denotes correction in the wrapper.', ['Offset', 'm',
                      'Raw$_F$', 'Gyro$_F$', 'Raw$_H$', 'Gyro$_H$'], lever_rows)
    for tag, tex in [('engine', engine_tex), ('roundtrip', transition_tex), ('mount', lever_tex)]:
        (ARTICLE / f'degradation-{tag}.tex-part').write_text('% Generated from the two case CSV bundles.\n' + tex)

    configure_svg(matplotlib)
    matplotlib.rcParams.update({'pgf.texsystem': 'pdflatex', 'pgf.rcfonts': False, 'font.size': 9})
    fig, axes = plt.subplots(1, 3, figsize=(8.3, 2.85), layout='constrained')
    colors = {'following': '#12667a', 'head': '#b05226'}
    for case, (engines, levers, segments) in data.items():
        color = colors[case]
        label = 'Following' if case == 'following' else 'Head'
        axes[0].plot(engines.index, engines.disp_3d_rms_m, 'o-', color=color, label=label)
        axes[1].plot(phase_names, segments.disp_z_rms_m, 'o-', color=color)
        for i, mode in enumerate(('unmodeled', 'wrapper')):
            values = [levers.loc[(1.0, axis, mode), 'disp_3d_rms_m'] for axis in AXES]
            offset = (-0.15 if case == 'following' else 0.15) + 0.08 * (i - 0.5)
            axes[2].scatter(np.arange(3) + offset, values, color=color, marker='o' if i == 0 else 'x',
                            label=('Uncorrected' if i == 0 else 'Gyro') if case == 'following' else None)
    axes[0].set(xlabel='Engine speed (rpm)', ylabel='3-D RMS (m)', title='Engine vibration',
                xticks=[0, 800, 1600, 2400, 3200])
    axes[1].set(ylabel='Heave RMS (m)', title='Low–high–low seas')
    axes[1].tick_params(axis='x', labelrotation=55)
    axes[2].set(ylabel='3-D RMS / at-CoG RMS', title='IMU offset of 1 m',
                xticks=np.arange(3), xticklabels=list(AXES.values()))
    axes[2].axhline(1, color='0.45', lw=0.75, linestyle='--', zorder=0)
    axes[0].legend(frameon=False, fontsize=8)
    axes[2].legend(frameon=False, fontsize=8)
    for ax in axes:
        ax.grid(axis='y', alpha=0.25)
        ax.spines[['right', 'top']].set_visible(False)
    save_svg(fig, OUT / 'underway_degradation.svg')
    save_figure(fig, OUT / 'underway_degradation.pgf')
    plt.close(fig)


if __name__ == '__main__':
    main()
