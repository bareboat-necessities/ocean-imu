#!/usr/bin/env python3
"""Paired statistical report and bounded publication figures for TAN research."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'reports/results/tangent_shared_modes'
DOC = ROOT/'doc/kalman_tangent_shared_modes'
ORDER = ['QSM', 'TFG', 'TAN', 'QSM-G', 'TFG-G', 'TAN-G', 'CARMA']


def check_evidence(verify_current=True):
    manifest = json.loads((OUT/'manifest.json').read_text())
    rows = list(csv.DictReader((OUT/'scores.csv').open()))
    assert rows and manifest['profile'] == 'acceptance'
    assert not manifest.get('source_changed_during_run', False)
    tree = manifest.get('implementation_source_tree', manifest['implementation_commit'])
    assert re.fullmatch(r'[0-9a-f]{40}', tree)
    if not verify_current:
        # A historical article is bound to the actual producer tree, not
        # merely to a label saying the current source need not be checked.
        subprocess.run(['git', 'cat-file', '-e', tree+'^{tree}'], cwd=ROOT, check=True)
    for path, sha in manifest['source_sha256'].items():
        assert len(sha) == 64
        if verify_current:
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == sha, f'stale source: {path}'
        else:
            source = subprocess.check_output(['git', 'show', tree+':'+path], cwd=ROOT)
            assert hashlib.sha256(source).hexdigest() == sha, f'stale pinned source: {path}'
    pairs = {}
    for r in rows:
        key = (r['scenario'], r['seed'], r['noise_profile'])
        pairs.setdefault(key, set()).add(r['sensor_hash'])
    assert all(len(v) == 1 for v in pairs.values()), 'unpaired sensors'
    assert len(pairs) == len(manifest['commands'])
    assert set(manifest['seeds']) == {int(r['seed']) for r in rows}
    assert len(manifest['seeds']) >= 6
    grouped = {}
    for row in rows:
        key = (row['scenario'], int(row['seed']), row['noise_profile'])
        grouped.setdefault(key, []).append(row)
    expected = {(case, seed, 'native') for case in manifest['cases'] for seed in manifest['seeds']}
    if 'wave' in manifest['cases']:
        expected.update(('wave', seed, mode) for seed in manifest['seeds'] for mode in ('noise', 'population'))
    assert set(grouped) == expected, 'missing or unexpected condition/seed'
    for (case, seed, mode), group in grouped.items():
        filters = set(ORDER) | ({'QSM-LG', 'TFG-LG', 'TAN-LG'} if case == 'wave' else set())
        end = 7200 if case == 'long' else 900
        windows = {(150, end), (150, 300), (300, 600), (600, 900)}
        if case == 'long': windows.add((5400, 7200))
        actual = {(r['filter'], float(r['start']), float(r['end'])) for r in group}
        assert actual == {(f, b, e) for f in filters for b, e in windows}, (case, seed, mode)
        assert len(group) == len(actual), 'duplicate score rows'
    for name in ['tangent_accuracy', 'tangent_bias_tracking', 'tangent_motion',
                 'tangent_degradation', 'tangent_consistency', 'tangent_observability']:
        assert (DOC/(name+'.pgf')).exists()
    assert (DOC/'tangent-results-generated.tex-part').exists()
    publication = json.loads((OUT/'publication.json').read_text())
    for path, sha in publication['sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == sha, f'stale publication: {path}'
    binding = 'current source' if verify_current else 'pinned source receipt'
    print(f'evidence completeness/publication/{binding} PASS: {len(pairs)} histories, {len(rows)} rows')


if '--fetch-pinned-source' in sys.argv:
    producer = json.loads((OUT/'manifest.json').read_text())['implementation_commit']
    assert re.fullmatch(r'[0-9a-f]{40}', producer)
    subprocess.run(['git', 'fetch', '--no-tags', '--depth=1', 'origin', producer], cwd=ROOT, check=True)

if '--check-only' in sys.argv:
    check_evidence(verify_current='--allow-pinned-sources' not in sys.argv)
    raise SystemExit(0)

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0, str(ROOT/'plots'))
from vector_export import save_figure
from tangent_shared_modes_observability import analyze, structural_checks
from tangent_shared_modes_inputs import verify_raw_inputs

plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Liberation Serif', 'DejaVu Serif'],
                     'font.size': 11.5, 'axes.labelsize': 11.5, 'legend.fontsize': 9,
                     'pgf.texsystem': 'lualatex', 'pgf.rcfonts': False,
                     'pgf.preamble': r'\usepackage{newtxtext,newtxmath}',
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'svg.hashsalt': 'ocean-imu-tangent-shared-modes'})
COLORS = dict(zip(ORDER, ['#51719b', '#d38a34', '#278c70', '#99b0cb', '#dcb185', '#a24464', '#777777']))


def ci(values):
    v = np.asarray(values, dtype=float); v = v[np.isfinite(v)]
    if not len(v): return dict(n=0, mean=None, low=None, high=None)
    mean = float(v.mean())
    radius = float(stats.t.ppf(.975, len(v)-1)*stats.sem(v)) if len(v) > 1 else 0.
    return dict(n=len(v), mean=mean, low=mean-radius, high=mean+radius)


def paired(frame, before, after, column):
    b = frame[frame['filter'] == before].set_index('seed')[column]
    c = frame[frame['filter'] == after].set_index('seed')[column]
    d = pd.concat([b.rename('b'), c.rename('c')], axis=1).dropna()
    return dict(before=before, after=after, metric=column,
                difference=ci(d.c-d.b), improvement_percent=ci(100*(d.b-d.c)/d.b))


def number(x, digits=2):
    return '--' if x is None or not np.isfinite(x) else f'{x:.{digits}f}'


def interval(x, digits=2):
    return f"{number(x['mean'], digits)} [{number(x['low'], digits)}, {number(x['high'], digits)}]"


def table(caption, columns, rows, wide=True):
    env = 'table*' if wide else 'table'
    return '\n'.join([f'\\begin{{{env}}}[t]\\centering\\footnotesize',
                      '\\caption{'+caption+'}', '\\begin{adjustbox}{max width=\\linewidth}',
                      '\\begin{tabular}{l'+'r'*(len(columns)-1)+'}\\toprule',
                      '&'.join(columns)+r'\\\midrule']+
                     ['&'.join(row)+r'\\' for row in rows]+
                     [r'\bottomrule\end{tabular}\end{adjustbox}', f'\\end{{{env}}}'])


def write_macro(name, body):
    return '\\newcommand{\\'+name+'}{%\n'+body+'%\n}\n'


def save(fig, name):
    fig.tight_layout()
    save_figure(fig, DOC/(name+'.svg'), metadata={'Date': None})
    save_figure(fig, DOC/(name+'.pgf'))
    for suffix in ('.svg', '.pgf'):
        path = DOC/(name+suffix)
        path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')
    plt.close(fig)


def trace(raw, case='wave', *, seed, model='native'):
    return pd.read_csv(raw/f'{case}_{seed}_{model}_trace.csv')


def ljung_box(values, lag=20):
    y = np.asarray(values, dtype=float)
    y = y[np.isfinite(y)]
    result = dict(lag=lag, samples=len(y), statistic=None, p=None)
    if len(y) <= lag or np.all(y == y[0]):
        return result
    y = y-y.mean()
    energy = np.dot(y, y)
    if not np.isfinite(energy) or energy <= 0:
        return result
    acf = np.array([np.dot(y[k:], y[:-k])/energy for k in range(1, lag+1)])
    q = len(y)*(len(y)+2)*np.sum(acf**2/(len(y)-np.arange(1, lag+1)))
    if np.isfinite(q):
        result.update(statistic=float(q), p=float(stats.chi2.sf(q, lag)))
    return result


def available_statistics(value):
    """Keep missing endpoint/failed-filter statistics unavailable in strict JSON."""
    if isinstance(value, dict):
        return {key: available_statistics(item) for key, item in value.items()}
    if isinstance(value, list):
        return [available_statistics(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--skip-information', action='store_true')
    parser.add_argument('--allow-pinned-sources', action='store_true',
                        help='Regenerate a historical cohort against its recorded source receipt.')
    a = parser.parse_args()
    manifest = json.loads((OUT/'manifest.json').read_text())
    verify_raw_inputs(a.raw, manifest, OUT/'scores.csv')
    representative = manifest['seeds'][0]
    OUT.mkdir(parents=True, exist_ok=True); DOC.mkdir(parents=True, exist_ok=True)
    # Canonical publication line endings; the original raw score receipts
    # remain untouched. This changes no score value.
    score_path = OUT/'scores.csv'
    score_path.write_bytes(score_path.read_bytes().replace(b'\r\n', b'\n'))
    df = pd.read_csv(OUT/'scores.csv')
    for prefix in ('ba', 'bg', 'v'):
        df[prefix+'_h'] = np.hypot(df[prefix+'_x'], df[prefix+'_y'])
        df[prefix+'_3d'] = np.sqrt(sum(df[prefix+'_'+axis]**2 for axis in 'xyz'))
    primary = df[(df.scenario == 'wave') & (df.noise_profile == 'native') & (df.start == 150) & (df.end == 900)]
    cohort = df[(df.noise_profile == 'native') & (df.start == 150) & (df.end == 900)]
    comparisons = [('QSM', 'TFG'), ('TFG', 'TAN'), ('QSM', 'QSM-G'), ('TFG', 'TFG-G'),
                   ('TAN', 'TAN-G'), ('TFG-G', 'TAN-G'), ('TAN', 'TAN-LG')]
    metrics = ['att_rms', 'heave_rms', 'v_h', 'v_z', 'ba_x', 'ba_y', 'ba_z', 'ba_h', 'bg_3d']
    pairs = [paired(primary, b, c, metric) for b, c in comparisons for metric in metrics]
    summary = dict(primary={f:{m:ci(primary[primary['filter'] == f][m]) for m in metrics} for f in ORDER},
                   paired=pairs, failures=int(df.groupby(['scenario', 'seed', 'noise_profile', 'filter']).failed.max().sum()),
                   failed_processes=sum(r['exit_code'] != 0 for r in manifest['commands']),
                   incomplete_experiments=manifest.get('incomplete_experiments', 0),
                   structural_observability=structural_checks())
    rows = []
    for f in ORDER:
        d = primary[primary['filter'] == f]
        rows.append([f, number(d.att_rms.mean(), 3), number(d.att_max.max(), 3)]+[number(d[c].mean()*s)
                         for c, s in [('heave_rms', 1000), ('v_h', 1000), ('v_z', 1000), ('ba_h', 1000), ('horizontal_drift', 1)]])
    main_table = table('Native acquired-reference waves, 150--900 s. Six holdout seeds; entries are means except worst attitude, which is the cohort maximum.',
                       ['Method', 'Att. [deg]', 'Worst att. [deg]', 'Heave [mm]', '$v_h$ [mm/s]', '$v_z$ [mm/s]', '$b_{a,h}$ [mm/s$^2$]', 'Wave drift [m]'], rows)
    rows = []
    for f in ORDER:
        d = primary[primary['filter'] == f]
        rows.append([f]+[number(d[p+'_'+axis].mean()*1000) for p in ('ba', 'bg') for axis in 'xyz'])
    bias_table = table('Per-axis bias RMS on the same waves. Accelerometer: mm/s$^2$; gyro: mrad/s.',
                       ['Method', '$b_{ax}$', '$b_{ay}$', '$b_{az}$', '$b_{gx}$', '$b_{gy}$', '$b_{gz}$'], rows)
    rows = []
    for b, c in comparisons:
        digits = 4 if (b, c) in [('TFG', 'TAN'), ('TFG-G', 'TAN-G')] else 2
        rows.append([b+r'$\to$'+c]+[interval(paired(primary, b, c, m)['improvement_percent'], digits) for m in ('att_rms', 'heave_rms', 'ba_h')])
    paired_table = table(r'Paired percentage improvement (positive is lower error), mean [95\% Student interval], six independent seeds.',
                         ['Comparison', r'Attitude [\%]', r'Heave [\%]', r'Horizontal bias [\%]'], rows)

    td = trace(a.raw, seed=representative)
    endpoint = pd.concat([trace(a.raw, seed=s).query('t == 899') for s in manifest['seeds']], ignore_index=True)
    bounds = stats.chi2.ppf([.025, .975], 15*len(manifest['seeds']))/(15*len(manifest['seeds']))
    bounds3 = stats.chi2.ppf([.025, .975], 3*len(manifest['seeds']))/(3*len(manifest['seeds']))
    consistency_rows = []
    for f in ORDER:
        d = primary[primary['filter'] == f]; e = endpoint[endpoint['filter'] == f]
        consistency_rows.append([f, number(e.nees15.mean()/15), number(e.bias_nees.mean()/3),
                                 number(e.acc_nis.mean()/3), number(e.mag_nis.mean()/3),
                                 number(d.acc_nis.mean()/3), number(d.mag_nis.mean()/3),
                                 number(d[['innovation_acf1_x', 'innovation_acf1_y', 'innovation_acf1_z']].abs().max(axis=1).mean(), 3)])
    consistency_table = table('Normalized consistency: NEES and packet NIS at fixed 899 s, plus accepted-packet NIS window means and lag-one innovations at 20 Hz. Nominal mean is one. Fixed-time 95\\% ensemble bounds: NEES$_{15}$ ['+number(bounds[0])+', '+number(bounds[1])+']; bias NEES/NIS$_3$ ['+number(bounds3[0])+', '+number(bounds3[1])+'].',
                              ['Method', 'NEES$_{15}$/15', 'Bias NEES/3', 'Acc. NIS/3', 'Mag. NIS/3', 'Mean acc./3', 'Mean mag./3', 'Max. $|\\rho_1|$'], consistency_rows)
    summary['nees15_ensemble_95_bounds_normalized'] = bounds.tolist()
    summary['nees3_nis3_ensemble_95_bounds_normalized'] = bounds3.tolist()
    summary['endpoint_consistency'] = {f:dict(nees15=float(endpoint[endpoint['filter'] == f].nees15.mean()/15),
                                             bias_nees=float(endpoint[endpoint['filter'] == f].bias_nees.mean()/3),
                                             acc_nis=float(endpoint[endpoint['filter'] == f].acc_nis.mean()/3),
                                             mag_nis=float(endpoint[endpoint['filter'] == f].mag_nis.mean()/3)) for f in ORDER}

    recovery_rows = []
    for case in ('rocking', 'late_bias', 'late_attitude', 'transition', 'vibration', 'long'):
        if case == 'long': d = df[(df.scenario == case) & (df.start == 5400)]
        else: d = cohort[cohort.scenario == case]
        for f in ('TFG', 'TAN', 'TAN-G'):
            x = d[d['filter'] == f]; times = x.bias_convergence[x.bias_convergence >= 0]
            recovery_rows.append([case.replace('_', r'\_'), f, number(x.ba_h.mean()*1000),
                                  f'{len(times)}/{len(x)}', number(times.median(), 1)])
    recovery_table = table('Degradation and bias convergence. Long uses 5400--7200 s. Other RMS windows are 150--900 s. Convergence time is absolute time since construction.',
                           ['Condition', 'Method', '$b_{a,h}$ [mm/s$^2$]', 'Converged', 'Median time [s]'], recovery_rows)
    robustness_rows = []
    for case in ('high', 'inclination87', 'yaw'):
        for f in ('TFG', 'TAN', 'TFG-G'):
            x = cohort[(cohort.scenario == case) & (cohort['filter'] == f)]
            robustness_rows.append([case, f, number(x.att_rms.mean(), 2),
                                    number(x.att_max.max(), 2), f'{int((x.att_max > 10).sum())}/{len(x)}',
                                    number(x.bias_enabled_fraction.min(), 3),
                                    number(x.mag_applied_fraction.min(), 3)])
    robustness_table = table('Accuracy excursions despite finite, positive covariance. Attitude columns are mean RMS and worst peak across six seeds. The 10-degree count is descriptive and changes no replay gate. Last columns show the worst learning/accepted-magnetometer fractions.',
                             ['Condition', 'Method', 'Att. RMS [deg]', 'Peak [deg]', 'Peak $>10$ deg', 'Bias enabled', 'Mag. accepted'], robustness_rows)
    summary['robustness'] = []
    for (case, f), x in cohort.groupby(['scenario', 'filter']):
        summary['robustness'].append(dict(scenario=case, filter=f,
            peak_attitude_deg=float(x.att_max.max()), peak_over_10_deg=int((x.att_max > 10).sum()),
            trials=len(x), minimum_bias_enabled_fraction=float(x.bias_enabled_fraction.min()),
            minimum_mag_applied_fraction=float(x.mag_applied_fraction.min())))
    cost_rows = []
    for f in ORDER:
        d = primary[primary['filter'] == f]
        cost_rows.append([f, str(int(d.states.iloc[0])), str(int(d.memory_bytes.iloc[0])),
                          number(d.step_us.median(), 2), number(d.step_us.min(), 2)+'--'+number(d.step_us.max(), 2)])
    cost_table = table('Host update CPU time after warmup, per IMU step; memory includes the heap-owned core. No ESP32 board timing is claimed.',
                       ['Method', 'Errors', 'Host bytes', r'Median [$\mu$s]', r'Seed range [$\mu$s]'], cost_rows)
    timing_rows = []
    for f in ORDER:
        d = primary[primary['filter'] == f]
        att = d.att_convergence[d.att_convergence >= 0]
        bias = d.bias_convergence[d.bias_convergence >= 0]
        timing_rows.append([f, number(d.live_time.median(), 1),
                            f'{len(att)}/{len(d)}', number(att.median(), 1),
                            f'{len(bias)}/{len(d)}', number(bias.median(), 1)])
    timing_table = table('Primary-wave startup and sustained convergence. Times are seconds since construction; convergence requires 30 s of evaluated samples below the declared threshold. Missing convergence is counted.',
                         ['Method', 'Median Live [s]', 'Att. converged', 'Median att. [s]', 'Bias converged', 'Median bias [s]'], timing_rows)

    ablation_rows = []
    ablations = [('wave', 'native'), ('exact_reference', 'native'), ('reference_error', 'native'),
                 ('wave', 'noise'), ('wave', 'population'), ('initial_bias', 'native'),
                 ('initial_attitude', 'native'), ('high', 'native')]
    for case, mode in ablations:
        d = df[(df.scenario == case) & (df.noise_profile == mode) & (df.start == 150) & (df.end == 900)]
        for f in ('TFG', 'TAN', 'TFG-G'):
            x = d[d['filter'] == f]
            ablation_rows.append([case.replace('_', r'\_')+'/'+mode, f,
                                  number(x.att_rms.mean(), 3), number(x.ba_h.mean()*1000),
                                  number(x.heave_rms.mean()*1000), number(x.nees15.mean()/15)])
    ablation_table = table('Reference, likelihood/prior and live-initialization ablations. Noise matches packet noise; population additionally matches the turn-on bias prior. NEES is a descriptive serial mean.',
                           ['Control', 'Method', 'Att. [deg]', '$b_{a,h}$ [mm/s$^2$]', 'Heave [mm]', 'NEES/15'], ablation_rows)

    startup = []
    for seed in manifest['seeds']:
        data = trace(a.raw, seed=seed)
        left = data[data['filter'] == 'TFG'].set_index('t')
        right = data[data['filter'] == 'TAN'].set_index('t')
        held = (left.bias_enabled == 0) & (right.bias_enabled == 0)
        cols = ['att', 'px', 'py', 'pz', 'vx', 'vy', 'vz', 'bax', 'bay', 'baz', 'bgx', 'bgy', 'bgz']
        startup.append(dict(seed=seed, held_trace_nodes=int(held.sum()),
                            maximum_physical_difference=float((left.loc[held, cols]-right.loc[held, cols]).abs().max().max()),
                            qsm_live_s=float(primary[(primary['filter'] == 'QSM') & (primary.seed == seed)].live_time.iloc[0]),
                            tfg_live_s=float(primary[(primary['filter'] == 'TFG') & (primary.seed == seed)].live_time.iloc[0]),
                            tangent_live_s=float(primary[(primary['filter'] == 'TAN') & (primary.seed == seed)].live_time.iloc[0])))
    summary['startup_parity'] = startup

    # Accuracy by condition, same seed-level evaluation window.
    cases = ['stationary', 'rocking', 'low', 'wave', 'high', 'yaw', 'transition', 'vibration', 'colored', 'inclination87']
    fig, axes = plt.subplots(2, 1, figsize=(10.6, 5.5), sharex=True)
    for f, offset in zip(('QSM', 'TFG', 'TAN', 'TAN-G'), (-.24, -.08, .08, .24)):
        for ax, metric, unit, label in zip(axes, ('ba_h', 'heave_rms'), (1000, 1000), ('Horizontal bias RMS [mm/s²]', 'Heave RMS [mm]')):
            means, err = [], []
            for case in cases:
                values = cohort[(cohort.scenario == case) & (cohort['filter'] == f)][metric]*unit
                c = ci(values); means.append(c['mean']); err.append((c['high']-c['low'])/2)
            ax.errorbar(np.arange(len(cases))+offset, means, yerr=err, marker='o', markersize=3, capsize=2, linewidth=.8, label=f, color=COLORS[f])
            ax.set_ylabel(label); ax.grid(alpha=.2)
    axes[0].legend(ncol=4); axes[-1].set_xticks(range(len(cases)), cases, rotation=25, ha='right')
    save(fig, 'tangent_accuracy')

    chosen = ['QSM', 'TFG', 'TAN', 'TAN-G']
    fig, axes = plt.subplots(3, 2, figsize=(10.6, 5.8), sharex=True)
    truth = td[td['filter'] == 'QSM']
    for i, axis in enumerate('xyz'):
        for j, prefix in enumerate(('ba', 'bg')):
            drift = trace(a.raw, 'bias_drift' if prefix == 'ba' else 'gyro_drift', seed=representative)
            truth_drift = drift[drift['filter'] == 'QSM']
            ax = axes[i, j]
            for f in chosen:
                d = drift[drift['filter'] == f]; ax.plot(d.t, d[prefix+axis]*1000, label=f, color=COLORS[f], lw=.8)
            ax.plot(truth_drift.t, truth_drift['true_'+prefix+axis]*1000, color='black', ls='--', lw=1, label='Truth')
            ax.set_ylabel(f'{prefix} {axis} '+('[mm/s²]' if j == 0 else '[mrad/s]')); ax.grid(alpha=.2)
    axes[0, 0].legend(ncol=3); axes[-1, 0].set_xlabel('Time [s]'); axes[-1, 1].set_xlabel('Time [s]')
    save(fig, 'tangent_bias_tracking')

    fig, axes = plt.subplots(3, 2, figsize=(10.6, 5.8))
    for f in chosen:
        d = td[td['filter'] == f]
        axes[0, 0].plot(d.t, d.att, label=f, color=COLORS[f], lw=.8)
        axes[1, 0].plot(d.t, -d.pz, color=COLORS[f], lw=.8)
        last = d[d.t >= 870]
        axes[1, 1].plot(last.t, -last.pz, color=COLORS[f], lw=.9)
        for ax, column in zip(axes[2], ('px', 'py')): ax.plot(last.t, last[column], color=COLORS[f], lw=.9)
    for axis, color in zip('xyz', ['#456990', '#cf7a35', '#498f6d']):
        d = td[td['filter'] == 'TAN']; axes[0, 1].plot(d.t, d['phi_'+axis], color=color, label=axis.upper(), lw=.8)
    for ax, column, mask, sign in [(axes[1, 0], 'true_pz', truth.t >= 0, -1), (axes[1, 1], 'true_pz', truth.t >= 870, -1),
                                   (axes[2, 0], 'true_px', truth.t >= 870, 1), (axes[2, 1], 'true_py', truth.t >= 870, 1)]:
        ax.plot(truth.t[mask], sign*truth[column][mask], color='black', ls='--', lw=1)
    for ax, label in zip(axes.flat, ['Attitude error [deg]', 'TAN rotation-error components [deg]', 'Heave [m]', 'Heave, final 30 s [m]', 'X, final 30 s [m]', 'Y, final 30 s [m]']):
        ax.set_ylabel(label); ax.set_xlabel('Time [s]'); ax.grid(alpha=.2)
    axes[0, 0].legend(ncol=2); axes[0, 1].legend(ncol=3)
    save(fig, 'tangent_motion')

    fig, axes = plt.subplots(2, 2, figsize=(10.6, 4.8), sharex=True)
    for j, case in enumerate(('transition', 'vibration')):
        data = trace(a.raw, case, seed=representative)
        for f in chosen:
            d = data[data['filter'] == f]
            axes[0, j].plot(d.t, d.true_pz-d.pz, label=f, color=COLORS[f], lw=.8)
            axes[1, j].plot(d.t, 1000*np.hypot(d.bax-d.true_bax, d.bay-d.true_bay), color=COLORS[f], lw=.8)
        axes[0, j].set_title(case); axes[0, j].set_ylabel('Heave error [m]'); axes[1, j].set_ylabel('Horizontal bias error [mm/s²]'); axes[1, j].set_xlabel('Time [s]')
        for ax in axes[:, j]:
            ax.grid(alpha=.2)
            if case == 'transition':
                for t in (300, 600): ax.axvline(t, color='grey', ls=':', lw=.8)
    axes[0, 0].legend(ncol=4)
    save(fig, 'tangent_degradation')

    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.4))
    x = np.arange(len(ORDER))
    axes[0].bar(x, [summary['endpoint_consistency'][f]['nees15'] for f in ORDER], color=[COLORS[f] for f in ORDER])
    axes[0].axhspan(*bounds, color='#8fbe93', alpha=.35); axes[0].axhline(1, color='black', lw=.8); axes[0].set_yscale('log'); axes[0].set_ylabel('NEES15 / 15 at 899 s')
    for metric, offset in [('acc_nis', -.17), ('mag_nis', .17)]:
        axes[1].bar(x+offset, [primary[primary['filter'] == f][metric].mean()/3 for f in ORDER],
                    width=.32, label='Accelerometer' if metric == 'acc_nis' else 'Magnetometer')
    axes[1].axhline(1, color='black', lw=.8); axes[1].set_ylabel('Mean NIS / 3'); axes[1].legend()
    axes[2].bar(x, [primary[primary['filter'] == f][['innovation_acf1_x', 'innovation_acf1_y', 'innovation_acf1_z']].abs().max(axis=1).mean() for f in ORDER], color=[COLORS[f] for f in ORDER]); axes[2].set_ylabel('Maximum |innovation lag 1| at 20 Hz')
    for ax in axes: ax.set_xticks(x, ORDER, rotation=50, ha='right'); ax.grid(axis='y', alpha=.2)
    save(fig, 'tangent_consistency')

    information_path = OUT/'observability.json'
    if not a.skip_information:
        information = {}
        for case in ('rocking', 'wave', 'yaw', 'inclination87'):
            for sd, label in [(0., 'known'), (2., 'uncertain'), (float('inf'), 'free')]:
                print('information', case, label, flush=True)
                information[case+'_'+label] = analyze(a.raw/f'{case}_{representative}_native', field_sd=sd)
        for case in ('rocking', 'wave'):
            information[case+'_motion4'] = analyze(a.raw/f'{case}_{representative}_native', prior_scale=4.)
        information_path.write_text(json.dumps(information, indent=2, allow_nan=False)+'\n')
    information = json.loads(information_path.read_text())
    for key, receipt in information.items():
        prefix = a.raw/(key.rsplit('_', 1)[0]+f'_{representative}_native')
        assert receipt['analysis_sha256'] == hashlib.sha256((ROOT/'tools/tangent_shared_modes_observability.py').read_bytes()).hexdigest()
        for suffix, sha in receipt['input_sha256'].items():
            assert hashlib.sha256(Path(str(prefix)+suffix).read_bytes()).hexdigest() == sha, f'stale information input: {prefix}'
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.4))
    for j, label in enumerate(('known', 'uncertain', 'free')):
        for case, marker in zip(('rocking', 'wave', 'yaw', 'inclination87'), ('o', 's', '^', 'd')):
            values = np.array(information[case+'_'+label]['bias_only']['eigenvalues'])
            axes[j].plot([1, 2, 3], np.maximum(values, 1e-12), marker=marker, label=case)
        axes[j].set_yscale('log'); axes[j].set_title('Magnetic reference: '+label); axes[j].set_xlabel('Bias-information eigenvalue index'); axes[j].grid(alpha=.2)
    axes[0].set_ylabel('Nuisance-eliminated information'); axes[0].legend()
    save(fig, 'tangent_observability')
    summary['window_solver_converged'] = all(r[k]['converged'] for r in information.values() for k in ('carma', 'sensor_only', 'bias_only'))
    summary['window_profiles'] = {k: dict(bias_eigenvalues=v['bias_only']['eigenvalues'], field=v['field_assumption']) for k, v in information.items()}

    # One-Hz innovations: descriptive Ljung--Box diagnostics, not extra
    # independent Monte Carlo samples. Report all rejections rather than
    # selecting a favorable axis or unreported multiple-test correction.
    whiteness = []
    for seed in manifest['seeds']:
        data = trace(a.raw, seed=seed)
        for f in ORDER:
            d = data[(data['filter'] == f) & (data.t >= 150)]
            for axis in 'xyz':
                whiteness.append(dict(seed=seed, filter=f, axis=axis,
                                      **ljung_box(d['white_'+axis].to_numpy())))
    summary['whiteness'] = whiteness
    summary['whiteness_available_tests'] = sum(r['statistic'] is not None for r in whiteness)
    summary['whiteness_rejections_at_05'] = sum(r['p'] is not None and r['p'] < .05 for r in whiteness)
    (OUT/'summary.json').write_text(json.dumps(available_statistics(summary), indent=2, allow_nan=False)+'\n')

    tangent_gain = paired(primary, 'TFG', 'TAN', 'ba_h')['improvement_percent']
    headline = ('On acquired-reference waves, the paired TAN versus TFG horizontal-bias RMS change is '
                +number(tangent_gain['mean'], 4)+'\\% improvement [95\\% interval '
                +number(tangent_gain['low'], 4)+', '+number(tangent_gain['high'], 4)+'\\%].')
    reliable = (tangent_gain['low'] is not None and tangent_gain['low'] > 0 and summary['failures'] == 0
                and summary['failed_processes'] == 0 and summary['incomplete_experiments'] == 0
                and summary['endpoint_consistency']['TAN']['nees15'] <= bounds[1])
    recommendation = ('The measured cohort warrants further evaluation of the tangent-bias candidate, subject to independent replay and target qualification.' if reliable else
                      'Do not promote the tangent-bias geometry or the conditional gain veto as a new production default from this evidence. Retain the experimental implementation and transfer only independently validated shared corrections. Geometry and a conservative gain policy have not jointly established the required accuracy, consistency, robustness and target qualification.')
    interpretation = ('The paired tables isolate a finite-retraction/reset effect from the bias-gain effect. '
                      'Noise/population-prior controls and the condition-specific scores are retained in the CSV, including numerical failures and missing convergence. '
                      'The conditional excitation monitor and the full nuisance-eliminated window information use different assumptions; agreement is not presumed. '
                      'Fixed-time normalized NEES ensemble bounds are ['+number(bounds[0])+', '+number(bounds[1])+']. '
                      'All covariance statements remain conditional on the declared CARMA, sensor-noise, calibration and bias priors. '
                      'The one-Hz Ljung--Box diagnostics and all p-values are retained in the JSON; they are descriptive serial tests, not new independent trials.')
    native_tfg = primary[primary['filter'] == 'TFG']
    exact_tfg = cohort[(cohort.scenario == 'exact_reference') & (cohort['filter'] == 'TFG')]
    matched_tfg = df[(df.scenario == 'wave') & (df.noise_profile == 'noise') &
                     (df['filter'] == 'TFG') & (df.start == 150) & (df.end == 900)]
    interpretation += (' Supplying the exact field changes TFG horizontal-bias RMS from '+number(native_tfg.ba_h.mean()*1000)+
                       ' to '+number(exact_tfg.ba_h.mean()*1000)+' mm/s$^2$; it does not remove the remaining error. '
                       'Matching the packet-noise model gives '+number(matched_tfg.ba_h.mean()*1000)+
                       ' mm/s$^2$ and normalized serial NEES '+number(matched_tfg.nees15.mean()/15)+
                       ', so packet-noise mismatch alone does not explain the weak bias estimate. '
                       'Across the recorded bias-held startup nodes the maximum TFG/TAN physical difference is '+
                       number(max(s['maximum_physical_difference'] for s in startup), 6)+
                       '; startup is shared rather than a source of the geometric comparison. ')
    interpretation += f' At the unadjusted 5\\% level, {summary["whiteness_rejections_at_05"]} of {summary["whiteness_available_tests"]} available axis/seed/method whiteness tests reject; these correlated multiple tests are not independent acceptance votes.'
    macros = [write_macro('StudyHeadline', headline), write_macro('StudyCohort',
              f'The cohort contains {len(manifest["commands"])} paired sensor histories, with '+str(len(manifest['seeds']))+
              ' independent holdout seeds per condition. Conditions share the declared seed set. The scores retain every declared variant and evaluation phase; failed runs are counted explicitly. '+
              f'Numerically failed filter/history pairs: {summary["failures"]}; failed native processes: {summary["failed_processes"]}; incomplete experiments: {summary["incomplete_experiments"]}.'),
              write_macro('PrimaryResultsTable', main_table), write_macro('BiasResultsTable', bias_table),
              write_macro('PairedResultsTable', paired_table), write_macro('ConsistencyResultsTable', consistency_table),
              write_macro('RecoveryResultsTable', recovery_table), write_macro('CostResultsTable', cost_table),
              write_macro('RobustnessResultsTable', robustness_table),
              write_macro('TimingResultsTable', timing_table),
              write_macro('AblationResultsTable', ablation_table),
              write_macro('StudyInterpretation', interpretation), write_macro('StudyRecommendation', recommendation),
              write_macro('ValidationResults', r'Validation receipts, untouched-main parity and pinned-replay outcomes are recorded in \texttt{validation.json}; they must be read alongside the statistical cohort.')]
    (DOC/'tangent-results-generated.tex-part').write_text('% Generated from source-bound paired scores.\n'+''.join(macros))
    publication_files = [OUT/name for name in ('scores.csv', 'manifest.json', 'observability.json', 'summary.json')]
    publication_files += sorted(DOC.glob('tangent_*.pgf')) + sorted(DOC.glob('tangent_*.svg'))
    publication_files += [DOC/'tangent-results-generated.tex-part']
    (OUT/'publication.json').write_text(json.dumps(dict(
        report_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in publication_files}), indent=2)+'\n')
    check_evidence(verify_current=not a.allow_pinned_sources)
    print(headline)


if __name__ == '__main__': main()
