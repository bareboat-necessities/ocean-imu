#!/usr/bin/env python3
"""Paired statistical report and bounded publication figures for TAN research."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'reports/results/tangent_shared_modes'
DOC = ROOT/'doc/kalman_tangent_shared_modes'
ORDER = ['QSM', 'TFG', 'TAN', 'QSM-G', 'TFG-G', 'TAN-G', 'CARMA']


def check_evidence():
    manifest = json.loads((OUT/'manifest.json').read_text())
    rows = list(csv.DictReader((OUT/'scores.csv').open()))
    assert rows and manifest['profile'] == 'acceptance'
    for path, sha in manifest['source_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == sha, f'stale source: {path}'
    pairs = {}
    for r in rows:
        key = (r['scenario'], r['seed'], r['noise_profile'])
        pairs.setdefault(key, set()).add(r['sensor_hash'])
    assert all(len(v) == 1 for v in pairs.values()), 'unpaired sensors'
    assert len(pairs) == len(manifest['commands'])
    assert set(manifest['seeds']) == {int(r['seed']) for r in rows}
    assert len(manifest['seeds']) >= 6
    for name in ['tangent_accuracy', 'tangent_bias_tracking', 'tangent_motion',
                 'tangent_degradation', 'tangent_consistency', 'tangent_observability']:
        assert (DOC/(name+'.pgf')).exists()
    assert (DOC/'tangent-results-generated.tex-part').exists()
    print(f'evidence structure/source binding PASS: {len(pairs)} histories, {len(rows)} rows')


if '--check-only' in sys.argv:
    check_evidence()
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

plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Liberation Serif', 'DejaVu Serif'],
                     'font.size': 11.5, 'axes.labelsize': 11.5, 'legend.fontsize': 9,
                     'pgf.texsystem': 'lualatex', 'pgf.rcfonts': False,
                     'pgf.preamble': r'\usepackage{newtxtext,newtxmath}',
                     'axes.spines.top': False, 'axes.spines.right': False})
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


def interval(x):
    return f"{number(x['mean'])} [{number(x['low'])}, {number(x['high'])}]"


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
    save_figure(fig, DOC/(name+'.svg'))
    save_figure(fig, DOC/(name+'.pgf'))
    plt.close(fig)


def trace(raw, case='wave', seed=4101, model='native'):
    return pd.read_csv(raw/f'{case}_{seed}_{model}_trace.csv')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--skip-information', action='store_true')
    a = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True); DOC.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(OUT/'scores.csv')
    for prefix in ('ba', 'bg', 'v'):
        df[prefix+'_h'] = np.hypot(df[prefix+'_x'], df[prefix+'_y'])
        df[prefix+'_3d'] = np.sqrt(sum(df[prefix+'_'+axis]**2 for axis in 'xyz'))
    manifest = json.loads((OUT/'manifest.json').read_text())
    primary = df[(df.scenario == 'wave') & (df.noise_profile == 'native') & (df.start == 150) & (df.end == 900)]
    cohort = df[(df.noise_profile == 'native') & (df.start == 150) & (df.end == 900)]
    comparisons = [('QSM', 'TFG'), ('TFG', 'TAN'), ('QSM', 'QSM-G'), ('TFG', 'TFG-G'),
                   ('TAN', 'TAN-G'), ('TFG-G', 'TAN-G'), ('TAN', 'TAN-LG')]
    metrics = ['att_rms', 'heave_rms', 'v_h', 'v_z', 'ba_x', 'ba_y', 'ba_z', 'ba_h', 'bg_3d']
    pairs = [paired(primary, b, c, metric) for b, c in comparisons for metric in metrics]
    summary = dict(primary={f:{m:ci(primary[primary['filter'] == f][m]) for m in metrics} for f in ORDER},
                   paired=pairs, failures=int(df.groupby(['scenario', 'seed', 'noise_profile', 'filter']).failed.max().sum()),
                   structural_observability=structural_checks())
    rows = []
    for f in ORDER:
        d = primary[primary['filter'] == f]
        rows.append([f]+[number(d[c].mean()*s, 3 if c == 'att_rms' else 2)
                         for c, s in [('att_rms', 1), ('heave_rms', 1000), ('v_h', 1000), ('v_z', 1000), ('ba_h', 1000), ('horizontal_drift', 1)]])
    main_table = table('Native acquired-reference waves, 150--900 s. Six holdout seeds; entries are means.',
                       ['Method', 'Att. [deg]', 'Heave [mm]', '$v_h$ [mm/s]', '$v_z$ [mm/s]', '$b_{a,h}$ [mm/s$^2$]', 'Wave drift [m]'], rows)
    rows = []
    for f in ORDER:
        d = primary[primary['filter'] == f]
        rows.append([f]+[number(d[p+'_'+axis].mean()*1000) for p in ('ba', 'bg') for axis in 'xyz'])
    bias_table = table('Per-axis bias RMS on the same waves. Accelerometer: mm/s$^2$; gyro: mrad/s.',
                       ['Method', '$b_{ax}$', '$b_{ay}$', '$b_{az}$', '$b_{gx}$', '$b_{gy}$', '$b_{gz}$'], rows)
    rows = []
    for b, c in comparisons:
        rows.append([b+r'$\to$'+c]+[interval(paired(primary, b, c, m)['improvement_percent']) for m in ('att_rms', 'heave_rms', 'ba_h')])
    paired_table = table(r'Paired percentage improvement (positive is lower error), mean [95\% Student interval], six independent seeds.',
                         ['Comparison', r'Attitude [\%]', r'Heave [\%]', r'Horizontal bias [\%]'], rows)

    td = trace(a.raw)
    endpoint = pd.concat([trace(a.raw, seed=s).query('t == 899') for s in manifest['seeds']], ignore_index=True)
    bounds = stats.chi2.ppf([.025, .975], 15*len(manifest['seeds']))/(15*len(manifest['seeds']))
    consistency_rows = []
    for f in ORDER:
        d = primary[primary['filter'] == f]; e = endpoint[endpoint['filter'] == f]
        consistency_rows.append([f, number(e.nees15.mean()/15), number(e.bias_nees.mean()/3),
                                 number(d.acc_nis.mean()/3), number(d.mag_nis.mean()/3),
                                 number(d[['innovation_acf1_x', 'innovation_acf1_y', 'innovation_acf1_z']].abs().max(axis=1).mean(), 3)])
    consistency_table = table('Normalized consistency: NEES at fixed 899 s; NIS window means; lag-one innovations at 20 Hz. Nominal normalized mean is one.',
                              ['Method', 'NEES$_{15}$/15', 'Bias NEES/3', 'Acc. NIS/3', 'Mag. NIS/3', 'Max. $|\rho_1|$'], consistency_rows)
    summary['nees15_ensemble_95_bounds_normalized'] = bounds.tolist()
    summary['endpoint_consistency'] = {f:dict(nees15=float(endpoint[endpoint['filter'] == f].nees15.mean()/15),
                                             bias_nees=float(endpoint[endpoint['filter'] == f].bias_nees.mean()/3)) for f in ORDER}

    recovery_rows = []
    for case in ('rocking', 'late_bias', 'late_attitude', 'transition', 'vibration', 'long'):
        if case == 'long': d = df[(df.scenario == case) & (df.start == 5400)]
        else: d = cohort[cohort.scenario == case]
        for f in ('TFG', 'TAN', 'TAN-G'):
            x = d[d['filter'] == f]; times = x.bias_convergence[x.bias_convergence >= 0]
            recovery_rows.append([case, f, number(x.ba_h.mean()*1000),
                                  f'{len(times)}/{len(x)}', number(times.median(), 1)])
    recovery_table = table('Degradation and bias convergence. Long uses 5400--7200 s. Other RMS windows are 150--900 s. Convergence time is absolute time since construction.',
                           ['Condition', 'Method', '$b_{a,h}$ [mm/s$^2$]', 'Converged', 'Median time [s]'], recovery_rows)
    cost_rows = []
    for f in ORDER:
        d = primary[primary['filter'] == f]
        cost_rows.append([f, str(int(d.states.iloc[0])), str(int(d.memory_bytes.iloc[0])),
                          number(d.step_us.median(), 2), number(d.step_us.min(), 2)+'--'+number(d.step_us.max(), 2)])
    cost_table = table('Host update CPU time after warmup, per IMU step; memory includes the heap-owned core. No ESP32 board timing is claimed.',
                       ['Method', 'Errors', 'Host bytes', r'Median [$\mu$s]', r'Seed range [$\mu$s]'], cost_rows)

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
    fig, axes = plt.subplots(3, 2, figsize=(10.6, 6.3), sharex=True)
    truth = td[td['filter'] == 'QSM']
    for i, axis in enumerate('xyz'):
        for j, prefix in enumerate(('ba', 'bg')):
            ax = axes[i, j]
            for f in chosen:
                d = td[td['filter'] == f]; ax.plot(d.t, d[prefix+axis]*1000, label=f, color=COLORS[f], lw=.8)
            ax.plot(truth.t, truth['true_'+prefix+axis]*1000, color='black', ls='--', lw=1, label='Truth')
            ax.set_ylabel(f'{prefix} {axis} '+('[mm/s²]' if j == 0 else '[mrad/s]')); ax.grid(alpha=.2)
    axes[0, 0].legend(ncol=3); axes[-1, 0].set_xlabel('Time [s]'); axes[-1, 1].set_xlabel('Time [s]')
    save(fig, 'tangent_bias_tracking')

    fig, axes = plt.subplots(3, 2, figsize=(10.6, 6.3))
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

    fig, axes = plt.subplots(2, 2, figsize=(10.6, 5.5), sharex=True)
    for j, case in enumerate(('transition', 'vibration')):
        data = trace(a.raw, case)
        for f in chosen:
            d = data[data['filter'] == f]
            axes[0, j].plot(d.t, d.pz-d.true_pz, label=f, color=COLORS[f], lw=.8)
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
        axes[1].bar(x+offset, [primary[primary['filter'] == f][metric].mean()/3 for f in ORDER], width=.32, label=metric)
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
                information[case+'_'+label] = analyze(a.raw/f'{case}_4101_native', field_sd=sd)
        for case in ('rocking', 'wave'):
            information[case+'_motion4'] = analyze(a.raw/f'{case}_4101_native', prior_scale=4.)
        information_path.write_text(json.dumps(information, indent=2, allow_nan=False)+'\n')
    information = json.loads(information_path.read_text())
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.4))
    for j, label in enumerate(('known', 'uncertain', 'free')):
        for case, marker in zip(('rocking', 'wave', 'yaw', 'inclination87'), ('o', 's', '^', 'd')):
            values = np.array(information[case+'_'+label]['bias_only']['eigenvalues'])
            axes[j].plot([1, 2, 3], np.maximum(values, 1e-12), marker=marker, label=case)
        axes[j].set_yscale('log'); axes[j].set_title('Magnetic reference: '+label); axes[j].set_xlabel('Scaled bias information eigenvalue'); axes[j].grid(alpha=.2)
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
                y = d['white_'+axis].to_numpy(); y = y[np.isfinite(y)]; y -= y.mean()
                acf = np.array([np.dot(y[k:], y[:-k])/np.dot(y, y) for k in range(1, 21)])
                q = len(y)*(len(y)+2)*np.sum(acf**2/(len(y)-np.arange(1, 21)))
                whiteness.append(dict(seed=seed, filter=f, axis=axis, lag=20, statistic=float(q), p=float(stats.chi2.sf(q, 20))))
    summary['whiteness'] = whiteness
    (OUT/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')

    tangent_gain = paired(primary, 'TFG', 'TAN', 'ba_h')['improvement_percent']
    headline = ('On acquired-reference waves, the paired TAN versus TFG horizontal-bias RMS change is '
                +number(tangent_gain['mean'])+'\\% improvement [95\\% interval '
                +number(tangent_gain['low'])+', '+number(tangent_gain['high'])+'\\%].')
    reliable = tangent_gain['low'] > 0 and summary['failures'] == 0 and summary['endpoint_consistency']['TAN']['nees15'] <= bounds[1]
    recommendation = ('The measured cohort warrants further evaluation of the tangent-bias candidate, subject to independent replay and target qualification.' if reliable else
                      'Do not promote the tangent-bias geometry or the conditional gain veto as a new production default from this evidence. Retain the experimental implementation and transfer only independently validated shared corrections. Geometry and a conservative gain policy have not jointly established the required accuracy, consistency, robustness and target qualification.')
    interpretation = ('The paired tables isolate a finite-retraction/reset effect from the bias-gain effect. '
                      'Noise/population-prior controls and the condition-specific scores are retained in the CSV, including numerical failures and missing convergence. '
                      'The conditional excitation monitor and the full nuisance-eliminated window information use different assumptions; agreement is not presumed. '
                      'Fixed-time normalized NEES ensemble bounds are ['+number(bounds[0])+', '+number(bounds[1])+']. '
                      'All covariance statements remain conditional on the declared CARMA, sensor-noise, calibration and bias priors. '
                      'The one-Hz Ljung--Box diagnostics and all p-values are retained in the JSON; they are descriptive serial tests, not new independent trials.')
    macros = [write_macro('StudyHeadline', headline), write_macro('StudyCohort',
              f'The cohort contains {len(manifest["commands"])} independent paired sensor histories, with '+str(len(manifest['seeds']))+
              ' holdout seeds. The scores retain every declared variant and evaluation phase; failed runs are counted explicitly.'),
              write_macro('PrimaryResultsTable', main_table), write_macro('BiasResultsTable', bias_table),
              write_macro('PairedResultsTable', paired_table), write_macro('ConsistencyResultsTable', consistency_table),
              write_macro('RecoveryResultsTable', recovery_table), write_macro('CostResultsTable', cost_table),
              write_macro('StudyInterpretation', interpretation), write_macro('StudyRecommendation', recommendation),
              write_macro('ValidationResults', r'Validation receipts, untouched-main parity and pinned-replay outcomes are recorded in \texttt{validation.json}; they must be read alongside the statistical cohort.')]
    (DOC/'tangent-results-generated.tex-part').write_text('% Generated from source-bound paired scores.\n'+''.join(macros))
    check_evidence()
    print(headline)


if __name__ == '__main__': main()
