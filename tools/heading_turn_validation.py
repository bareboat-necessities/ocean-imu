#!/usr/bin/env python3
"""Summarize turn-model replays without replacing the existing acceptance bars.

First run heading_kinematics_diagnosis.py with --variants released, all eight
records and all eight families. This tool verifies its receipts and compares
the resulting full-rate scores with the committed pre-fix measurements.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

from heading_kinematics_diagnosis import FAMILIES, cached_trace_valid
from sim_dataset import ROOT, CASES, REFERENCE_CSV_SHA256_BY_CASE, digest


def summarize(work: Path) -> dict:
    baseline = {case: json.loads((ROOT/'reports/results/sim_v1.4.0'/case/
                                 'measurements.json').read_text()) for case in CASES}
    attribution = {}
    for name in ('measurements', 'medium'):
        report = json.loads((ROOT/'reports/results/heading_kinematics'/f'{name}.json').read_text())
        attribution.update({(r['case'], r['input'], r['family']): r['metrics']
                            for r in report['runs'] if r['variant'] == 'released'})
    rows = []
    for case in CASES:
        for name in REFERENCE_CSV_SHA256_BY_CASE[case]:
            for family in FAMILIES:
                directory = work/'runs'/case/'released'/family/Path(name).stem
                receipt = json.loads((directory/'receipt.json').read_text())
                if not cached_trace_valid(directory, receipt):
                    raise ValueError(f'{directory}: missing or changed trace; rerun the replay')
                if digest(directory/'run.log') != receipt['log_sha256']:
                    raise ValueError(f'{directory}: changed replay log')
                if receipt['input_sha256'] != REFERENCE_CSV_SHA256_BY_CASE[case][name]:
                    raise ValueError(f'{directory}: not the released input')
                if any(receipt[k] != v for k, v in
                       dict(case=case, input=name, family=family, variant='released').items()):
                    raise ValueError(f'{directory}: wrong receipt identity')
                before = next(r for r in baseline[case]['records'][family]['default']
                              if r['input'] == name)
                log = (directory/'run.log').read_text().splitlines()
                measured = [s for s in log if s.startswith('VALIDATION_METRICS ')]
                if len(measured) != 1:
                    raise ValueError(f'{directory}: expected one complete score')
                after = dict(item.split('=', 1) for item in measured[0].split()[1:])
                if (after['samples'], after['start_s'], after['window_s']) != ('180000', '300', '900'):
                    raise ValueError(f'{directory}: incomplete scoring window')
                errors = [s for s in log if s.startswith('ERROR:')]
                if receipt['exit_code'] not in (0, 1) or bool(errors) != bool(receipt['exit_code']):
                    raise ValueError(f'{directory}: exit status disagrees with the quality gates')
                pair = dict(case=case, input=name, family=family,
                            input_sha256=receipt['input_sha256'],
                            binary_sha256=receipt['binary_sha256'],
                            trace_sha256=receipt['trace_sha256'], log_sha256=receipt['log_sha256'],
                            quality_gate_pass=receipt['exit_code'] == 0, gate_failures=errors)
                for metric in ('disp_x_rms_m', 'disp_y_rms_m', 'disp_z_rms_m',
                               'roll_rms_deg', 'pitch_rms_deg', 'yaw_rms_deg',
                               'accel_bias_3d_rms_mps2', 'gyro_bias_3d_rms_radps'):
                    pair[metric+'_before'] = before[metric]
                    pair[metric+'_after'] = float(after[metric])
                pair['y_reduction_percent'] = 100*(1-float(after['disp_y_rms_m'])/before['disp_y_rms_m'])
                pair.update(receipt['metrics'])
                old = attribution.get((case, name, family))
                if old is not None:
                    old_amp = old['acc_bias_y_10s_amplitude']
                    pair['bias_y_10s_amplitude_before'] = old_amp
                    pair['bias_y_10s_reduction_percent'] = 100*(1-pair['acc_bias_y_10s_amplitude']/old_amp)
                rows.append(pair)
    return dict(
        role='Runtime comparison; existing quality-gate failures are retained',
        baseline_revision='5426da8a394b33f5f473f70faa05247f832f85fa',
        score='Final 900 s, all 180000 samples at 200 Hz; unchanged references and sensor draws',
        frames=dict(displacement='geographic NWU', bias='sensor body after (x,y,z) -> (y,x,-z)'),
        runs=rows,
        all_y_improved=all(r['y_reduction_percent'] > 0 for r in rows),
        quality_gate_passes=sum(r['quality_gate_pass'] for r in rows),
        quality_gate_total=len(rows))


def plot_tails(work: Path, baseline_traces: Path, destination: Path) -> None:
    """Display the unsmoothed bias offset as well as the steering harmonic."""
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT/'plots'))
    from vector_export import save_figure
    name = 'wave_data_jonswap_H0.270_L14.047_A30.00_P60.00.csv'
    old_report = json.loads((ROOT/'reports/results/heading_kinematics/measurements.json').read_text())
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 5.8), sharex=True, sharey='row')
    keep = ['time', 'disp_ref_y', 'disp_est_y', 'acc_bias_y', 'acc_bias_est_y']
    for col, case in enumerate(CASES):
        old = baseline_traces/case/'released'/'ou3'/Path(name).stem
        new = work/'runs'/case/'released'/'ou3'/Path(name).stem
        old_receipt = next(r for r in old_report['runs'] if
                           (r['case'], r['input'], r['family'], r['variant']) ==
                           (case, name, 'ou3', 'released'))
        traces = []
        for directory, receipt in ((old, old_receipt),
                                   (new, json.loads((new/'receipt.json').read_text()))):
            if not cached_trace_valid(directory, receipt):
                raise ValueError(f'{directory}: changed plot source')
            path = next(directory.glob('w3d*.csv'))
            with path.open() as stream:
                columns = stream.readline().strip().split(',')
            values = np.loadtxt(path, delimiter=',', skiprows=1,
                                usecols=[columns.index(k) for k in keep])
            traces.append(values[values[:, 0] >= 1140][::10])
        data = np.column_stack((traces[1][:, :2], traces[0][:, 2], traces[1][:, 2],
                                traces[1][:, 3], traces[0][:, 4], traces[1][:, 4]))
        np.savetxt(destination/f'runtime-{case}-small-ou3-tail.csv', data,
                   delimiter=',', comments='', fmt='%.9g',
                   header='time,disp_true_y,disp_before_y,disp_after_y,bias_true_y,bias_before_y,bias_after_y')
        for row, truth, estimate in ((0, 1, 2), (1, 3, 4)):
            ax = axes[row, col]
            ax.plot(traces[1][:, 0], traces[1][:, truth], lw=1.3, color='#263445', label='True')
            ax.plot(traces[0][:, 0], traces[0][:, estimate], lw=1, color='#bc652a', label='Before')
            ax.plot(traces[1][:, 0], traces[1][:, estimate], lw=1.2, color='#167daf', label='Turn model')
            ax.spines[['top', 'right']].set_visible(False)
            ax.grid(alpha=.18)
        axes[0, col].set_title(case.capitalize()+' seas')
        axes[1, col].set_xlabel('Time [s]')
    axes[0, 0].set_ylabel('Wave displacement Y [m]')
    axes[1, 0].set_ylabel('Body accelerometer bias Y [m/s²]')
    axes[0, 0].legend(frameon=False, ncol=3, loc='upper left')
    fig.suptitle('OU-III · 0.27 m incident waves · same IMU and noise\n'
                 'Steering response removed; slow bias offset remains', fontsize=12)
    fig.tight_layout()
    figure = destination/'runtime-small-ou3.svg'
    with matplotlib.rc_context({'svg.hashsalt': 'heading-turn-runtime'}):
        save_figure(fig, figure, metadata={'Date': None})
    figure.write_text('\n'.join(line.rstrip() for line in figure.read_text().splitlines())+'\n')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, required=True)
    parser.add_argument('--runtime-revision', required=True,
                        help='Source revision used to build the measured runtime; not inferred from HEAD')
    parser.add_argument('--build-description', required=True,
                        help='Compiler and effective flags used for the comparison executables')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--baseline-traces', type=Path,
                        help='Optional pre-fix dense traces, used only for the before/after figure')
    args = parser.parse_args()
    report = summarize(args.work_dir.resolve())
    report['runtime_revision'] = args.runtime_revision
    report['build_description'] = args.build_description
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    fields = ['case', 'input', 'family', 'quality_gate_pass', 'y_reduction_percent',
              'disp_y_rms_m_before', 'disp_y_rms_m_after', 'disp_z_rms_m_before',
              'disp_z_rms_m_after', 'accel_bias_3d_rms_mps2_before',
              'accel_bias_3d_rms_mps2_after', 'bias_y_10s_amplitude_before',
              'acc_bias_y_10s_amplitude', 'bias_y_10s_reduction_percent']
    with args.output.with_suffix('.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore', lineterminator='\n')
        writer.writeheader()
        writer.writerows(report['runs'])
    if args.baseline_traces:
        plot_tails(args.work_dir.resolve(), args.baseline_traces.resolve(), args.output.parent)
    print(f"Y improved in {sum(r['y_reduction_percent'] > 0 for r in report['runs'])}/{len(report['runs'])} cells; "
          f"existing quality gates pass in {report['quality_gate_passes']}/{report['quality_gate_total']}.")


if __name__ == '__main__':
    main()
