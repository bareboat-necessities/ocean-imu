#!/usr/bin/env python3
"""Run the external-record studies for one pinned direction in isolated outputs.

CI crosses STUDIES with following/head. Synthetic-only Lie-group, CARMA prior,
calibration and stability experiments retain their own physical scenarios.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from sim_dataset import CASE, CASES, RELEASE, ROOT, REFERENCE_CSV_SHA256, input_provenance, digest

OU = [('kalman_ou_ii', 'kalman_ou_ii-sim'), ('kalman_ou_iii', 'kalman_ou_iii-sim')]
TFG = [('kalman_tfg', 'kalman_tfg-sim')]
CARMA = [('kalman_carma', 'carma3-sim')]
CSM = [('kalman_carma_shared_modes', 'kalman_carma_shared_modes-sim')]
TCSM = [('kalman_tfg_shared_modes', 'kalman_tfg_shared_modes-sim')]
OBS = [('pii_observer', 'pii_observer-adaptive'), ('nlo', 'nlo-sim')]
# name: executables required by the existing study implementation
STUDIES = {
    'ou_validation': OU,
    'ou_robustness': OU,
    'tfg_comparison': OU + TFG,
    'model_mismatch_ablation': OU + TFG,
    'engine_noise_degradation': OU + TFG,
    'ou3_engine_noise_mitigation': OU[1:],
    'ou_roundtrip_transition': OU,
    'startup_ablation': OU + TFG,
    'startup_transient_study': OU + TFG + OBS,
    'direction_rao_ablation': OU,
    'ou3_lever_arm_study': OU[1:],
    'carma_lever_arm_study': [('kalman_carma', 'carma_lever_arm-sim')],
    'carma3_benchmark': OU[1:] + CARMA,
    'tfg_carma3_benchmark': OU[1:] + CARMA,
    'carma_excitation_study': CARMA,
    'carma_shared_modes_study': CSM,
    'carma_shared_modes_validation': OU + TFG + CSM + OBS,
    'tfg_shared_modes_dataset': CARMA + CSM + TCSM,
    'live_bias_csv_replay': CSM + TCSM,
    'live_bias_csv_timing': [],
    'heave_only_comparison': OU + TFG + [('heave_only', 'heave_only-sim')],
}
DATASET_STUDIES = {'carma3_benchmark': 'carma3_benchmark',
                   'tfg_carma3_benchmark': 'tfg_carma3',
                   'carma_excitation_study': 'carma_excitation',
                   'carma_shared_modes_study': 'kalman_carma_shared_modes'}


def run_commands(commands, env, out):
    """Preserve diagnostics when a legacy study replaces its results directory."""
    if not commands:
        raise ValueError('at least one study command is required')
    # A sibling log survives output cleanup and remains available if CI
    # interrupts the job before the child returns.
    out.parent.mkdir(parents=True, exist_ok=True)
    log_path = out.parent/f'{out.name}.log'
    with log_path.open('w', encoding='utf-8') as log:
        for command in commands:
            try:
                result = subprocess.run(command, env=env, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
            except OSError as error:
                log.write(f'Unable to launch {command!r}: {error}\n')
                result = subprocess.CompletedProcess(command, 127)
            if result.returncode:
                break
        out.mkdir(parents=True, exist_ok=True)
        log.flush()
        shutil.copyfile(log_path, out/'study.log')
    return result


def command_for(study, data, out, jobs, mode):
    command = [sys.executable, '-u', str(ROOT/'tools'/f'{study}.py')]
    common = ['--data-dir', str(data), '--jobs', str(jobs)]
    if study == 'live_bias_csv_replay':
        return command + ['--all-records', '--jobs', str(jobs), '--output-dir', str(out)]
    if study == 'live_bias_csv_timing':
        return command + ['--output-dir', str(out)]
    if study in ('ou_validation', 'ou_robustness', 'tfg_comparison'):
        return command + common + ['--mode', mode, '--output-dir', str(out), '--skip-build']
    if study in DATASET_STUDIES:
        return command + ['dataset', '--jobs', str(jobs)]
    if study == 'carma_shared_modes_validation':
        return command + ['all', *common, '--results', str(out), '--outdir', str(out/'article'),
                          '--seeds', '10' if mode == 'full' else '1']
    if study == 'tfg_shared_modes_dataset':
        return command + common + ['--results', str(out), '--seeds', '8' if mode == 'full' else '1']
    if study == 'ou_roundtrip_transition':
        return command + ['--output-dir', str(out)]
    if study == 'startup_transient_study':
        return command + common + ['--output-dir', str(out), '--work-dir', str(out/'work')]
    if study in ('ou3_lever_arm_study', 'carma_lever_arm_study'):
        return command + common + ['--mode', mode, '--output-dir', str(out)]
    return command + common + ['--output-dir', str(out)]


def dataset_chart(out):
    """A dataset-only report must not import synthetic results from another run."""
    sys.path.insert(0, str(ROOT/'plots'))
    from vector_export import save_figure
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    replay = out/'csv_rollout_all.json'
    if replay.exists():
        rows = []
        for run in json.loads(replay.read_text()):
            for line in run['metrics']:
                if line.startswith('VALIDATION_METRICS '):
                    fields = dict(pair.split('=', 1) for pair in line.split()[1:])
                    rows.append(dict(input=run['input'], family=run['filter']+'/'+run['contract'],
                                     disp_z_rms_m=float(fields['disp_z_rms_m'])))
        if len(rows) != 48:
            raise ValueError('expected 48 complete paired bias/timing replay rows')
        pd.DataFrame(rows).to_csv(out/'dataset_runs.csv', index=False)
    timing = out/'csv_timing.json'
    if timing.exists():
        data = pd.DataFrame(json.loads(timing.read_text()))
        data.to_csv(out/'csv_timing.csv', index=False)
        columns = ['gyro_endpoint_increment_rms', 'gyro_centered_increment_rms']
        ax = data.set_index('file')[columns].plot.bar(figsize=(13, 6), ylabel='Gyro increment residual RMS (rad/s)',
                                                     title=f'{RELEASE} · {CASES[CASE]["description"]}')
        ax.set_xticklabels([str(n).replace('wave_data_', '').replace('.csv', '') for n in data.file], rotation=30, ha='right', fontsize=8)
        ax.figure.tight_layout()
        save_figure(ax.figure, out/'csv_timing.svg')
        plt.close(ax.figure)
    files = list(out.rglob('dataset_runs.csv')) + list(out.rglob('dataset_summary.csv'))
    for path in files:
        data = pd.read_csv(path)
        group = next((n for n in ('family', 'variant', 'method') if n in data), None)
        sea = next((n for n in ('input', 'record') if n in data), None)
        if not group or not sea or 'disp_z_rms_m' not in data:
            continue
        table = data.pivot_table(index=sea, columns=group, values='disp_z_rms_m', aggfunc='mean')
        table.to_csv(path.with_name('heave_mean_m.csv'))
        ax = table.plot.bar(figsize=(13, 6), ylabel='Heave RMS (m)',
                            title=f'{RELEASE} · {CASES[CASE]["description"]}')
        ax.set_xticklabels([str(n).replace('wave_data_', '').replace('.csv', '') for n in table.index], rotation=30, ha='right', fontsize=8)
        ax.figure.tight_layout()
        save_figure(ax.figure, path.with_name('heave_mean_m.svg'))
        plt.close(ax.figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study', choices=STUDIES)
    parser.add_argument('--jobs', type=int, default=2)
    parser.add_argument('--mode', choices=('smoke', 'full'), default='full')
    parser.add_argument('--skip-build', action='store_true')
    parser.add_argument('--output-root', type=Path, default=ROOT/'reports/results/sim_v1.4.0')
    args = parser.parse_args()
    out = args.output_root.resolve()/CASE/'studies'/args.study
    out.mkdir(parents=True, exist_ok=True)
    data = ROOT/'.sim-data/cases'/CASE
    targets = STUDIES[args.study]
    destinations = [data]
    if args.study != 'heave_only_comparison':
        destinations += [ROOT/'plots/kalman_ou_ii'] + [ROOT/'tests'/d for d,_ in targets]
    subprocess.run([sys.executable, str(ROOT/'tools/sim_dataset.py'), '--case', CASE,
                    '--dest', *map(str, destinations)], check=True)
    if not args.skip_build:
        for directory, binary in targets:
            subprocess.run(['make', '-C', str(ROOT/'tests'/directory), binary], check=True)
    env = {**os.environ, 'OCEAN_IMU_SIM_CASE': CASE,
           'OCEAN_IMU_STUDY_ROOT': str(out), 'OPENBLAS_NUM_THREADS': '1'}
    command = command_for(args.study, data, out, args.jobs, args.mode)
    receipt = dict(study=args.study, case=CASE, mode=args.mode,
                   inputs=input_provenance(data/n for n in REFERENCE_CSV_SHA256),
                   source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                   binary_sha256={f'{d}/{b}': digest(ROOT/'tests'/d/b) for d,b in targets},
                   command=command, status='running')
    receipt_path = out/'case-receipt.json'
    receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
    begin = time.monotonic()
    commands = [command]
    def tool(name, *arguments):
        return [sys.executable, str(ROOT/'tools'/f'{name}.py'), *map(str, arguments)]
    if args.study == 'startup_ablation':
        commands += [tool('startup_publication_sync', '--results', out, '--doc-dir', out/'article'),
                     tool('ou2_dual_mse_coefficients', '--exact', '--data-dir', data,
                          '--operating-points', out/'startup_runs.csv', '--output-dir', out/'ou2_dual_mse'),
                     tool('ou3_reduced_mse_study', '--data-dir', data,
                          '--operating-points', out/'startup_runs.csv', '--output-dir', out/'ou3_reduced_mse',
                          '--doc-dir', out/'article')]
    elif args.study == 'direction_rao_ablation':
        commands.append(tool('direction_rao_report', '--results', out, '--doc-dir', out/'article'))
    elif args.study == 'engine_noise_degradation':
        commands.append(tool('engine_noise_guard_families', '--data-dir', data, '--jobs', args.jobs,
                             '--degradation-summary', out/'engine_noise_summary.csv', '--output-dir', out))
    receipt['commands'] = commands
    result = run_commands(commands, env, out)
    receipt.update(exit_code=result.returncode, seconds=time.monotonic()-begin,
                   status='postprocessing' if result.returncode == 0 else 'failed', log_sha256=digest(out/'study.log'))
    receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
    if result.returncode:
        print((out/'study.log').read_text()[-6000:], file=sys.stderr)
        return result.returncode
    try:
        if args.study in DATASET_STUDIES or args.study in ('tfg_shared_modes_dataset', 'live_bias_csv_replay', 'live_bias_csv_timing'):
            dataset_chart(out)
        after = {f'{d}/{b}': digest(ROOT/'tests'/d/b) for d,b in targets}
        if after != receipt['binary_sha256']:
            raise ValueError('study executable changed during the run')
    except Exception as error:
        receipt.update(status='failed', exit_code=1, error=str(error))
        receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
        raise
    receipt.update(status='completed', seconds=time.monotonic()-begin)
    receipt_path.write_text(json.dumps(receipt, indent=2)+'\n')
    print(f'{CASE}/{args.study}: completed; outputs in {out}', flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
