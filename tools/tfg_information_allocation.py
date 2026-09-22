#!/usr/bin/env python3
"""Passive TFG accelerometer information-allocation replay after core refit.

This is a diagnostic, not another tuning campaign. It compares retained defaults,
the sealed challenger and the already-rejected interior point on the same eight
records and one previously used refinement seed. It never opens new holdout data
or selects/promotes coefficients. H_j K_j and H_j K_j r share residual units.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import re

from model_mismatch_ablation import FAMILIES, RECORDS
from tfg_rao_core_refit import KEYS, arm
from tfg_rao_refit_replay import execute, write_json

ROOT = Path(__file__).resolve().parents[1]


def analyse(directory: Path) -> dict:
    rows = json.loads((directory / 'replay/runs.json').read_text())
    records = {r.filename: r for r in RECORDS}
    output = []
    for row in rows:
        record = records[row['input']]
        path = directory / 'replay' / f"{row['config']}_{record.spectrum}_{record.hs_m}_{row['seed']}.log"
        text = path.read_text()
        lines = re.findall(r'^TFG_ALLOCATION_SUMMARY (.*)$', text, re.M)
        if len(lines) != 1:
            raise ValueError(f'expect exactly one passive summary: {path}')
        values = {key: float(value) for token in lines[0].split() for key, value in [token.split('=', 1)]}
        if not all(math.isfinite(value) for value in values.values()):
            raise ValueError(f'non-finite allocation: {path}')
        if values['invalid_samples'] != 0 or not 17999 <= values['samples'] <= 18001:
            raise ValueError(f'incomplete 900-second diagnostic: {path}')
        if values['live'] != 1 or values['accepted'] != 1:
            raise ValueError(f'allocation denominator requires accepted Live samples: {path}')
        total_trace = sum(values[block + '_trace'] for block in ('theta', 'aw', 'ba'))
        component_energy = sum(values[block + '_energy'] for block in ('theta', 'aw', 'ba'))
        reconstructed = component_energy + 2 * sum(values[key] for key in ('theta_aw_cross', 'theta_ba_cross', 'aw_ba_cross'))
        if abs(reconstructed - values['net_energy']) > 1e-11 * max(1.0, abs(reconstructed)):
            raise ValueError('component/cross-term energy identity failed')
        if values['closure_sq'] > 1e-20:
            raise ValueError('H K block-decomposition closure failed')
        output.append(dict(config=row['config'], input=row['input'], seed=row['seed'],
            coefficients=row['env'], physical_metrics=row['metrics'], violations=row['violations'],
            log_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), allocation=values,
            signed_trace_shares={block: values[block + '_trace'] / total_trace for block in ('theta', 'aw', 'ba')},
            component_energy_shares={block: values[block + '_energy'] / component_energy for block in ('theta', 'aw', 'ba')},
            net_to_component_energy=values['net_energy'] / component_energy))
    result = dict(records=8, seed='31', arms=3, completed_replays=len(rows),
        sample_rate_hz=20, scoring_window_s=900, rows=output,
        interpretation={
            'operator': 'G_j = H_j K_j for j = theta, a_w, b_a from the retained, actual accelerometer update.',
            'energy': 'u_j = G_j r in m/s^2. Report individual energies AND signed cross terms; individual energies do not generally add to the net correction energy.',
            'signed_trace': 'Signed shares of trace(H K), not probabilities and not percentages of Fisher information.',
            'rank': 'The instantaneous 3-row accelerometer relation cannot independently identify all nine theta/a_w/b_a coordinates. This is NOT a claim of temporal unobservability.',
            'scope': 'Passive accelerometer allocation only. No new claim about magnetic or integral-update information, reset defects, or pre/post-reset covariance loss is made.',
            'selection': 'No coefficient is selected from this diagnostic. The interior arm was rejected before the sealed holdout.'})
    write_json(directory / 'allocation.json', result)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'reports/results/tfg_attitude_bias_diagnostics')
    parser.add_argument('--core-dir', type=Path, default=ROOT / 'reports/results/tfg_rao_core_refit')
    parser.add_argument('--binary', type=Path, default=FAMILIES['TFG'])
    parser.add_argument('--jobs', type=int, default=4)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    directory = args.output_dir.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    if args.execute:
        verdict_path = args.core_dir / 'holdout-verdict.json'
        if not verdict_path.exists():
            raise ValueError('complete the core refit before allocation diagnostics')
        verdict = json.loads(verdict_path.read_text())
        if verdict['holdout_passed']:
            raise ValueError('this fixed diagnostic design assumes the completed nonpromotion decision')
        frozen = json.loads((args.core_dir / 'candidate.json').read_text())
        configs = [arm('baseline'), arm('sealed_challenger', *(frozen['coefficients'][key] for key in KEYS)),
                   arm('rejected_interior', 1.03, 1.10, .75, 1.30)]
        for config in configs:
            config['env']['TFG_ALLOCATION_TRACE'] = '1'
        write_json(directory / 'design.json', dict(configs=configs, seeds=['31'], record_indices=list(range(8)),
            coefficient_search=False, binary_path=str(args.binary),
            core_verdict_sha256=hashlib.sha256(verdict_path.read_bytes()).hexdigest()))
        FAMILIES['TFG'] = args.binary.resolve()
        execute(configs, 'TFG', list(range(8)), ['31'], directory / 'replay', args.jobs, directory / 'cache')
        control = execute([arm('baseline')], 'TFG', [0], ['31'], directory / 'parity-off', 1, directory / 'cache')[0]
        observed = next(row for row in json.loads((directory / 'replay/runs.json').read_text())
                        if row['config'] == 'baseline' and row['input'] == control['input'])
        if control['metrics'] != observed['metrics'] or control['violations'] != observed['violations']:
            raise ValueError('passive logging changed physical metrics or gate results')
        write_json(directory / 'passivity.json', dict(same_binary_trace_on_off=True,
            physical_metrics_identical=True, quality_gate_results_identical=True,
            input=control['input'], seed='31', binary_sha256=control['binary_sha256']))
    result = analyse(directory)
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=2))
    for row in result['rows']:
        print(row['config'], row['input'], row['signed_trace_shares'], 'rS floor=', row['allocation']['rs_floor'])


if __name__ == '__main__':
    main()
