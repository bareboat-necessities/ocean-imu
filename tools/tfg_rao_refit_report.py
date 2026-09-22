#!/usr/bin/env python3
"""Report completed RAO phases and enforce the predeclared sealed holdout.

No selection occurs here. Pooled RMS, mean case RMS, balanced geometric changes,
per-record results and per-seed consistency remain distinct. Failed cases remain
visible and make a candidate ineligible rather than disappearing from a score.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import numpy as np
from tfg_rao_core_refit import PROTOCOL, KEYS, LEGACY, cases
from tfg_rao_refit_analysis import analyse, METRICS, WEIGHTS
from tfg_rao_refit_replay import write_json

ROOT = Path(__file__).resolve().parents[1]
ALL_METRICS = ('disp_z_rms_m', *METRICS, 'accel_3d_rms_mps2')


def compare(rows, candidate, baseline):
    groups = cases(rows)
    current, reference = groups[candidate], groups[baseline]
    if set(current) != set(reference):
        raise ValueError('unpaired comparison histories')
    keys = sorted(reference)
    inputs = {k[0] for k in keys}
    seeds = {k[1] for k in keys}
    if len(inputs) != 8 or len(keys) != len(inputs) * len(seeds):
        raise ValueError('comparison requires the complete eight-record cross product')
    output = dict(candidate=candidate, reference=baseline, cases=len(keys), metrics={})
    invalid = [k for k in keys if current[k].get('error') or reference[k].get('error') or
               any(not isinstance(r['metrics'].get(m), (int, float)) or
                   not math.isfinite(r['metrics'][m]) or r['metrics'][m] <= 0
                   for r in (current[k], reference[k]) for m in ALL_METRICS)]
    output.update(invalid_cases=invalid,
        same_input_hashes=all(current[k].get('input_sha256') is not None and
                             current[k]['input_sha256'] == reference[k].get('input_sha256') for k in keys),
        new_gate_failure_cases=[k for k in keys if current[k]['exit_code'] != 0 and reference[k]['exit_code'] == 0],
        extra_gate_violations=sum(max(0, len(current[k]['violations'])-len(reference[k]['violations'])) for k in keys),
        candidate_gate_failed_cases=sum(current[k]['exit_code'] != 0 for k in keys),
        reference_gate_failed_cases=sum(reference[k]['exit_code'] != 0 for k in keys),
        candidate_gate_violations=sum(len(current[k]['violations']) for k in keys),
        reference_gate_violations=sum(len(reference[k]['violations']) for k in keys))
    def gate_names(row):
        return {re.sub(r'\([^\n]*', '', v['message']).strip() for v in row['violations']}
    output['new_individual_gate_types'] = [dict(input=k[0], seed=k[1],
        gates=sorted(gate_names(current[k])-gate_names(reference[k]))) for k in keys
        if gate_names(current[k])-gate_names(reference[k])]
    if invalid:
        return output
    logratios = []
    for metric in ALL_METRICS:
        a = np.array([current[k]['metrics'][metric] for k in keys])
        b = np.array([reference[k]['metrics'][metric] for k in keys])
        logs = np.log(a/b)
        def split(part, pooled):
            result = {}
            for value in sorted({k[part] for k in keys}):
                index = [i for i, k in enumerate(keys) if k[part] == value]
                aa, bb = a[index], b[index]
                ratio = np.sqrt(np.mean(aa*aa))/np.sqrt(np.mean(bb*bb)) if pooled else aa.mean()/bb.mean()
                result[value] = float(100*(ratio-1))
            return result
        records, means = split(0, True), split(0, False)
        output['metrics'][metric] = dict(candidate_mean=float(a.mean()), reference_mean=float(b.mean()),
            mean_rms_change_pct=float(100*(a.mean()/b.mean()-1)),
            candidate_pooled_rms=float(np.sqrt(np.mean(a*a))), reference_pooled_rms=float(np.sqrt(np.mean(b*b))),
            pooled_change_pct=float(100*(np.sqrt(np.mean(a*a))/np.sqrt(np.mean(b*b))-1)),
            balanced_geometric_change_pct=float(100*np.expm1(logs.mean())),
            case_wins=int((a < b).sum()), record_wins=sum(v < -1e-7 for v in records.values()),
            per_record_change_pct=records, per_record_mean_rms_change_pct=means,
            mean_rms_record_wins=sum(v < -1e-7 for v in means.values()), per_seed_change_pct=split(1, True),
            worst_paired_change_pct=float(100*(a/b-1).max()))
        if metric in METRICS:
            logratios.append(logs)
    composite = np.array(logratios).T @ WEIGHTS / WEIGHTS.sum()
    output['balanced_composite_change_pct'] = float(100*np.expm1(composite.mean()))
    output['composite_case_wins'] = int((composite < 0).sum())
    output['composite_per_record_change_pct'] = {name: float(100*np.expm1(composite[
        [i for i, k in enumerate(keys) if k[0] == name]].mean())) for name in sorted(inputs)}
    output['composite_per_seed_change_pct'] = {seed: float(100*np.expm1(composite[
        [i for i, k in enumerate(keys) if k[1] == seed]].mean())) for seed in sorted(seeds)}
    return output


def holdout_verdict(rows, candidate):
    data = compare(rows, 'candidate', 'baseline')
    frozen = candidate['coefficients']
    digest = hashlib.sha256(json.dumps(frozen, sort_keys=True).encode()).hexdigest()
    checks = dict(frozen_coefficient_hash=digest == candidate['coefficient_sha256'],
        coefficients_unchanged=all(tuple(float(r['env'][k]) for k in KEYS) ==
            (tuple(frozen[k] for k in KEYS) if r['config'] == 'candidate' else LEGACY)
            for r in rows if r['config'] in ('baseline', 'candidate')),
        gyro_bias_process_noise_retained=all(float(r['env']['SF_GYRO_BIAS_RW_VAR']) == PROTOCOL['fixed_gyro_bias_rw_var'] for r in rows),
        exact_sealed_holdout_seeds={str(r['seed']) for r in rows} == set(PROTOCOL['sealed_holdout_seeds']),
        all_cases_complete_and_finite=not data['invalid_cases'] and data['cases'] == 8*len(PROTOCOL['sealed_holdout_seeds']),
        identical_pinned_inputs=data['same_input_hashes'])
    limits = PROTOCOL['promotion_checks']
    checks['no_new_gate_failure_cases'] = len(data['new_gate_failure_cases']) <= limits['new_gate_failure_cases_vs_tfg_baseline']
    checks['no_additional_gate_violations'] = data['extra_gate_violations'] <= limits['additional_gate_violations_vs_tfg_baseline']
    for metric in limits['primary_metrics']:
        result = data['metrics'].get(metric, {})
        for field, label in [('pooled_change_pct', 'pooled'), ('mean_rms_change_pct', 'mean_rms')]:
            checks[metric+'_'+label+'_nonregression'] = result.get(field, math.inf) <= limits['maximum_pooled_primary_metric_regression_pct']
        for field, label in [('record_wins', 'record'), ('mean_rms_record_wins', 'mean_rms_record')]:
            checks[metric+'_'+label+'_consistency'] = result.get(field, 0) >= limits['minimum_record_wins_for_each_primary_metric']
    for metric in limits['secondary_metrics']:
        for field, label in [('pooled_change_pct', 'pooled'), ('mean_rms_change_pct', 'mean_rms')]:
            checks[metric+'_'+label+'_guard'] = data['metrics'].get(metric, {}).get(field, math.inf) <= limits['maximum_pooled_secondary_metric_regression_pct']
    checks['minimum_balanced_improvement'] = data.get('balanced_composite_change_pct', math.inf) <= -limits['minimum_balanced_composite_improvement_pct']
    return dict(holdout_passed=all(checks.values()), checks=checks,
                failed_checks=[k for k, v in checks.items() if not v], comparison=data,
                shipping_default_promotion_requires='Additionally pass the unchanged compiled-default deterministic suite.')


def metric_table(data):
    lines = ['| Metric | Reference pooled RMS | Candidate pooled RMS | Pooled change | Mean-RMS change | Balanced change | Record wins | Case wins |',
             '|---|---:|---:|---:|---:|---:|---:|---:|']
    for metric, d in data['metrics'].items():
        lines.append(f"| {metric} | {d['reference_pooled_rms']:.8g} | {d['candidate_pooled_rms']:.8g} | {d['pooled_change_pct']:+.3f}% | {d['mean_rms_change_pct']:+.3f}% | {d['balanced_geometric_change_pct']:+.3f}% | {d['record_wins']}/8 | {d['case_wins']}/{data['cases']} |")
    return lines + ['']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=ROOT/'reports/results/tfg_rao_core_refit')
    root = parser.parse_args().output_dir
    document = ['# TFG RAO coefficient-refit evidence', '',
        'All eight pinned 28-ft vessel-RAO records; Q_bg = 3.75e-10. Negative percentage changes mean lower error. '
        'Pooled RMS is sqrt(mean(case_RMS^2)) over equal-duration histories. Mean case RMS and equal-history geometric changes are also reported. '
        'No failed candidate is discarded and no quality gate is loosened.', '',
        'Training seeds: '+', '.join(PROTOCOL['training_seeds'])+'. Refinement seeds: '+', '.join(PROTOCOL['refinement_seeds'])+'. '
        'Sealed holdout seeds: '+', '.join(PROTOCOL['sealed_holdout_seeds'])+'.', '']
    report = {'protocol': PROTOCOL, 'phases': {}}
    for phase in ('broad', 'extend', 'refine', 'fine', 'regime-check', 'interior-check', 'training-confirmation', 'feasibility', 'constraint-check'):
        path = root/phase/'runs.json'
        if not path.exists():
            continue
        rows = json.loads(path.read_text())
        if not any(r['config'] == 'baseline' for r in rows):
            continue
        result = analyse(rows)
        write_json(root/phase/'analysis.json', result)
        ranked = result['ranking']
        report['phases'][phase] = dict(cases=len(rows), configurations=len(ranked), winner=ranked[0],
            bounds=result['bounds'], boundary_contacts=result['top_boundary_contacts'],
            coarse_joint_face_contacts=result['coarse_joint_face_contacts'],
            regime_winners=result['regime_winners'], surface=result['surface_diagnostic'])
        document += [f'## {phase}: {len(ranked)} configurations, {len(rows)} case results (including reused controls)', '',
            '| Configuration | c_tau | c_sigma | k_Sx | k_Sy | Physical score | Paired gate penalty | Gate failures |',
            '|---|---:|---:|---:|---:|---:|---:|---:|']
        for row in ranked[:10]:
            if row['admissible']:
                document.append('| '+row['config']+' | '+' | '.join(f'{x:.6g}' for x in row['parameters'])+
                    f" | {row['raw_score']:+.6f} | {row['penalty']:.3f} | {row['failed_cases']} |")
            else:
                document.append('| '+row['config']+' | ineligible; retained in full evidence |')
        document += ['', f"Full per-metric response surface: `{phase}/response_surface.csv`; every case and failure: `{phase}/runs.json`.", '',
            'Exploratory quadratic surface R²: '+f"{result['surface_diagnostic']['r_squared']:.5f}. This fit summarizes interactions; it does not certify global optimality.", '',
            'Low-/high-sea physical-score winners: '+str(result['regime_winners'])+'.', '']
        iso = [r for r in ranked if r['admissible'] and math.isclose(r['parameters'][2], r['parameters'][3])]
        ani = [r for r in ranked if r['admissible'] and not math.isclose(r['parameters'][2], r['parameters'][3])]
        if iso and ani:
            document += [f"Best symmetric score: {iso[0]['score']:+.6f}; best asymmetric score: {ani[0]['score']:+.6f} ({ani[0]['config']}).", '']
        with (root/phase/'response_surface.csv').open('w', newline='') as handle:
            columns = ['config', *KEYS, 'score', 'raw_score', 'extra_gate_violations', *METRICS]
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            for row in ranked:
                output = dict(config=row['config'], **dict(zip(KEYS, row['parameters'])),
                    score=row.get('score'), raw_score=row.get('raw_score'), extra_gate_violations=row.get('extra_gate_violations'))
                output.update({m: row.get('metrics', {}).get(m, {}).get('balanced_geomean_change_pct') for m in METRICS})
                writer.writerow(output)
    frozen, holdout = root/'candidate.json', root/'holdout'/'runs.json'
    if frozen.exists():
        candidate = json.loads(frozen.read_text())
        report['candidate'] = candidate
        document += ['## Frozen candidate', '', json.dumps(candidate['coefficients'], sort_keys=True), '',
            'Coefficient SHA-256: `'+candidate['coefficient_sha256']+'`.', '']
    if frozen.exists() and holdout.exists():
        rows = json.loads(holdout.read_text())
        verdict = holdout_verdict(rows, candidate)
        report['holdout'] = verdict
        write_json(root/'holdout-verdict.json', verdict)
        document += ['## Sealed holdout', '',
            '**'+('PASS' if verdict['holdout_passed'] else 'REJECT FOR PROMOTION')+'**. Failed checks: '+', '.join(verdict['failed_checks'])+'.', '']
        document += metric_table(verdict['comparison'])
        document += [f"Gate-failing cases: baseline {verdict['comparison']['reference_gate_failed_cases']}, candidate {verdict['comparison']['candidate_gate_failed_cases']}; additional paired violations {verdict['comparison']['extra_gate_violations']}.", '']
        other = root/'holdout-ou3'/'runs.json'
        if other.exists():
            all_rows = rows + json.loads(other.read_text())
            for name in ('candidate', 'baseline'):
                comparison = compare(all_rows, name, 'ou3_shipping')
                report[name+'_vs_ou3'] = comparison
                document += ['### '+name+' versus unmodified shipping OU-III', ''] + metric_table(comparison)
        document += ['### Per-record and per-seed consistency', '',
            'Complete per-record and per-seed changes, individual gate changes and paired win counts are in `holdout-verdict.json` and `report.json`. '
            'These are eight fixed RAO motion histories with independent sensor/initialization seeds, not validation on independent vessels or independently generated wave histories.', '']
    else:
        document += ['## Holdout', '', 'Not evaluated. No new shipping coefficients are validated by this report.', '']
    write_json(root/'report.json', report)
    (root/'README.md').write_text('\n'.join(document)+'\n')


if __name__ == '__main__':
    main()
