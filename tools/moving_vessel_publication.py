#!/usr/bin/env python3
"""Publish the complete, passing underway replay report without rerunning it."""
import csv
import json

from sim_dataset import ROOT, ARCHIVE, SHA256, REFERENCE_CSV_SHA256
from filter_sim_gates import measured_worst
from carma_sim_table import format_table
from ou_sim_table import RECORDS

FAMILIES = {
    'ou2': 'OU-II', 'ou3': 'OU-III', 'tfg': 'TFG', 'carma': 'CARMA',
    'csm': 'Shared QMEKF', 'tcarma': 'TFG-CARMA', 'tcsm': 'Shared TFG',
    'tfg_legacy': 'TFG cubic',
}
TABLES = {
    'carma': ('kalman_carma', 'carma-sim-results-generated.tex-part', 'tab:carma_sim_results'),
    'csm': ('kalman_carma_shared_modes', 'csm-sim-results-generated.tex-part', 'tab:csm_sim_results'),
    'tcarma': ('kalman_tfg_carma', 'tfg-carma-sim-results-generated.tex-part', 'tab:tfg_carma_sim_results'),
    'tcsm': ('kalman_tfg_shared_modes', 'tfg-shared-sim-results-generated.tex-part', 'tab:tfg_shared_sim_results'),
}


def publish():
    report = json.loads((ROOT / 'reports/results/filter_sim_gates/measurements.json').read_text())
    if report.get('acceptance') != 'all replay gates passed':
        raise ValueError('publication requires a complete passing replay report')
    if report['inputs']['archive'] != ARCHIVE or report['inputs']['archive_sha256'] != SHA256:
        raise ValueError('publication report uses another dataset')
    if report['inputs']['input_sha256'] != REFERENCE_CSV_SHA256:
        raise ValueError('publication report uses different motion bytes')
    summary = []
    for family, label in FAMILIES.items():
        profiles = report['records'][family]
        for profile, rows in profiles.items():
            worst = measured_worst(rows)
            summary.append({'filter': label, 'profile': profile,
                'heave_jonswap_pct_hs': worst['err_limit_percent_z_jonswap'],
                'heave_pmstokes_pct_hs': worst['err_limit_percent_z_pmstokes'],
                'disp_3d_rms_m': max(r['disp_3d_rms_m'] for r in rows),
                'roll_rms_deg': worst['err_limit_roll_deg'],
                'pitch_rms_deg': worst['err_limit_pitch_deg'],
                'yaw_rms_deg': worst['err_limit_yaw_deg'],
                'accel_bias_3d_rms_mps2': max(r['accel_bias_3d_rms_mps2'] for r in rows),
                'gyro_bias_3d_rms_radps': max(r['gyro_bias_3d_rms_radps'] for r in rows)})
        if family in TABLES:
            directory, filename, label_id = TABLES[family]
            by_name = {r['input']: r for r in profiles['default']}
            formatted = []
            for wave, hs, name in RECORDS:
                row = by_name[name]
                formatted.append((wave, hs, {
                    'pct_hs': tuple(100 * row[f'disp_{a}_rms_m'] / hs for a in 'xyz'),
                    'd3': (row['disp_3d_rms_m'], row['disp_3d_pct_refmax']),
                    'angles': tuple(row[f'{a}_rms_deg'] for a in ('roll', 'pitch', 'yaw'))}))
            table = {'tool': 'tools/moving_vessel_publication.py', 'variant': FAMILIES[family],
                     'caption': FAMILIES[family] + r' on the 5 kn underway records (last \SI{900}{s} RMS window).',
                     'label': label_id}
            (ROOT / 'doc' / directory / filename).write_text(format_table(formatted, table))
    out = ROOT / 'reports/results/filter_sim_gates/underway_summary.csv'
    with out.open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    prose = (
        'The current external-data baseline uses oceanography-waves-lib v1.3.0: '
        'a 28 ft estimated hull RAO moving at 5 kn with heading oscillating '
        'plus/minus 3 degrees about north every 120 s. Eight fixed '
        'JONSWAP and PM-Stokes records run for 1200 s at 200 Hz; the last 900 s '
        'are scored with default sensor and initialization seeds. Displacement '
        'and velocity truth are relative to the declared cruise track. '
        'All recorded specific force, turning acceleration, changing attitude '
        'and body rates reach the estimator unchanged. The 3D denominator '
        'uses wave displacement, not forward travel. These deterministic '
        'regression measurements are not confidence intervals or sea trials. '
        'A dataset change is not an estimator improvement.'
    )
    md = ['# Underway replay results', '', prose, '',
          '| Filter | JONSWAP heave (% Hs) | PM-Stokes heave (% Hs) | 3D RMS (m) | Yaw RMS (deg) | Accel bias RMS (mm/s²) |',
          '| --- | ---: | ---: | ---: | ---: | ---: |']
    tex = [r'\section{Current Underway Replay Baseline}', prose.replace('3D', '3-D'), '',
           r'\begin{center}\small\setlength{\tabcolsep}{4pt}', r'\begin{tabular}{lrrrr}',
           r'\hline Filter & JONSWAP & PM--Stokes & 3-D & Yaw \\',
           r' & \%$H_s$ & \%$H_s$ & m & deg \\ \hline']
    for row in summary:
        if row['profile'] != 'default':
            continue
        md.append(f"| {row['filter']} | {row['heave_jonswap_pct_hs']:.3f} | {row['heave_pmstokes_pct_hs']:.3f} | {row['disp_3d_rms_m']:.4f} | {row['yaw_rms_deg']:.3f} | {1000*row['accel_bias_3d_rms_mps2']:.2f} |")
        tex.append(f"{row['filter']} & {row['heave_jonswap_pct_hs']:.2f} & {row['heave_pmstokes_pct_hs']:.2f} & {row['disp_3d_rms_m']:.3f} & {row['yaw_rms_deg']:.2f}" + r' \\')
    pii = report['other_replays']['pii']
    pii_j = max(r['disp_z_pct_hs'] for r in pii if 'jonswap' in r['input'])
    pii_p = max(r['disp_z_pct_hs'] for r in pii if 'pmstokes' in r['input'])
    pii_yaw = max(r['yaw_rms_deg'] for r in pii)
    md.append(f'| Adaptive PII | {pii_j:.3f} | {pii_p:.3f} | — | {pii_yaw:.3f} | — |')
    tex.append(f'Adaptive PII & {pii_j:.2f} & {pii_p:.2f} & -- & {pii_yaw:.2f}' + r' \\')
    tex += [r'\hline\end{tabular}\end{center}',
            'Each column reports its worst record; maxima need not come from the same sea. '
            r'Production and interval-rate limits are qualified separately with a 0.5\% '
            'margin and upward rounding. The analytic hull uses encounter corrections, '
            'without hull-specific forward-speed hydrodynamics or steering dynamics.', '']
    md += ['', 'Each column is the worst of eight records; maxima need not occur together.',
           'The compact per-record measurements and both sampling profiles are in '
           '`reports/results/filter_sim_gates/measurements.json` and `underway_summary.csv`.',
           'See [the scoring and surrogate convention](moving-vessel-replays.md).', '']
    (ROOT / 'docs/moving-vessel-results.md').write_text('\n'.join(md))
    for directory in ('kalman_ou_iii', 'kalman_tfg', 'kalman_carma',
                      'kalman_carma_shared_modes', 'kalman_tfg_carma', 'kalman_tfg_shared_modes'):
        (ROOT / 'doc' / directory / 'underway-baseline.tex-part').write_text('\n'.join(tex))


if __name__ == '__main__':
    publish()
