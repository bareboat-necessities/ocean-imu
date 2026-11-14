#!/usr/bin/env python3
"""Publish the complete, passing underway replay report without rerunning it."""
import csv
import json

from sim_dataset import ROOT, ARCHIVE, SHA256, REFERENCE_CSV_SHA256
from filter_sim_gates import measured_worst
from carma_sim_table import format_table
from ou_sim_table import RECORDS
from sim_dataset import RELEASE

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


def publish_comparison(report):
    before = json.loads((ROOT / 'reports/results/filter_sim_gates/stationary_control.json').read_text())
    if before['release'] != RELEASE:
        raise ValueError('zero-speed and underway controls must use the same clock-corrected release')
    families, lines, angles = {}, [], []
    for name, label in {**FAMILIES, 'pii': 'Adaptive PII', 'nlo': 'NLO'}.items():
        control = before['families'][name]
        if control['binary_sha256'] != report['binary_sha256'][name]:
            raise ValueError(f'{name}: stationary and underway executable hashes differ')
        old = control['records']
        if name == 'pii':
            new = report['other_replays']['pii']
        elif name == 'nlo':
            new = [dict(input=n, disp_z_pct_hs=v) for n,v in zip(sorted(REFERENCE_CSV_SHA256), report['other_replays']['nlo_z_percent_hs'])]
        else:
            new = report['records'][name]['default']
        if {r['input'] for r in old} != set(REFERENCE_CSV_SHA256) or {r['input'] for r in new} != set(REFERENCE_CSV_SHA256):
            raise ValueError(f'{name}: unpaired input names')
        worst = lambda rr,k: max(r[k] for r in rr)
        b,a = worst(old,'disp_z_pct_hs'),worst(new,'disp_z_pct_hs')
        row = dict(heave_before_pct_hs=b,heave_after_pct_hs=a,heave_change_percent=100*(a/b-1),
                   after_worst_heave_input=max(new,key=lambda r:r['disp_z_pct_hs'])['input'])
        cells = [label,f'{b:.2f} → {a:.2f}',f'{row["heave_change_percent"]:+.0f}%']
        for key, scale, precision in [('disp_3d_rms_m',1,3),('accel_bias_3d_rms_mps2',1000,2)]:
            if key in old[0] and key in new[0]:
                row[key+'_before'],row[key+'_after'] = worst(old,key),worst(new,key)
                cells.append(f'{scale*row[key+"_before"]:.{precision}f} → {scale*row[key+"_after"]:.{precision}f}')
            else:
                cells.append('—')
        lines.append('| '+' | '.join(cells)+' |')
        if name != 'nlo':
            angle_cells = [label]
            for angle in ('roll','pitch','yaw'):
                key=angle+'_rms_deg'
                row[key+'_before'],row[key+'_after']=worst(old,key),worst(new,key)
                angle_cells.append(f'{row[key+"_before"]:.3f} → {row[key+"_after"]:.3f}')
            angles.append('| '+' | '.join(angle_cells)+' |')
        families[name]=row
    ahrs = before['families']['ahrs']
    if ahrs['binary_sha256'] != report['binary_sha256']['ahrs']:
        raise ValueError('AHRS: stationary and underway executable hashes differ')
    angle_cells, row = ['AHRS (last 60 s)'], {}
    for angle in ('roll','pitch','yaw'):
        b = max(r[angle] for r in ahrs['records'])
        a = max(r[angle] for r in report['other_replays']['ahrs'])
        row[angle+'_rms_deg_before'],row[angle+'_rms_deg_after']=b,a
        angle_cells.append(f'{b:.3f} → {a:.3f}')
    angles.append('| '+' | '.join(angle_cells)+' |')
    families['ahrs']=row
    receipt = dict(before=before['archive'],after=report['inputs'],before_inputs=before['input_sha256'],
                   release=RELEASE,scope='same executable and default draws in both arms; each column is its own worst record; physical input change, not an estimator improvement',
                   binary_sha256=report['binary_sha256'],families=families)
    (ROOT/'reports/results/filter_sim_gates/dataset_comparison.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    text = ['# Stationary versus moving-boat measurements','',
            f'Both arms are fresh replays of oceanography-waves-lib {RELEASE}, with exact 200 Hz clocks. Before is the zero-speed 28 ft vessel archive; after is the 5 kn archive with ±3° heading oscillation every 120 s. Each pair uses the identical rebuilt executable, default sensor and initialization draws, and trailing 900 s scoring. No estimator retuning is applied. Hashes and per-record controls are retained in `stationary_control.json`, `measurements.json`, and `dataset_comparison.json`.','',
            '| Filter | Worst heave before → after (% Hs) | Relative increase | Worst 3D RMS before → after (m) | Accel bias RMS before → after (mm/s²) |',
            '| --- | ---: | ---: | ---: | ---: |',*lines,'',
            'Each column reports its own worst record; these maxima need not occur in the same sea. These deterministic comparisons do not provide population confidence intervals.','',
            '| Filter | Roll before → after (deg) | Pitch before → after (deg) | Yaw before → after (deg) |',
            '| --- | ---: | ---: | ---: |',*angles,'',
            'NLO has the largest relative heave deterioration, followed by OU-III and cubic TFG. OU-II and the CARMA wrappers do not let still-water detection attenuate wave variance, which a 0.27 m following sea at 5 kn triggers. TFG has a substantial accelerometer-bias and attitude penalty; the shared-mode CARMA variants retain nearly unchanged worst absolute 3D RMS despite worse normalized heave in the smallest sea. These paired inputs change encounter frequency and steering together. The [regression diagnosis](moving-vessel-regression-diagnosis.md) explains the stillness fix and the limits of attributing the remaining differences.','',
            'The [engine-noise study](../reports/results/engine_noise_degradation/engine_noise_report.md) separately tests machinery stress on the moving records. At 800 rpm, OU-III with the guard plus covariance inflation reaches 2.417 times its engine-off 3D RMS; the old zero-speed 1.7-times claim does not carry over. Its pitch offset also worsens when covariance inflation is added. This is a measured operating limitation, not a change to deployed tuning.','',
            'Position and velocity truth are expressed about the declared cruise track. All specific force, turning acceleration and angular rates remain intact. See the [coordinate convention](moving-vessel-replays.md) and [current results](moving-vessel-results.md).','']
    (ROOT/'docs/moving-vessel-comparison.md').write_text('\n'.join(text))


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
        'The current external-data baseline uses oceanography-waves-lib v1.3.1: '
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
    publish_comparison(report)


if __name__ == '__main__':
    publish()
