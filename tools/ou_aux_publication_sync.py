#!/usr/bin/env python3
"""Refresh auxiliary OU article numbers from the retained full-rate studies."""
import csv
import math
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / 'doc/kalman_ou_iii'
RESULTS = ROOT / 'reports/results'


def rows(study, name):
    with (RESULTS / study / name).open() as stream:
        return list(csv.DictReader(stream))


def replace(text, pattern, value):
    result, count = re.subn(pattern, lambda _: value, text, flags=re.M)
    if count != 1:
        raise ValueError(f'expected one article field, found {count}: {pattern}')
    return result


def num(row, key):
    return float(row[key])


def main():
    mismatch = {r['family']: r for r in rows('model_mismatch_ablation', 'model_mismatch_summary.csv')}
    path = DOC / 'w3d-model-mismatch-ablation.tex-part'
    text = path.read_text().replace('equal-length stationary records', 'equal-length v1.3.1 underway records').replace('over eight stationary\n  records', 'over eight underway\n  records')
    for family, row in mismatch.items():
        label = family.replace('OU-', 'OU--')
        values = [f'{num(row,k):.3f}' for k in ('disp_x_rms_m','disp_y_rms_m','disp_z_rms_m','disp_3d_rms_m')]
        values.append(f"{num(row,'disp_z_pct_refrms'):.2f}")
        text = replace(text, '^    ' + re.escape(label) + r' & .*$', '    ' + label + ' & ' + ' & '.join(values) + r' \\')
    values = [(num(mismatch[f], 'disp_z_rms_m'), num(mismatch[f], 'disp_3d_rms_m')) for f in ('OU-III','OU-II','TFG')]
    prose = (f"OU--III vertical and 3-D RMS are ${values[0][0]:.3f}$ and ${values[0][1]:.3f}$~m, versus ${values[1][0]:.3f}$ and\n"
             f"${values[1][1]:.3f}$~m for OU--II and ${values[2][0]:.3f}$ and ${values[2][1]:.3f}$~m for TFG\n")
    text = replace(text, r'OU--III vertical and 3-D RMS are [\s\S]*?(?=\(Fig\.)', prose)
    path.write_text(text)

    engine = {(r['family'],r['label']): r for r in rows('engine_noise_degradation','engine_noise_summary.csv')}
    guard = {(r['condition'],r['guard']): r for r in rows('engine_noise_mitigation','guard_summary.csv')}
    path = DOC / 'w3d-engine-noise-degradation.tex-part'
    text = path.read_text().replace('Archived Zero-Speed Engine-Noise Degradation','Underway Engine-Noise Degradation').replace('eight stationary v1.2.1 vessel-response records','eight v1.3.1 5 kn vessel-response records')
    text = text.replace('unchanged; encounter-frequency and propulsion-induced motion are outside\nthis experiment.', 'unchanged within each paired arm. The reference records include the prescribed\nencounter and steering kinematics; hull-specific propulsion hydrodynamics remain outside this experiment.')
    fractions = []
    for family in ('OU-II','OU-III','TFG'):
        label = family.replace('OU-','OU--')
        for key, arm in [('engine off','off'),('speed 2400 rpm','on')]:
            r = engine[family,key]
            offset = math.sqrt(sum(num(r,k)**2 for k in ('disp_x_mean_m','disp_y_mean_m','disp_z_mean_m')))
            line = f"{label} & {arm} & {num(r,'disp_3d_rms_m'):.3f} & {offset:.3f} & {num(r,'pitch_rms_deg'):.3f} & {num(r,'yaw_rms_deg'):.2f}" + r' \\'
            text = replace(text,'^'+re.escape(label+' & '+arm)+r' & .*$',line)
            if arm == 'on':
                fractions.append(label + rf": \pct{{{100*offset/num(r,'disp_3d_rms_m'):.1f}}}")
    off,on = engine['OU-III','engine off'],engine['OU-III','speed 2400 rpm']
    text = replace(text,r"OU--III's 3-D RMS rises from .*",rf"OU--III's 3-D RMS rises from \SI{{{num(off,'disp_3d_rms_m'):.3f}}}{{m}} to \SI{{{num(on,'disp_3d_rms_m'):.3f}}}{{m}}, a factor of \num{{{num(on,'disp_3d_rms_m')/num(off,'disp_3d_rms_m'):.2f}}}.")
    text = replace(text,r'The 3-D offset fractions are .*','The 3-D offset fractions are '+', '.join(fractions)+'.')
    text = replace(text,r'Pitch RMS and pitch offset are .*',rf"Pitch RMS and pitch offset are \SI{{{num(on,'pitch_rms_deg'):.3f}}}{{\degree}} and \SI{{{abs(num(on,'pitch_mean_deg')):.3f}}}{{\degree}}.")
    delta = max(100*abs(num(engine[f,'accelerometer only'],'disp_3d_rms_m')/num(engine[f,'speed 2400 rpm'],'disp_3d_rms_m')-1) for f in ('OU-II','OU-III','TFG'))
    text = replace(text,r'at most \\pct\{[^}]+\} across the families\.',rf'at most \pct{{{delta:.3f}}} across the families.')
    contrasts = ', '.join(f"{f.replace('OU-','OU--')}: {100*(num(engine[f,'accelerometer only'],'disp_3d_rms_m')/num(engine[f,'speed 2400 rpm'],'disp_3d_rms_m')-1):+.2f}\\%" for f in ('OU-II','OU-III','TFG'))
    text = replace(text, r'Removing the modeled gyroscope perturbation[\s\S]*?(?=\\subsubsection\{Speed)',
                   'Removing the modeled gyroscope perturbation changes nominal 3-D RMS by\n'
                   + contrasts + '. These paired accelerometer-only contrasts change the nonlinear\n'
                   'estimator trajectory; they do not establish that the gyro path is negligible\n'
                   'or provide an additive causal decomposition of the total error.\n\n')
    for f in ('OU-II','OU-III'):
        spreads = []
        for prefix in ('bandwidth','matched'):
            v = [num(engine[f,k],'disp_3d_rms_m') for k in (prefix+' 20 Hz',prefix+' 40 Hz','speed 2400 rpm',prefix+' 160 Hz')]
            spreads.append(max(v)/min(v))
        label = f.replace('OU-','OU--')
        text = replace(text,'^ For '+label+r', the maximum/minimum.*',rf' For {label}, the maximum/minimum 3-D spread falls from \num{{{spreads[0]:.2f}}} to \num{{{spreads[1]:.2f}}}.')
    base = num(guard['engine off','off'],'disp_3d_rms_m')
    labels = {'engine off':'Engine off','800 rpm':r'\SI{800}{rpm}','1600 rpm':r'\SI{1600}{rpm}','2400 rpm':r'\SI{2400}{rpm}','3200 rpm':r'\SI{3200}{rpm}','2400 rpm, quiet mount':'Quiet mount','2400 rpm, engine bed':'Engine bed','2400 rpm, wide sensor':'Wide sensor'}
    for condition,label in labels.items():
        arms = [guard[condition,a] for a in ('off','guard','guard+R')]
        v = [num(r,'disp_3d_rms_m') for r in arms]+[abs(num(r,'pitch_mean_deg')) for r in arms]+[num(arms[-1],'disp_3d_rms_m')/base]
        text = replace(text,'^'+re.escape(label)+r' & .*',label+' & '+' & '.join(f'{x:.3f}' for x in v)+r' \\')
    g,b = guard['2400 rpm','guard'],guard['2400 rpm','guard+R']
    text = replace(text,r'At nominal cruise, the guard reduces .*',rf"At nominal cruise, the guard reduces OU--III 3-D RMS to \SI{{{num(g,'disp_3d_rms_m'):.3f}}}{{m}}; adding covariance inflation reduces it to \SI{{{num(b,'disp_3d_rms_m'):.3f}}}{{m}}, still {100*(num(b,'disp_3d_rms_m')/base-1):.1f}\% above the \SI{{{base:.3f}}}{{m}} engine-off baseline. Improvement is not uniform:")
    g,b = guard['800 rpm','guard'],guard['800 rpm','guard+R']
    text = replace(text,r'At \\SI\{800\}\{rpm\}, pitch offset rises .*',rf"At \SI{{800}}{{rpm}}, pitch offset rises from \SI{{{abs(num(g,'pitch_mean_deg')):.3f}}}{{\degree}} to \SI{{{abs(num(b,'pitch_mean_deg')):.3f}}}{{\degree}} with covariance inflation. Its 3-D RMS is \num{{{num(b,'disp_3d_rms_m')/base:.3f}}} times the engine-off baseline.")
    text = text.replace('These fixed-heading, zero-speed RAO replays are sensor-path\nstress studies','These prescribed moving-vessel RAO replays are sensor-path\nstress studies')
    path.write_text(text)
    for study in ('model_mismatch_ablation','engine_noise_degradation','engine_noise_mitigation','ou_rs_law','ou3_lever_arm_study'):
        for svg in (RESULTS/study).glob('*.svg'):
            shutil.copy2(svg,DOC/svg.name)
    shutil.copy2(RESULTS/'ou_rs_law/ou_rs_roundtrip_scores.tex',DOC/'w3d-roundtrip-transition-scores-generated.tex-part')


if __name__ == '__main__':
    main()
