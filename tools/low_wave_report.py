#!/usr/bin/env python3
"""Render the retained low-wave study without selecting parameters or changing gates."""
import gzip
import json
import math
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'reports/results/low_wave_noise'


def main():
    with gzip.open(OUT / 'runs.jsonl.gz', 'rt') as stream:
        rows = [json.loads(line) for line in stream]
    lines = ['# Low-wave filter study', '',
        'Full 1200-second records, final 900 seconds; unchanged sensor injection and quality gates. '
        'Fresh validation uses paired IMU/initialization seeds 31013, 37003, 41011 and 47017 on all eight pinned records. '
        'These are new sensor draws, not independent stochastic seas. See the manifest for binary/source provenance and invalid arms.', '',
        '| Family | Weighting | Violations / 32 | Low X acceleration RMS, m/s² | Low Y acceleration RMS, m/s² | Roll, deg | Pitch, deg | Yaw, deg | 3-D displacement, m |',
        '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for family, study in [('OU-II', 'ou2-snr-fresh'), ('OU-III', 'ou3-snr-fresh')]:
        for config in ['current', 'snr2']:
            rr = [r for r in rows if r['study'] == study and r['config'] == config]
            low = [r for r in rr if 'H0.270' in r['input']]
            values = [sum(len(r['violations']) for r in rr)]
            values += [statistics.fmean(r['metrics'][k] for r in low) for k in
                       ['accel_x_rms_mps2', 'accel_y_rms_mps2']]
            values += [statistics.fmean(r['metrics'][k] for r in rr) for k in
                       ['roll_rms_deg', 'pitch_rms_deg', 'yaw_rms_deg', 'disp_3d_rms_m']]
            lines.append('| ' + ' | '.join([family, 'off' if config == 'current' else 'SNR-dependent',
                         str(values[0]), *[f'{v:.6f}' for v in values[1:]]]) + ' |')
    lines += ['', 'Acceleration columns average only the two lowest-wave records across four draws; '
        'attitude/displacement columns average all 32 cases. OU-II still fails 30/32 cases and OU-III 32/32. '
        'Small pitch/yaw tradeoffs remain. Reduced acceleration noise is not a claim that low-wave attitude or direction is solved.', '',
        'The separate heading/angle validation uses seed 53017: eight rigid world rotations and eight recomputed hull responses '
        'at incident angles −60° and 75°, retaining the prescribed harmonic phases. Weighting leaves violation counts unchanged '
        '(OU-II 73/16 cases; OU-III 65/16), while small attitude tradeoffs persist.', '',
        'Direction bias correction is paired within a single executable. It does not alter upstream attitude/displacement/tuning. '
        'Zero axes are unavailable; conditional RMS cannot be interpreted without availability. At the default low-wave draw, '
        'corrected axis availability is below 1% and travel direction is unresolved more than 99% of the time.', '',
        'Rejected defaults include direct RAO scaling of the OU prior, longer tau alone, cadence-matched tau changes, '
        'globally enlarged accelerometer covariance, asymmetric R_S factors, increased bias priors, and magnetic-weight changes. '
        'OU-III R_S factors 0.5 improve training displacement but add a fresh-draw violation and worsen mean pitch. '
        'Global OU-II horizontal sigma ×4 adds a fresh-draw bias violation. '
        'The SNR policy preserves the original weighting in stronger seas.', '',
        'Yaw screens and the magnetic-error diagnostic show sensitivity to both hard-iron offset and soft-iron distortion; '
        'single-heading low-wave motion does not identify a full magnetic calibration. No new yaw default is retained.', '',
        'Applied OU-III tau/sigma/R_S telemetry now reports the active model. Earlier OU-III CSV fields with those names '
        'reported staged commands; they cannot locate actual covariance synchronization events. '
        'Tuning-runner checks reject variables absent from the executable and record the starting producer/source/binary hashes.', '']
    (OUT / 'study.md').write_text('\n'.join(lines))


if __name__ == '__main__':
    main()
