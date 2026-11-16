#!/usr/bin/env python3
"""Gate full-rate heave comparisons and retain compact, case-specific traces."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd

from sim_dataset import CASE, DESCRIPTION, ROOT, REFERENCE_CSV_SHA256, input_provenance
from filter_sim_gates import REPLAYS, profiles_for, parse_w3d, parse_heave

FAMILIES = ('ou2', 'ou3', 'tfg', 'heave')
COLUMNS = ('time', 'disp_ref_z', 'disp_est_z', 'vel_ref_z', 'vel_est_z')


def retain_tail(source, destination):
    """Keep the final 60 s at original rate; plotted decimation happens later."""
    tail = pd.DataFrame(columns=COLUMNS)
    count = 0
    with pd.read_csv(source, usecols=COLUMNS, chunksize=24000) as chunks:
        for chunk in chunks:
            if not np.isfinite(chunk.to_numpy()).all():
                raise ValueError(f'nonfinite heave trace: {source}')
            count += len(chunk)
            tail = chunk if tail.empty else pd.concat((tail, chunk), ignore_index=True)
            tail = tail.tail(12002)
    if count != 240000:
        raise ValueError(f'incomplete heave trace: {source}: {count} rows')
    tail = tail[tail.time >= tail.time.iloc[-1]-60]
    destination.parent.mkdir(parents=True, exist_ok=True)
    tail.to_csv(destination, index=False, lineterminator='\n')
    source.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    data, out = args.data_dir.resolve(), args.output_dir.resolve()
    inputs = input_provenance(data/n for n in REFERENCE_CSV_SHA256)
    for name in ('work', 'traces', 'charts'):
        if (out/name).exists():
            shutil.rmtree(out/name)
    env = {k:v for k,v in os.environ.items()
           if not k.startswith(('W3D_', 'SF_', 'OU_', 'HO_', 'NLO_', 'TFG_'))}
    env.update(W3D_WRITE_TIMESERIES='1', W3D_VALIDATION_WINDOW_SEC='900', W3D_COLLECT_ALL_GATES='1')

    def run(family):
        folder, binary, _ = REPLAYS[family]
        work = out/'work'/folder
        work.mkdir(parents=True, exist_ok=True)
        profiles = profiles_for('heave') if family == 'heave' else [('default', [])]
        rows = []
        for profile, flags in profiles:
            log = out/f'{family}-{profile}.log'
            with log.open('w') as stream:
                for record in REFERENCE_CSV_SHA256:
                    # One record at a time bounds temporary CSV storage.
                    link = work/record
                    link.symlink_to(data/record)
                    done = subprocess.run([str(ROOT/'tests'/folder/binary), *flags,
                                           '--input', record], cwd=work, env=env,
                                          stdout=stream, stderr=subprocess.STDOUT)
                    link.unlink()
                    if done.returncode:
                        raise RuntimeError(f'{family}/{profile}/{record}: exit {done.returncode}; see {log}')
                    traces = list(work.glob('w3d_*.csv'))
                    if not traces:
                        raise ValueError(f'{family}/{profile}/{record}: no trace')
                    for trace in traces:
                        retain_tail(trace, out/'traces'/folder/trace.name)
            text = log.read_text()
            if 'QUALITY_GATE: PASS=0' in text:
                raise ValueError(f'{family}/{profile}: failed accuracy gate')
            if family == 'heave':
                rows.extend(dict(family='heave_'+r['method'], profile=profile,
                                 input=Path(r['record']).name, disp_z_pct_hs=float(r['z_pct_hs']))
                            for r in parse_heave(text, profile))
            else:
                rows.extend(dict(family=family, profile=profile, **r) for r in parse_w3d(text))
        return rows

    with ThreadPoolExecutor(args.jobs) as pool:
        rows = [r for group in pool.map(run, FAMILIES) for r in group]
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (out/'heave_comparison_metrics.csv').open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    for frontend in ('mahony', 'proxy', 'truth'):
        subprocess.run([sys.executable, str(ROOT/'plots/heave_only/heave_only-plots.py'),
                        '--tests-dir', str(out/'traces'), '--output-dir', str(out/'charts'),
                        '--frontend', frontend, '--case-label', DESCRIPTION], env=env, check=True)
    charts = list((out/'charts').glob('*.svg'))
    if len(rows) != 136 or len(charts) != 102:
        raise ValueError(f'incomplete heave comparison: {len(rows)} rows, {len(charts)} charts')
    (out/'comparison.json').write_text(json.dumps(dict(case=CASE, inputs=inputs,
        acceptance='all 136 full-rate record gates passed', charts=len(charts),
        trace_window_s=60, trace_rate_hz=200, plotted_rate_hz=50), indent=2)+'\n')


if __name__ == '__main__':
    main()
