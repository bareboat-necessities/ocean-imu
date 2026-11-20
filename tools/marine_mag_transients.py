#!/usr/bin/env python3
"""Score clock-aligned yaw traces; truth is read only by this evaluator."""
import argparse
import os
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from marine_mag_study import BINARY, OUT, ROOT, write


def wrap(x):
    return (x+180)%360-180


def case_scores(case):
    source=next((ROOT/'.sim-data/cases'/case).glob('wave_data_jonswap_H1.500_*.csv'))
    truth=pd.read_csv(source,usecols=['time','roll_deg','pitch_deg'])
    t=truth.time.to_numpy()
    window=t>=300
    signals=[np.sin(2*np.pi*t/20),np.cos(2*np.pi*t/20),truth.roll_deg,truth.pitch_deg]
    x=np.column_stack([np.ones(len(t)),t/1200,*signals])[window]
    rows=[]
    with tempfile.TemporaryDirectory() as tmp:
        def run(key):
            v,s=key
            file=Path(tmp)/f'{v}-{s}.bin'
            env={k:vv for k,vv in os.environ.items() if not k.startswith(('W3D_','OCEAN_IMU_SIM_CASE'))}
            env.update(OCEAN_IMU_SIM_CASE=case,MARINE_MAG_YAW_TRACE=str(file))
            subprocess.run([BINARY,source,str(v),s,'66.5'],env=env,stdout=subprocess.DEVNULL,check=True)
            y=np.fromfile(file,dtype=np.float32).astype(float)
            assert len(y)==len(t) and np.all(np.isfinite(y))
            return key,y
        with ThreadPoolExecutor(2) as pool:
            traces=dict(pool.map(run,[(v,s) for v in [-1,0,6,7,8] for s in ['combined','burst']]))
        for (v,s),y in traces.items():
            b=np.linalg.lstsq(x,y[window],rcond=None)[0]
            r=y[window]-x@b
            difference=wrap(y-traces[-1,s])
            burst_delta=wrap(y-traces[v,'combined'])
            after=abs(burst_delta[(t>=610)&(t<1200)]).reshape(-1,200).max(axis=1)
            good=np.convolve((after<=.1).astype(int),np.ones(10,dtype=int),'valid')
            hits=np.flatnonzero(good==10)
            rows.append(dict(case=case,scenario=s,variant=v,early_mean_deg=y[(t>=300)&(t<400)].mean(),
                late_mean_deg=y[t>=1100].mean(),trend_deg_per_min=b[1]/20,
                steering_amplitude_deg=np.hypot(b[2],b[3]),roll_gain=b[4],pitch_gain=b[5],
                partial_wave_rms_deg=np.std(x[:,4:]@b[4:]),residual_centered_deg=np.std(r),
                paired_max_increment_deg=max(abs(wrap(np.diff(difference[window])))),
                burst_peak_deg=max(abs(burst_delta[(t>=600)&(t<610)])),
                recovery_s=int(hits[0]) if len(hits) else -1))
    return rows


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=OUT)
    args=p.parse_args()
    write(args.out/'transients.csv',[r for case in ('following','head') for r in case_scores(case)])
