#!/usr/bin/env python3
"""Audit the actual pinned CSV clock and gyro sampling representation."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.spatial.transform import Rotation
from sim_dataset import ROOT,SHA256,REFERENCE_CSV_SHA256,digest

def main():
 rows=[];base=ROOT/'.sim-data'/SHA256
 for name,expected_hash in sorted(REFERENCE_CSV_SHA256.items()):
  path=base/name
  if digest(path)!=expected_hash:raise ValueError('input hash mismatch: '+name)
  columns=['time','gyro_x','gyro_y','gyro_z','q_wb_zu_w','q_wb_zu_x','q_wb_zu_y','q_wb_zu_z']
  d=pd.read_csv(path,usecols=columns);time=d.time.to_numpy();rate=d[['gyro_x','gyro_y','gyro_z']].to_numpy()
  wb=Rotation.from_quat(d[['q_wb_zu_x','q_wb_zu_y','q_wb_zu_z','q_wb_zu_w']].to_numpy())
  # CSV world->body Z-up convention. The repository's proper basis map to
  # NED preserves the norm of this body angular increment and gyro residual.
  interval=(wb[:-1]*wb[1:].inv()).as_rotvec()/.005
  rows.append(dict(file=name,sha256=expected_hash,n=len(time),grid_max_error=float(np.max(abs(time-np.arange(len(time))*.005))),min_interval=float(np.diff(time).min()),max_interval=float(np.diff(time).max()),gyro_endpoint_increment_rms=float(np.sqrt(np.mean(np.sum((rate[1:]-interval)**2,axis=1)))),gyro_centered_increment_rms=float(np.sqrt(np.mean(np.sum((.5*(rate[1:]+rate[:-1])-interval)**2,axis=1))))))
 out=ROOT/'reports/results/live_bias_consistency/csv_timing.json';out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(rows,indent=2)+'\n')
 print('Audited',len(rows),'pinned CSV records; gyro residuals are rad/s.')
if __name__=='__main__':main()
