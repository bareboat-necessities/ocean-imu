#!/usr/bin/env python3
"""Pair legacy and explicit instantaneous contracts on an existing replay."""
import argparse,concurrent.futures,json,os,subprocess
from pathlib import Path
from sim_dataset import ROOT,SHA256,REFERENCE_CSV_SHA256,digest

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--jobs',type=int,default=2);args=parser.parse_args()
 out=ROOT/'reports/results/live_bias_consistency';dest=out/'runs/csv';dest.mkdir(parents=True,exist_ok=True)
 name='wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv';pinned=ROOT/'.sim-data'/SHA256/name
 if digest(pinned)!=REFERENCE_CSV_SHA256[name]:raise ValueError('pinned input hash mismatch')
 def run(case):
  geometry,contract=case
  executable=ROOT/'tests'/('kalman_carma_shared_modes/kalman_carma_shared_modes-sim' if geometry=='qmekf' else 'kalman_tfg_shared_modes/kalman_tfg_shared_modes-sim')
  env=os.environ.copy();env['W3D_WRITE_TIMESERIES']='0'
  if contract=='instantaneous':env['W3D_CSM_GYRO_INPUT']='instantaneous'
  else:env.pop('W3D_CSM_GYRO_INPUT',None)
  result=subprocess.run([str(executable),'--input',str(pinned)],env=env,capture_output=True,text=True)
  (dest/(geometry+'_'+contract+'.log')).write_text(result.stdout+result.stderr)
  if result.returncode:raise RuntimeError(result.stderr)
  required='REGRESSION_PROTOCOL: GATED='+('1' if contract=='interval' else '0')
  if required not in result.stdout:raise RuntimeError('nondefault replay environment: '+required)
  keys=('VALIDATION_METRICS','REPLAY_COMPLETE','REGRESSION_PROTOCOL','CSM_GYRO_INPUT')
  return dict(filter=geometry,contract=contract,exit_code=result.returncode,metrics=[line for line in result.stdout.splitlines() if any(k in line for k in keys)])
 with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
  rows=list(pool.map(run,[(f,c) for f in ('qmekf','tfg') for c in ('interval','instantaneous')]))
 (out/'csv_replay.json').write_text(json.dumps(rows,indent=2)+'\n');print('Completed four paired existing-infrastructure replays.')
if __name__=='__main__':main()
