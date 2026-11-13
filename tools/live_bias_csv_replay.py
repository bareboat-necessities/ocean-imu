#!/usr/bin/env python3
"""Pair legacy and explicit instantaneous contracts on an existing replay."""
import argparse,concurrent.futures,hashlib,json,os,subprocess
from sim_dataset import ROOT,SHA256,REFERENCE_CSV_SHA256,digest

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--jobs',type=int,default=2);parser.add_argument('--all-records',action='store_true');parser.add_argument('--localize',action='store_true');args=parser.parse_args()
 out=ROOT/'reports/results/live_bias_consistency';dest=out/'runs/csv';dest.mkdir(parents=True,exist_ok=True)
 names=list(REFERENCE_CSV_SHA256) if args.all_records else ['wave_data_jonswap_H1.500_L50.710_A-30.00_P120.00.csv']
 for name in names:
  if digest(ROOT/'.sim-data'/SHA256/name)!=REFERENCE_CSV_SHA256[name]:raise ValueError('pinned input hash mismatch')
 def run(case):
  geometry,contract,name=case;pinned=ROOT/'.sim-data'/SHA256/name
  executable=ROOT/'tests'/('kalman_carma_shared_modes/kalman_carma_shared_modes-sim' if geometry=='qmekf' else 'kalman_tfg_shared_modes/kalman_tfg_shared_modes-sim')
  env=os.environ.copy();env['W3D_WRITE_TIMESERIES']='0'
  if contract not in ('interval','defaults'):env['W3D_CSM_GYRO_INPUT']=contract
  else:env.pop('W3D_CSM_GYRO_INPUT',None)
  command=[str(executable),'--input',str(pinned)]
  if contract=='interval':command+=['--regression']
  result=subprocess.run(command,env=env,capture_output=True,text=True)
  (dest/(geometry+'_'+contract+'_'+pinned.stem+'.log')).write_text(result.stdout+result.stderr)
  if result.returncode not in (0,1):raise RuntimeError(result.stdout+result.stderr)
  required='REGRESSION_PROTOCOL: GATED='+('1' if contract in ('interval','defaults') else '0')
  if required not in result.stdout:raise RuntimeError('nondefault replay environment: '+required)
  keys=('VALIDATION_METRICS','REPLAY_COMPLETE','REGRESSION_PROTOCOL','CSM_GYRO_INPUT','ERROR:','QUALITY_GATE:','Bias error RMS','Max TRUE bias')
  return dict(filter=geometry,contract=contract,input=name,input_sha256=REFERENCE_CSV_SHA256[name],binary_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),exit_code=result.returncode,metrics=[line for line in result.stdout.splitlines() if any(k in line for k in keys)])
 with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
  contracts=['interval','instantaneous','defaults']+(['core','frontend','both_schmidt','core_schmidt'] if args.localize else [])
  rows=list(pool.map(run,[(f,c,n) for f in ('qmekf','tfg') for c in contracts for n in names]))
 (out/('csv_rollout_all.json' if args.all_records else 'csv_replay.json')).write_text(json.dumps(rows,indent=2)+'\n')
 failed=sum(r['exit_code']!=0 for r in rows);print(f'Completed {len(rows)} paired existing-infrastructure replays; {failed} numeric failures.')
 if failed:raise SystemExit(1)
if __name__=='__main__':main()
