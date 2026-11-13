#!/usr/bin/env python3
"""Paired continuous histories for the OU-II/OU-III/TFG default rollout."""
import argparse,concurrent.futures,hashlib,json,shutil,subprocess
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'reports/results/live_bias_consistency'
def main():
 p=argparse.ArgumentParser();p.add_argument('--jobs',type=int,default=3);p.add_argument('--duration',type=float,default=7200);p.add_argument('--seeds',type=int,nargs='+',default=list(range(1001,1007)));a=p.parse_args()
 runroot=OUT/'runs/main_rollout';runroot.mkdir(parents=True,exist_ok=True)
 binary=runroot/'estimator';shutil.copy2(ROOT/'tests/live_bias_consistency/main-gyro-rollout',binary)
 digest=hashlib.sha256(binary.read_bytes()).hexdigest();cases=[(f,s) for s in a.seeds for f in ['OU2','OU3','TFG']]
 def run(case):
  family,seed=case;prefix=runroot/f'{family}_{seed}_{a.duration:g}'
  r=subprocess.run([str(binary),family,str(seed),str(a.duration),str(prefix)],capture_output=True,text=True)
  prefix.with_suffix('.log').write_text(r.stdout+r.stderr)
  if r.returncode:raise RuntimeError(r.stdout+r.stderr)
  rows=pd.read_csv(str(prefix)+'_scores.csv',dtype={'sensor_hash':str});assert len(rows)==2 and rows.sensor_hash.nunique()==1
  rows['family']=family;return rows
 rows=[]
 with concurrent.futures.ThreadPoolExecutor(a.jobs) as pool:
  for i,x in enumerate(pool.map(run,cases)):
   rows.append(x);pd.concat(rows,ignore_index=True).to_csv(OUT/'main_rollout_scores.csv',index=False,float_format='%.12g');print(i+1,len(cases),flush=True)
 d=pd.concat(rows,ignore_index=True)
 assert (d.groupby('seed').sensor_hash.nunique()==1).all(),'different physical histories'
 manifest=dict(binary_sha256=digest,pairs=len(cases),filter_runs=2*len(cases),seeds=a.seeds,duration=a.duration,paired=True,ordinary_learned_magnetic_reference=True,fixed_zero_hard_iron=True,command=' '.join(__import__('sys').argv))
 (OUT/'main_rollout_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':main()
