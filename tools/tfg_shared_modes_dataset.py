#!/usr/bin/env python3
"""Paired sensor/initialization-seed replay on the eight pinned vessel records.

The wave realization is fixed by each published record. Unlike the primary
synthetic study, these intervals are conditional on that fixed wave history.
"""
import argparse,json,subprocess,tempfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import pandas as pd
from sim_dataset import input_provenance,REFERENCE_CSV_SHA256
from tfg_shared_modes_study import ROOT,OUT,clean_env,ci,sha
import sys
sys.path.insert(0,str(ROOT/'plots'));import ou_validation as ov
BIN={
 'Q independent':ROOT/'tests/kalman_carma/carma3-sim',
 'T independent':ROOT/'tests/kalman_carma/carma3-sim',
 'Q shared':ROOT/'tests/kalman_carma_shared_modes/kalman_carma_shared_modes-sim',
 'T shared':ROOT/'tests/kalman_tfg_shared_modes/kalman_tfg_shared_modes-sim'}
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--data-dir',type=Path,default=ROOT/'tests/kalman_carma_shared_modes');ap.add_argument('--seeds',type=int,default=8);ap.add_argument('--jobs',type=int,default=2)
ap.add_argument('--results',type=Path,default=OUT);args=ap.parse_args();OUT=args.results;OUT.mkdir(parents=True,exist_ok=True)
inputs=[args.data_dir/n for n in REFERENCE_CSV_SHA256]
protocol=dict(inputs=input_provenance(inputs),seeds=list(range(1001,1001+args.seeds)),
 binaries={k:sha(v) for k,v in BIN.items()},runner_sha256=sha(__file__),sample_rate_hz=200,score_window_s='last 900 seconds',wave_excitation='Evidence on all four wrappers',role='fixed-record conditional paired sensor and initialization resampling')
p=OUT/'dataset_protocol.json'
if p.exists() and json.loads(p.read_text())!=json.loads(json.dumps(protocol)):raise RuntimeError('Dataset replay protocol changed')
p.write_text(json.dumps(protocol,indent=2)+'\n');cache=OUT/'dataset_cache';cache.mkdir(exist_ok=True)
def run(task):
 """Replay one pinned wave record with a common sensor/initialization seed."""
 path,seed,method=task;key=f'{path.stem}_{seed}_{method.replace(" ","_")}'
 saved=cache/(key+'.json')
 if saved.exists():return json.loads(saved.read_text())
 env=clean_env();env.update(W3D_IMU_SEED=str(seed),W3D_INIT_SEED=str(seed+10000),W3D_WRITE_TIMESERIES='0',W3D_VALIDATION_WINDOW_SEC='900',W3D_COLLECT_ALL_GATES='1',W3D_CARMA3_EXCITATION='C')
 if method=='T independent':env['W3D_CARMA3_VARIANT']='T41'
 with tempfile.TemporaryDirectory(prefix='tfg-shared-replay-') as td:
  cmd=[str(BIN[method]),'--input',str(path.resolve())]
  done=subprocess.run(cmd,cwd=td,env=env,text=True,capture_output=True)
  if done.returncode not in (0,1):raise RuntimeError(done.stdout[-1500:]+done.stderr[-1500:])
  metrics=ov.parse_validation_metrics(done.stdout)
  rec=dict(record=path.name,seed=seed,method=method,**metrics)
  saved.write_text(json.dumps(rec)+'\n');return rec
tasks=[(p,s,m) for p in inputs for s in protocol['seeds'] for m in BIN]
rows=[]
with ThreadPoolExecutor(max_workers=args.jobs) as pool:
 futures=[pool.submit(run,t) for t in tasks]
 for i,f in enumerate(as_completed(futures),1):rows.append(f.result());print(f'{i}/{len(tasks)}',flush=True)
df=pd.DataFrame(rows).sort_values(['record','seed','method']).reset_index(drop=True)
df.to_csv(OUT/'dataset_runs.csv',index=False,float_format='%.9g');paired=[]
for record,g in df.groupby('record'):
 a=g[g.method=='T shared'].set_index('seed');b=g[g.method=='Q shared'].set_index('seed')
 for metric in ['disp_x_rms_m','disp_y_rms_m','disp_z_rms_m','disp_3d_rms_m','roll_rms_deg','pitch_rms_deg','yaw_rms_deg','accel_bias_3d_rms_mps2','gyro_bias_3d_rms_radps']:
  n,d,lo,hi,pv=ci(100*(a[metric]/b[metric]-1));paired.append(dict(record=record,metric=metric,n=n,percent=d,ci_low=lo,ci_high=hi,p=pv))
pd.DataFrame(paired).to_csv(OUT/'dataset_paired.csv',index=False,float_format='%.9g')
print('Pinned dataset replay complete',flush=True)
